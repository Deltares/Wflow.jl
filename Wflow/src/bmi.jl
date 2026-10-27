# Basic Model Interface (BMI) implementation based on
# https://github.com/Deltares/BasicModelInterface.jl

# Mapping of grid identifier to a key, to get the active indices of the model domain.
# See also function active_indices(network, key::AbstractString).
const GRIDS = Dict{Int, String}(
    0 => "reservoir",
    1 => "drain",
    2 => "river",
    3 => "land",
    4 => "land",
    5 => "land",
)

"""
    BMI.initialize(::Type{<:Wflow.Model}, config_file)::Model

Initialize the model. Reads the input settings and data as defined in the Config object
generated from the configuration file `config_file`. Will return a Model that is ready to
run.
"""
function BMI.initialize(::Type{<:Model}, config_file::AbstractString)::Model
    config = Config(config_file)
    model_type = config.model.type

    @info "Initialize model variables for model type `$model_type`."

    type = if model_type == ModelType.sbm
        SbmModel()
    elseif model_type == ModelType.sbm_gwf
        SbmGwfModel()
    elseif model_type == ModelType.sediment
        SedimentModel()
    end
    model = Model(config, type)
    load_fixed_forcing!(model)

    return model
end

"""
    BMI.update(model::Model)::Nothing

Update the model for a single timestep.
"""
function BMI.update(model::Model)::Nothing
    run_timestep!(model)
    return nothing
end

function BMI.update_until(model::Model, time::Float64)::Nothing
    t = BMI.get_current_time(model)
    _div, _rem = divrem(time - t, model.clock.dt.value)
    steps = Int(_div)
    if steps < 0
        error("The current model timestamp $t is larger than provided `time` $time")
    elseif abs(_rem) > eps()
        error_message = string(
            "Provided `time` $time minus the current model timestamp $t",
            " is not an integer multiple of model time step $(model.clock.dt.value)",
        )
        error(error_message)
    end
    for _ in 1:steps
        run_timestep!(model)
    end
    return nothing
end

"Write state output to netCDF and close files."
function BMI.finalize(model::Model)::Nothing
    (; config, writer) = model
    (; endstate_writer) = writer
    # it is possible that the state dataset has been closed by `save_state`
    if !isnothing(endstate_writer.output_dataset) && isopen(endstate_writer.output_dataset)
        write_netcdf_timestep(model, writer.endstate_writer)
    end
    reset_clock!(model.clock, config)
    close_files(model; delete_output = false)
    return nothing
end

function BMI.get_component_name(model::Model)::String
    return string(model.config.model.type)
end

function BMI.get_input_item_count(model::Model)::Int
    return length(BMI.get_input_var_names(model))
end

function BMI.get_output_item_count(model::Model)::Int
    return length(BMI.get_output_var_names(model))
end

"""
    BMI.get_input_var_names(model::Model)::Vector{String}

Returns model input variables, based on the `API` section in the model configuration file.
This `API` sections contains a list of `Model` components for which variables can be
exchanged.
"""
function BMI.get_input_var_names(model::Model)::Vector{String}
    (; config, land) = model
    if do_api(config)
        var_names = config.API.variables
        idx = Int[]
        for (i, var) in enumerate(var_names)
            var = get_standard_name(var)
            if isnothing(get_metadata(var, land; model))
                push!(idx, i)
                @warn(
                    "$var is not listed as variable for BMI exchange and removed from list"
                )
            end
        end
        deleteat!(var_names, idx)
        return var_names
    else
        @warn("TOML file does not contain section [API] to extract model var names")
        return String[]
    end
end

"Returns input variables from `BMI.get_input_var_names(model::Model)`, there is no
distinction between input - and output variables."
function BMI.get_output_var_names(model::Model)::Vector{String}
    return BMI.get_input_var_names(model)
end

function BMI.get_var_grid(::Model, name::String)::Int
    return if occursin("reservoir", name)
        0
    elseif occursin("drain", name)
        1
    elseif occursin("river", name) || occursin("floodplain", name)
        2
    elseif occursin("land_surface_water__x_component", name) # LocalInertialOverlandFlowModel
        3
    elseif occursin("land_surface_water__y_component", name) # LocalInertialOverlandFlowModel
        4
    else
        5
    end
end

function BMI.get_var_type(model::Model, name::String)::String
    value = BMI.get_value_ptr(model, name)
    return repr(eltype(first(value)))
end

function BMI.get_var_units(model::Model, name::String)::String
    (; land) = model
    name = get_standard_name(name)
    metadata = get_metadata(name, land; model)
    return to_string(to_SI(metadata.unit); BMI_standard = true)
end

function BMI.get_var_itemsize(model::Model, name::String)::Int
    value = BMI.get_value_ptr(model, name)
    return sizeof(eltype(first(value)))
end

function BMI.get_var_nbytes(model::Model, name::String)::Int
    return sizeof(BMI.get_value_ptr(model, name))
end

function BMI.get_var_location(model::Model, name::String)::String
    name = get_standard_name(name)
    (; lens) = get_metadata(name; model)
    element_type = grid_element_type(model, lens)
    return element_type
end

function BMI.get_current_time(model::Model)::Float64
    (; config, clock) = model
    (; starttime, calendar) = config.time
    starttime = cftime(starttime, calendar)
    # ms => s
    return to_SI(Float64(Dates.value(clock.time - starttime)), MS)
end

function BMI.get_start_time(::Model)::Float64
    return 0.0
end

function BMI.get_end_time(model::Model)::Float64
    (; starttime, endtime, calendar) = model.config.time
    starttime_ = cftime(starttime, calendar)
    endtime_ = cftime(endtime, calendar)
    # ms => s
    return to_SI(Float64(Dates.value(endtime_ - starttime_)), MS)
end

function BMI.get_time_units(model::Model)::String
    return "s"
end

function BMI.get_time_step(model::Model)::Float64
    return model.config.time.timestepsecs
end

function BMI.get_value(model::Model, name::String, dest::Vector{Float64})::Vector{Float64}
    dest .= copy(BMI.get_value_ptr(model, name))
    return dest
end

"""
    BMI.get_value_ptr(model::Model, name::String)::SubArray

Return a reference to the values of model variable `name`, limited to the active elements
of its grid.

Unlike the `DenseVector` return type in the BasicModelInterface.jl specification, this
always returns a view (`SubArray`) into the model's internal array. For most variables the
view is contiguous and covers the whole internal array, but it is required for:
- variables with an internal array that is longer than the grid size, e.g. the land x- and
  y-component of flow for local inertial routing, that include a boundary point;
- a single layer of a layered soil variable, stored as a vector of `SVector`s, which is a
  strided (non-contiguous) view.

Changes to the returned values modify the model state.
"""
function BMI.get_value_ptr(model::Model, name::String)::SubArray
    (; domain) = model
    n = length(active_indices(domain, name))

    if startswith(name, "soil_layer_") && occursin(r"soil_layer_\d+_", name)
        name_2d, ind = soil_layer_standard_name(name)
        model_vals, _ = get_field_in_model(model, name_2d)
        el_type = eltype(first(model_vals))
        dim = length(first(model_vals))
        value = reshape(reinterpret(el_type, model_vals), dim, :)
        return @view value[ind, 1:n]
    else
        (; lens) = get_metadata(name; model)
        if isnothing(lens)
            error("Accessing '$name' is not supported.")
        else
            vec = lens(model)
            return @view(vec[1:n])
        end
    end
end

function BMI.get_value_at_indices(
        model::Model,
        name::String,
        dest::Vector{Float64},
        inds::Vector{Int},
    )::Vector{Float64}
    dest .= BMI.get_value_ptr(model, name)[inds]
    return dest
end

"""
    BMI.set_value(model::Model, name::String, src::Vector{Float64})::Nothing

Set a model variable `name` to the values in vector `src`, overwriting the current contents.
The type and size of `src` must match the model's internal array.
"""
function BMI.set_value(model::Model, name::String, src::Vector{Float64})::Nothing
    BMI.get_value_ptr(model, name) .= src
    return nothing
end

"""
    BMI.set_value_at_indices(model::Model, name::String, inds::Vector{Int}, src::Vector{Float64})::Nothing

Set a model variable `name` to the values in vector `src`, at indices `inds`.
"""
function BMI.set_value_at_indices(
        model::Model,
        name::String,
        inds::Vector{Int},
        src::Vector{Float64},
    )::Nothing
    BMI.get_value_ptr(model, name)[inds] .= src
    return nothing
end

function BMI.get_grid_type(model::Model, grid::Int)::String
    if grid in 0:1
        return "points"
    elseif grid in 2:5
        return "unstructured"
    else
        error("unknown grid type $grid")
    end
end

function BMI.get_grid_rank(model::Model, grid::Int)::Int
    if grid in 0:5
        return 2
    else
        error("unknown grid type $grid")
    end
end

function BMI.get_grid_x(model::Model, grid::Int, x::Vector{Float64})::Vector{Float64}
    (; reader, domain) = model
    (; dataset) = reader
    sel = active_indices(domain, GRIDS[grid])
    inds = [sel[i][1] for i in eachindex(sel)]
    x_nc = read_x_axis(dataset)
    x .= x_nc[inds]
    return x
end

function BMI.get_grid_y(model::Model, grid::Int, y::Vector{Float64})::Vector{Float64}
    (; reader, domain) = model
    (; dataset) = reader
    sel = active_indices(domain, GRIDS[grid])
    inds = [sel[i][2] for i in eachindex(sel)]
    y_nc = read_y_axis(dataset)
    y .= y_nc[inds]
    return y
end

function BMI.get_grid_node_count(model::Model, grid::Int)::Int
    return length(active_indices(model.domain, GRIDS[grid]))
end

function BMI.get_grid_size(model::Model, grid::Int)::Int
    return length(active_indices(model.domain, GRIDS[grid]))
end

function BMI.get_grid_edge_count(model::Model, grid::Int)::Int
    (; domain) = model
    if grid == 2
        return ne(domain.river.network.graph)
    elseif grid == 3
        return length(domain.land.network.edge_indices.ind_x_up)
    elseif grid == 4
        return length(domain.land.network.edge_indices.ind_y_up)
    elseif grid in (0, 1, 5)
        @warn("edges are not provided for grid type $grid (variables are located at nodes)")
        return 0
    else
        error("unknown grid type $grid")
    end
end

function BMI.get_grid_edge_nodes(model::Model, grid::Int, edge_nodes::Vector{Int})::Vector{Int}
    (; domain) = model
    n = length(edge_nodes)
    m = div(n, 2)
    # inactive nodes (boundary/ghost points) are set at -999
    if grid == 2
        nodes_at_edge = adjacent_nodes_at_edge(domain.river.network.graph)
        nodes_at_edge.dst[nodes_at_edge.dst .== m + 1] .= -999
        edge_nodes[range(1, n; step = 2)] = nodes_at_edge.src
        edge_nodes[range(2, n; step = 2)] = nodes_at_edge.dst
        return edge_nodes
    elseif grid == 3
        ind_x_up = copy(domain.land.network.edge_indices.ind_x_up)
        edge_nodes[range(1, n; step = 2)] = 1:m
        ind_x_up[ind_x_up .== m + 1] .= -999
        edge_nodes[range(2, n; step = 2)] = ind_x_up
        return edge_nodes
    elseif grid == 4
        ind_y_up = copy(domain.land.network.edge_indices.ind_y_up)
        edge_nodes[range(1, n; step = 2)] = 1:m
        ind_y_up[ind_y_up .== m + 1] .= -999
        edge_nodes[range(2, n; step = 2)] = ind_y_up
        return edge_nodes
    elseif grid in (0, 1, 5)
        @warn("edges are not provided for grid type $grid (variables are located at nodes)")
        return edge_nodes
    else
        error("unknown grid type $grid")
    end
end

# Extension of BMI functions (state handling and start time), required for OpenDA coupling.
# May also be useful for other external software packages.
function load_state(model::Model)::Nothing
    set_states!(model)
    return nothing
end

function save_state(model::Model)::Nothing
    (; endstate_writer) = model.writer
    (; output_path, output_dataset) = endstate_writer
    if !isnothing(output_path)
        @info "Write output states to netCDF file `$output_path`."
    end
    write_netcdf_timestep(model, endstate_writer)
    close(output_dataset)
    return nothing
end

function get_start_unix_time(model::Model)::Float64
    return datetime2unix(model.config.time.starttime)
end

# BMI helper functions.
"""
Return the standard name `name_layered` representing a layered soil model variable (vector
of svectors) and the layer index `layer_index` based on a standard `name` representing a
layer of a layered soil model variable.
"""
function soil_layer_standard_name(name::AbstractString)
    # Parse new naming format: soil_layer_N_water_... where N is the layer number
    parts = split(name, "_")
    if length(parts) >= 3 && parts[1] == "soil" && parts[2] == "layer"
        layer_number = parts[3]
        layer_index = tryparse(Int, layer_number)
        if !isnothing(layer_index)
            # Remove the layer number to get the base name
            name_layered = join([parts[1], parts[2], parts[4:end]...], "_")
            return name_layered, layer_index
        end
    end
    # Fallback for unexpected format
    @warn "Unable to parse layer standard name: $name"
    return name, nothing
end

"""
    grid_element_type(model, lens::ComposedFunction)
    grid_element_type(::T, var::PropertyLens)
    grid_element_type(model, var::PropertyLens)

Return the grid element type of a model variable (PropertyLens `var`) based on a `lens`. A
`lens` allows access to a nested model variable.
"""
function grid_element_type(
        ::T,
        var::PropertyLens,
    ) where {T <: Union{RiverFlowModel{<:LocalInertial}, OverlandFlowModel{<:LocalInertial}}}
    vars = (PropertyLens(x) for x in (:q, :q_average, :qx, :qy))
    element_type = if var in vars
        "edge"
    else
        "node"
    end
    return element_type
end

grid_element_type(model, var::PropertyLens) = "node"

function grid_element_type(model::Model, lens::ComposedFunction)
    lens_components = decompose(lens)
    var = lens_components[1]
    element_type = if PropertyLens(:river_flow) in lens_components
        grid_element_type(model.routing.river_flow, var)
    elseif PropertyLens(:overland_flow) in lens_components
        grid_element_type(model.routing.overland_flow, var)
    else
        grid_element_type(model, var)
    end
    return element_type
end

function get_standard_name(name::String)::String
    if startswith(name, "soil_layer_") && occursin(r"soil_layer_\d+_", name)
        # map to standard name for layered soil model variable (not available per layer)
        name, _ = soil_layer_standard_name(name)
    end
    return name
end
