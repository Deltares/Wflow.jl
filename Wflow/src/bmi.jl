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
    BMI.initialize(::Type{<:Wflow.Model}, config_file)

Initialize the model. Reads the input settings and data as defined in the Config object
generated from the configuration file `config_file`. Will return a Model that is ready to
run.
"""
function BMI.initialize(::Type{<:Model}, config_file::AbstractString)
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
    BMI.update(model::Model)

Update the model for a single timestep.
"""
function BMI.update(model::Model)
    run_timestep!(model)
    return nothing
end

function BMI.update_until(model::Model, time::Float64)
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
function BMI.finalize(model::Model)
    (; config, writer) = model
    # it is possible that the state dataset has been closed by `save_state`
    if !isnothing(writer.state_dataset) && isopen(writer.state_dataset)
        write_netcdf_timestep(model, writer.state_dataset, writer.state_parameters)
    end
    reset_clock!(model.clock, config)
    close_files(model; delete_output = false)
    return nothing
end

function BMI.get_component_name(model::Model)
    return string(model.config.model.type)
end

function BMI.get_input_item_count(model::Model)
    return length(BMI.get_input_var_names(model))
end

function BMI.get_output_item_count(model::Model)
    return length(BMI.get_output_var_names(model))
end

"""
    BMI.get_input_var_names(model::Model)

Returns model input variables, based on the `API` section in the model configuration file.
This `API` sections contains a list of `Model` components for which variables can be
exchanged.
"""
function BMI.get_input_var_names(model::Model)
    (; config, land) = model
    if do_api(config)
        var_names = config.API.variables
        idx = []
        for (i, var) in enumerate(var_names)
            if startswith(var, "soil_layer_") && occursin(r"soil_layer_\d+_", var)
                # map to standard name for layered soil model variable (not available per layer)
                var, _ = soil_layer_standard_name(var)
            end
            if !haskey(standard_name_map(land), var)
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
        return []
    end
end

"Returns input variables from `BMI.get_input_var_names(model::Model)`, there is no
distinction between input - and output variables."
function BMI.get_output_var_names(model::Model)
    return BMI.get_input_var_names(model)
end

function BMI.get_var_grid(::Model, name::String)
    return if occursin("reservoir", name)
        0
    elseif occursin("drain", name)
        1
    elseif occursin("river", name) || occursin("floodplain", name)
        2
    elseif occursin("land_surface_water__x_component", name) # LocalInertialOverlandFlow
        3
    elseif occursin("land_surface_water__y_component", name) # LocalInertialOverlandFlow
        4
    else
        5
    end
end

function BMI.get_var_type(model::Model, name::String)
    value = BMI.get_value_ptr(model, name)
    return repr(eltype(first(value)))
end

function BMI.get_var_units(model::Model, name::String)
    (; land) = model
    (; unit) = standard_name_map(land)[name]
    return to_string(to_SI(unit); BMI_standard = true)
end

function BMI.get_var_itemsize(model::Model, name::String)
    value = BMI.get_value_ptr(model, name)
    return sizeof(eltype(first(value)))
end

function BMI.get_var_nbytes(model::Model, name::String)
    return sizeof(BMI.get_value_ptr(model, name))
end

function BMI.get_var_location(model::Model, name::String)
    (; land) = model
    lens = standard_name_map(land)[name].lens
    element_type = grid_element_type(model, lens)
    return element_type
end

function BMI.get_current_time(model::Model)
    (; config, clock) = model
    (; starttime, calendar) = config.time
    starttime = cftime(starttime, calendar)
    return Dates.value(clock.time - starttime)
end

function BMI.get_start_time(::Model)
    return 0.0
end

function BMI.get_end_time(model::Model)
    (; starttime, endtime, calendar) = model.config.time
    starttime_ = cftime(starttime, calendar)
    endtime_ = cftime(endtime, calendar)
    return Dates.value(endtime_ - starttime_)
end

function BMI.get_time_units(model::Model)
    return "s"
end

function BMI.get_time_step(model::Model)
    return model.config.time.timestepsecs
end

function BMI.get_value(model::Model, name::String, dest::Vector{Float64})
    dest .= copy(BMI.get_value_ptr(model, name))
    return dest
end

function BMI.get_value_ptr(model::Model, name::String)
    (; land, domain) = model
    n = length(active_indices(domain, name))

    if startswith(name, "soil_layer_") && occursin(r"soil_layer_\d+_", name)
        name_2d, ind = soil_layer_standard_name(name)
        lens = standard_name_map(land)[name_2d].lens
        model_vals = lens(model)
        el_type = eltype(first(model_vals))
        dim = length(first(model_vals))
        value = reshape(reinterpret(el_type, model_vals), dim, :)
        return @view value[ind, 1:n]
    else
        lens = standard_name_map(land)[name].lens
        return @view(lens(model)[1:n])
    end
end

function BMI.get_value_at_indices(
    model::Model,
    name::String,
    dest::Vector{Float64},
    inds::Vector{Int},
)
    dest .= BMI.get_value_ptr(model, name)[inds]
    return dest
end

"""
    BMI.set_value(model::Model, name::String, src::Vector{Float64})

Set a model variable `name` to the values in vector `src`, overwriting the current contents.
The type and size of `src` must match the model's internal array.
"""
function BMI.set_value(model::Model, name::String, src::Vector{Float64})
    return BMI.get_value_ptr(model, name) .= src
end

"""
    BMI.set_value_at_indices(model::Model, name::String, inds::Vector{Int}, src::Vector{Float64})

    Set a model variable `name` to the values in vector `src`, at indices `inds`.
"""
function BMI.set_value_at_indices(
    model::Model,
    name::String,
    inds::Vector{Int},
    src::Vector{Float64},
)
    return BMI.get_value_ptr(model, name)[inds] .= src
end

function BMI.get_grid_type(model::Model, grid::Int)
    if grid in 0:2
        "points"
    elseif grid in 3:6
        "unstructured"
    else
        error("unknown grid type $grid")
    end
end

function BMI.get_grid_rank(model::Model, grid::Int)
    if grid in 0:6
        2
    else
        error("unknown grid type $grid")
    end
end

function BMI.get_grid_x(model::Model, grid::Int, x::Vector{Float64})
    (; reader, domain) = model
    (; dataset) = reader
    sel = active_indices(domain, GRIDS[grid])
    inds = [sel[i][1] for i in eachindex(sel)]
    x_nc = read_x_axis(dataset)
    x .= x_nc[inds]
    return x
end

function BMI.get_grid_y(model::Model, grid::Int, y::Vector{Float64})
    (; reader, domain) = model
    (; dataset) = reader
    sel = active_indices(domain, GRIDS[grid])
    inds = [sel[i][2] for i in eachindex(sel)]
    y_nc = read_y_axis(dataset)
    y .= y_nc[inds]
    return y
end

function BMI.get_grid_node_count(model::Model, grid::Int)
    return length(active_indices(model.domain, GRIDS[grid]))
end

function BMI.get_grid_size(model::Model, grid::Int)
    return length(active_indices(model.domain, GRIDS[grid]))
end

function BMI.get_grid_edge_count(model::Model, grid::Int)
    (; domain) = model
    if grid == 3
        return ne(domain.river.network.graph)
    elseif grid == 4
        return length(domain.land.network.edge_indices.xu)
    elseif grid == 5
        return length(domain.land.network.edge_indices.yu)
    elseif grid in 0:2 || grid == 6
        @warn("edges are not provided for grid type $grid (variables are located at nodes)")
    else
        error("unknown grid type $grid")
    end
end

function BMI.get_grid_edge_nodes(model::Model, grid::Int, edge_nodes::Vector{Int})
    (; domain) = model
    n = length(edge_nodes)
    m = div(n, 2)
    # inactive nodes (boundary/ghost points) are set at -999
    if grid == 3
        nodes_at_edge = adjacent_nodes_at_edge(domain.river.network.graph)
        nodes_at_edge.dst[nodes_at_edge.dst .== m + 1] .= -999
        edge_nodes[range(1, n; step = 2)] = nodes_at_edge.src
        edge_nodes[range(2, n; step = 2)] = nodes_at_edge.dst
        return edge_nodes
    elseif grid == 4
        xu = domain.land.network.edge_indices.xu
        edge_nodes[range(1, n; step = 2)] = 1:m
        xu[xu .== m + 1] .= -999
        edge_nodes[range(2, n; step = 2)] = xu
        return edge_nodes
    elseif grid == 5
        yu = domain.land.network.edge_indices.yu
        edge_nodes[range(1, n; step = 2)] = 1:m
        yu[yu .== m + 1] .= -999
        edge_nodes[range(2, n; step = 2)] = yu
        return edge_nodes
    elseif grid in 0:2 || grid == 6
        @warn("edges are not provided for grid type $grid (variables are located at nodes)")
    else
        error("unknown grid type $grid")
    end
end

# Extension of BMI functions (state handling and start time), required for OpenDA coupling.
# May also be useful for other external software packages.
function load_state(model::Model)
    set_states!(model)
    return nothing
end

function save_state(model::Model)
    (; writer) = model
    if !isnothing(writer.state_nc_path)
        @info "Write output states to netCDF file `$(writer.state_nc_path)`."
    end
    write_netcdf_timestep(model, writer.state_dataset, writer.state_parameters)
    close(writer.state_dataset)
    return nothing
end

function get_start_unix_time(model::Model)
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
) where {T <: Union{LocalInertialRiverFlow, LocalInertialOverlandFlow}}
    vars = (PropertyLens(x) for x in (:q, :q_av, :qx, :qy))
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
