# Basic Model Interface (BMI) implementation based on
# https://github.com/Deltares/BasicModelInterface.jl

# Mapping of grid identifier to a key, to get the active indices of the model domain.
# See also function active_indices(network, key::Tuple).
grids = Dict{Int,Tuple}(
    1 => (:reservoir,),
    2 => (:lake,),
    3 => (:drain,),
    4 => (:river,),
    5 => (:land,),
    6 => (:land,),
    7 => (:land,),
)

"""
    BMI.initialize(::Type{<:Wflow.Model}, config_file)

Initialize the model. Reads the input settings and data as defined in the Config object
generated from the configuration file `config_file`. Will return a Model that is ready to
run.
"""
function BMI.initialize(::Type{<:Model}, config_file)
    config = Config(config_file)
    modeltype = config.model.type
    model = if modeltype == "sbm"
        initialize_sbm_model(config)
    elseif modeltype == "sbm_gwf"
        initialize_sbm_gwf_model(config)
    elseif modeltype == "hbv"
        initialize_hbv_model(config)
    elseif modeltype == "sediment"
        initialize_sediment_model(config)
    elseif modeltype == "flextopo"
        initialize_flextopo_model(config)
    else
        error("unknown model type")
    end
    load_fixed_forcing(model)
    return model
end

"""
    BMI.update(model::Model; run = nothing)

Update the model for a single timestep.
# Arguments
- `run = nothing`: to update a model partially.
"""
function BMI.update(model::Model; run = nothing)
    @unpack clock, network, config = model
    if isnothing(run)
        model = run_timestep(model)
    elseif run == "sbm_until_recharge"
        model = run_timestep(
            model,
            update_func = update_until_recharge,
            write_model_output = false,
        )
    elseif run == "sbm_after_subsurfaceflow"
        model = run_timestep(model, update_func = update_after_subsurfaceflow)
    end
    return model
end

function BMI.update_until(model::Model, time::Float64)
    @unpack clock, network, config = model
    curtime = BMI.get_current_time(model)
    n = Int(max(0, (time - curtime) / model.clock.Δt.value))
    for _ = 1:n
        model = run_timestep(model)
    end
    return model
end

"Write state output to netCDF and close files."
function BMI.finalize(model::Model)
    @unpack config, writer, clock = model
    write_netcdf_timestep(model, writer.state_dataset, writer.state_parameters)
    reset_clock!(model.clock, config)
    close_files(model, delete_output = false)
end

function BMI.get_component_name(model::Model)
    @unpack config = model
    return config.model.type
end

function BMI.get_input_item_count(model::Model)
    length(BMI.get_input_var_names(model))
end

function BMI.get_output_item_count(model::Model)
    length(BMI.get_output_var_names(model))
end

"""
    BMI.get_input_var_names(model::Model)

Returns model input variables, based on the `API` section in the model configuration file.
This `API` sections contains a list of `Model` components for which variables can be
exchanged.
"""
function BMI.get_input_var_names(model::Model)
    @unpack config = model
    if haskey(config, "API")
        var_names = Vector{String}()
        for c in config.API.components
            type = typeof(param(model, c))
            inds = findall(x -> x != 0, exchange(type))
            field_names = fieldnames(type)[inds]

            for name in field_names
                var = string(c, ".", name)
                model_var = param(model, var)
                if eltype(model_var) <: SVector
                    for i = 1:length(first(model_var))
                        push!(var_names, string(var, "[", i, "]"))
                    end
                elseif ndims(model_var) > 1
                    for i = 1:length(first(model_var))
                        push!(var_names, string(var, "[", i, "]"))
                    end
                else
                    push!(var_names, var)
                end
            end
        end
        return var_names
    else
        @warn("TOML file does not contain section [API] to extract model var names")
        return []
    end
end

"Returns input variables from `BMI.get_input_var_names(model::Model)`, there is no
distinction between input - and output variables."
function BMI.get_output_var_names(model::Model)
    BMI.get_input_var_names(model)
end

function BMI.get_var_grid(model::Model, name::String)
    s = split(name, "[")
    key = symbols(first(s))
    if exchange(param(model, key[1:end-1]), key[end]) == 1
        gridtype = grid_type(param(model, key))
        type = typeof(param(model, key))
        if gridtype == "scalar"
            0
        elseif :reservoir in key
            1
        elseif :lake in key
            2
        elseif :drain in key
            3
        elseif :river in key
            4
        elseif type <: ShallowWaterLand && occursin("x", s[end])
            5
        elseif type <: ShallowWaterLand && occursin("y", s[end])
            6
        else
            7
        end
    else
        error("$name not listed as variable for BMI exchange")
    end
end

function BMI.get_var_type(model::Model, name::String)
    value = BMI.get_value_ptr(model, name)
    repr(eltype(first(value)))
end

function BMI.get_var_units(model::Model, name::String)
    key = symbols(first(split(name, "[")))
    if exchange(param(model, key[1:end-1]), key[end]) == 1
        get_units(param(model, key[1:end-1]), key[end])
    else
        error("$name not listed as variable for BMI exchange")
    end
end

function BMI.get_var_itemsize(model::Model, name::String)
    value = BMI.get_value_ptr(model, name)
    sizeof(eltype(first(value)))
end

function BMI.get_var_nbytes(model::Model, name::String)
    sizeof(BMI.get_value_ptr(model, name))
end

function BMI.get_var_location(model::Model, name::String)
    key = symbols(first(split(name, "[")))
    if exchange(param(model, key[1:end-1]), key[end]) == 1
        return grid_location(model, key)
    else
        error("$name not listed as variable for BMI exchange")
    end
end

function BMI.get_current_time(model::Model)
    0.001 * Dates.value(model.clock.time - model.config.starttime)
end

function BMI.get_start_time(model::Model)
    0.0
end

function BMI.get_end_time(model::Model)
    0.001 * Dates.value(model.config.endtime - model.config.starttime)
end

function BMI.get_time_units(model::Model)
    "s"
end

function BMI.get_time_step(model::Model)
    Float64(model.config.timestepsecs)
end

function BMI.get_value(model::Model, name::String, dest::Vector{T}) where {T<:AbstractFloat}
    dest .= copy(BMI.get_value_ptr(model, name))
    return dest
end

function BMI.get_value_ptr(model::Model, name::String)
    @unpack network = model
    s = split(name, "[")
    key = symbols(first(s))
    if exchange(param(model, key[1:end-1]), key[end]) == 1
        n = length(active_indices(network, key))
        if occursin("[", name)
            ind = tryparse(Int, split(s[end], "]")[1])
            if eltype(param(model, key)) <: SVector
                model_vals = param(model, key)
                el_type = eltype(first(model_vals))
                dim = length(first(model_vals))
                value = reshape(reinterpret(el_type, model_vals), dim, :)
                return @view value[ind, 1:n]
            else
                value = @view param(model, key)[ind, 1:n]
                return value
            end
        else
            return @view(param(model, key)[1:n])
        end
    else
        error("$name not listed as variable for BMI exchange")
    end
end

function BMI.get_value_at_indices(
    model::Model,
    name::String,
    dest::Vector{T},
    inds::Vector{Int},
) where {T<:AbstractFloat}
    dest .= BMI.get_value_ptr(model, name)[inds]
    return dest
end

"""
    BMI.set_value(model::Model, name::String, src::Vector{T}) where T<:AbstractFloat

Set a model variable `name` to the values in vector `src`, overwriting the current contents.
The type and size of `src` must match the model’s internal array.
"""
function BMI.set_value(model::Model, name::String, src::Vector{T}) where {T<:AbstractFloat}
    BMI.get_value_ptr(model, name) .= src
end

"""
    BMI.set_value_at_indices(model::Model, name::String, inds::Vector{Int}, src::Vector{T})
    where T<:AbstractFloat

    Set a model variable `name` to the values in vector `src`, at indices `inds`.
"""
function BMI.set_value_at_indices(
    model::Model,
    name::String,
    inds::Vector{Int},
    src::Vector{T},
) where {T<:AbstractFloat}
    BMI.get_value_ptr(model, name)[inds] .= src
end

function BMI.get_grid_type(model::Model, grid::Int)
    if grid == 0
        "scalar"
    elseif grid in range(1, 3)
        "points"
    elseif grid in range(4, 7)
        "unstructured"
    else
        error("unknown grid type $grid")
    end
end

function BMI.get_grid_rank(model::Model, grid::Int)
    if grid == 0
        0
    elseif grid in range(1, 7)
        2
    else
        error("unknown grid type $grid")
    end
end

function BMI.get_grid_x(model::Model, grid::Int, x::Vector{T}) where {T<:AbstractFloat}
    @unpack reader, config = model
    @unpack dataset = reader
    sel = active_indices(model.network, grids[grid])
    inds = [sel[i][1] for i in eachindex(sel)]
    x_nc = read_x_axis(dataset)
    x .= x_nc[inds]
    return x
end

function BMI.get_grid_y(model::Model, grid::Int, y::Vector{T}) where {T<:AbstractFloat}
    @unpack reader, config = model
    @unpack dataset = reader
    sel = active_indices(model.network, grids[grid])
    inds = [sel[i][2] for i in eachindex(sel)]
    y_nc = read_y_axis(dataset)
    y .= y_nc[inds]
    return y
end

function BMI.get_grid_node_count(model::Model, grid::Int)
    return length(active_indices(model.network, grids[grid]))
end

function BMI.get_grid_edge_count(model::Model, grid::Int)
    @unpack network = model
    if grid == 4
        return ne(network.river.graph)
    elseif grid == 5
        return length(network.land.staggered_indices.xu)
    elseif grid == 6
        return length(network.land.staggered_indices.yu)
    elseif grid in range(0, 3) || grid == 7
        warn("edges are not provided for grid type $grid (variables are located at nodes)")
    else
        error("unknown grid type $grid")
    end
end

function BMI.get_grid_edge_nodes(model::Model, grid::Int, edge_nodes::Vector{Int})
    @unpack network = model
    n = length(edge_nodes)
    m = div(n, 2)
    # inactive nodes (boundary/ghost points) are set at -999
    if grid == 4
        nodes_at_edge = adjacent_nodes_at_link(network.river.graph)
        nodes_at_edge.dst[nodes_at_edge.dst.==m+1] .= -999
        edge_nodes[range(1, n, step = 2)] = nodes_at_edge.src
        edge_nodes[range(2, n, step = 2)] = nodes_at_edge.dst
        return edge_nodes
    elseif grid == 5
        xu = network.staggered_indices.xu
        edge_nodes[range(1, n, step = 2)] = range(1, m)
        xu[xu.==m+1] .= -999
        edge_nodes[range(2, n, step = 2)] = xu
        return edge_nodes
    elseif grid == 6
        yu = network.staggered_indices.yu
        edge_nodes[range(1, n, step = 2)] = range(1, m)
        yu[yu.==m+1] .= -999
        edge_nodes[range(2, n, step = 2)] = yu
        return edge_nodes
    elseif grid in range(0, 3) || grid == 7
        warn("edges are not provided for grid type $grid (variables are located at nodes)")
    else
        error("unknown grid type $grid")
    end
end
