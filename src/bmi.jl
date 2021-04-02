# Basic Model Interface (BMI) implementation based on
# https://github.com/Deltares/BasicModelInterface.jl

# BMI grid type based on grid identifier
const gridtype = Dict{Int,String}(
    0 => "unstructured",
    1 => "unstructured",
    2 => "unstructured",
    3 => "unstructured",
    4 => "unstructured",
)

# Mapping of grid identifier to a key, to get the active indices of the model domain.
# See also function active_indices(network, key::Tuple).
const grids =
    Dict{Int,String}(0 => "reservoir", 1 => "lake", 2 => "river", 3 => "drain", 4 => "land")

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
    else
        error("unknown model type")
    end
    return model
end

"""
    BMI.update(model::Model; run = nothing)

Update the model for a single timestep. 
# Arguments
- `run = nothing`: to update a model partially.
"""
function BMI.update(model::Model; run = nothing)
    @unpack network, config = model
    if isnothing(run)
        update_func = config.model.type == "sbm_gwf" ? update_sbm_gwf : update
    elseif run == "sbm_until_recharge"
        update_func = update_until_recharge
    elseif run == "sbm_after_subsurfaceflow"
        update_func = update_after_subsurfaceflow
    end
    return update_func(model)
end

function BMI.update_until(model::Model, time::Float64)
    @unpack network, config = model
    update_func = config.model.type == "sbm_gwf" ? update_sbm_gwf : update
    curtime = BMI.get_current_time(model)
    n_iter = Int(max(0, (time - curtime) / model.clock.Δt.value))
    end_time = curtime + n_iter * config.timestepsecs
    @info("update model until $end_time")
    for i = 1:n_iter
        update_func(model)
    end
    return model
end

"Write state output to netCDF and close files."
function BMI.finalize(model::Model)
    @unpack config, writer, clock = model
    rewind!(clock)
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
            append!(
                var_names,
                collect(string.(c, ".", fieldnames(typeof(param(model, c))))),
            )
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
    if name in BMI.get_input_var_names(model)
        if occursin("reservoir", name)
            0
        elseif occursin("lake", name)
            1
        elseif occursin("river", name)
            2
        elseif occursin("drain", name)
            3
        else
            4
        end
    else
        error("Model variable $name not listed as input of output variable")
    end
end

function BMI.get_var_type(model::Model, name::String)
    string(typeof(param(model, name)))
end

function BMI.get_var_units(model::Model, name::String)
    s = split(name, ".")
    get_units(param(model, join(s[1:end-1], ".")), Symbol(s[end]))
end

function BMI.get_var_itemsize(model::Model, name::String)
    sizeof(param(model, name)[1])
end

function BMI.get_var_nbytes(model::Model, name::String)
    sizeof(param(model, name))
end

function BMI.get_var_location(model::Model, name::String)
    if name in BMI.get_input_var_names(model)
        return "node"
    else
        error("$name not in get_input_var_names")
    end
end

function BMI.get_current_time(model::Model)
    datetime2unix(reinterpret(DateTime, model.clock.time))
end

function BMI.get_start_time(model::Model)
    datetime2unix(DateTime(model.config.starttime))
end

function BMI.get_end_time(model::Model)
    datetime2unix(DateTime(model.config.endtime))
end

function BMI.get_time_units(::Type{<:Model})
    string("seconds since ", unix2datetime(0))
end

function BMI.get_time_step(model::Model)
    Float64(model.config.timestepsecs)
end

function BMI.get_value(model::Model, name::String)
    copy(BMI.get_value_ptr(model, name))
end

function BMI.get_value_ptr(model::Model, name::String)
    if name in BMI.get_input_var_names(model)
        return param(model, name)
    else
        error("$name not defined as an output BMI variable")
    end
end

"""
    BMI.get_value_at_indices(model::Model, name::String, inds::Vector{Int})

Returns values of a model variable `name` at indices `inds`.
"""
function BMI.get_value_at_indices(model::Model, name::String, inds::Vector{Int})
    BMI.get_value_ptr(model, name)[inds]
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
    gridtype[grid]
end

function BMI.get_grid_rank(model::Model, grid::Int)
    if grid in keys(gridtype)
        return 2
    else
        error("unknown grid type $grid")
    end
end

function BMI.get_grid_shape(model::Model, grid::Int)
    size(active_indices(model.network, symbols(grids[grid])), 1)
end

function BMI.get_grid_size(model::Model, grid::Int)
    BMI.get_grid_shape(model, grid)
end

function BMI.get_grid_x(model::Model, grid::Int)
    @unpack reader, config = model
    @unpack dataset = reader
    sel = active_indices(model.network, symbols(grids[grid]))
    inds = [sel[i][1] for i in eachindex(sel)]
    x_nc = "x" in keys(dataset.dim) ? ncread(dataset, "x") : ncread(dataset, "lon")
    return x_nc[inds]
end

function BMI.get_grid_y(model::Model, grid::Int)
    @unpack reader, config = model
    @unpack dataset = reader
    sel = active_indices(model.network, symbols(grids[grid]))
    inds = [sel[i][2] for i in eachindex(sel)]
    y_nc = "y" in keys(dataset.dim) ? ncread(dataset, "y") : ncread(dataset, "lat")
    return y_nc[inds]
end

function BMI.get_grid_node_count(model::Model, grid::Int)
    BMI.get_grid_size(model, grid)
end
