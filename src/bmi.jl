#=
For the Basic Model Interface (BMI) implementation of wflow all grids are defined as 
unstructured grids:
    https://bmi-spec.readthedocs.io/en/latest/model_grids.html#unstructured-grids

While the input (forcing and model parameters) is structured (uniform rectilinear), 
internally wflow works with one dimensional arrays based on the active grid cells of the
2D model domain.

The variables that wflow can exchange through BMI should be listed under the section API in 
the config TOML file of the model type. Below an example for the sbm model type:

[API]
components = [
  "vertical",
  "lateral.subsurface",
  "lateral.land",
  "lateral.river",
  "lateral.river.reservoir"
]

See also the function BMI.get_input_var_names(model::Model).
=#

# BMI grid type based on grid identifier
gridtype = Dict(0 => "unstructured",
                1 => "unstructured",
                2 => "unstructured",
                3 => "unstructured",
                4 => "unstructured")

#=
Mapping of grid identifier to a key, to get the active indices of the model domain. See also
function active_indices(network, key::Tuple).
=#
grids = Dict(0 => "reservoir",
             1 => "lake",
             2 => "river",
             3 => "drain",
             4 => "land")

"""
    BMI.initialize(config_file)

Initialize the model. Reads the input settings and data as defined in the Config object
generated from the configuration file `config_file`. Will return a Model that is ready to
run. 
"""           
function BMI.initialize(config_file)
    config = Config(config_file)
    modeltype = config.model.type
    model = if modeltype == "sbm"
        initialize_sbm_model(config)
    elseif modeltype == "sbm_gwf"
        initialize_sbm_gwf_model(config)
    elseif modeltype == "hbv"
        initialize_hbv_model(config)
    else
        error("unknown model type")
    end
    return model
end

"""
    BMI.update(model::Model)

Advance model state by one time step.
"""
function BMI.update(model::Model)
    @unpack network, config = model
    update_func = config.model.type == "sbm_gwf" ? update_sbm_gwf : update
    return update_func(model)
end

"""
    BMI.update_until(model::Model, time::Float64)

Advance model state until and including the given `time`.

Time in units and epoch returned by the `BMI.get_time_units`.
"""
function BMI.update_until(model::Model, time::Float64)
    @unpack network, config = model
    update_func = config.model.type == "sbm_gwf" ? update_sbm_gwf : update
    curtime = BMI.get_current_time(model)
    n_iter = Int(max(0, (time - curtime) / model.clock.Δt.value))
    end_time = curtime + n_iter * config.timestepsecs
    @info("update model until $end_time")
    for i in 1:n_iter
        update_func(model)
    end
    return model
end

"""
    BMI.finalize(model::Model)

Write state output to netCDF and close files.
"""
function BMI.finalize(model::Model)
    @unpack config, writer, clock = model
    rewind!(clock)
    write_netcdf_timestep(model, writer.state_dataset, writer.state_parameters)
    reset_clock!(model.clock, config)
    close_files(model, delete_output = false)
end

"Returns the type of the model from the model configuration file."
function BMI.get_component_name(model::Model)
    @unpack config = model
    return config.model.type
end

"Returns the number of input variables of the model that can be exchanged."
function BMI.get_input_item_count(model::Model)
    length(BMI.get_input_var_names(model))
end

"Returns the number of output variables of the model that can be exchanged."
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
    if haskey(config,"API")
        var_names = Vector{String}()
        for c in config.API.components
            append!(var_names, collect(string.(c,".", fieldnames(typeof(param(model, c))))))
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

"""
    BMI.get_var_grid(model::Model, name::String)

Returns the grid identifier integer based on the model variable `name`.  
"""
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
        @warn("Model variable $name not listed as input of output variable")
        return -1
    end
end

"""
    BMI.get_var_type(model::Model, name::String)
Returns the type (as String) of a model variable `name`.
"""
function BMI.get_var_type(model::Model, name::String)
    string(typeof(param(model,name)))
end

""" 
    BMI.get_var_units(model::Model, name::String)
Returns the units of a model variable `name`.
"""
function BMI.get_var_units(model::Model, name::String)
    #TODO implement properly
    "mm"
end

"""
    BMI.get_var_itemsize(model::Model, name::String)
Returns the size (bytes) of one element of a model variable `name`.
"""
function BMI.get_var_itemsize(model::Model, name::String)
    sizeof(param(model, name)[1])
end

"""
    BMI.get_var_nbytes(model::Model, name::String)
Returns the total size (bytes) of a model variable `name`.
"""
function BMI.get_var_nbytes(model::Model, name::String)
    sizeof(param(model, name))
end

""" 
    BMI.get_var_location(model::Model, name::String)
Returns grid element type of a model variable `name`.
"""
function BMI.get_var_location(model::Model, name::String)
    if name in BMI.get_input_var_names(model)
        return "node"
    end
end

""" 
    BMI.get_current_time(model::Model)
Returns current time of the model.
"""
function BMI.get_current_time(model::Model)
    datetime2unix(model.clock.time)
end

""" 
    BMI.get_start_time(model::Model)
Returns start time of the model.
"""
function BMI.get_start_time(model::Model)
    datetime2unix(model.config.starttime)
end

"""
    BMI.get_end_time(model::Model)
Returns end time of the model.
"""
function BMI.get_end_time(model::Model)
    datetime2unix(model.config.endtime)
end

"""
    BMI.get_time_units()
Returns time units as reported by the BMI.
"""
function BMI.get_time_units()
    string("seconds since ", unix2datetime(0))
end

"""
    BMI.get_time_step(model::Model)
Returns time units used in the model.
"""
function BMI.get_time_step(model::Model)
    Float64(model.config.timestepsecs)
end

"""
    BMI.get_value(model::Model, name::String)
Returns copied values of a model variable `name`.
"""
function BMI.get_value(model::Model, name::String)
    copy(BMI.get_value_ptr(model, name))
end

"""
    BMI.get_value_ptr(model::Model, name::String)
Returns a reference to a model variable `name`.
"""
function BMI.get_value_ptr(model::Model, name::String)
    if name in BMI.get_input_var_names(model)
        param(model, name)
    else
        error("$name not defined as an output BMI variable")
    end        
end

"""
    BMI.BMI.get_value_at_indices(model::Model, name::String, inds::Vector{Int})
Returns values of a model variable `name` at indices `inds`.
"""
function BMI.get_value_at_indices(model::Model, name::String, inds::Vector{Int})
    BMI.get_value_ptr(model, name)[inds]
end

"""
    BMI.BMI.set_value(model::Model, name::String, src::Vector{T}) where T<:AbstractFloat
Set a model variable `name` to the values in vector `src`, overwriting the current contents.
The type and size of `src` must match the model’s internal array.
"""
function BMI.set_value(model::Model, name::String, src::Vector{T}) where T<:AbstractFloat
    BMI.get_value_ptr(model, name) .= src
end

"""
    BMI.set_value_at_indices(model::Model, name::String, inds::Vector{Int}, src::Vector{T})
    where T<:AbstractFloat 
Set a model variable `name` to the values in vector `src`, at indices `inds`.
"""
function BMI.set_value_at_indices(model::Model, name::String, inds::Vector{Int}, src::Vector{T}) where T<:AbstractFloat
    BMI.get_value_ptr(model, name)[inds] .= src
end

"""
    function BMI.get_grid_type(grid::Int)
Returns the type of a grid based on the grid identifier `grid`.
"""
function BMI.get_grid_type(grid::Int)
    gridtype[grid]
end

"""
    BMI.get_grid_rank(grid::Int)
Returns the rank of a grid based on the grid identifier `grid`.
"""
function BMI.get_grid_rank(grid::Int)
    if grid in keys(gridtype)
        2
    else
        @warn("unknow grid type $grid")
    end
end

"""
    BMI.get_grid_shape(model::Model, grid::Int)
Returns the dimensions of the model grid based on the grid identifier `grid`. 
"""
function BMI.get_grid_shape(model::Model, grid::Int)
     n = size(active_indices(model.network, symbols(grids[grid])))[1]
end

"""
    function BMI.get_grid_size(model::Model, grid::Int)
Returns the total number of nodes of the model grid based on the grid identifier `grid`. 
"""
function BMI.get_grid_size(model::Model, grid::Int)
    BMI.get_grid_shape(model, grid)
end

"""
    function BMI.get_grid_x(model::Model, grid::Int)
Returns the locations of the grid nodes in the first coordinate direction. 
"""
function BMI.get_grid_x(model::Model, grid::Int)
    @unpack reader, config = model
    @unpack dataset = reader
    sel = active_indices(model.network, symbols(grids[grid]))
    inds = [sel[i][1] for i in eachindex(sel)]
    x_nc = "x" in keys(dataset.dim) ? ncread(dataset, "x") : ncread(dataset, "lon")
    x_nc[inds]
end

"""
    function BMI.get_grid_y(model::Model, grid::Int)
Returns the locations of the grid nodes in the second coordinate direction. 
"""
function BMI.get_grid_y(model::Model, grid::Int)
    @unpack reader, config = model
    @unpack dataset = reader
    sel = active_indices(model.network, symbols(grids[grid]))
    inds = [sel[i][2] for i in eachindex(sel)]
    y_nc = "y" in keys(dataset.dim) ? ncread(dataset, "y") : ncread(dataset, "lat")
    y_nc[inds]
end

"""
    function BMI.get_grid_node_count(model::Model, grid::Int)
Returns the number of nodes in the model grid based on the grid identifier `grid`. 
"""
function BMI.get_grid_node_count(model::Model, grid::Int)
    BMI.get_grid_size(model, grid)
end