#=
Code to support input and output of data and configuration.
Input data can be loaded from NetCDF files.
Output data can be written to NetCDF or CSV files.
For configuration files we use TOML.
=#

"""Turn "a.aa.aaa" into (:a, :aa, :aaa)"""
symbols(s) = Tuple(Symbol(x) for x in split(s, '.'))

"""Turn symbols"a.aa.aaa" into (:a, :aa, :aaa)"""
macro symbols_str(s)
    Tuple(Symbol(x) for x in split(s, '.'))
end

"Get a nested field using a tuple of Symbols"
param(obj, fields) = foldl(getproperty, fields; init = obj)
param(obj, fields::AbstractString) = param(obj, symbols(fields))
function param(obj, fields, default)
    try
        return param(obj, fields)
    catch
        return default
    end
end

"""
    Config(path::AbstractString)
    Config(dict::AbstractDict)
    Config(dict::Dict{String,Any}, path::Union{String,Nothing})

Struct that contains the parsed TOML configuration, as well as a reference to the TOML path,
if it exists. It behaves largely like a distionary, but it overloads `getproperty` and
`setproperty` to support syntax like `config.model.reinit = false`.
"""
struct Config
    dict::Dict{String,Any}  # nested key value mapping of all settings
    path::Union{String,Nothing}  # path to the TOML file, or nothing
end

Config(path::AbstractString) = Config(TOML.parsefile(path), path)
Config(dict::AbstractDict) = Config(dict, nothing)

# allows using getproperty, e.g. config.input.time instead of config["input"]["time"]
function Base.getproperty(config::Config, f::Symbol)
    dict = Dict(config)
    path = pathof(config)
    a = dict[String(f)]
    # if it is a Dict, wrap the result in Config to keep the getproperty behavior
    return a isa AbstractDict ? Config(a, path) : a
end

function Base.setproperty!(config::Config, f::Symbol, x)
    dict = Dict(config)
    return dict[String(f)] = x
end

"Get a value from the Config with either the key or an alias of the key."
function get_alias(config::Config, key, alias, default)
    alias_or_default = param(config, alias, default)
    return param(config, key, alias_or_default)
end

"Get a value from the Config with a key, throwing an error if it is not one of the options."
function get_options(config::Config, key, options, default)
    dict = Dict(config)
    a = get(dict, key, default)
    if !(a in options)
        error("TOML key `$key` is set to $(repr(a)). Valid options are: $options.")
    end
    return a
end

# also used in autocomplete
Base.propertynames(config::Config) = collect(keys(Dict(config)))
Base.haskey(config::Config, key) = haskey(Dict(config), key)
Base.keys(config::Config) = keys(Dict(config))
Base.values(config::Config) = values(Dict(config))
Base.pairs(config::Config) = pairs(Dict(config))
Base.get(config::Config, key, default) = get(Dict(config), key, default)
Base.getindex(config::Config, key) = getindex(Dict(config), key)
Base.setindex!(config::Config, val, key) = setindex!(Dict(config), val, key)
Base.pop!(config::Config, key) = pop!(Dict(config), key)
Base.pop!(config::Config, key, default) = pop!(Dict(config), key, default)
Base.Dict(config::Config) = getfield(config, :dict)
Base.pathof(config::Config) = getfield(config, :path)
Base.iterate(config::Config) = iterate(Dict(config))
Base.iterate(config::Config, state) = iterate(Dict(config), state)

function Base.dirname(config::Config)
    path = pathof(config)
    return path === nothing ? nothing : dirname(path)
end

function combined_path(config::Config, dir::AbstractString, path::AbstractString)
    tomldir = dirname(config)
    return normpath(tomldir, dir, path)
end

"Construct a path relative to both the TOML directory and the optional `dir_input`"
function input_path(config::Config, path::AbstractString)
    dir = get(config, "dir_input", ".")
    return combined_path(config, dir, path)
end

"Construct a path relative to both the TOML directory and the optional `dir_output`"
function output_path(config::Config, path::AbstractString)
    dir = get(config, "dir_output", ".")
    return combined_path(config, dir, path)
end

"Extract NetCDF variable name `ncname` from `var` (type `String` or `Config`). If `var` has 
type `Config`, either `scale`, `offset` and an optional `index` are expected (with `ncname`) 
or a `value` (uniform value), these are stored as part of `NamedTuple` `modifier`."
function ncvar_name_modifier(var; config = nothing)
    ncname = nothing
    modifier = (scale = 1.0, offset = 0.0, value = nothing, index = nothing)
    if isa(var, Config)
        if haskey(var, "netcdf") &&
           haskey(var.netcdf, "variable") &&
           haskey(var.netcdf.variable, "name")
            ncname = var.netcdf.variable.name
            scale = param(var, "scale", 1.0)
            offset = param(var, "offset", 0.0)
            if haskey(var, "layer") || haskey(var, "class")
                haskey(var, "layer") ? dim_name = "layer" : dim_name = "class"
                if length(var[dim_name])>1
                    # if modifier is provided as a list for each dim item
                    indices = []
                    for i = 1:length(var[dim_name])
                        index = get_index_dimension(var, config, var[dim_name][i])
                        @info "NetCDF parameter `$ncname` is modified with scale `$(scale[i])` and offset `$(offset[i])` at index `$index`."
                        push!(indices, index)
                    end
                    modifier = (scale = scale, offset = offset, value = nothing, index = indices)
                else 
                    index = get_index_dimension(var, config, var[dim_name])
                    modifier = (scale = scale, offset = offset, value = nothing, index = index)
                    @info "NetCDF parameter `$ncname` is modified with scale `$scale` and offset `$offset` at index `$index`."
                end
            else 
                modifier = (scale = scale, offset = offset, value = nothing, index = nothing)
                @info "NetCDF parameter `$ncname` is modified with scale `$scale` and offset `$offset`."
            end
        elseif haskey(var, "value")
            modifier = (scale = 1.0, offset = 0.0, value = param(var, "value"), index = nothing)
        else
            error("Unrecognized modifier $(Dict(var))")
        end
    elseif isa(var, String)
        ncname = var
    else
        error("Unknown type")
    end
    return ncname, modifier
end

"Extract a NetCDF variable at a given time"
function get_at(
    ds::CFDataset,
    varname::AbstractString,
    times::AbstractVector{<:TimeType},
    t::TimeType,
)
    # this behaves like a backward fill interpolation
    i = findfirst(>=(t), times)
    t < first(times) && throw(DomainError("time $t before dataset begin $(first(times))"))
    i === nothing && throw(DomainError("time $t after dataset end $(last(times))"))
    return get_at(ds, varname, i)
end

function get_at(ds::CFDataset, varname::AbstractString, i)
    return read_standardized(ds, varname, (x = :, y = :, time = i))
end

function get_param_res(model)
    Dict(
        symbols"vertical.precipitation" => model.lateral.river.reservoir.precipitation,
        symbols"vertical.potential_evaporation" =>
            model.lateral.river.reservoir.evaporation,
    )
end

function get_param_lake(model)
    Dict(
        symbols"vertical.precipitation" => model.lateral.river.lake.precipitation,
        symbols"vertical.potential_evaporation" => model.lateral.river.lake.evaporation,
    )
end

mover_params = (symbols"vertical.precipitation", symbols"vertical.potential_evaporation")

function load_fixed_forcing(model)
    @unpack reader, network, config = model
    @unpack forcing_parameters = reader

    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool

    mover_params =
        (symbols"vertical.precipitation", symbols"vertical.potential_evaporation")
    reverse_indices = network.land.reverse_indices
    if do_reservoirs
        sel_reservoirs = network.reservoir.indices_coverage
        param_res = get_param_res(model)
    end
    if do_lakes
        sel_lakes = network.lake.indices_coverage
        param_lake = get_param_lake(model)
    end

    for (par, ncvar) in forcing_parameters
        if ncvar.name === nothing
            val = ncvar.value * ncvar.scale + ncvar.offset
            param_vector = param(model, par)
            param_vector .= val
            # set fixed precipitation and evaporation over the lakes and reservoirs and put
            # these into the lakes and reservoirs structs and set the precipitation and
            # evaporation to 0 in the vertical model
            if par in mover_params
                if do_reservoirs
                    for (i, sel_reservoir) in enumerate(sel_reservoirs)
                        param_vector[reverse_indices[sel_reservoir]] .= 0
                        param_res[par][i] = val
                    end
                end
                if do_lakes
                    for (i, sel_lake) in enumerate(sel_lakes)
                        param_vector[reverse_indices[sel_lake]] .= 0
                        param_lake[par][i] = val
                    end
                end
            end
        end
    end
end

"Get dynamic NetCDF input for the given time"
function update_forcing!(model)
    @unpack vertical, clock, reader, network, config = model
    @unpack dataset, forcing_parameters = reader
    nctimes = dataset["time"][:]

    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool

    if do_reservoirs
        sel_reservoirs = network.reservoir.indices_coverage
        param_res = get_param_res(model)
    end
    if do_lakes
        sel_lakes = network.lake.indices_coverage
        param_lake = get_param_lake(model)
    end

    # Wflow expects `right` labeling of the forcing time interval, e.g. daily precipitation
    # at 01-02-2000 00:00:00 is the accumulated total precipitation between 01-01-2000
    # 00:00:00 and 01-02-2000 00:00:00.

    # load from NetCDF into the model according to the mapping
    for (par, ncvar) in forcing_parameters
        # no need to update fixed values
        ncvar.name === nothing && continue

        time = convert(eltype(nctimes), clock.time)
        data = get_at(dataset, ncvar.name, nctimes, time)

        if ncvar.scale != 1.0 || ncvar.offset != 0.0
            data .= data .* ncvar.scale .+ ncvar.offset
        end

        # calculate the mean precipitation and evaporation over the lakes and reservoirs
        # and put these into the lakes and reservoirs structs
        # and set the precipitation and evaporation to 0 in the vertical model
        if par in mover_params
            if do_reservoirs
                for (i, sel_reservoir) in enumerate(sel_reservoirs)
                    avg = mean(data[sel_reservoir])
                    data[sel_reservoir] .= 0
                    param_res[par][i] = avg
                end
            end
            if do_lakes
                for (i, sel_lake) in enumerate(sel_lakes)
                    avg = mean(data[sel_lake])
                    data[sel_lake] .= 0
                    param_lake[par][i] = avg
                end
            end
        end

        param_vector = param(model, par)
        sel = active_indices(network, par)
        data_sel = data[sel]
        if any(ismissing, data_sel)
            print(par)
            msg = "Forcing data has missing values on active model cells for $(ncvar.name)"
            throw(ArgumentError(msg))
        end
        param_vector .= data_sel
    end

    return model
end

"""
    monthday_passed(curr, avail)

Given two monthday tuples such as (12, 31) and (12, 15), return true if the first argument
falls after or on the same day as the second argument, assuming the same year. The tuples
generally come from `Dates.monthday`.

# Examples
```julia-repl
julia> monthday_passed((12, 31), (12, 15))
true
```
"""
monthday_passed(curr, avail) = (curr[1] >= avail[1]) && (curr[2] >= avail[2])

"Get dynamic and cyclic NetCDF input"
function load_dynamic_input!(model)
    update_forcing!(model)
    if haskey(model.config.input, "cyclic")
        update_cyclic!(model)
    end
end

"Get cyclic NetCDF input for the given time"
function update_cyclic!(model)
    @unpack vertical, clock, reader, network, config = model
    @unpack cyclic_dataset, cyclic_times, cyclic_parameters = reader
    #sel = network.land.indices

    # pick up the data that is valid for the past 24 hours
    month_day = monthday(clock.time - Day(1))

    is_first_timestep = clock.iteration == 1
    if is_first_timestep || (month_day in cyclic_times)
        # time for an update of the cyclic forcing
        i = findlast(t -> monthday_passed(month_day, t), cyclic_times)
        isnothing(i) && error("Could not find applicable cyclic timestep for $month_day")

        # load from NetCDF into the model according to the mapping
        for (par, ncvar) in cyclic_parameters
            data = get_at(cyclic_dataset, ncvar.name, i)
            param_vector = param(model, par)
            sel = active_indices(network, par)
            param_vector .= data[sel]
            if ncvar.scale != 1.0 || ncvar.offset != 0.0
                param_vector .= param_vector .* ncvar.scale .+ ncvar.offset
            end
        end
    end
end

"""
    nc_handles::Dict{String, NCDataset{Nothing}}

For each NetCDF file that will be opened for writing, store an entry in this Dict from the
absolute path of the file to the NCDataset. This allows us to close the NCDataset if we try
to create them twice in the same session, and thus providing a workaround for this issue:
https://github.com/Alexander-Barth/NCDatasets.jl/issues/106

Note that using this will prevent automatic garbage collection and thus closure of the
NCDataset.
"""
const nc_handles = Dict{String,NCDataset{Nothing}}()

"Safely create a NetCDF file, even if it has already been opened for creation"
function create_tracked_netcdf(path)
    abs_path = abspath(path)
    # close existing NCDataset if it exists
    if haskey(nc_handles, abs_path)
        # fine if it was already closed
        close(nc_handles[abs_path])
    end
    # create directory if needed
    mkpath(dirname(path))
    ds = NCDataset(path, "c")
    nc_handles[abs_path] = ds
    return ds
end

"prepare an output dataset for scalar data"
function setup_scalar_netcdf(
    path,
    ncvars,
    modelmap,
    calendar,
    time_units,
    extra_dim,
    config,
    float_type = Float32,
)
    ds = create_tracked_netcdf(path)
    defDim(ds, "time", Inf)  # unlimited
    defVar(
        ds,
        "time",
        Float64,
        ("time",),
        attrib = ["units" => time_units, "calendar" => calendar],
    )
    set_extradim_netcdf(ds, extra_dim)
    for (nc, netcdfvars) in zip(ncvars, config.netcdf.variable)
        # Delft-FEWS requires the attribute :cf_role = "timeseries_id" when a NetCDF file 
        # contains more than one location list 
        defVar(
            ds,
            nc.location_dim,
            nc.locations,
            (nc.location_dim,),
            attrib = ["cf_role" => "timeseries_id"],
        )
        v = param(modelmap, nc.par)
        if eltype(v) <: AbstractFloat
            defVar(
                ds,
                nc.var,
                float_type,
                (nc.location_dim, "time"),
                attrib = ["_FillValue" => float_type(NaN)],
            )
        elseif eltype(v) <: SVector
            if haskey(netcdfvars, extra_dim.name)
                # `extra_dim.name` as specified in the TOML file is used to index
                defVar(
                    ds,
                    nc.var,
                    float_type,
                    (nc.location_dim, "time"),
                    attrib = ["_FillValue" => float_type(NaN)],
                )
            else
                defVar(
                    ds,
                    nc.var,
                    float_type,
                    (nc.location_dim, extra_dim.name, "time"),
                    attrib = ["_FillValue" => float_type(NaN)],
                )
            end
        else
            error("Unsupported output type: ", typeof(v))
        end
    end
    return ds
end

"set extra dimension in output NetCDF file"
function set_extradim_netcdf(
    ds,
    extra_dim::NamedTuple{(:name, :value),Tuple{String,Vector{T}}} where T <: Union{String, Float64},
)
    # the axis attribute `Z` is required to import this type of 3D data by Delft-FEWS the
    # values of this dimension `extra_dim.value` should be of type Float64
    if extra_dim.name == "layer"
        attributes =
            ["long_name" => "layer_index", "standard_name" => "layer_index", "axis" => "Z"]
    elseif extra_dim.name == "classes"
        attributes =
            ["long_name" => extra_dim.name, "standard_name" => extra_dim.name, "axis" => "Z"]
    end
    defVar(ds, extra_dim.name, extra_dim.value, (extra_dim.name,), attrib = attributes)
    return nothing
end

set_extradim_netcdf(ds, extra_dim::Nothing) = nothing

"prepare an output dataset for grid data"
function setup_grid_netcdf(
    path,
    ncx,
    ncy,
    parameters,
    calendar,
    time_units,
    extra_dim,
    sizeinmetres;
    float_type = Float32,
    deflatelevel = 0,
)

    ds = create_tracked_netcdf(path)
    defDim(ds, "time", Inf)  # unlimited
    if sizeinmetres
        defVar(
            ds,
            "x",
            ncx,
            ("x",),
            attrib = [
                "_FillValue" => NaN,
                "long_name" => "x coordinate of projection",
                "standard_name" => "projection_x_coordinate",
                "axis" => "X",
                "units" => "m",
            ],
            deflatelevel = deflatelevel,
        )
        defVar(
            ds,
            "y",
            ncy,
            ("y",),
            attrib = [
                "_FillValue" => NaN,
                "long_name" => "y coordinate of projection",
                "standard_name" => "projection_y_coordinate",
                "axis" => "Y",
                "units" => "m",
            ],
            deflatelevel = deflatelevel,
        )

    else
        defVar(
            ds,
            "lon",
            ncx,
            ("lon",),
            attrib = [
                "_FillValue" => NaN,
                "long_name" => "longitude",
                "standard_name" => "longitude",
                "axis" => "X",
                "units" => "degrees_east",
            ],
        )
        defVar(
            ds,
            "lat",
            ncy,
            ("lat",),
            attrib = [
                "_FillValue" => NaN,
                "long_name" => "latitude",
                "standard_name" => "latitude",
                "axis" => "Y",
                "units" => "degrees_north",
            ],
            deflatelevel = deflatelevel,
        )
    end
    set_extradim_netcdf(ds, extra_dim)
    defVar(
        ds,
        "time",
        Float64,
        ("time",),
        attrib = ["units" => time_units, "calendar" => calendar],
        deflatelevel = deflatelevel,
    )
    if sizeinmetres
        for (key, val) in pairs(parameters)
            if eltype(val.vector) <: AbstractFloat
                # all floats are saved as Float32
                defVar(
                    ds,
                    key,
                    float_type,
                    ("x", "y", "time"),
                    attrib = ["_FillValue" => float_type(NaN)],
                    deflatelevel = deflatelevel,
                )
            elseif eltype(val.vector) <: SVector
                # for SVectors an additional dimension (`extra_dim`) is required
                defVar(
                    ds,
                    key,
                    float_type,
                    ("x", "y", extra_dim.name, "time"),
                    attrib = ["_FillValue" => float_type(NaN)],
                    deflatelevel = deflatelevel,
                )
            else
                error("Unsupported output type: ", typeof(val.vector))
            end
        end
    else
        for (key, val) in pairs(parameters)
            if eltype(val.vector) <: AbstractFloat
                # all floats are saved as Float32
                defVar(
                    ds,
                    key,
                    float_type,
                    ("lon", "lat", "time"),
                    attrib = ["_FillValue" => float_type(NaN)],
                    deflatelevel = deflatelevel,
                )
            elseif eltype(val.vector) <: SVector
                # for SVectors an additional dimension (`extra_dim`) is required
                defVar(
                    ds,
                    key,
                    float_type,
                    ("lon", "lat", extra_dim.name, "time"),
                    attrib = ["_FillValue" => float_type(NaN)],
                    deflatelevel = deflatelevel,
                )
            else
                error("Unsupported output type: ", typeof(val.vector))
            end
        end
    end
    return ds
end

"Add a new time to the unlimited time dimension, and return the index"
function add_time(ds, time)
    i = length(ds["time"]) + 1
    ds["time"][i] = time
    return i
end

struct NCReader
    dataset::CFDataset
    cyclic_dataset::Union{NCDataset,Nothing}
    cyclic_times::Vector{Tuple{Int,Int}}
    forcing_parameters::Dict{Tuple{Symbol,Vararg{Symbol}},NamedTuple}
    cyclic_parameters::Dict{Tuple{Symbol,Vararg{Symbol}},NamedTuple}
end

struct Writer
    dataset::Union{NCDataset,Nothing}           # dataset (NetCDF) for grid data
    parameters::Dict{String,Any}                # mapping of NetCDF variable names to model parameters (arrays)
    nc_path::Union{String,Nothing}              # path NetCDF file (grid data)
    csv_path::Union{String,Nothing}             # path of CSV file
    csv_cols::Vector                            # model parameter (arrays) and associated reducer function for CSV output
    csv_io::IO                                  # file handle to CSV file
    state_dataset::Union{NCDataset,Nothing}     # dataset with model states (NetCDF)
    state_parameters::Dict{String,Any}          # mapping of NetCDF variable names to model states (arrays)
    state_nc_path::Union{String,Nothing}        # path NetCDF file with states
    dataset_scalar::Union{NCDataset,Nothing}    # dataset(NetCDF) for scalar data
    nc_scalar::Vector                           # model parameter (arrays) and associated reducer function for NetCDF scalar output
    ncvars_dims::Vector                         # model parameter (String) and associated NetCDF variable, location dimension and location name for scalar data
    nc_scalar_path::Union{String,Nothing}       # path NetCDF file (scalar data)
    extra_dim::Union{NamedTuple,Nothing}        # name and values for extra dimension (to store SVectors) 
end

function prepare_reader(config)
    path_forcing = config.input.path_forcing
    abspath_forcing = input_path(config, path_forcing)
    cyclic_path = input_path(config, config.input.path_static)

    @info "Cyclic parameters are provided by `$cyclic_path`."

    # absolute paths are not supported, see Glob.jl#2
    # the path separator in a glob pattern is always /
    if isabspath(path_forcing)
        parts = splitpath(path_forcing)
        # use the root/drive as the dir, to support * in directory names as well
        glob_dir = parts[1]
        glob_path = join(parts[2:end], '/')
    else
        tomldir = dirname(config)
        dir_input = get(config, "dir_input", ".")
        glob_dir = normpath(tomldir, dir_input)
        glob_path = replace(path_forcing, '\\' => '/')
    end
    @info "Forcing parameters are provided by `$abspath_forcing`."

    dynamic_paths = glob(glob_path, glob_dir)  # expand "data/forcing-year-*.nc"
    if isempty(dynamic_paths)
        error("No files found with name '$glob_path' in '$glob_dir'")
    end
    dataset = NCDataset(dynamic_paths, aggdim = "time", deferopen = false)

    # check for cyclic parameters
    do_cyclic = haskey(config.input, "cyclic")

    # TODO:include leaf_area_index climatology in update() vertical SBM model
    # we currently assume the same dimension ordering as the forcing
    if do_cyclic == true
        cyclic_dataset = NCDataset(cyclic_path)
        cyclic_nc_times = collect(cyclic_dataset["time"])
        cyclic_times = timecycles(cyclic_nc_times)
    else
        cyclic_dataset = nothing
        cyclic_times = Tuple{Int,Int}[]
    end

    # create map from internal location to NetCDF variable name for forcing parameters
    forcing_parameters = Dict{Tuple{Symbol,Vararg{Symbol}},NamedTuple}()
    for par in config.input.forcing
        fields = symbols(par)
        ncname, mod = ncvar_name_modifier(param(config.input, fields))
        forcing_parameters[fields] =
            (name = ncname, scale = mod.scale, offset = mod.offset, value = mod.value)

        @info "Set `$par` using NetCDF variable `$ncname` as forcing parameter."
    end

    # create map from internal location to NetCDF variable name for cyclic parameters
    if do_cyclic == true
        cyclic_parameters = Dict{Tuple{Symbol,Vararg{Symbol}},NamedTuple}()
        for par in config.input.cyclic
            fields = symbols(par)
            ncname, mod = ncvar_name_modifier(param(config.input, fields))
            cyclic_parameters[fields] =
                (name = ncname, scale = mod.scale, offset = mod.offset)

            @info "Set `$par` using NetCDF variable `$ncname` as cyclic parameter."
        end
    else
        cyclic_parameters = Dict{Tuple{Symbol,Vararg{Symbol}},NamedTuple}()
    end

    # check if there is overlap
    if do_cyclic == true
        forcing_keys = keys(forcing_parameters)
        cyclic_keys = keys(cyclic_parameters)
        for forcing_key in forcing_keys
            if forcing_key in cyclic_keys
                error("parameter specified in both forcing and cyclic")
            end
        end
    end

    return NCReader(
        dataset,
        cyclic_dataset,
        cyclic_times,
        forcing_parameters,
        cyclic_parameters,
    )
end

"Get a Vector of all unique location ids from a 2D map"
function locations_map(ds, mapname, config)
    map_2d = ncread(
        ds,
        config,
        mapname;
        optional = false,
        type = Union{Int,Missing},
        allow_missing = true,
    )
    ids = unique(skipmissing(map_2d))
    return ids
end

"Get a Vector{Tuple} with model parameter and associated NetCDF variable name, dimension and location names for scalar data"
function nc_variables_dims(nc_variables, dataset, config)
    ncvars_dims = []
    for nc_var in nc_variables
        var = nc_var["name"]::String
        par = nc_var["parameter"]::String
        if haskey(nc_var, "map")
            mapname = nc_var["map"]
            ids = string.(locations_map(dataset, mapname, config))
            location_dim = string(var, '_', nc_var["map"])
            push!(
                ncvars_dims,
                (par = par, var = var, location_dim = location_dim, locations = ids),
            )
        else
            push!(
                ncvars_dims,
                (
                    par = par,
                    var = var,
                    location_dim = nc_var["location"],
                    locations = [nc_var["location"]],
                ),
            )
        end
    end
    return ncvars_dims
end

"Get a Vector{String} of all columns names for the CSV header, except the first, time"
function csv_header(cols, dataset, config)
    header = [col["header"] for col in cols]
    header = String[]
    for col in cols
        h = col["header"]::String
        if haskey(col, "map")
            mapname = col["map"]
            ids = locations_map(dataset, mapname, config)
            hvec = [string(h, '_', id) for id in ids]
            append!(header, hvec)
        else
            push!(header, h)
        end
    end
    return header
end

"""
Flatten a nested dictionary, keeping track of the full address of the keys.
Useful for converting TOML of the format:

    [a.b]
    field_of_b = "name_in_output"

    [a.d.c]
    field_of_c = "other_name_in_output"

to a non-nested format:

Dict(
    symbols"a.b.field_of_b" => "name_in_output,
    symbols"a.d.c.field_of_c" => "other_name_in_output,
)
"""
function flat!(d, path, el::Dict)
    for (k, v) in pairs(el)
        flat!(d, string(path, '.', k), v)
    end
    return d
end

function flat!(d, path, el)
    k = symbols(path)
    d[k] = el
    return d
end

"""
    ncnames(dict)

Create a flat mapping from internal parameter locations to NetCDF variable names.

Ignores top level values in the Dict. This function is used to convert a TOML such as:

```toml
[output]
path = "path/to/file.nc"

[output.vertical]
canopystorage = "my_canopystorage"

[output.lateral.river]
q = "my_q"
```

To a dictionary of the flattened parameter locations and NetCDF names. The top level
values are ignored since the output path is not a NetCDF name.

```julia
Dict(
    (:vertical, :canopystorage) => "my_canopystorage,
    (:lateral, :river, :q) => "my_q,
)
```
"""
function ncnames(dict)
    ncnames_dict = Dict{Tuple{Symbol,Vararg{Symbol}},String}()
    for (k, v) in dict
        if v isa Dict  # ignore top level values (e.g. output.path)
            flat!(ncnames_dict, k, v)
        end
    end
    return ncnames_dict
end

"""
    out_map(ncnames_dict, modelmap)

Create a Dict that maps parameter NetCDF names to arrays in the Model.
"""
function out_map(ncnames_dict, modelmap)
    output_map = Dict{String,Any}()
    for (par, ncname) in ncnames_dict
        A = param(modelmap, par)
        output_map[ncname] = (par = par, vector = A)
    end
    return output_map
end


function get_reducer_func(col, rev_inds, args...)
    parameter = col["parameter"]
    if occursin("reservoir", parameter)
        reducer_func = reducer(col, rev_inds.reservoir, args...)
    elseif occursin("lake", parameter)
        reducer_func = reducer(col, rev_inds.lake, args...)
    elseif occursin("river", parameter)
        reducer_func = reducer(col, rev_inds.river, args...)
    elseif occursin("drain", parameter)
        reducer_func = reducer(col, rev_inds.drain, args...)
    else
        reducer_func = reducer(col, rev_inds.land, args...)
    end
end

function prepare_writer(
    config,
    modelmap,
    rev_inds,
    x_nc,
    y_nc,
    nc_static;
    extra_dim = nothing,
)
    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool

    calendar = get(config, "calendar", "standard")::String
    time_units = get(config, "time_units", CFTime.DEFAULT_TIME_UNITS)

    # create an output NetCDF that will hold all timesteps of selected parameters for grid
    # data but only if config.output.path has been set
    if haskey(config, "output") && haskey(config.output, "path")
        nc_path = output_path(config, config.output.path)
        deflatelevel = get(config.output, "compressionlevel", 0)::Int
        @info "Create an output NetCDF file `$nc_path` for grid data, using compression level `$deflatelevel`."
        # create a flat mapping from internal parameter locations to NetCDF variable names
        output_ncnames = ncnames(config.output)
        # fill the output_map by mapping parameter NetCDF names to arrays
        output_map = out_map(output_ncnames, modelmap)
        ds = setup_grid_netcdf(
            nc_path,
            x_nc,
            y_nc,
            output_map,
            calendar,
            time_units,
            extra_dim,
            sizeinmetres,
            deflatelevel = deflatelevel
        )
    else
        nc_path = nothing
        output_map = Dict{String,Any}()
        ds = nothing
    end

    # create a separate state output NetCDF that will hold the last timestep of all states
    # but only if config.state.path_output has been set
    if haskey(config, "state") && haskey(config.state, "path_output")
        state_ncnames = ncnames(config.state)
        state_map = out_map(state_ncnames, modelmap)
        nc_state_path = output_path(config, config.state.path_output)
        @info "Create a state output NetCDF file `$nc_state_path`."
        ds_outstate = setup_grid_netcdf(
            nc_state_path,
            x_nc,
            y_nc,
            state_map,
            calendar,
            time_units,
            extra_dim,
            sizeinmetres;
            float_type = Float,
        )
    else
        ds_outstate = nothing
        state_map = Dict{String,Any}()
        nc_state_path = nothing
    end

    # create an output NetCDF that will hold all timesteps of selected parameters for scalar
    # data, but only if config.netcdf.variable has been set. 
    if haskey(config, "netcdf") && haskey(config.netcdf, "variable")
        nc_scalar_path = output_path(config, config.netcdf.path)
        @info "Create an output NetCDF file `$nc_state_path` for scalar data."
        # get NetCDF info for scalar data (variable name, locationset (dim) and
        # location ids)
        ncvars_dims = nc_variables_dims(config.netcdf.variable, nc_static, config)
        ds_scalar = setup_scalar_netcdf(
            nc_scalar_path,
            ncvars_dims,
            modelmap,
            calendar,
            time_units,
            extra_dim,
            config,
        )
        # create a vector of (parameter, reducer) named tuples which will be used to
        # retrieve and reduce data during a model run
        nc_scalar = []
        for var in config.netcdf.variable
            parameter = var["parameter"]
            reducer_func =
                get_reducer_func(var, rev_inds, x_nc, y_nc, config, nc_static, "NetCDF")
            push!(nc_scalar, (parameter = parameter, reducer = reducer_func))
        end
    else
        ds_scalar = nothing
        nc_scalar = []
        ncvars_dims = []
        nc_scalar_path = nothing
    end

    if haskey(config, "csv") && haskey(config.csv, "column")
        # open CSV file and write header
        csv_path = output_path(config, config.csv.path)
        @info "Create an output CSV file `$csv_path` for scalar data."
        # create directory if needed
        mkpath(dirname(csv_path))
        csv_io = open(csv_path, "w")
        print(csv_io, "time,")
        header = csv_header(config.csv.column, nc_static, config)
        println(csv_io, join(header, ','))
        flush(csv_io)

        # create a vector of (parameter, reducer) named tuples which will be used to
        # retrieve and reduce the CSV data during a model run
        csv_cols = []
        for col in config.csv.column
            parameter = col["parameter"]
            reducer_func =
                get_reducer_func(col, rev_inds, x_nc, y_nc, config, nc_static, "CSV")
            push!(csv_cols, (parameter = parameter, reducer = reducer_func))
        end
    else
        # no CSV file is checked by isnothing(csv_path)
        csv_path = nothing
        csv_cols = []
        csv_io = IOBuffer()
    end


    return Writer(
        ds,
        output_map,
        nc_path,
        csv_path,
        csv_cols,
        csv_io,
        ds_outstate,
        state_map,
        nc_state_path,
        ds_scalar,
        nc_scalar,
        ncvars_dims,
        nc_scalar_path,
        extra_dim,
    )
end

"Write a new timestep with scalar data to a NetCDF file"
function write_netcdf_timestep(model, dataset)
    @unpack writer, clock, config = model

    time_index = add_time(dataset, clock.time)
    for (nt, nc) in zip(writer.nc_scalar, config.netcdf.variable)
        A = param(model, nt.parameter)
        elemtype = eltype(A)
        # could be a value, or a vector in case of map
        if elemtype <: AbstractFloat
            v = nt.reducer(A)
            dataset[nc["name"]][:, time_index] .= v
        elseif elemtype <: SVector
            # check if an extra dimension and index is specified in the TOML file 
            if haskey(nc, writer.extra_dim.name)
                i = get_index_dimension(nc, model)
                v = nt.reducer(getindex.(A, i))
                dataset[nc["name"]][:, time_index] .= v
            else
                nlayer = length(first(A))
                for i = 1:nlayer
                    v = nt.reducer(getindex.(A, i))
                    dataset[nc["name"]][:, i, time_index] .= v
                end
            end
        else
            error("Unsupported output type: ", elemtype)
        end
    end
    return model
end

"Write a new timestep with grid data to a NetCDF file"
function write_netcdf_timestep(model, dataset, parameters)
    @unpack vertical, clock, reader, network = model

    time_index = add_time(dataset, clock.time)

    buffer = zeros(Union{Float,Missing}, size(model.network.land.reverse_indices))
    for (key, val) in parameters
        @unpack par, vector = val
        sel = active_indices(network, par)
        # write the active cells vector to the 2d buffer matrix
        elemtype = eltype(vector)
        if elemtype <: AbstractFloat
            # ensure no other information is written
            fill!(buffer, missing)
            # cut off possible boundary conditions/ ghost points with [1:length(sel)]
            buffer[sel] .= vector[1:length(sel)]
            dataset[key][:, :, time_index] = buffer
        elseif elemtype <: SVector
            nlayer = length(first(vector))
            for i = 1:nlayer
                # ensure no other information is written
                fill!(buffer, missing)
                buffer[sel] .= getindex.(vector, i)
                dataset[key][:, :, i, time_index] = buffer
            end
        else
            error("Unsupported output type: ", elemtype)
        end
    end

    return model
end

# don't do anything for no dataset, used if no output NetCDF is needed
write_netcdf_timestep(model, dataset::Nothing, parameters) = model
write_netcdf_timestep(model, dataset::Nothing) = model

"Write model output"
function write_output(model)
    @unpack vertical, clock, reader, network, writer = model
    @unpack dataset, dataset_scalar, parameters = writer

    write_csv_row(model)
    write_netcdf_timestep(model, dataset, parameters)
    write_netcdf_timestep(model, dataset_scalar)

    return model
end

"""
    timecycles(times)

Given a vector of times, return a tuple of (month, day) for each time entry, to use as a
cyclic time series that repeats every year. By using `monthday` rather than `dayofyear`,
leap year offsets are avoided.

It can generate such a series from eiher TimeTypes given that the year is constant, or
it will interpret integers as either months or days of year if possible.
"""
function timecycles(times)
    if eltype(times) <: TimeType
        # all timestamps are from the same year
        year1 = year(first(times))
        if !all(==(year1), year.(times))
            error("unsupported cyclic timeseries")
        end
        # returns a (month, day) tuple for each date
        return monthday.(times)
    elseif eltype(times) <: Integer
        if length(times) == 12
            months = Date(2000, 1, 1):Month(1):Date(2000, 12, 31)
            return monthday.(months)
        elseif length(times) == 365
            days = Date(2001, 1, 1):Day(1):Date(2001, 12, 31)
            return monthday.(days)
        elseif length(times) == 366
            days = Date(2000, 1, 1):Day(1):Date(2000, 12, 31)
            return monthday.(days)
        else
            error("unsupported cyclic timeseries")
        end
    else
        error("unsupported cyclic timeseries")
    end
end

"Close input and output datasets that are opened on model initialization"
function close_files(model; delete_output::Bool = false)
    @unpack reader, writer, config = model

    close(reader.dataset)
    if haskey(config.input, "cyclic")
        close(reader.cyclic_dataset)
    end
    writer.dataset === nothing || close(writer.dataset)
    writer.dataset_scalar === nothing || close(writer.dataset_scalar)
    close(writer.csv_io)  # can be an IOBuffer
    writer.state_dataset === nothing || close(writer.state_dataset)

    if delete_output
        if writer.nc_path !== nothing
            isfile(writer.nc_path) && rm(writer.nc_path)
        end
        if writer.csv_path !== nothing
            isfile(writer.csv_path) && rm(writer.csv_path)
        end
        if writer.state_nc_path !== nothing
            isfile(writer.state_nc_path) && rm(writer.state_nc_path)
        end
        if writer.nc_scalar_path !== nothing
            isfile(writer.nc_scalar_path) && rm(writer.nc_scalar_path)
        end
    end
    return nothing
end

"Mapping from reducer strings in the TOML to functions"
function reducerfunction(reducer::AbstractString)
    functionmap = Dict{String,Function}(
        "maximum" => maximum,
        "minimum" => minimum,
        "mean" => mean,
        "median" => median,
        "first" => first,
        "last" => last,
        "only" => only,
    )
    f = get(functionmap, reducer, nothing)
    isnothing(f) && error("unknown reducer")
    return f
end

"Get a reducer function based on output settings for scalar data defined in a dictionary"
function reducer(col, rev_inds, x_nc, y_nc, config, dataset, fileformat)
    param = col["parameter"]
    if haskey(col, "map")
        # assumes the parameter in "map" has a 2D input map, with
        # integers indicating the points or zones that are to be aggregated
        mapname = col["map"]
        # if no reducer is given, pick "only", this is the only safe reducer,
        # and makes sense in the case of a gauge map
        reducer_name = get(col, "reducer", "only")
        f = reducerfunction(reducer_name)
        map_2d = ncread(
            dataset,
            config,
            mapname;
            optional = false,
            type = Union{Int,Missing},
            allow_missing = true,
        )
        @info "Adding scalar output for a map with a reducer function." fileformat param mapname reducer_name
        ids = unique(skipmissing(map_2d))
        # from id to list of internal indices
        inds = Dict{Int,Vector{Int}}(id => Vector{Int}() for id in ids)
        for i in eachindex(map_2d)
            v = map_2d[i]
            ismissing(v) && continue
            v::Int
            vector = inds[v]
            ind = rev_inds[i]
            if iszero(ind)
                error("""inactive cell found in requested scalar output
                    map `$mapname` value $v for parameter $param""")
            end
            push!(vector, ind)
        end
        return A -> (f(A[inds[id]]) for id in ids)
    elseif haskey(col, "reducer")
        # reduce over all active cells
        # needs to be behind the map if statement, because it also can use a reducer
        reducer_name = col["reducer"]
        @info "Adding scalar output of all active cells with reducer function." fileformat param reducer_name
        return reducerfunction(reducer_name)
    elseif haskey(col, "index")
        index = col["index"]
        if index isa Int
            # linear index into the internal vector of active cells
            # this one mostly makes sense for debugging, or for vectors of only a few elements
            @info "Adding scalar output for linear index." fileformat param index
            return x -> getindex(x, index)
        elseif index isa Dict
            # index into the 2D input/output arrays
            # the first always corresponds to the x dimension, then the y dimension
            # this is 1-based
            i = index["x"]::Int
            j = index["y"]::Int
            ind = rev_inds[i, j]
            @info "Adding scalar output for 2D index." fileformat param index
            iszero(ind) && error("inactive loc specified for output")
            return A -> getindex(A, ind)
        else
            error("unknown index used")
        end
    elseif haskey(col, "coordinate")
        x = col["coordinate"]["x"]::Float64
        y = col["coordinate"]["y"]::Float64
        # find the closest cell center index
        _, iy = findmin(abs.(y_nc .- y))
        _, ix = findmin(abs.(x_nc .- x))
        I = CartesianIndex(ix, iy)
        i = rev_inds[I]
        @info "Adding scalar output for coordinate." fileformat param x y
        iszero(i) && error("inactive coordinate specified for output")
        return A -> getindex(A, i)
    else
        error("unknown reducer")
    end
end

function write_csv_row(model)
    @unpack writer, clock, config = model
    isnothing(writer.csv_path) && return nothing
    io = writer.csv_io
    print(io, string(clock.time))
    for (nt, col) in zip(writer.csv_cols, config.csv.column)
        A = param(model, nt.parameter)
        # v could be a value, or a vector in case of map
        if eltype(A) <: SVector
            # indexing is required in case of a SVector and CSV output
            i = get_index_dimension(col, model)
            v = nt.reducer(getindex.(A, i))
        else
            v = nt.reducer(A)
        end
        # numbers are also iterable
        for el in v
            print(io, ',', el)
        end
    end
    println(io)
end

"From a time and a calendar, create the right CFTime DateTimeX type"
function cftime(time, calendar)
    timetype = CFTime.timetype(calendar)
    # invalid Gregorian dates like 2020-02-30 cannot be made into a DateTime
    # even though they may exist in the target calendar
    # though this constructor does offer support for string formatted dates
    dt = DateTime(time)
    cal_time = reinterpret(timetype, dt)
    return cal_time
end

function reset_clock!(clock::Clock, config)
    new_clock = Clock(config)
    # we want this method to be mutating
    clock.time = new_clock.time
    clock.iteration = new_clock.iteration
    clock.Δt = new_clock.Δt
    return clock
end

function advance!(clock)
    clock.iteration += 1
    clock.time += clock.Δt
    return clock
end

function rewind!(clock)
    clock.iteration -= 1
    clock.time -= clock.Δt
    return clock
end

"Read a rating curve from CSV into a NamedTuple of vectors"
function read_sh_csv(path)
    data, header = readdlm(path, ',', Float, header = true)
    names = vec(uppercase.(header))
    idx_h = findfirst(==("H"), names)
    idx_s = findfirst(==("S"), names)

    if isnothing(idx_h) || isnothing(idx_s)
        error("$path needs to provide H and S columns, got $names")
    end

    return (H = data[:, idx_h], S = data[:, idx_s])
end

"Read a specific storage curve from CSV into a NamedTuple of vectors"
function read_hq_csv(path)
    data = readdlm(path, ',', Float, skipstart = 1)
    # Q is a matrix with 365 columns, one for each day in the year
    return (H = data[:, 1], Q = data[:, 2:end])
end

# these represent the type of the rating curve and specific storage data
const SH = NamedTuple{(:H, :S),Tuple{Vector{Float},Vector{Float}}}
const HQ = NamedTuple{(:H, :Q),Tuple{Vector{Float},Matrix{Float}}}

is_increasing(v) = last(v) > first(v)

"""
    nc_dim_name(dims::Vector{Symbol}, name::Symbol)
    nc_dim_name(ds::CFDataset, name::Symbol)

Given a NetCDF dataset or list of dimensions, and an internal dimension name, return the
corresponding NetCDF dimension name. Certain common alternatives are supported, e.g. :lon or
:longitude instead of :x.
"""
function nc_dim_name(dims::Vector{Symbol}, name::Symbol)
    # direct naming
    name in dims && return name
    # list of common alternative names
    if name == :x
        for candidate in (:lon, :longitude)
            candidate in dims && return candidate
        end
        error("No x dimension candidate found")
    elseif name == :y
        for candidate in (:lat, :latitude)
            candidate in dims && return candidate
        end
        error("No y dimension candidate found")
    else
        error("Unknown dimension $name")
    end
end

nc_dim_name(ds::CFDataset, name::Symbol) = nc_dim_name(Symbol.(keys(ds.dim)), name)

"""
    nc_dim(ds::CFDataset, name::Symbol)

Return the dimension coordinate, based on the internal name (:x, :y, :`extra_dim.name`,
:time), `extra_dim` depends on the model type, which will map to the correct NetCDF name
using `nc_dim_name`.
"""
nc_dim(ds::CFDataset, name) = ds[nc_dim_name(ds, name)]


"""
    internal_dim_name(name::Symbol)

Given a NetCDF dimension name string, return the corresponding internal dimension name.
"""
function internal_dim_name(name::Symbol)
    if name in (:x, :lon, :longitude)
        return :x
    elseif name in (:y, :lat, :latitude)
        return :y
    elseif name in (:time, :layer, :classes)
        return name
    else
        error("Unknown dimension $name")
    end
end

"""
    read_dims(A::NCDatasets.CFVariable, dim_sel)

Return the data of a NetCDF data variable as an Array. Only dimensions in `dim_sel`, a
NamedTuple like (x=:, y=:, time=1). Other dimensions that may be present need to be size 1,
otherwise an error is thrown.

`dim_sel` keys should be the internal dimension names; :x, :y, :time, :`extra_dim.name`.
"""
function read_dims(A::NCDatasets.CFVariable, dim_sel::NamedTuple)
    dimsizes = dimsize(A)
    indexer = []
    data_dim_order = Symbol[]
    # need to iterate in this order, to get the indices in order
    for (dim_name, dim_size) in pairs(dimsizes)
        if dim_size == 1
            push!(indexer, 1)
        else
            dim = internal_dim_name(dim_name)
            # each dim should either be in dim_sel or be size 1
            if dim in keys(dim_sel)
                idx = dim_sel[dim]
                push!(indexer, idx)
                # Would be nice to generalize this to anything that drops the dimension
                # when used to index. Not sure how we could support `begin` or `end` here.
                if !(idx isa Int)
                    push!(data_dim_order, dim)
                end
            else
                throw(ArgumentError("""NetCDF dimension $dim_name has length $dim_size.
                    Only extra dimensions of length 1 are supported."""))
            end
        end
    end
    data = A[indexer...]
    order = Tuple(data_dim_order)
    @assert ndims(data) == length(order) "dimensions not properly recorded"
    return data, order
end

"""
    dim_directions(ds, dim_names)

Given a NCDataset and a list of internal dimension names, return a NamedTuple, which maps
the dimension name to `true` if it is increasing, or `false` otherwise.

For the layer dimension, we allow coordinate arrays to be missing, in which case we
consider it increasing, going from the top layer (1) to deeper layers. This is to keep
accepting data that we have accepted before.
"""
function dim_directions(ds::CFDataset, dim_names)
    pairs = Pair{Symbol,Bool}[]
    for d in dim_names
        if d == :layer && !(haskey(ds, "layer"))
            inc = true
        else
            inc = is_increasing(nc_dim(ds, d))
        end
        push!(pairs, d => inc)
    end
    return NamedTuple(pairs)
end

"""
    permute_data(data, dim_names)

Given an Array of data, and a list of its dimension names, return a permuted version of the
data such that the dimension order will follow (:x, :y, `extra_dim.name`). No permutation is
done if this is not needed.
"""
function permute_data(data, dim_names)
    @assert ndims(data) == length(dim_names)
    if :layer in dim_names
        desired_order = (:x, :y, :layer)
    elseif :classes in dim_names
        desired_order = (:x, :y, :classes)
    else
        desired_order = (:x, :y)
    end
    @assert dim_names ⊆ desired_order
    if length(dim_names) == 2
        if first(dim_names) == :x
            return data, dim_names
        else
            return permutedims(data), reverse(dim_names)
        end
    elseif length(dim_names) == 3
        permutation = [findfirst(==(d), dim_names) for d in desired_order]
        if permutation == (1, 2, 3)
            return data, dim_names
        else
            return permutedims(data, permutation), dim_names[permutation]
        end
    else
        error("Unsupported number of dimensions")
    end
end

"""
    reverse_data!(data, dims_increasing)

Reverse the data as needed, such that it is increasing along all dimensions.
"""
function reverse_data!(data, dims_increasing)
    # for the reverse call it is important that the dims_increasing tuple is ordered in the
    # desired internal ordering, just like the data is after permutation
    if length(dims_increasing) == 2
        dims_increasing_ordered = (x = dims_increasing.x, y = dims_increasing.y)
    elseif length(dims_increasing) == 3 && haskey(dims_increasing, :layer)
        dims_increasing_ordered =
            (x = dims_increasing.x, y = dims_increasing.y, layer = dims_increasing.layer)
    elseif length(dims_increasing) == 3 && haskey(dims_increasing, :classes)
        dims_increasing_ordered = (
            x = dims_increasing.x,
            y = dims_increasing.y,
            classes = dims_increasing.classes,
        )
    else
        error("Unsupported number of dimensions")
    end
    # the non increasing dimensions should be reversed, to make them all increasing
    dims = Tuple(findall(.!values(dims_increasing_ordered)))
    return reverse!(data; dims)
end

"""
    read_standardized(ds::CFDataset, varname::AbstractString, dim_names)

Read the dimensions listed in dim_names from a variable with name `varname` from a NetCDF
dataset `ds`. `dim_sel` should be a NamedTuple like (x=:, y=:, time=1), which will return
a 2 dimensional array with x and y axes, representing the first index in the time dimension.
"""
function read_standardized(ds::CFDataset, varname::AbstractString, dim_sel::NamedTuple)
    data, data_dim_order = read_dims(ds[varname], dim_sel)
    data, new_dim_order = permute_data(data, data_dim_order)
    dims_increasing = dim_directions(ds, new_dim_order)
    data = reverse_data!(data, dims_increasing)
    return data
end

"""
    read_x_axis(ds::CFDataset)

Return the x coordinate Vector{Float64}, whether it is called x, lon or longitude.
Also sorts the vector to be increasing, to match `read_standardized`.
"""
function read_x_axis(ds::CFDataset)::Vector{Float64}
    candidates = ("x", "lon", "longitude")
    for candidate in candidates
        if haskey(ds, candidate)
            return sort!(Float64.(ds[candidate][:]))
        end
    end
    error("no x axis found in $(path(ds))")
end

"""
    read_y_axis(ds::CFDataset)

Return the y coordinate Vector{Float64}, whether it is called y, lat or latitude.
Also sorts the vector to be increasing, to match `read_standardized`.
"""
function read_y_axis(ds::CFDataset)::Vector{Float64}
    candidates = ("y", "lat", "latitude")
    for candidate in candidates
        if haskey(ds, candidate)
            return sort!(Float64.(ds[candidate][:]))
        end
    end
    error("no y axis found in $(path(ds))")
end

"Get `index` for dimension name `layer` or `classes` based on `model`"
function get_index_dimension(var, model)::Int
    @unpack vertical = model
    if haskey(var, "layer")
        inds = collect(1:vertical.maxlayers)
        index = inds[var["layer"]]
    elseif haskey(var,"class")
        index = findfirst(x->x==var["class"], vertical.classes)
    else
        error("Unrecognized or missing dimension name to index $(var)")
    end
    return index
end

"Get `index` for dimension name `layer` or `classes` based on `config` (TOML file)"
function get_index_dimension(var, config::Config, dim_value)::Int
    if haskey(var, "layer")
        v = get(config.model, "thicknesslayers", Float[])
        inds = collect(1:length(v)+1)
        index = inds[dim_value]
    elseif haskey(var,"class")
        classes = get(config.model, "classes", "")
        index = findfirst(x->x==dim_value, classes)
    else
        error("Unrecognized or missing dimension name to index $(var)")
    end
    return index
end