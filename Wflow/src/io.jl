"""Turn "a.aa.aaa" into (:a, :aa, :aaa)"""
symbols(s::AbstractString) = Tuple(Symbol(x) for x in split(s, '.'))

"Get a nested field using a tuple of Symbols"
param(obj, fields::Tuple{Vararg{Symbol}}) = foldl(getproperty, fields; init = obj)
param(obj, fields::AbstractString) = param(obj, symbols(fields))
function param(obj, fields, default)
    try
        return param(obj, fields)
    catch
        return default
    end
end

"Extract a netCDF variable at a given time"
function get_at(
    ds::CFDataset,
    var::InputEntry,
    times::AbstractVector{<:TimeType},
    t::TimeType,
)
    # this behaves like a backward fill interpolation
    i = findfirst(>=(t), times)
    t < first(times) && throw(DomainError("time $t before dataset begin $(first(times))"))
    i === nothing && throw(DomainError("time $t after dataset end $(last(times))"))
    return get_at(ds, var, i)
end

function get_at(ds::CFDataset, var::InputEntry, i)
    (; scale, offset) = var
    data = read_standardized(ds, variable_name(var), (x = :, y = :, time = i))
    if scale != 1.0 || offset != 0.0
        data .= data .* scale .+ offset
    end
    return data
end

function get_param_res(model)
    return Dict(
        "atmosphere_water__precipitation_volume_flux" =>
            model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation,
        "land_surface_water__potential_evaporation_volume_flux" =>
            model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.evaporation,
    )
end

const mover_params = (
    "atmosphere_water__precipitation_volume_flux",
    "land_surface_water__potential_evaporation_volume_flux",
)

function get_param(model, parameter::AbstractString)
    (; land) = model
    lens = get_lens(parameter, land)
    param = lens(model)
    return param
end

function routing_with_reservoirs(model)
    (; config) = model
    return config.model.reservoir__flag
end

function routing_with_reservoirs(model::AbstractModel{<:SedimentModel})
    return false
end

"""
    load_fixed_forcing!(model)
Get fixed netCDF forcing input."
"""
function load_fixed_forcing!(model)
    (; reader, domain) = model
    (; forcing_parameters) = reader

    do_reservoirs = routing_with_reservoirs(model)

    reverse_indices = domain.land.network.reverse_indices
    if do_reservoirs
        sel_reservoirs = domain.reservoir.network.indices_coverage
        param_res = get_param_res(model)
    end

    for (par, ncvar) in forcing_parameters
        if variable_name(ncvar) === nothing
            val = only(ncvar.value) * only(ncvar.scale) + only(ncvar.offset)
            param = get_param(model, par)
            param .= val
            # set fixed precipitation and evaporation over the reservoirs and put these into
            # the reservoir structs and set the precipitation and evaporation to 0 in the
            # land model
            if par in mover_params
                if do_reservoirs
                    for (i, sel_reservoir) in enumerate(sel_reservoirs)
                        param[reverse_indices[sel_reservoir]] .= 0
                        param_res[par][i] = val
                    end
                end
            end
        end
    end
    return nothing
end

"""
    update_forcing!(model::Model)

Get dynamic netCDF input for the given time. Wflow expects `right` labeling of the forcing
time interval, e.g. daily precipitation at 01-02-2000 00:00:00 is the accumulated total
precipitation between 01-01-2000 00:00:00 and 01-02-2000 00:00:00.
"""
function update_forcing!(model)
    (; clock, reader, domain) = model
    (; dataset, dataset_times, forcing_parameters) = reader

    do_reservoirs = routing_with_reservoirs(model)

    if do_reservoirs
        sel_reservoirs = domain.reservoir.network.indices_coverage
        param_res = get_param_res(model)
    end

    # load from netCDF into the model according to the mapping
    for (par, ncvar) in forcing_parameters
        # no need to update fixed values
        variable_name(ncvar) === nothing && continue

        time = convert(eltype(dataset_times), clock.time)
        data = get_at(dataset, ncvar, dataset_times, time)

        # calculate the mean precipitation and evaporation over reservoirs and put these
        # into the reservoir structs and set the precipitation and evaporation to 0 in the
        # land model
        if par in mover_params
            if do_reservoirs
                for (i, sel_reservoir) in enumerate(sel_reservoirs)
                    avg = mean(data[sel_reservoir])
                    data[sel_reservoir] .= 0
                    param_res[par][i] = avg
                end
            end
        end
        sel = active_indices(domain, par)
        # missing data for observed reservoir outflow is allowed at reservoir location(s)
        if par == "reservoir_water__outgoing_observed_volume_flow_rate"
            data_sel = nomissing(data[sel], MISSING_VALUE)
        else
            data_sel = data[sel]
        end
        if any(ismissing, data_sel)
            msg = "Forcing data at $time has missing values on active model cells for $(variable_name(ncvar))"
            throw(ArgumentError(msg))
        end
        param = get_param(model, par)
        param .= data_sel
    end
    return nothing
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

"Get dynamic and cyclic netCDF input"
function load_dynamic_input!(model)
    update_forcing!(model)
    if do_cyclic(model.config)
        update_cyclic!(model)
    end
    return nothing
end

"Get cyclic netCDF input for the given time"
function update_cyclic!(model)
    (; clock, reader, domain) = model
    (; cyclic_dataset, cyclic_times, cyclic_parameters) = reader

    # pick up the data that is valid for the past model time step
    month_day = monthday(clock.time - clock.dt)

    is_first_timestep = clock.iteration == 1

    for (par, ncvar) in cyclic_parameters
        if is_first_timestep || (month_day in cyclic_times[par])
            # time for an update of the cyclic forcing
            i = findlast(t -> monthday_passed(month_day, t), cyclic_times[par])
            isnothing(i) &&
                error("Could not find applicable cyclic timestep for $month_day")
            # load from netCDF into the model according to the mapping
            data = get_at(cyclic_dataset, ncvar, i)
            sel = active_indices(domain, par)
            # missing data for observed reservoir outflow is allowed at reservoir
            # location(s)
            if par == "reservoir_water__outgoing_observed_volume_flow_rate"
                data_sel = nomissing(data[sel], MISSING_VALUE)
            else
                data_sel = data[sel]
            end
            if any(ismissing, data_sel)
                msg = "Cyclic data at month $(month_day[1]) and day $(month_day[2]) has missing values on active model cells for $(variable_name(ncvar))"
                throw(ArgumentError(msg))
            end
            param = get_param(model, par)
            param .= data_sel
        end
    end
    return nothing
end

"""
    NC_HANDLES::Dict{String, NCDataset{Nothing}}

For each netCDF file that will be opened for writing, store an entry in this Dict from the
absolute path of the file to the NCDataset. This allows us to close the NCDataset if we try
to create them twice in the same session, and thus providing a workaround for this issue:
https://github.com/Alexander-Barth/NCDatasets.jl/issues/106

Note that using this will prevent automatic garbage collection and thus closure of the
NCDataset.
"""
const NC_HANDLES = Dict{String, NCDataset{Nothing}}()

"Safely create a netCDF file, even if it has already been opened for creation"
function create_tracked_netcdf(path)
    abs_path = abspath(path)
    # close existing NCDataset if it exists
    if haskey(NC_HANDLES, abs_path)
        # fine if it was already closed
        close(NC_HANDLES[abs_path])
    end
    # create directory if needed
    mkpath(dirname(path))
    ds = NCDataset(path, "c")
    NC_HANDLES[abs_path] = ds
    return ds
end

"prepare an output dataset for scalar data"
function setup_scalar_netcdf(
    path,
    dataset,
    modelmap,
    calendar,
    time_units,
    extra_dim,
    config,
    float_type = Float32,
)
    (; land) = modelmap
    ds = create_tracked_netcdf(path)
    defDim(ds, "time", Inf)  # unlimited
    defVar(
        ds,
        "time",
        Float64,
        ("time",);
        attrib = ["units" => time_units, "calendar" => convert(String, calendar)],
    )
    set_extradim_netcdf(ds, extra_dim)
    for scalar_variable in config.output.netcdf_scalar.variable
        (; map, _location_dim, location, parameter, name) = scalar_variable
        # Delft-FEWS requires the attribute :cf_role = "timeseries_id" when a netCDF file
        # contains more than one location list
        if _location_dim ∉ keys(ds.dim)
            locations =
                isnothing(map) ? [location] : string.(locations_map(dataset, map, config))
            defVar(
                ds,
                _location_dim,
                locations,
                (_location_dim,);
                attrib = ["cf_role" => "timeseries_id"],
            )
        end
        v = if haskey(standard_name_map(land), parameter)
            lens = get_lens(parameter, land)
            lens(modelmap)
        else
            param(modelmap, parameter)
        end
        if eltype(v) <: AbstractFloat
            defVar(
                ds,
                name,
                float_type,
                (_location_dim, "time");
                attrib = ["_FillValue" => float_type(NaN)],
            )
        elseif eltype(v) <: SVector
            if haskey(netcdfvars, extra_dim.name)
                # `extra_dim.name` as specified in the TOML file is used to index
                defVar(
                    ds,
                    name,
                    float_type,
                    (_location_dim, "time");
                    attrib = ["_FillValue" => float_type(NaN)],
                )
            else
                defVar(
                    ds,
                    name,
                    float_type,
                    (_location_dim, extra_dim.name, "time");
                    attrib = ["_FillValue" => float_type(NaN)],
                )
            end
        else
            error("Unsupported output type: ", typeof(v))
        end
    end
    return ds
end

"set extra dimension in output netCDF file"
function set_extradim_netcdf(
    ds,
    extra_dim::NamedTuple{
        (:name, :value),
        Tuple{String, Vector{T}},
    } where {T <: Union{String, Float64}},
)
    # the axis attribute `Z` is required to import this type of 3D data by Delft-FEWS the
    # values of this dimension `extra_dim.value` should be of type Float64
    if extra_dim.name == "layer"
        attributes =
            ["long_name" => "layer_index", "standard_name" => "layer_index", "axis" => "Z"]
    end
    defVar(ds, extra_dim.name, extra_dim.value, (extra_dim.name,); attrib = attributes)
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
    cell_length_in_meter;
    float_type = Float32,
    deflatelevel = 0,
)
    ds = create_tracked_netcdf(path)
    defDim(ds, "time", Inf)  # unlimited
    if cell_length_in_meter
        defVar(
            ds,
            "x",
            ncx,
            ("x",);
            attrib = [
                "long_name" => "x coordinate of projection",
                "standard_name" => "projection_x_coordinate",
                "axis" => "X",
                "units" => "m",
            ],
            deflatelevel,
        )
        defVar(
            ds,
            "y",
            ncy,
            ("y",);
            attrib = [
                "long_name" => "y coordinate of projection",
                "standard_name" => "projection_y_coordinate",
                "axis" => "Y",
                "units" => "m",
            ],
            deflatelevel,
        )

    else
        defVar(
            ds,
            "lon",
            ncx,
            ("lon",);
            attrib = [
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
            ("lat",);
            attrib = [
                "long_name" => "latitude",
                "standard_name" => "latitude",
                "axis" => "Y",
                "units" => "degrees_north",
            ],
            deflatelevel,
        )
    end
    set_extradim_netcdf(ds, extra_dim)
    defVar(
        ds,
        "time",
        Float64,
        ("time",);
        attrib = ["units" => time_units, "calendar" => convert(String, calendar)],
        deflatelevel,
    )
    if cell_length_in_meter
        for (key, val) in pairs(parameters)
            if eltype(val.vector) <: AbstractFloat
                # all floats are saved as Float32
                defVar(
                    ds,
                    key,
                    float_type,
                    ("x", "y", "time");
                    attrib = ["_FillValue" => float_type(NaN)],
                    deflatelevel,
                )
            elseif eltype(val.vector) <: SVector
                # for SVectors an additional dimension (`extra_dim`) is required
                defVar(
                    ds,
                    key,
                    float_type,
                    ("x", "y", extra_dim.name, "time");
                    attrib = ["_FillValue" => float_type(NaN)],
                    deflatelevel,
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
                    ("lon", "lat", "time");
                    attrib = ["_FillValue" => float_type(NaN)],
                    deflatelevel,
                )
            elseif eltype(val.vector) <: SVector
                # for SVectors an additional dimension (`extra_dim`) is required
                defVar(
                    ds,
                    key,
                    float_type,
                    ("lon", "lat", extra_dim.name, "time");
                    attrib = ["_FillValue" => float_type(NaN)],
                    deflatelevel,
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

struct NCReader{T}
    dataset::CFDataset
    dataset_times::Vector{T}
    cyclic_dataset::Union{NCDataset, Nothing}
    cyclic_times::Dict{String, Vector{Tuple{Int, Int}}}
    forcing_parameters::InputEntries
    cyclic_parameters::InputEntries
end

struct Writer
    dataset::Union{NCDataset, Nothing}           # dataset (netCDF) for grid data
    parameters::Dict{String, Any}                # mapping of netCDF variable names to model parameters (arrays)
    nc_path::Union{String, Nothing}              # path netCDF file (grid data)
    csv_path::Union{String, Nothing}             # path of CSV file
    csv_cols::Vector{CSVColumn}                  # model parameter (arrays) for CSV output
    csv_io::IO                                   # file handle to CSV file
    state_dataset::Union{NCDataset, Nothing}     # dataset with model states (netCDF)
    state_parameters::Dict{String, Any}          # mapping of netCDF variable names to model states (arrays)
    state_nc_path::Union{String, Nothing}        # path netCDF file with states
    dataset_scalar::Union{NCDataset, Nothing}    # dataset (netCDF) for scalar data
    nc_scalar::Vector{NetCDFScalarVariable}      # model parameter (arrays) and for netCDF scalar output                        # model parameter (String) and associated netCDF variable, location dimension and location name for scalar data
    nc_scalar_path::Union{String, Nothing}       # path netCDF file (scalar data)
    extra_dim::Union{NamedTuple, Nothing}        # name and values for extra dimension (to store SVectors)
    reducer::Dict{Union{CSVColumn, NetCDFScalarVariable}, Function} # The reducer associated with output variables
end

function NCReader(config)
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
        glob_dir = normpath(tomldir, config.dir_input)
        glob_path = replace(path_forcing, '\\' => '/')
    end
    @info "Forcing parameters are provided by `$abspath_forcing`."

    dynamic_paths = glob(glob_path, glob_dir)  # expand "data/forcing-year-*.nc"
    if isempty(dynamic_paths)
        error("No files found with name '$glob_path' in '$glob_dir'")
    end
    dataset = NCDataset(dynamic_paths; aggdim = "time", deferopen = false)

    if haskey(dataset["time"].attrib, "_FillValue")
        @warn "Time dimension contains `_FillValue` attribute, this is not in line with CF conventions."
        nctimes = dataset["time"][:]
        times_dropped = collect(skipmissing(nctimes))
        # check if lenght has changed (missings in time dimension are not allowed), and throw
        # an error if the lenghts are different
        if length(times_dropped) != length(nctimes)
            error("Time dimension in `$abspath_forcing` contains missing values")
        else
            nctimes = times_dropped
            nctimes_type = eltype(nctimes)
        end
    else
        nctimes = dataset["time"][:]
        nctimes_type = eltype(nctimes)
    end

    for (par, var) in config.input.forcing
        ncname = variable_name(var)
        variable_info(var)
        @info "Set `$par` using netCDF variable `$ncname` as forcing parameter."
    end

    # create map from internal location to netCDF variable name for cyclic parameters and
    # store cyclic times for each internal location (duplicate cyclic times are possible
    # this way, however it seems not worth to keep track of unique cyclic times for now
    # (memory usage))
    if do_cyclic(config)
        cyclic_dataset = NCDataset(cyclic_path)
        cyclic_times = Dict{String, Vector{Tuple{Int, Int}}}()
        for (par, var) in config.input.cyclic
            ncname = variable_name(var)
            i = findfirst(x -> startswith(x, "time"), dimnames(cyclic_dataset[ncname]))
            dimname = dimnames(cyclic_dataset[ncname])[i]
            cyclic_nc_times = collect(cyclic_dataset[dimname])
            cyclic_times[par] = timecycles(cyclic_nc_times)
            variable_info(var)
            @info "Set `$par` using netCDF variable `$ncname` as cyclic parameter, with `$(length(cyclic_nc_times))` timesteps."
        end
    else
        cyclic_dataset = nothing
        cyclic_times = Dict{String, Vector{Tuple{Int, Int}}}()
    end

    return NCReader{nctimes_type}(
        dataset,
        nctimes,
        cyclic_dataset,
        cyclic_times,
        config.input.forcing,
        config.input.cyclic,
    )
end

"Get a Vector of all unique location ids from a 2D map"
function locations_map(ds, mapname, config)
    map_2d = ncread(
        ds,
        config,
        mapname;
        optional = false,
        type = Union{Int, Missing},
        allow_missing = true,
    )
    ids = unique(skipmissing(map_2d))
    return ids
end

"Get a Vector{String} of all columns names for the CSV header, except the first, time"
function csv_header(cols, dataset, config)
    out = String[]
    for col in cols
        (; header, map) = col
        if !isnothing(map)
            ids = locations_map(dataset, map, config)
            hvec = [string(header, '_', id) for id in ids]
            append!(out, hvec)
        else
            push!(out, header)
        end
    end
    return out
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
function flat!(d, path, el::AbstractDict)
    for (k, v) in pairs(el)
        flat!(d, string(path, '.', k), v)
    end
    return nothing
end

function flat!(d, path, el)
    k = path
    d[k] = el
    return nothing
end

"""
    ncnames(dict)

Create a flat mapping from internal parameter locations to netCDF variable names.

Ignores top level values in the Dict. This function is used to convert a TOML such as:

```toml
[output]
path = "path/to/file.nc"

[output.land]
canopystorage = "my_canopystorage"

[output.routing.river_flow]
q = "my_q"
```

To a dictionary of the flattened parameter locations and netCDF names. The top level
values are ignored since the output path is not a netCDF name.

```julia
Dict(
    (:land, :canopystorage) => "my_canopystorage,
    (:routing, :river_flow, :q) => "my_q,
)
```
"""
function ncnames(dict)
    ncnames_dict = Dict{String, String}()
    for (k, v) in dict
        flat!(ncnames_dict, k, v)
    end
    return ncnames_dict
end

"""
    out_map(ncnames_dict, modelmap)

Create a Dict that maps parameter netCDF names to arrays in the Model.
"""
function out_map(ncnames_dict, modelmap)
    output_map = Dict{String, Any}()
    (; land) = modelmap
    for (par, ncname) in ncnames_dict
        A = if haskey(standard_name_map(land), par)
            lens = get_lens(par, land)
            lens(modelmap)
        else
            param(modelmap, par)
        end
        output_map[ncname] = (par = par, vector = A)
    end
    return output_map
end

function get_reducer_func(col, domain, args...)
    (; parameter) = col
    if occursin("reservoir", parameter)
        reducer_func = reducer(col, domain.reservoir.network.reverse_indices, args...)
    elseif occursin("river", parameter)
        reducer_func = reducer(col, domain.river.network.reverse_indices, args...)
    elseif occursin("drain", parameter)
        reducer_func = reducer(col, domain.drain.network.reverse_indices, args...)
    else
        reducer_func = reducer(col, domain.land.network.reverse_indices, args...)
    end
end

function Writer(config, modelmap, domain, nc_static; extra_dim = nothing)
    x_coords = read_x_axis(nc_static)
    y_coords = read_y_axis(nc_static)

    reducer = Dict{Union{CSVColumn, NetCDFScalarVariable}, Function}()

    # create an output netCDF that will hold all timesteps of selected parameters for grid
    # data
    if do_netcdf_grid(config)
        nc_path = output_path(config, config.output.netcdf_grid.path)
        deflatelevel = config.output.netcdf_grid.compressionlevel
        @info "Create an output netCDF file `$nc_path` for grid data, using compression level `$deflatelevel`."
        # create a flat mapping from internal parameter locations to netCDF variable names
        output_ncnames = ncnames(config.output.netcdf_grid.variables)
        # fill the output_map by mapping parameter netCDF names to arrays
        output_map = out_map(output_ncnames, modelmap)
        ds = setup_grid_netcdf(
            nc_path,
            x_coords,
            y_coords,
            output_map,
            config.time.calendar,
            config.time.time_units,
            extra_dim,
            config.model.cell_length_in_meter__flag;
            deflatelevel,
        )
    else
        nc_path = nothing
        output_map = Dict{String, Any}()
        ds = nothing
    end

    # create a separate state output netCDF that will hold the last timestep of all states
    # but only if config.state.path_output has been set
    if !isnothing(config.state.path_output)
        state_ncnames = check_states(config)
        state_map = out_map(state_ncnames, modelmap)
        nc_state_path = output_path(config, config.state.path_output)
        @info "Create a state output netCDF file `$nc_state_path`."
        ds_outstate = setup_grid_netcdf(
            nc_state_path,
            x_coords,
            y_coords,
            state_map,
            config.time.calendar,
            config.time.time_units,
            extra_dim,
            config.model.cell_length_in_meter__flag;
            float_type = Float64,
        )
    else
        ds_outstate = nothing
        state_map = Dict{String, Any}()
        nc_state_path = nothing
    end

    # create an output netCDF that will hold all timesteps of selected parameters for scalar
    # data, but only if config.netcdf.path and config.netcdf.variable have been set.
    if do_netcdf_scalar(config)
        nc_scalar_path = output_path(config, config.output.netcdf_scalar.path)
        @info "Create an output netCDF file `$nc_scalar_path` for scalar data."
        # get netCDF info for scalar data (variable name, locationset (dim) and
        # location ids)
        ds_scalar = setup_scalar_netcdf(
            nc_scalar_path,
            nc_static,
            modelmap,
            config.time.calendar,
            config.time.time_units,
            extra_dim,
            config,
        )

        for var in config.output.netcdf_scalar.variable
            reducer[var] =
                get_reducer_func(var, domain, x_coords, y_coords, config, nc_static)
        end
    else
        ds_scalar = nothing
        nc_scalar_path = nothing
    end

    if do_csv(config)
        # open CSV file and write header
        csv_path = output_path(config, config.output.csv.path)
        @info "Create an output CSV file `$csv_path` for scalar data."
        # create directory if needed
        mkpath(dirname(csv_path))
        csv_io = open(csv_path, "w")
        print(csv_io, "time,")
        header = csv_header(config.output.csv.column, nc_static, config)
        println(csv_io, join(header, ','))
        flush(csv_io)

        for col in config.output.csv.column
            reducer[col] =
                get_reducer_func(col, domain, x_coords, y_coords, config, nc_static)
        end
    else
        # no CSV file is checked by isnothing(csv_path)
        csv_path = nothing
        csv_io = IOBuffer()
    end

    return Writer(
        ds,
        output_map,
        nc_path,
        csv_path,
        config.output.csv.column,
        csv_io,
        ds_outstate,
        state_map,
        nc_state_path,
        ds_scalar,
        config.output.netcdf_scalar.variable,
        nc_scalar_path,
        extra_dim,
        reducer,
    )
end

"Write a new timestep with scalar data to a netCDF file"
function write_netcdf_timestep(model, dataset)
    (; writer, land, clock) = model

    time_index = add_time(dataset, clock.time)
    for var in writer.nc_scalar
        (; name, parameter) = var
        reducer = writer.reducer[var]
        A = if haskey(standard_name_map(land), parameter)
            lens = get_lens(parameter, land)
            lens(model)
        else
            param(model, parameter)
        end
        elemtype = eltype(A)
        # could be a value, or a vector in case of map
        if elemtype <: AbstractFloat
            v = reducer(A)
            dataset[name][:, time_index] .= v
        elseif elemtype <: SVector
            # check if an extra dimension and index is specified in the TOML file
            if haskey(nc, writer.extra_dim.name)
                i = only(nc.layer)
                v = nt.reducer(getindex.(A, i))
                dataset[name][:, time_index] .= v
            else
                nlayer = length(first(A))
                for i in 1:nlayer
                    v = nt.reducer(getindex.(A, i))
                    dataset[name][:, i, time_index] .= v
                end
            end
        else
            error("Unsupported output type: ", elemtype)
        end
    end
    return model
end

"Write a new timestep with grid data to a netCDF file"
function write_netcdf_timestep(model, dataset, parameters)
    (; clock, domain) = model

    time_index = add_time(dataset, clock.time)

    buffer = zeros(Union{Float64, Missing}, domain.land.network.modelsize)
    for (key, val) in parameters
        (; par, vector) = val
        sel = active_indices(domain, par)
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
            for i in 1:nlayer
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

# don't do anything for no dataset, used if no output netCDF is needed
write_netcdf_timestep(model, dataset::Nothing, parameters) = model
write_netcdf_timestep(model, dataset::Nothing) = model

"Write model output"
function write_output(model)
    (; writer) = model
    (; dataset, dataset_scalar, parameters) = writer

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

It can generate such a series from either TimeTypes given that the year is constant, or
it will interpret integers as either months or days of year if possible.
"""
function timecycles(times)
    if eltype(times) <: TimeType
        # all timestamps are from the same year
        year1 = year(first(times))
        if !all(==(year1), year.(times))
            error("unsupported cyclic timeseries")
        end
        # sub-daily time steps are not allowed
        min_tstep = Second(minimum(diff(times)))
        if min_tstep < Second(Day(1))
            error("unsupported cyclic timeseries")
        else
            # returns a (month, day) tuple for each date
            return monthday.(times)
        end
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
    (; reader, writer, config) = model

    close(reader.dataset)
    if do_cyclic(config)
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

const function_map = Dict{ReducerType.T, Function}(
    ReducerType.maximum => maximum,
    ReducerType.minimum => minimum,
    ReducerType.mean => mean,
    ReducerType.median => median,
    ReducerType.first => first,
    ReducerType.last => last,
    ReducerType.only => only,
    ReducerType.sum => sum,
)

reducer_func(::Nothing) = only
reducer_func(reducer_type::ReducerType.T) = function_map[reducer_type]

"Get a reducer function based on output settings for scalar data defined in a dictionary"
function reducer(col, rev_inds, x_nc, y_nc, config, dataset)
    (; parameter, map, reducer, index, coordinate) = col
    fileformat = col isa CSVColumn ? "CSV" : "NetCDF"
    f = reducer_func(reducer)
    if !isnothing(map)
        # assumes the parameter in "map" has a 2D input map, with
        # integers indicating the points or zones that are to be aggregated
        # if no reducer is given, pick "only", this is the only safe reducer,
        # and makes sense in the case of a gauge map
        map_2d = ncread(
            dataset,
            config,
            map;
            type = Union{Int, Missing},
            allow_missing = true,
            logging = false,
        )
        @info "Adding scalar output for a map with a reducer function." fileformat param =
            parameter mapname = map reducer_name = String(nameof(f))
        ids = unique(skipmissing(map_2d))
        # from id to list of internal indices
        inds = Dict{Int, Vector{Int}}(id => Vector{Int}() for id in ids)
        for i in eachindex(map_2d)
            v = map_2d[i]
            ismissing(v) && continue
            v::Int
            vector = inds[v]
            ind = rev_inds[i]
            if iszero(ind)
                error("""inactive cell found in requested scalar output
                    map `$map` value $v for parameter $param""")
            end
            push!(vector, ind)
        end
        return A -> (f(A[inds[id]]) for id in ids)
    elseif !isnothing(reducer)
        # reduce over all active cells
        # needs to be behind the map if statement, because it also can use a reducer
        @info "Adding scalar output of all active cells with reducer function." fileformat param =
            parameter reducer_name = String(nameof(f))
        return function_map[reducer]
    elseif do_index(index)
        if !isnothing(index.i)
            # linear index into the internal vector of active cells
            # this one mostly makes sense for debugging, or for vectors of only a few elements
            @info "Adding scalar output for linear index." fileformat param = parameter index
            return x -> getindex(x, index.i)
        else
            # index into the 2D input/output arrays
            # the first always corresponds to the x dimension, then the y dimension
            # this is 1-based
            ind = rev_inds[index.x, index.y]
            @info "Adding scalar output for linear index." fileformat param = parameter index
            iszero(ind) && error("inactive loc specified for output")
            return A -> getindex(A, ind)
        end
    elseif !isnothing(coordinate)
        (; x, y) = coordinate
        # find the closest cell center index
        _, iy = findmin(abs.(y_nc .- y))
        _, ix = findmin(abs.(x_nc .- x))
        I = CartesianIndex(ix, iy)
        i = rev_inds[I]
        @info "Adding scalar output for coordinate." fileformat param = parameter x y
        iszero(i) && error("inactive coordinate specified for output")
        return A -> getindex(A, i)
    else
        error("unknown reducer")
    end
end

function write_csv_row(model)
    (; writer, land, clock) = model
    (; csv_path, csv_io, csv_cols) = writer
    isnothing(csv_path) && return nothing
    print(csv_io, string(clock.time))
    for col in csv_cols
        (; parameter) = col
        reducer = writer.reducer[col]
        A = if haskey(standard_name_map(land), parameter)
            lens = get_lens(parameter, land)
            lens(model)
        else
            param(model, parameter)
        end
        # v could be a value, or a vector in case of map
        if eltype(A) <: SVector
            # indexing is required in case of a SVector and CSV output
            i = only(col.layer)
            v = reducer(getindex.(A, i))
        else
            v = reducer(A)
        end
        # numbers are also iterable
        for el in v
            print(csv_io, ',', el)
        end
    end
    return println(csv_io)
end

"From a time and a calendar, create the right CFTime DateTimeX type"
function cftime(time, calendar)
    timetype = CFTime.timetype(convert(String, calendar))
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
    clock.dt = new_clock.dt
    return nothing
end

function advance!(clock)
    clock.iteration += 1
    clock.time += clock.dt
    return nothing
end

function rewind!(clock)
    clock.iteration -= 1
    clock.time -= clock.dt
    return nothing
end

"Read a rating curve from CSV into a NamedTuple of vectors"
function read_sh_csv(path)
    data, header = readdlm(path, ',', Float64; header = true)
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
    data = readdlm(path, ',', Float64; skipstart = 1)
    # Q is a matrix with 365 columns, one for each day in the year
    return (H = data[:, 1], Q = data[:, 2:366])
end

# these represent the type of the rating curve and specific storage data
const SH = NamedTuple{(:H, :S), Tuple{Vector{Float64}, Vector{Float64}}}
const HQ = NamedTuple{(:H, :Q), Tuple{Vector{Float64}, Matrix{Float64}}}

is_increasing(v) = last(v) > first(v)

"""
    nc_dim_name(dims::Vector{Symbol}, name::Symbol)
    nc_dim_name(ds::CFDataset, name::Symbol)

Given a netCDF dataset or list of dimensions, and an internal dimension name, return the
corresponding netCDF dimension name. Certain common alternatives are supported, e.g. :lon or
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
:time), `extra_dim` depends on the model type, which will map to the correct netCDF name
using `nc_dim_name`.
"""
nc_dim(ds::CFDataset, name) = ds[nc_dim_name(ds, name)]

"""
    internal_dim_name(name::Symbol)

Given a netCDF dimension name string, return the corresponding internal dimension name.
"""
function internal_dim_name(name::Symbol)
    if name in (:x, :lon, :longitude)
        return :x
    elseif name in (:y, :lat, :latitude)
        return :y
    elseif name in (:time, :layer, :flood_depth)
        return name
    elseif startswith(string(name), "time")
        return :time
    else
        error("Unknown dimension $name")
    end
end

"""
    read_dims(A::CFVariable_MF, dim_sel)

Return the data of a netCDF data variable as an Array. Only dimensions in `dim_sel`, a
NamedTuple like (x=:, y=:, time=1). Other dimensions that may be present need to be size 1,
otherwise an error is thrown.

`dim_sel` keys should be the internal dimension names; :x, :y, :time, :`extra_dim.name`.
"""
function read_dims(A::CFVariable_MF, dim_sel::NamedTuple)
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
    pairs = Pair{Symbol, Bool}[]
    for d in dim_names
        if d == :layer && !(haskey(ds, "layer"))
            inc = true
        elseif d == :flood_depth && !(haskey(ds, "flood_depth"))
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
    elseif :flood_depth in dim_names
        desired_order = (:x, :y, :flood_depth)
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
    elseif length(dims_increasing) == 3 && haskey(dims_increasing, :flood_depth)
        dims_increasing_ordered = (
            x = dims_increasing.x,
            y = dims_increasing.y,
            flood_depth = dims_increasing.flood_depth,
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

Read the dimensions listed in dim_names from a variable with name `varname` from a netCDF
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
        if haskey(ds.dim, candidate)
            return sort!(Float64.(ds[candidate][:]))
        end
    end
    return error("no x axis found in $(path(ds))")
end

"""
    read_y_axis(ds::CFDataset)

Return the y coordinate Vector{Float64}, whether it is called y, lat or latitude.
Also sorts the vector to be increasing, to match `read_standardized`.
"""
function read_y_axis(ds::CFDataset)::Vector{Float64}
    candidates = ("y", "lat", "latitude")
    for candidate in candidates
        if haskey(ds.dim, candidate)
            return sort!(Float64.(ds[candidate][:]))
        end
    end
    return error("no y axis found in $(path(ds))")
end
