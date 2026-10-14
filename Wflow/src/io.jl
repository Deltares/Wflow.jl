"NetCDF forcing files with only the currently selected file open."
mutable struct RollingDataset
    paths::Vector{String}
    file_end_indices::Vector{Int}
    dataset::NC.Dataset
    file_index::Int
end

function RollingDataset(paths::Vector{String})::Tuple{RollingDataset, Vector{DateTime}}
    times_per_file = Vector{DateTime}[]
    for path in paths
        dataset = NC.Dataset(path)
        push!(times_per_file, NC.read_time(dataset, "time"))
        close(dataset)
    end
    # glob sorts paths lexicographically, which is not necessarily chronologically
    order = sortperm(times_per_file; by = first)
    paths = paths[order]
    times_per_file = times_per_file[order]

    file_end_indices = cumsum(length.(times_per_file))
    times = reduce(vcat, times_per_file)
    dataset = NC.Dataset(first(paths))
    return RollingDataset(paths, file_end_indices, dataset, 1), times
end

"Select the forcing dataset and local index for a global time index."
function dataset_index!(dataset::RollingDataset, index::Int)::Tuple{NC.Dataset, Int}
    file_index = searchsortedfirst(dataset.file_end_indices, index)
    checkbounds(dataset.paths, file_index)
    if file_index != dataset.file_index
        close(dataset.dataset)
        dataset.file_index = 0
        dataset.dataset = NC.Dataset(dataset.paths[file_index])
        dataset.file_index = file_index
    end
    previous_end =
        file_index == firstindex(dataset.file_end_indices) ? 0 :
        dataset.file_end_indices[file_index - 1]
    return dataset.dataset, index - previous_end
end

Base.close(dataset::RollingDataset) = close(dataset.dataset)

"""Turn "a.aa.aaa" into (:a, :aa, :aaa)"""
symbols(s::AbstractString) = Tuple(Symbol(x) for x in split(s, '.'))

"Get a nested field using a tuple of Symbols"
param(obj, fields::Tuple{Vararg{Symbol}}) = foldl(getproperty, fields; init = obj)
param(obj, fields::AbstractString) = param(obj, symbols(fields))

"Extract a netCDF variable at a given time"
function get_at(
        ds::RollingDataset,
        var::InputEntry,
        unit::Unit,
        times::Vector{DateTime},
        t::DateTime,
        dt::Float64,
    )::Matrix{Union{Missing, Float64}}
    # this behaves like a backward fill interpolation
    i = findfirst(>=(t), times)
    if t < first(times)
        msg = "time $(format_datetime(t)) before dataset begin $(format_datetime(first(times)))"
        throw(DomainError(msg))
    end
    if i === nothing
        msg = "time $(format_datetime(t)) after dataset end $(format_datetime(last(times)))"
        throw(DomainError(msg))
    end
    return get_at(ds, var, unit, i, dt)
end

function get_at(
        ds::RollingDataset,
        var::InputEntry,
        unit::Unit,
        i::Int,
        dt::Float64,
    )::Matrix{Union{Missing, Float64}}
    dataset, local_index = dataset_index!(ds, i)
    return get_at(dataset, var, unit, local_index, dt)
end

function get_at(
        ds::NC.Dataset,
        var::InputEntry,
        unit::Unit,
        i::Int,
        dt::Float64,
    )::Matrix{Union{Missing, Float64}}
    data = read_standardized(ds, something(variable_name(var)), (x = :, y = :, time = i))
    apply_affine_transform!(data, var)
    # dynamic input is always floating point, only a unit conversion is needed
    return to_SI!(data, unit; dt_val = dt)
end

function get_param_res(model)
    (; boundary_conditions) = model.routing.river_flow.boundary_conditions.reservoir
    return Dict{String, Vector{Float64}}(
        "atmosphere_water__precipitation_volume_flux" =>
            boundary_conditions.precipitation,
        "land_surface_water__potential_evaporation_volume_flux" =>
            boundary_conditions.evaporation,
    )
end

const mover_params = (
    "atmosphere_water__precipitation_volume_flux",
    "land_surface_water__potential_evaporation_volume_flux",
)

# Map domain name to regular expression for matching with standard name or model path
const domain_parameter_map = Dict{String, Regex}(
    "reservoir" => r"reservoir|routing.river_flow.boundary_conditions.reservoir",
    "river" => r"river|floodplain|routing.river_flow",
    "drain" => r"land_drain|routing.subsurface_flow.boundaries.drain",
)

function routing_with_reservoirs(model)
    (; config) = model
    return config.model.reservoir__flag
end

function routing_with_reservoirs(model::AbstractModel{<:SedimentModel})
    return false
end

"Set a model field from input data, for fields that are a vector of numbers."
function set_input!(field::AbstractVector{<:Real}, data)::Nothing
    field .= data
    return nothing
end

set_input!(field, data) =
    error("Setting input is only supported for vectors of numbers, not $(typeof(field)).")

"""
    set_dynamic_input!(model, par::String, data)

Set the model field corresponding to standard name `par` from forcing or cyclic `data`.
"""
function set_dynamic_input!(model, par::String, data)::Nothing
    with_model_field(model, par) do field, metadata
        if !metadata.allow_dynamic_input
            error(
                "Tried to set '$par' dynamically via cyclic/forcing input, which is not allowed.",
            )
        end
        set_input!(field, data)
    end
    return nothing
end

"Zero the input over reservoir coverage cells, as it is added to the reservoirs instead."
function zero_reservoir_coverage!(
        field::AbstractVector{Float64},
        reverse_indices::Matrix{Int},
        sel_reservoirs::Vector{Vector{CartesianIndex{2}}},
    )::Nothing
    for sel_reservoir in sel_reservoirs
        field[reverse_indices[sel_reservoir]] .= 0
    end
    return nothing
end

zero_reservoir_coverage!(field, reverse_indices, sel_reservoirs) =
    error("Unsupported field type $(typeof(field)) for reservoir coverage.")

"""
    load_fixed_forcing!(model)
Get fixed netCDF forcing input."
"""
function load_fixed_forcing!(model)
    (; reader, domain, clock) = model
    (; forcing_parameters) = reader
    dt = tosecond(clock.dt)

    do_reservoirs = routing_with_reservoirs(model)

    reverse_indices = domain.land.network.reverse_indices

    for (par, ncvar) in forcing_parameters
        if variable_name(ncvar) === nothing
            value = something(ncvar.value)
            val = only(value) * only(ncvar.scale) + only(ncvar.offset)

            with_model_field(model, par) do param, metadata
                val_SI = to_SI(val, metadata.unit; dt_val = dt)
                set_input!(param, val_SI)
                # set fixed precipitation and evaporation over the reservoirs and put these into
                # the reservoir structs and set the precipitation and evaporation to 0 in the
                # land model
                if par in mover_params && do_reservoirs
                    sel_reservoirs = domain.reservoir.network.indices_coverage
                    zero_reservoir_coverage!(param, reverse_indices, sel_reservoirs)
                    get_param_res(model)[par] .= val_SI
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
    (; clock, reader, domain, land) = model
    (; dataset, dataset_times, forcing_parameters) = reader

    dt = tosecond(clock.dt)

    do_reservoirs = routing_with_reservoirs(model)

    # load from netCDF into the model according to the mapping
    for (par, ncvar) in forcing_parameters
        # no need to update fixed values
        variable_name(ncvar) === nothing && continue

        time = clock.time
        metadata = get_metadata(par, typeof(land))
        metadata === nothing && error("No metadata found for forcing parameter '$par'.")
        if !metadata.allow_dynamic_input
            error("Tried to update '$par' forcing from input, which is not allowed.")
        end
        data = get_at(dataset, ncvar, metadata.unit, dataset_times, time, dt)

        # calculate the mean precipitation and evaporation over reservoirs and put these
        # into the reservoir structs and set the precipitation and evaporation to 0 in the
        # land model
        if par in mover_params
            if do_reservoirs
                sel_reservoirs = domain.reservoir.network.indices_coverage
                param_res = get_param_res(model)
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
            msg = "Forcing data at $(format_datetime(time)) has missing values on active model cells for $(variable_name(ncvar))"
            throw(ArgumentError(msg))
        end
        set_dynamic_input!(model, par, data_sel)
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
    (; clock, reader, domain, land) = model
    (; cyclic_times, cyclic_parameters) = reader
    cyclic_dataset = reader.cyclic_dataset::NC.Dataset

    dt = tosecond(clock.dt)

    # pick up the data that is valid for the past model time step
    month_day = monthday(clock.time - clock.dt)

    is_first_timestep = clock.iteration == 1

    for (par, ncvar) in cyclic_parameters
        metadata = get_metadata(par, typeof(land))
        metadata === nothing && error("No metadata found for cyclic parameter '$par'.")
        if !metadata.allow_dynamic_input
            error("Tried to update '$par' from cyclic input, which is not allowed.")
        end
        if is_first_timestep || (month_day in cyclic_times[par])
            # time for an update of the cyclic forcing
            i = findlast(t -> monthday_passed(month_day, t), cyclic_times[par])
            isnothing(i) &&
                error("Could not find applicable cyclic timestep for $month_day")
            # load from netCDF into the model according to the mapping
            data = get_at(cyclic_dataset, ncvar, metadata.unit, i, dt)
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
            set_dynamic_input!(model, par, data_sel)
        end
    end
    return nothing
end

"""
    NC_HANDLES::Dict{String, NC.Dataset}

For each netCDF file that will be opened for writing, store an entry in this Dict from the
absolute path of the file to the dataset. This allows us to close the dataset if we try
to create them twice in the same session.
"""
const NC_HANDLES = Dict{String, NC.Dataset}()

"NetCDF global attributes following the CF conventions."
const CF_GLOBAL_ATTRIB = OrderedDict{String, String}(
    "Conventions" => "CF-1.12",
    "references" => "https://deltares.github.io/Wflow.jl/",
    "wflow_version" => string(VERSION),
)

"Safely create a netCDF file, even if it has already been opened for creation"
function create_tracked_netcdf(path::AbstractString)::NC.Dataset
    abs_path = abspath(path)
    # close existing dataset if it exists
    if haskey(NC_HANDLES, abs_path)
        # fine if it was already closed
        close(NC_HANDLES[abs_path])
    end
    # create directory if needed
    mkpath(dirname(abs_path))
    ds = NC.create(abs_path)
    for (name, value) in CF_GLOBAL_ATTRIB
        NC.put_global_attribute(ds, name, value)
    end
    NC_HANDLES[abs_path] = ds
    return ds
end

"""
Reduction of model output to scalar values, for scalar netCDF or CSV output. This is a
concrete type rather than a closure, to avoid dynamic dispatch.

- `map`: reduce over cells grouped by zone ids (`ids`, with cell indices `inds`)
- `all`: reduce over all active cells
- `index`: select a single cell `index`, `0` for an inactive cell (output NaN)
"""
@enumx ReducerKind map all index

struct Reducer
    kind::ReducerKind.T
    reducer_type::ReducerType.T
    index::Int
    ids::Vector{Int}
    inds::Vector{Vector{Int}}
end

Reducer(kind::ReducerKind.T, reducer_type::ReducerType.T) =
    Reducer(kind, reducer_type, 0, Int[], Vector{Int}[])
Reducer(index::Int) =
    Reducer(ReducerKind.index, ReducerType.only, index, Int[], Vector{Int}[])

"Name and coordinate values of the extra dimension used to store SVectors, e.g. layers."
struct ExtraDim
    name::String
    value::Vector{Float64}
end

ExtraDim() = ExtraDim("", Float64[])
has_extra_dim(extra_dim::ExtraDim) = !isempty(extra_dim.name)

"Model output variable, identified by standard name `par`."
@kwdef struct OutputData
    par::String
    unit::Unit
    # whether the model field is a vector of SVectors, written with an extra dimension
    layered::Bool = false
end

@kwdef struct NCWriter
    # Path to the NetCDF file
    output_path::Union{String, Nothing} = nothing
    # NetCDF dataset
    output_dataset::Union{NC.Dataset, Nothing} = nothing
    # mapping of netCDF variable names to model parameters
    output_map::OrderedDict{String, OutputData} = OrderedDict{String, OutputData}()
    # The reducer associated with the output variables (scalar output only)
    reducer::Vector{Tuple{NetCDFScalarVariable, Reducer}} =
        Tuple{NetCDFScalarVariable, Reducer}[]
    # CF time units of the time dimension
    time_units::String = CFTime.DEFAULT_TIME_UNITS
end

@kwdef struct CSVWriter
    # Path to the CSV file
    output_path::Union{String, Nothing} = nothing
    # File handle to CSV file
    output_io::IOStream = devnull_stream()
    # Mapping of CSV variable names to model parameters
    output_map::Dict{String, OutputData} = Dict{String, OutputData}()
    # The reducer associated with the output variables
    reducer::Vector{Tuple{CSVColumn, Reducer}} = Tuple{CSVColumn, Reducer}[]
end

"An IOStream placeholder for when no CSV output is requested."
devnull_stream() = IOStream("<no output>")

@kwdef struct Writer
    # Writer for transient grid output (no reducer)
    grid_writer::NCWriter
    # Writer for transient scalar output (with reducer)
    scalar_writer::NCWriter
    # Writer for for simulation end state output (no reducer)
    endstate_writer::NCWriter
    # Writer for CSV output
    csv_writer::CSVWriter
    # Name and values for extra dimension (to store SVectors)
    extra_dim::ExtraDim
end

"Define a CF time dimension and variable."
function def_time_netcdf(ds::NC.Dataset, calendar, time_units::String; deflatelevel = 0)
    NC.def_dim(ds, "time", NC.NC_UNLIMITED)
    id = NC.def_var(ds, "time", Float64, ["time"]; deflatelevel)
    NC.put_attribute(ds, id, "units", time_units)
    NC.put_attribute(ds, id, "calendar", convert(String, calendar))
    NC.put_attribute(ds, id, "standard_name", "time")
    NC.put_attribute(ds, id, "axis", "T")
    return nothing
end

"Define a 1D coordinate variable with values and attributes."
function def_coordinate_netcdf(
        ds::NC.Dataset,
        name::String,
        values::Vector{Float64},
        attributes::Vector{Pair{String, String}},
    )::Nothing
    NC.def_dim(ds, name, length(values))
    id = NC.def_var(ds, name, Float64, [name])
    for (key, value) in attributes
        NC.put_attribute(ds, id, key, value)
    end
    NC.write_hyperslab(ds, name, values, [1], [length(values)])
    return nothing
end

"Define a data variable of type `float_type` with a NaN fill value and units."
function def_data_netcdf(
        ds::NC.Dataset,
        name::String,
        float_type::Type{F},
        dims::Vector{String},
        unit::Unit;
        deflatelevel::Int = 0,
    )::Nothing where {F <: AbstractFloat}
    id = NC.def_var(ds, name, float_type, dims; deflatelevel)
    NC.put_attribute(ds, id, "_FillValue", float_type(NaN))
    NC.put_attribute(ds, id, "units", to_string(unit; BMI_standard = true))
    return nothing
end

"prepare an output dataset for scalar data"
function setup_scalar_netcdf(
        path,
        dataset,
        output_map::OrderedDict{String, OutputData},
        calendar,
        time_units,
        extra_dim::ExtraDim,
        config,
        indices;
        float_type::Type{F} = Float32,
    ) where {F <: AbstractFloat}
    ds = create_tracked_netcdf(path)
    def_time_netcdf(ds, calendar, time_units)
    set_extradim_netcdf(ds, extra_dim)
    for scalar_variable in config.output.netcdf_scalar.variable
        (; map, _location_dim, location, name) = scalar_variable
        (; unit, layered) = output_map[name]
        # Delft-FEWS requires the attribute :cf_role = "timeseries_id" when a netCDF file
        # contains more than one location list
        if !NC.hasdim(ds, _location_dim)
            locations = if isnothing(map)
                [something(location)]
            else
                string.(locations_map(dataset, map, config, indices))
            end
            NC.def_dim(ds, _location_dim, length(locations))
            id = NC.def_var(ds, _location_dim, String, [_location_dim])
            NC.put_attribute(ds, id, "cf_role", "timeseries_id")
            NC.write_strings(ds, _location_dim, locations)
        end
        dims = if layered && isnothing(scalar_variable.layer)
            [_location_dim, extra_dim.name, "time"]
        else
            [_location_dim, "time"]
        end
        def_data_netcdf(ds, name, float_type, dims, unit)
    end
    return ds
end

"set extra dimension in output netCDF file"
function set_extradim_netcdf(ds::NC.Dataset, extra_dim::ExtraDim)::Nothing
    has_extra_dim(extra_dim) || return nothing
    # the axis attribute `Z` is required to import this type of 3D data by Delft-FEWS the
    # values of this dimension `extra_dim.value` should be of type Float64
    attributes = if extra_dim.name == "layer"
        ["long_name" => "layer_index", "standard_name" => "layer_index", "axis" => "Z"]
    else
        Pair{String, String}[]
    end
    def_coordinate_netcdf(ds, extra_dim.name, extra_dim.value, attributes)
    return nothing
end

"prepare an output dataset for grid data"
function setup_grid_netcdf(
        path,
        ncx::Vector{Float64},
        ncy::Vector{Float64},
        output_map::OrderedDict{String, OutputData},
        calendar,
        time_units::String,
        extra_dim::ExtraDim,
        cell_length_in_meter::Bool,
        float_type::Type{F} = Float32;
        deflatelevel::Int = 0,
    ) where {F <: AbstractFloat}
    base_dims, attrib_x, attrib_y = if cell_length_in_meter
        ("x", "y", "time"),
            ("x coordinate of projection", "projection_x_coordinate", "m"),
            ("y coordinate of projection", "projection_y_coordinate", "m")
    else
        ("lon", "lat", "time"),
            ("longitude", "longitude", "degrees_east"),
            ("latitude", "latitude", "degrees_north")
    end
    ds = create_tracked_netcdf(path)
    def_coordinate_netcdf(
        ds,
        base_dims[1],
        ncx,
        [
            "long_name" => attrib_x[1],
            "standard_name" => attrib_x[2],
            "axis" => "X",
            "units" => attrib_x[3],
        ],
    )
    def_coordinate_netcdf(
        ds,
        base_dims[2],
        ncy,
        [
            "long_name" => attrib_y[1],
            "standard_name" => attrib_y[2],
            "axis" => "Y",
            "units" => attrib_y[3],
        ],
    )
    set_extradim_netcdf(ds, extra_dim)
    def_time_netcdf(ds, calendar, time_units; deflatelevel)
    for (name, output_data) in output_map
        (; unit, layered) = output_data
        dims = if layered
            # for SVectors an additional dimension (`extra_dim`) is required
            [base_dims[1], base_dims[2], extra_dim.name, base_dims[3]]
        else
            [base_dims[1], base_dims[2], base_dims[3]]
        end
        def_data_netcdf(ds, name, float_type, dims, unit; deflatelevel)
    end
    return ds
end

"Add a new time to the unlimited time dimension, and return the index"
function add_time(ds::NC.Dataset, time::DateTime, time_units::String)::Int
    i = NC.dimlength(ds, "time") + 1
    NC.write_hyperslab(ds, "time", [NC.encode_time(time, time_units)], [i], [1])
    return i
end

struct NCReader
    dataset::RollingDataset
    dataset_times::Vector{DateTime}
    cyclic_dataset::Union{NC.Dataset, Nothing}
    cyclic_times::Dict{String, Vector{Tuple{Int, Int}}}
    forcing_parameters::InputEntries
    cyclic_parameters::InputEntries
end

"Match a file name against a glob pattern supporting `*` and `?`."
function glob_match(pattern::AbstractString, name::AbstractString)::Bool
    p = collect(pattern)
    s = collect(name)
    # dynamic programming over pattern and name positions
    matches = falses(length(p) + 1, length(s) + 1)
    matches[1, 1] = true
    for i in 1:length(p)
        if p[i] == '*'
            matches[i + 1, 1] = matches[i, 1]
        end
    end
    for i in 1:length(p), j in 1:length(s)
        matches[i + 1, j + 1] = if p[i] == '*'
            matches[i, j + 1] || matches[i + 1, j]
        elseif p[i] == '?' || p[i] == s[j]
            matches[i, j]
        else
            false
        end
    end
    return matches[end, end]
end

"""
    expand_glob(glob_path, glob_dir)::Vector{String}

Expand a glob pattern relative to `glob_dir`. Only the last path component may contain
wildcards (`*` and `?`). Results are sorted lexicographically.
"""
function expand_glob(glob_path::AbstractString, glob_dir::AbstractString)::Vector{String}
    dir_part, file_pattern = splitdir(glob_path)
    if any(c -> c in ('*', '?', '['), dir_part)
        error("Wildcards are only supported in the file name of `$glob_path`.")
    end
    dir = normpath(glob_dir, dir_part)
    if !any(c -> c in ('*', '?'), file_pattern)
        path = joinpath(dir, file_pattern)
        return isfile(path) ? [path] : String[]
    end
    isdir(dir) || return String[]
    return sort!([joinpath(dir, name) for name in readdir(dir) if glob_match(file_pattern, name)])
end

function NCReader(config)
    path_forcing = config.input.path_forcing
    abspath_forcing = input_path(config, path_forcing)
    cyclic_path = input_path(config, config.input.path_static)

    @info "Cyclic parameters are provided by `$cyclic_path`."

    # the path separator in a glob pattern is always /
    if isabspath(path_forcing)
        glob_dir = ""
        glob_path = path_forcing
    else
        tomldir = dirname(config)
        glob_dir = normpath(tomldir, config.dir_input)
        glob_path = replace(path_forcing, '\\' => '/')
    end
    @info "Forcing parameters are provided by `$abspath_forcing`."

    dynamic_paths = expand_glob(glob_path, glob_dir)  # expand "data/forcing-year-*.nc"
    if isempty(dynamic_paths)
        error("No files found with name '$glob_path' in '$glob_dir'")
    end
    dataset, nctimes = RollingDataset(dynamic_paths)

    land_type = config.model.type == ModelType.sediment ? SoilLossModel : LandHydrologySBM
    for (par, var) in config.input.forcing
        ncname = variable_name(var)
        variable_info(var)
        metadata = get_metadata(par, land_type)
        unit_str = isnothing(metadata) ? "<unit not found>" : string(metadata.unit)
        @info "Set `$par [$unit_str]` using netCDF variable `$ncname` as forcing parameter."
    end

    # create map from internal location to netCDF variable name for cyclic parameters and
    # store cyclic times for each internal location (duplicate cyclic times are possible
    # this way, however it seems not worth to keep track of unique cyclic times for now
    # (memory usage))
    cyclic_times = Dict{String, Vector{Tuple{Int, Int}}}()
    if do_cyclic(config)
        cyclic_dataset = NC.Dataset(cyclic_path)
        for (par, var) in config.input.cyclic
            ncname = something(variable_name(var))
            dim_names = NC.var_dimnames(cyclic_dataset, ncname)
            i = findfirst(x -> startswith(x, "time"), dim_names)
            isnothing(i) && error("No time dimension found for cyclic parameter `$ncname`.")
            dimname = dim_names[i]
            cyclic_times[par] = if NC.is_time_variable(cyclic_dataset, dimname)
                timecycles(NC.read_time(cyclic_dataset, dimname))
            else
                timecycles(
                    round.(Int, NC.read_vector_nomissing(cyclic_dataset, dimname)),
                )
            end
            variable_info(var)
            metadata = get_metadata(par, land_type)
            unit_str = isnothing(metadata) ? "<unit not found>" : string(metadata.unit)
            @info "Set `$par [$unit_str]` using netCDF variable `$ncname` as cyclic parameter, with `$(length(cyclic_times[par]))` timesteps."
        end
    else
        cyclic_dataset = nothing
    end

    return NCReader(
        dataset,
        nctimes,
        cyclic_dataset,
        cyclic_times,
        config.input.forcing,
        config.input.cyclic,
    )
end

"Get a Vector of all unique location ids from a 2D map"
function locations_map(ds, mapname, config, indices)::Vector{Int}
    map_2d = ncread(
        ds,
        config,
        mapname,
        Writer;
        sel = indices,
        metadata = ParameterMetadata(; type = Int, allow_missing = true),
    )::Vector{Union{Missing, Int}}
    ids = unique(skipmissing(map_2d))
    return ids
end

"Get a Vector{String} of all columns names for the CSV header, except the first, time"
function csv_header(cols, dataset, config, indices)
    out = String[]
    for col in cols
        (; header, map) = col
        if !isnothing(map)
            ids = locations_map(dataset, map, config, indices)
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
    "a.b.field_of_b" => "name_in_output,
    "a.d.c.field_of_c" => "other_name_in_output,
)
"""
function flat!(d::Dict{String, String}, path::String, @nospecialize(el))::Nothing
    # use an explicit stack rather than recursion, since the values are of unknown type
    stack = Tuple{String, Any}[(path, el)]
    while !isempty(stack)
        key, value = pop!(stack)
        if value isa Dict{String, Any}
            for (k, v) in value
                push!(stack, (string(key, '.', k), v))
            end
        elseif value isa String
            d[key] = value
        else
            error("Expected a string value for `$key` in the TOML.")
        end
    end
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
"""
function ncnames(dict)::Dict{String, String}
    ncnames_dict = Dict{String, String}()
    for (k, v) in dict
        flat!(ncnames_dict, String(k), v)
    end
    return ncnames_dict
end

"Whether a model field is layered (a vector of SVectors), and should be written as such."
is_layered_field(::AbstractVector{<:SVector}) = true
is_layered_field(::Any) = false

"""
    out_map(output_names, modelmap)

Create a Dict that maps output names to the standard names and units in the Model.
"""
function out_map(output_names, modelmap)::OrderedDict{String, OutputData}
    output_map = OrderedDict{String, OutputData}()
    for (par, output_name) in output_names
        unit, layered = with_model_field(modelmap, par) do field, metadata
            metadata.unit, is_layered_field(field)
        end
        output_map[output_name] = OutputData(; par, unit, layered)
    end
    return output_map
end

function get_reducer_func(col, domain, x_nc, y_nc, config, dataset)::Reducer
    (; parameter) = col
    (; indices) = domain.land.network
    rev_inds = if startswith(parameter, domain_parameter_map["reservoir"])
        domain.reservoir.network.reverse_indices
    elseif startswith(parameter, domain_parameter_map["river"])
        domain.river.network.reverse_indices
    elseif startswith(parameter, domain_parameter_map["drain"])
        domain.drain.network.reverse_indices
    else
        domain.land.network.reverse_indices
    end
    return reducer(col, rev_inds, indices, x_nc, y_nc, config, dataset)
end

function Writer(
        config::Config,
        modelmap::NamedTuple,
        domain,
        nc_static;
        extra_dim::ExtraDim = ExtraDim(),
    )
    x_coords = read_x_axis(nc_static)
    y_coords = read_y_axis(nc_static)
    time_units = config.time.time_units

    # create an output netCDF that will hold all timesteps of selected parameters for grid
    # data
    grid_writer = if do_netcdf_grid(config)
        output_path_grid = output_path(config, config.output.netcdf_grid.path)
        # create a flat mapping from internal parameter locations to netCDF variable names
        output_ncnames = ncnames(config.output.netcdf_grid.variables)
        # fill the output_map by mapping parameter netCDF names to arrays
        output_map = out_map(output_ncnames, modelmap)
        deflatelevel = config.output.netcdf_grid.compressionlevel
        output_dataset = setup_grid_netcdf(
            output_path_grid,
            x_coords,
            y_coords,
            output_map,
            config.time.calendar,
            time_units,
            extra_dim,
            config.model.cell_length_in_meter__flag;
            deflatelevel,
        )
        @info "Created an output netCDF file `$output_path_grid` for grid data, using compression level `$deflatelevel`."
        NCWriter(; output_path = output_path_grid, output_dataset, output_map, time_units)
    else
        NCWriter()
    end

    # create a separate state output netCDF that will hold the last timestep of all states
    # but only if config.state.path_output has been set
    endstate_writer = if !isnothing(config.state.path_output)
        output_path_endstate = output_path(config, config.state.path_output)
        output_ncnames = check_states(config)
        output_map = out_map(output_ncnames, modelmap)
        output_dataset = setup_grid_netcdf(
            output_path_endstate,
            x_coords,
            y_coords,
            output_map,
            config.time.calendar,
            time_units,
            extra_dim,
            config.model.cell_length_in_meter__flag,
            Float64,
        )
        @info "Created a state output netCDF file `$output_path_endstate`."
        NCWriter(;
            output_path = output_path_endstate,
            output_dataset,
            output_map,
            time_units,
        )
    else
        NCWriter()
    end

    # create an output netCDF that will hold all timesteps of selected parameters for scalar
    # data, but only if config.netcdf.path and config.netcdf.variable have been set.
    scalar_writer = if do_netcdf_scalar(config)
        output_path_scalar = output_path(config, config.output.netcdf_scalar.path)
        output_ncnames =
            [var.parameter => var.name for var in config.output.netcdf_scalar.variable]
        output_map = out_map(output_ncnames, modelmap)
        # get netCDF info for scalar data (variable name, locationset (dim) and
        # location ids)
        indices = domain.land.network.indices
        output_dataset = setup_scalar_netcdf(
            output_path_scalar,
            nc_static,
            output_map,
            config.time.calendar,
            time_units,
            extra_dim,
            config,
            indices,
        )
        reducer = Tuple{NetCDFScalarVariable, Reducer}[
            (var, get_reducer_func(var, domain, x_coords, y_coords, config, nc_static)) for var in config.output.netcdf_scalar.variable
        ]
        @info "Created an output netCDF file `$output_path_scalar` for scalar data."
        NCWriter(;
            output_path = output_path_scalar,
            output_dataset,
            output_map,
            reducer,
            time_units,
        )
    else
        NCWriter()
    end

    csv_writer = if do_csv(config)
        # open CSV file and write header
        output_path_csv = output_path(config, config.output.csv.path)
        # create directory if needed
        mkpath(dirname(output_path_csv))
        output_io = open(output_path_csv, "w")
        # Add header
        print(output_io, "time,")
        indices = domain.land.network.indices
        header = csv_header(config.output.csv.column, nc_static, config, indices)
        println(output_io, join(header, ','))
        flush(output_io)
        output_csvnames = [col.parameter => col.header for col in config.output.csv.column]
        output_map = Dict{String, OutputData}(out_map(output_csvnames, modelmap))
        reducer = Tuple{CSVColumn, Reducer}[
            (col, get_reducer_func(col, domain, x_coords, y_coords, config, nc_static)) for col in config.output.csv.column
        ]
        @info "Created an output CSV file `$output_path_csv` for scalar data."
        CSVWriter(; output_path = output_path_csv, output_io, output_map, reducer)
    else
        CSVWriter()
    end

    return Writer(; grid_writer, scalar_writer, endstate_writer, csv_writer, extra_dim)
end

"Collect the values of a model field as a Vector{Float64}, selecting a layer for SVectors."
field_values(field::AbstractVector{<:AbstractFloat}, layer)::Vector{Float64} =
    Vector{Float64}(field)
field_values(field::AbstractVector{<:SVector}, layer::Int)::Vector{Float64} =
    Float64[x[layer] for x in field]
field_values(field::AbstractVector{<:SVector}, ::Nothing)::Vector{Float64} =
    error("A layer index is required for output of layered variables.")
field_values(field::Real, layer)::Vector{Float64} = Float64[field]
field_values(field, layer)::Vector{Float64} =
    error("Unsupported output type: $(typeof(field))")

"Write a new timestep with scalar data to a netCDF file"
function write_netcdf_scalar_timestep(model::AbstractModel, writer::NCWriter)::Nothing
    (; clock) = model
    (; output_dataset, output_map, reducer, time_units) = writer
    isnothing(output_dataset) && return nothing
    (; extra_dim) = model.writer
    dt_val = tosecond(clock.dt)
    time_index = add_time(output_dataset, clock.time, time_units)
    for (var, reducer_) in reducer
        (; name, layer) = var
        (; par, unit, layered) = output_map[name]
        with_model_field(model, par) do field, _
            if layered && isnothing(layer)
                nlayer = length(extra_dim.value)
                for i in 1:nlayer
                    v = from_SI(apply_reducer(reducer_, field_values(field, i)), unit; dt_val)
                    NC.write_hyperslab(
                        output_dataset,
                        name,
                        v,
                        [1, i, time_index],
                        [length(v), 1, 1],
                    )
                end
            else
                v = from_SI(apply_reducer(reducer_, field_values(field, layer)), unit; dt_val)
                NC.write_hyperslab(output_dataset, name, v, [1, time_index], [length(v), 1])
            end
        end
    end
    return nothing
end

"Write a 2D grid buffer from model field values at the active cells `sel`."
function fill_grid_buffer!(
        buffer::Matrix{Float64},
        sel::Vector{CartesianIndex{2}},
        values::Vector{Float64},
        unit::Unit,
        dt_val::Float64,
    )::Nothing
    # ensure no other information is written
    fill!(buffer, NaN)
    # cut off possible boundary conditions/ ghost points with [1:length(sel)]
    for (i, index) in enumerate(sel)
        buffer[index] = values[i]
    end
    from_SI!(buffer, unit; dt_val)
    return nothing
end

"Write a new timestep with grid data to a netCDF file"
function write_netcdf_timestep(model::AbstractModel, writer::NCWriter)::Nothing
    (; output_dataset, output_map, time_units) = writer
    isnothing(output_dataset) && return nothing
    (; clock, domain) = model

    time_index = add_time(output_dataset, clock.time, time_units)
    dt = tosecond(clock.dt)

    modelsize = domain.land.network.modelsize
    buffer = zeros(Float64, modelsize[1], modelsize[2])
    nx, ny = size(buffer)
    for (key, val) in output_map
        (; par, unit, layered) = val
        sel = active_indices(domain, par)
        with_model_field(model, par) do field, _
            if layered
                nlayer = length(model.writer.extra_dim.value)
                for i in 1:nlayer
                    fill_grid_buffer!(buffer, sel, field_values(field, i), unit, dt)
                    NC.write_hyperslab(
                        output_dataset,
                        key,
                        buffer,
                        [1, 1, i, time_index],
                        [nx, ny, 1, 1],
                    )
                end
            else
                fill_grid_buffer!(buffer, sel, field_values(field, nothing), unit, dt)
                NC.write_hyperslab(
                    output_dataset,
                    key,
                    buffer,
                    [1, 1, time_index],
                    [nx, ny, 1],
                )
            end
        end
    end
    return nothing
end

"Write model output"
function write_output(model)
    (; writer) = model
    (; grid_writer, scalar_writer, csv_writer) = writer

    write_csv_row(model, csv_writer)
    write_netcdf_timestep(model, grid_writer)
    write_netcdf_scalar_timestep(model, scalar_writer)

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
function timecycles(times::Vector{DateTime})::Vector{Tuple{Int, Int}}
    # all timestamps are from the same year
    year1 = year(first(times))
    if !all(==(year1), year.(times))
        error("unsupported cyclic timeseries")
    end
    # sub-daily time steps are not allowed
    min_tstep = minimum(diff(times))
    if min_tstep < Dates.Millisecond(Day(1))
        error("unsupported cyclic timeseries")
    end
    # returns a (month, day) tuple for each date
    return monthday.(times)
end

function timecycles(times::Vector{Int})::Vector{Tuple{Int, Int}}
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
end

"Close input and output datasets that are opened on model initialization"
function close_files(model; delete_output::Bool = false)
    (; reader, writer, config) = model
    (; grid_writer, scalar_writer, endstate_writer, csv_writer) = writer

    # Input files
    close(reader.dataset)
    if !isnothing(reader.cyclic_dataset)
        close(reader.cyclic_dataset)
    end

    # Output NetCDF files
    for writer_ in (grid_writer, scalar_writer, endstate_writer)
        (; output_dataset, output_path) = writer_
        isnothing(output_dataset) && continue

        close(output_dataset)
        delete_output && rm(something(output_path))
    end

    # Output CSV file
    if !isnothing(csv_writer.output_path)
        close(csv_writer.output_io)
        delete_output && rm(csv_writer.output_path)
    end
    return nothing
end

"""
Median of a vector, NaN if any value is NaN. `Statistics.median` is not used since its
error path is not compatible with `juliac --trim`.
"""
function median_value(A::AbstractVector{Float64})::Float64
    isempty(A) && throw(ArgumentError("median of an empty array is undefined"))
    any(isnan, A) && return NaN
    sorted = sort(A)
    n = length(sorted)
    mid = n ÷ 2
    return isodd(n) ? sorted[mid + 1] : (sorted[mid] + sorted[mid + 1]) / 2
end

"Apply a reduction function to a vector."
function reduce_values(reducer_type::ReducerType.T, A::AbstractVector{Float64})::Float64
    return if reducer_type == ReducerType.maximum
        maximum(A)
    elseif reducer_type == ReducerType.minimum
        minimum(A)
    elseif reducer_type == ReducerType.mean
        mean(A)
    elseif reducer_type == ReducerType.median
        median_value(A)
    elseif reducer_type == ReducerType.first
        first(A)
    elseif reducer_type == ReducerType.last
        last(A)
    elseif reducer_type == ReducerType.only
        only(A)
    else
        sum(A)
    end
end

"Apply a reducer to model output values, returning one value per location."
function apply_reducer(reducer::Reducer, A::Vector{Float64})::Vector{Float64}
    (; kind, reducer_type, index, ids, inds) = reducer
    return if kind == ReducerKind.map
        Float64[reduce_values(reducer_type, A[inds[i]]) for i in eachindex(ids)]
    elseif kind == ReducerKind.all
        Float64[reduce_values(reducer_type, A)]
    else
        Float64[iszero(index) ? NaN : A[index]]
    end
end

"Get a reducer based on output settings for scalar data defined in a dictionary"
function reducer(col, rev_inds, indices, x_nc, y_nc, config, dataset)::Reducer
    (; parameter, map, index, coordinate) = col
    reducer_type = col.reducer
    fileformat = col isa CSVColumn ? "CSV" : "NetCDF"
    if !isnothing(map)
        # assumes the parameter in "map" has a 2D input map, with
        # integers indicating the points or zones that are to be aggregated
        # if no reducer is given, pick "only", this is the only safe reducer,
        # and makes sense in the case of a gauge map
        f = something(reducer_type, ReducerType.only)
        map_1d = ncread(
            dataset,
            config,
            map,
            Writer;
            sel = indices,
            logging = false,
            metadata = ParameterMetadata(; allow_missing = true, type = Int),
        )::Vector{Union{Missing, Int}}
        @info "Adding scalar output for a map with a reducer function." fileformat param =
            parameter mapname = map reducer_name = string(f)
        ids = unique(skipmissing(map_1d))
        inds = Dict{Int, Vector{Int}}(id => Vector{Int}() for id in ids)

        for i in eachindex(map_1d)
            v = map_1d[i]
            ismissing(v) && continue
            # translate from land-domain position -> target-domain position
            ind = rev_inds[indices[i]]
            if iszero(ind)
                error(
                    """inactive cell found in requested scalar output
                    map `$map` value $v for parameter $parameter"""
                )
            end
            push!(inds[v], ind)
        end
        return Reducer(ReducerKind.map, f, 0, ids, [inds[id] for id in ids])
    elseif !isnothing(reducer_type)
        # reduce over all active cells
        # needs to be behind the map if statement, because it also can use a reducer
        @info "Adding scalar output of all active cells with reducer function." fileformat param =
            parameter reducer_name = string(reducer_type)
        return Reducer(ReducerKind.all, reducer_type)
    elseif do_index(index)
        if !isnothing(index.i)
            # linear index into the internal vector of active cells
            # this one mostly makes sense for debugging, or for vectors of only a few elements
            @info "Adding scalar output for linear index." fileformat param = parameter index =
                index.i
            return Reducer(index.i)
        else
            # index into the 2D input/output arrays
            # the first always corresponds to the x dimension, then the y dimension
            # this is 1-based
            ind = rev_inds[something(index.x), something(index.y)]
            if iszero(ind)
                @warn "Inactive index specified for output, skipping" fileformat param =
                    parameter x = index.x y = index.y
            else
                @info "Adding scalar output for linear index." fileformat param = parameter x =
                    index.x y = index.y
            end
            return Reducer(ind)
        end
    elseif !isnothing(coordinate)
        (; x, y) = coordinate
        # find the closest cell center index
        _, iy = findmin(abs.(y_nc .- y))
        _, ix = findmin(abs.(x_nc .- x))
        I = CartesianIndex(ix, iy)
        i = rev_inds[I]
        if iszero(i)
            @warn "Inactive coordinate specified for output, skipping" fileformat param =
                parameter x y
        else
            @info "Adding scalar output for coordinate." fileformat param = parameter x y
        end
        return Reducer(i)
    else
        error("unknown reducer")
    end
end

function write_csv_row(model::AbstractModel, csv_writer::CSVWriter)::Nothing
    (; clock) = model
    (; output_path, output_io, output_map, reducer) = csv_writer
    isnothing(output_path) && return nothing
    print(output_io, format_datetime(clock.time))
    dt_val = tosecond(clock.dt)
    for (col, reducer_) in reducer
        (; header, layer) = col
        (; par, unit) = output_map[header]
        v = with_model_field(model, par) do field, _
            from_SI(apply_reducer(reducer_, field_values(field, layer)), unit; dt_val)
        end
        for el in v
            print(output_io, ',', el)
        end
    end
    println(output_io)
    return nothing
end

"Zero pad a non-negative integer to `width` digits."
function zero_pad(x::Integer, width::Int)::String
    s = string(x)
    while length(s) < width
        s = "0" * s
    end
    return s
end

"""
Format a DateTime as ISO 8601, like `string(::DateTime)`. The Dates formatting code is not
compatible with `juliac --trim`.
"""
function format_datetime(t::DateTime)::String
    s = string(
        zero_pad(year(t), 4),
        '-',
        zero_pad(Dates.month(t), 2),
        '-',
        zero_pad(Dates.day(t), 2),
        'T',
        zero_pad(Dates.hour(t), 2),
        ':',
        zero_pad(Dates.minute(t), 2),
        ':',
        zero_pad(Dates.second(t), 2),
    )
    ms = Dates.millisecond(t)
    return iszero(ms) ? s : string(s, '.', zero_pad(ms, 3))
end

"""
From a time and a calendar, create the model time. Only standard calendars are supported,
represented by `DateTime`.
"""
function cftime(time, calendar::CalendarType.T)::DateTime
    NC.check_calendar(convert(String, calendar))
    return DateTime(time)
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

"""
    read_csv_matrix(path)::Tuple{Matrix{Float64}, Vector{String}}

Read a comma separated file with a header line and numeric values. Returns the values as a
matrix (rows are lines in the file) and the column names. This replaces
`DelimitedFiles.readdlm`, which is not compatible with `juliac --trim`.
"""
function read_csv_matrix(path::AbstractString)::Tuple{Matrix{Float64}, Vector{String}}
    # `readlines` is not used, `EachLine` stores a `Function` field
    lines = filter!(!isempty, strip.(split(read(path, String), '\n')))
    isempty(lines) && error("No data found in CSV file `$path`.")
    header = String[strip(name) for name in split(first(lines), ',')]
    ncol = length(header)
    data = Matrix{Float64}(undef, length(lines) - 1, ncol)
    for (row, line) in enumerate(Iterators.drop(lines, 1))
        values = split(line, ',')
        length(values) == ncol ||
            error("Expected $ncol values on line $(row + 1) of CSV file `$path`.")
        for (col, value) in enumerate(values)
            data[row, col] = parse(Float64, value)
        end
    end
    return data, header
end

"Read a rating curve from CSV into a NamedTuple of vectors"
function read_sh_csv(path)
    data, header = read_csv_matrix(path)
    names = uppercase.(header)
    idx_h = findfirst(==("H"), names)
    idx_s = findfirst(==("S"), names)

    if isnothing(idx_h) || isnothing(idx_s)
        error("$path needs to provide H and S columns, got $(join(names, ", "))")
    end

    return (H = data[:, idx_h], S = data[:, idx_s])
end

"Read a specific storage curve from CSV into a NamedTuple of vectors"
function read_hq_csv(path)
    data, _ = read_csv_matrix(path)
    # Q is a matrix with 365 columns, one for each day in the year
    return (H = data[:, 1], Q = data[:, 2:366])
end

# these represent the type of the rating curve and specific storage data
const SH = NamedTuple{(:H, :S), Tuple{Vector{Float64}, Vector{Float64}}}
const HQ = NamedTuple{(:H, :Q), Tuple{Vector{Float64}, Matrix{Float64}}}

is_increasing(v) = last(v) > first(v)

"""
    nc_dim_name(dims::Vector{Symbol}, name::Symbol)
    nc_dim_name(ds::NC.Dataset, name::Symbol)

Given a netCDF dataset or list of dimensions, and an internal dimension name, return the
corresponding netCDF dimension name. Certain common alternatives are supported, e.g. :lon or
:longitude instead of :x.
"""
function nc_dim_name(dims::Vector{Symbol}, name::Symbol)::Symbol
    # direct naming
    name in dims && return name
    # list of common alternative names
    return if name == :x
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

nc_dim_name(ds::NC.Dataset, name::Symbol)::Symbol =
    nc_dim_name(Symbol.(NC.dimnames(ds)), name)

"""
    nc_dim(ds::NC.Dataset, name::Symbol)

Return the dimension coordinate, based on the internal name (:x, :y, :`extra_dim.name`,
:time), `extra_dim` depends on the model type, which will map to the correct netCDF name
using `nc_dim_name`.
"""
nc_dim(ds::NC.Dataset, name::Symbol)::Vector{Float64} =
    NC.read_vector_nomissing(ds, String(nc_dim_name(ds, name)))

"""
    internal_dim_name(name::Symbol)

Given a netCDF dimension name string, return the corresponding internal dimension name.
"""
function internal_dim_name(name::Symbol)::Symbol
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

"Number of dimensions selected with a `Colon` in `dim_sel`, as a `Val`."
@generated function selected_ndims(::NamedTuple{K, T}) where {K, T}
    return :(Val($(count(t -> t === Colon, T.parameters))))
end

"""
    read_dims(ds::NC.Dataset, varname, dim_sel)

Return the data of a netCDF data variable as an Array. Only dimensions in `dim_sel`, a
NamedTuple like (x=:, y=:, time=1). Other dimensions that may be present need to be size 1,
otherwise an error is thrown.

`dim_sel` keys should be the internal dimension names; :x, :y, :time, :`extra_dim.name`.
The number of dimensions of the returned array is the number of `Colon` entries in
`dim_sel`.
"""
function read_dims(
        ds::NC.Dataset,
        varname::AbstractString,
        dim_sel::NamedTuple,
    )::Tuple{Array{Union{Missing, Float64}}, Vector{Symbol}}
    return read_dims(ds, varname, dim_sel, selected_ndims(dim_sel))
end

function read_dims(
        ds::NC.Dataset,
        varname::AbstractString,
        dim_sel::NamedTuple,
        ::Val{N},
    )::Tuple{Array{Union{Missing, Float64}, N}, Vector{Symbol}} where {N}
    dim_names = NC.var_dimnames(ds, varname)
    dim_sizes = NC.var_size(ds, varname)
    start = Int[]
    count = Int[]
    data_dim_order = Symbol[]
    data_size = Int[]
    # need to iterate in this order, to get the indices in order
    for (dim_name, dim_size) in zip(dim_names, dim_sizes)
        if dim_size == 1
            push!(start, 1)
            push!(count, 1)
        else
            dim = internal_dim_name(Symbol(dim_name))
            # each dim should either be in dim_sel or be size 1
            if dim in keys(dim_sel)
                idx = dim_sel[dim]
                if idx isa Int
                    push!(start, idx)
                    push!(count, 1)
                else
                    push!(start, 1)
                    push!(count, dim_size)
                    push!(data_dim_order, dim)
                    push!(data_size, dim_size)
                end
            else
                throw(
                    ArgumentError(
                        """NetCDF dimension $dim_name has length $dim_size.
                        Only extra dimensions of length 1 are supported."""
                    )
                )
            end
        end
    end
    flat = NC.read_hyperslab(ds, varname, start, count)
    if length(data_size) != N
        error(
            "NetCDF variable `$varname` has $(length(data_size)) dimensions of length > 1, expected $N.",
        )
    end
    data = reshape(flat, ntuple(i -> data_size[i], Val(N)))
    return data, data_dim_order
end

"""
    dim_directions(ds, dim_names)

Given a netCDF dataset and a list of internal dimension names, return a Vector{Bool}, which
is `true` if the dimension is increasing, or `false` otherwise.

For the layer dimension, we allow coordinate arrays to be missing, in which case we
consider it increasing, going from the top layer (1) to deeper layers. This is to keep
accepting data that we have accepted before.
"""
function dim_directions(ds::NC.Dataset, dim_names::Vector{Symbol})::Vector{Bool}
    increasing = Bool[]
    for d in dim_names
        if d == :layer && !NC.hasvar(ds, "layer")
            inc = true
        elseif d == :flood_depth && !NC.hasvar(ds, "flood_depth")
            inc = true
        else
            inc = is_increasing(nc_dim(ds, d))
        end
        push!(increasing, inc)
    end
    return increasing
end

"""
    permute_data(data, dim_names)

Given an Array of data, and a list of its dimension names, return a permuted version of the
data such that the dimension order will follow (:x, :y, `extra_dim.name`). No permutation is
done if this is not needed.
"""
function permute_data(
        data::Array{T, N},
        dim_names::Vector{Symbol},
    )::Tuple{Array{T, N}, Vector{Symbol}} where {T, N}
    @assert N == length(dim_names)
    if :layer in dim_names
        desired_order = [:x, :y, :layer]
    elseif :flood_depth in dim_names
        desired_order = [:x, :y, :flood_depth]
    else
        desired_order = [:x, :y]
    end
    @assert dim_names ⊆ desired_order
    permutation = Int[findfirst(==(d), dim_names) for d in desired_order]
    if permutation == 1:N
        return data, dim_names
    else
        return permutedims(data, permutation), dim_names[permutation]
    end
end

"""
    reverse_data!(data, dims_increasing)

Reverse the data as needed, such that it is increasing along all dimensions.
"""
function reverse_data!(data::Array, dims_increasing::Vector{Bool})
    for (dim, increasing) in enumerate(dims_increasing)
        increasing || reverse!(data; dims = dim)
    end
    return data
end

"""
    read_standardized(ds::NC.Dataset, varname::AbstractString, dim_sel::NamedTuple)

Read the dimensions listed in dim_names from a variable with name `varname` from a netCDF
dataset `ds`. `dim_sel` should be a NamedTuple like (x=:, y=:, time=1), which will return
a 2 dimensional array with x and y axes, representing the first index in the time dimension.
"""
function read_standardized(ds::NC.Dataset, varname::AbstractString, dim_sel::NamedTuple)
    return read_standardized(ds, varname, dim_sel, selected_ndims(dim_sel))
end

function read_standardized(
        ds::NC.Dataset,
        varname::AbstractString,
        dim_sel::NamedTuple,
        n::Val{N},
    )::Array{Union{Missing, Float64}, N} where {N}
    data, data_dim_order = read_dims(ds, varname, dim_sel, n)
    data, new_dim_order = permute_data(data, data_dim_order)
    dims_increasing = dim_directions(ds, new_dim_order)
    return reverse_data!(data, dims_increasing)
end

"Return a sorted coordinate vector using the first candidate dimension name found."
function read_axis(ds::NC.Dataset, candidates::NTuple{3, String})::Vector{Float64}
    for candidate in candidates
        if NC.hasdim(ds, candidate)
            return sort!(NC.read_vector_nomissing(ds, candidate))
        end
    end
    return error("no $(first(candidates)) axis found in $(ds.path)")
end

"""
    read_x_axis(ds::NC.Dataset)

Return the x coordinate Vector{Float64}, whether it is called x, lon or longitude.
Also sorts the vector to be increasing, to match `read_standardized`.
"""
read_x_axis(ds::NC.Dataset)::Vector{Float64} = read_axis(ds, ("x", "lon", "longitude"))
read_x_axis(ds::RollingDataset)::Vector{Float64} = read_x_axis(ds.dataset)

"""
    read_y_axis(ds::NC.Dataset)

Return the y coordinate Vector{Float64}, whether it is called y, lat or latitude.
Also sorts the vector to be increasing, to match `read_standardized`.
"""
read_y_axis(ds::NC.Dataset)::Vector{Float64} = read_axis(ds, ("y", "lat", "latitude"))
read_y_axis(ds::RollingDataset)::Vector{Float64} = read_y_axis(ds.dataset)
