#=
Code to support input and output of data and configuration.
Input data can be loaded from NetCDF files.
Output data can be written to NetCDF or CSV files.
For configuration files we use TOML.
=#

"Parsed TOML configuration"
struct Config
    dict::Dict{String,Any}  # nested key value mapping of all settings
    path::Union{String, Nothing}  # path to the TOML file, or nothing
end

Config(path::AbstractString) = Config(parsefile(path), path)
Config(dict::AbstractDict) = Config(dict, nothing)

# allows using getproperty, e.g. config.input.time instead of config["input"]["time"]
function Base.getproperty(config::Config, f::Symbol)
    dict = Dict(config)
    path = pathof(config)
    a = dict[String(f)]
    # if it is a Dict, wrap the result in Config to keep the getproperty behavior
    return a isa AbstractDict ? Config(a, path) : a
end

# also used in autocomplete
Base.propertynames(config::Config) = collect(keys(Dict(config)))
Base.haskey(config::Config, key) = haskey(Dict(config), key)
Base.get(config::Config, key, default) = get(Dict(config), key, default)
Base.Dict(config::Config) = getfield(config, :dict)
Base.pathof(config::Config) = getfield(config, :path)
Base.dirname(config::Config) = dirname(pathof(config))

"Extract a NetCDF variable at a given time"
function get_at!(
    buffer,
    var::NCDatasets.CFVariable,
    times::AbstractVector{<:TimeType},
    t::TimeType,
)
    # this behaves like a forward fill interpolation
    i = findfirst(>=(t), times)
    i === nothing && throw(DomainError("time $t after dataset end $(last(times))"))
    return get_at!(buffer, var, i)
end

function get_at!(buffer, var::NCDatasets.CFVariable, i)
    # assumes the dataset has 12 time steps, from January to December
    dim = findfirst(==("time"), NCDatasets.dimnames(var))
    # load in place, using a lower level NCDatasets function
    # currently all indices must be of the same type, so create three ranges
    # https://github.com/Alexander-Barth/NCDatasets.jl/blob/fa742ee1b36c9e4029a40581751a21c140f01f84/src/variable.jl#L372
    spatialdim1 = 1:size(buffer, 1)
    spatialdim2 = 1:size(buffer, 2)

    if dim == 1
        NCDatasets.load!(var.var, buffer, i:i, spatialdim1, spatialdim2)
    elseif dim == 3
        NCDatasets.load!(var.var, buffer, spatialdim1, spatialdim2, i:i)
    else
        error("Time dimension expected at position 1 or 3")
    end
    return buffer
end

"Get dynamic NetCDF input for the given time"
function update_forcing!(model)
    @unpack vertical, clock, reader = model
    @unpack dataset, forcing_parameters, buffer, inds = reader
    nctimes = ncread(dataset, "time")

    # load from NetCDF into the model according to the mapping
    for (param, ncvarname) in forcing_parameters
        buffer = get_at!(buffer, dataset[ncvarname], nctimes, clock.time)
        param_vector = getproperty(vertical, param)
        param_vector .= buffer[inds]
    end

    return model
end

"Get cyclic NetCDF input for the given time"
function update_cyclic!(model)
    @unpack vertical, clock, reader = model
    @unpack cyclic_dataset, cyclic_times, cyclic_parameters, buffer, inds = reader

    month_day = monthday(clock.time)
    if monthday(clock.time) in cyclic_times
        # time for an update of the cyclic forcing
        i = findfirst(==(month_day), cyclic_times)

        # load from NetCDF into the model according to the mapping
        for (param, ncvarname) in cyclic_parameters
            buffer = get_at!(buffer, cyclic_dataset[ncvarname], i)
            param_vector = getproperty(vertical, param)
            param_vector .= buffer[inds]
        end
    end
end

"prepare an output dataset"
function setup_netcdf(
    output_path,
    nclon,
    nclat,
    parameters,
    calendar,
    time_units,
    row,
    maxlayers,
)
    ds = NCDataset(output_path, "c")
    defDim(ds, "time", Inf)  # unlimited
    defVar(
        ds,
        "lon",
        nclon,
        ("lon",),
        attrib = [
            "_FillValue" => NaN,
            "long_name" => "longitude",
            "units" => "degrees_east",
        ],
    )
    defVar(
        ds,
        "lat",
        nclat,
        ("lat",),
        attrib = [
            "_FillValue" => NaN,
            "long_name" => "latitude",
            "units" => "degrees_north",
        ],
    )
    defVar(ds, "layer", collect(1:maxlayers), ("layer",))
    defVar(
        ds,
        "time",
        Float64,
        ("time",),
        attrib = ["units" => time_units, "calendar" => calendar],
    )
    for parameter in parameters
        cell = getproperty(row, Symbol(parameter))[1]
        if cell isa AbstractFloat
            # all floats are saved as Float32
            defVar(
                ds,
                parameter,
                Float32,
                ("lon", "lat", "time"),
                attrib = ["_FillValue" => Float32(NaN)],
            )
        elseif cell isa SVector
            # SVectors are used to store layers
            defVar(
                ds,
                parameter,
                Float32,
                ("lon", "lat", "layer", "time"),
                attrib = ["_FillValue" => Float32(NaN)],
            )
        else
            error("Unsupported output type: ", typeof(cell))
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

function checkdims(dims)
    # TODO check if the x y ordering is equal to the staticmaps NetCDF
    @assert length(dims) == 3
    @assert "time" in dims
    @assert ("x" in dims) || ("lon" in dims)
    @assert ("y" in dims) || ("lat" in dims)
    @assert dims[2] != "time"
    return dims
end

struct NCReader{T}
    dataset::NCDataset
    cyclic_dataset::NCDataset
    cyclic_times::Vector{Tuple{Int,Int}}
    forcing_parameters::Dict{Symbol,String}
    cyclic_parameters::Dict{Symbol,String}
    buffer::Matrix{T}
    inds::Vector{CartesianIndex{2}}
end

struct NCWriter
    dataset::NCDataset
    parameters::Vector{String}
end

"Convert a piece of Config to a Dict{Symbol, String} used for parameter lookup"
parameter_lookup_table(config) = Dict(Symbol(k) => String(v) for (k, v) in Dict(config))

function prepare_reader(path, cyclic_path, inds, config)
    dataset = NCDataset(path)
    var = dataset[config.forcing_parameters.precipitation].var

    fillvalue = get(var.attrib, "_FillValue", nothing)
    scale_factor = get(var.attrib, "scale_factor", nothing)
    add_offset = get(var.attrib, "add_offset", nothing)
    # TODO support scale_factor and add_offset with in place loading
    # TODO check other forcing parameters as well
    @assert isnothing(fillvalue) || isnan(fillvalue)
    @assert isnothing(scale_factor) || isone(scale_factor)
    @assert isnothing(add_offset) || iszero(add_offset)

    T = eltype(var)
    dims = dimnames(var)
    checkdims(dims)
    timelast = last(dims) == "time"
    lateral_size = timelast ? size(var)[1:2] : size(var)[2:3]
    buffer = zeros(T, lateral_size)

    # TODO:include LAI climatology in update() vertical SBM model
    # we currently assume the same dimension ordering as the forcing
    cyclic_dataset = NCDataset(cyclic_path)
    cyclic_nc_times = collect(cyclic_dataset["time"])
    cyclic_times = Wflow.timecycles(cyclic_nc_times)

    forcing_parameters = parameter_lookup_table(config.forcing_parameters)
    cyclic_parameters = parameter_lookup_table(config.cyclic_parameters)
    forcing_keys = keys(forcing_parameters)
    cyclic_keys = keys(cyclic_parameters)
    for forcing_key in forcing_keys
        if forcing_key in cyclic_keys
            error("parameter specified in both forcing and cyclic")
        end
    end

    return NCReader(
        dataset,
        cyclic_dataset,
        cyclic_times,
        forcing_parameters,
        cyclic_parameters,
        buffer,
        inds,
    )
end

function prepare_writer(config, reader, output_path, row, maxlayers)
    # TODO remove random string from the filename
    # this makes it easier to develop for now, since we don't run into issues with open files
    base, ext = splitext(output_path)
    randomized_path = string(base, '_', randstring('a':'z', 4), ext)

    nclon = ncread(reader.dataset, "lon"; type = Float64)
    nclat = ncread(reader.dataset, "lat"; type = Float64)

    output_parameters = config.output.parameters
    calendar = get(config.input, "calendar", "proleptic_gregorian")
    time_units = get(config.input, "time_units", CFTime.DEFAULT_TIME_UNITS)
    ds = Wflow.setup_netcdf(
        randomized_path,
        nclon,
        nclat,
        output_parameters,
        calendar,
        time_units,
        row,
        maxlayers,
    )
    return NCWriter(ds, output_parameters)
end

"Write NetCDF output"
function write_output(model, writer::NCWriter)
    @unpack vertical, clock, reader = model
    @unpack buffer, inds = reader
    @unpack dataset, parameters = writer

    time_index = add_time(dataset, clock.time)

    for parameter in parameters
        # write the active cells vector to the 2d buffer matrix
        param = Symbol(parameter)
        vector = getproperty(vertical, param)

        elemtype = eltype(vector)
        if elemtype <: AbstractFloat
            # ensure no other information is written
            fill!(buffer, NaN)
            buffer[inds] .= vector
            dataset[parameter][:, :, time_index] = buffer
        elseif elemtype <: SVector
            nlayer = length(first(vector))
            for i = 1:nlayer
                # ensure no other information is written
                fill!(buffer, NaN)
                buffer[inds] .= getindex.(vector, i)
                dataset[parameter][:, :, i, time_index] = buffer
            end
        else
            error("Unsupported output type: ", elemtype)
        end
    end

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
function close_files(model)
    @unpack reader, writer = model

    close(reader.dataset)
    close(reader.cyclic_dataset)
    close(writer.dataset)
end
