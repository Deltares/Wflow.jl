#=
Code to support input and output of data and configuration.
Input data can be loaded from NetCDF files.
Output data can be written to NetCDF or CSV files.
For configuration files we use TOML.
=#

"Parsed TOML configuration"
struct Config
    dict::Dict{String,Any}  # nested key value mapping of all settings
    path::Union{String,Nothing}  # path to the TOML file, or nothing
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
Base.keys(config::Config) = keys(Dict(config))
Base.values(config::Config) = values(Dict(config))
Base.pairs(config::Config) = pairs(Dict(config))
Base.get(config::Config, key, default) = get(Dict(config), key, default)
Base.getindex(config::Config, i) = getindex(Dict(config), i)
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
    for (key, val) in pairs(parameters)
        if eltype(val) <: AbstractFloat
            # all floats are saved as Float32
            defVar(
                ds,
                key,
                Float32,
                ("lon", "lat", "time"),
                attrib = ["_FillValue" => Float32(NaN)],
            )
        elseif eltype(val) <: SVector
            # SVectors are used to store layers
            defVar(
                ds,
                key,
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
    inds_riv::Vector{CartesianIndex{2}}
end

struct Writer
    dataset::NCDataset
    parameters::Dict{String,Any}
    csv_path::Union{String,Nothing}
    csv_cols::Vector
    csv_io::IO
    statenames::Tuple{String,Vararg{String}}
end

"Convert a piece of Config to a Dict{Symbol, String} used for parameter lookup"
parameter_lookup_table(config) = Dict(Symbol(k) => String(v) for (k, v) in Dict(config))

function prepare_reader(path, cyclic_path, inds, inds_riv, config)
    dataset = NCDataset(path)
    var = dataset[config.dynamic.parameters.precipitation].var

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

    forcing_parameters = parameter_lookup_table(config.dynamic.parameters)
    cyclic_parameters = parameter_lookup_table(config.cyclic.parameters)
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
        inds_riv,
    )
end

function prepare_writer(config, reader, output_path, modelmap, maxlayers, statenames)

    nclon = ncread(reader.dataset, "lon"; type = Float64)
    nclat = ncread(reader.dataset, "lat"; type = Float64)

    # fill the output_map by mapping parameters to arrays
    output_map = Dict{String,Vector}()
    for (param, ncname) in pairs(config.output.parameters)
        parsym = Symbol(param)
        for (name, component) in pairs(modelmap)
            if parsym in propertynames(component)
                # this array will be reused
                # therefore the reference in this dict can be used for all timesteps.
                if haskey(output_map, param)
                    @warn "output parameter $param not unique, currently set to $name"
                end
                output_map[param] = getproperty(component, parsym)
            end
        end
        if !haskey(output_map, param)
            error("output parameter $param not found in any model component")
        end
    end

    calendar = get(config, "calendar", "proleptic_gregorian")
    time_units = get(config, "time_units", CFTime.DEFAULT_TIME_UNITS)
    ds =
        setup_netcdf(output_path, nclon, nclat, output_map, calendar, time_units, maxlayers)
    tomldir = dirname(config)

    if haskey(config, "csv") && haskey(config.csv, "column")
        # open CSV file and write header
        csv_path = joinpath(tomldir, config.csv.path)
        csv_io = open(csv_path, "w")
        header = join((col["header"] for col in config.csv.column), ',')
        println(csv_io, "time,", header)
        flush(csv_io)

        # cleate a vector of (lens, reducer) named tuples which will be used to
        # retrieve and reduce the CSV data during a model run
        csv_cols = []
        for col in config.csv.column
            lens = Wflow.paramap[col["parameter"]]
            reducer = Wflow.reducer(col)
            push!(csv_cols, (; lens, reducer))
        end
    else
        # no CSV file is checked by isnothing(csv_path)
        csv_path = nothing
        csv_cols = []
        csv_io = IOBuffer()
    end


    return Writer(ds, output_map, csv_path, csv_cols, csv_io, statenames)
end

"Write model output"
function write_output(model, writer::Writer)
    @unpack vertical, clock, reader = model
    @unpack buffer, inds, inds_riv = reader
    @unpack dataset, parameters = writer

    write_csv_row(model)

    time_index = add_time(dataset, clock.time)

    for (key, vector) in pairs(parameters)
        inds_used = Symbol(key) in propertynames(model.lateral.river) ? inds_riv : inds
        # write the active cells vector to the 2d buffer matrix
        elemtype = eltype(vector)
        if elemtype <: AbstractFloat
            # ensure no other information is written
            fill!(buffer, NaN)
            buffer[inds_used] .= vector
            dataset[key][:, :, time_index] = buffer
        elseif elemtype <: SVector
            nlayer = length(first(vector))
            for i = 1:nlayer
                # ensure no other information is written
                fill!(buffer, NaN)
                buffer[inds_used] .= getindex.(vector, i)
                dataset[key][:, :, i, time_index] = buffer
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
function close_files(model; delete_output::Bool = false)
    @unpack reader, writer = model

    output_nc_path = try
        path(writer.dataset)
    catch
        nothing
    end
    # TODO patch NCDatasets to not throw on calling close on closed dataset, like Base
    # and perhaps the same for the `path` function above
    try
        close(reader.dataset)
    catch
    end
    try
        close(reader.cyclic_dataset)
    catch
    end
    try
        close(writer.dataset)
    catch
    end
    close(writer.csv_io)

    if delete_output
        isnothing(output_nc_path) || rm(output_nc_path)
        isfile(writer.csv_path) && rm(writer.csv_path)
    end
    return nothing
end

"Get a reducer function based on CSV output settings defined in a dictionary"
function reducer(col)
    if haskey(col, "reducer")
        if col["reducer"] == "maximum"
            return maximum
        elseif col["reducer"] == "minimum"
            return minimum
        elseif col["reducer"] == "mean"
            return mean
        elseif col["reducer"] == "median"
            return median
        elseif col["reducer"] == "first"
            return first
        elseif col["reducer"] == "last"
            return last
        else
            error("unknown reducer")
        end
    elseif haskey(col, "index")
        return x -> getindex(x, col["index"])
    elseif haskey(col, "id")
        # id = true  # creat n columns, "volume_101"
        error("id not implemented")
    elseif haskey(col, "loc")
        # loc = [100,50]
        error("loc not implemented")
    elseif haskey(col, "coordinate")
        # coordinate = [53.2, 5.6]
        error("coordinate not implemented")
    else
        error("unknown reducer")
    end
end

function write_csv_row(model)
    @unpack writer, clock = model
    isnothing(writer.csv_path) && return nothing
    io = writer.csv_io
    print(io, clock.time)
    for nt in writer.csv_cols
        A = get(model, nt.lens)
        v = nt.reducer(A)
        print(io, ',', v)
    end
    println(io)
end

function reset_clock!(clock::Clock, config::Config)
    clock.time = config.starttime
    clock.iteration = 1
    clock.Δt = Second(config.timestepsecs)
end
