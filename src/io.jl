
"Extract a NetCDF variable at a given time"
function get_at!(
    buffer,
    var::NCDatasets.CFVariable,
    times::AbstractVector{<:TimeType},
    t::TimeType,
)
    dim = findfirst(==("time"), NCDatasets.dimnames(var))
    i = findfirst(>=(t), times)
    i === nothing && throw(DomainError("time $t after dataset end $(last(times))"))
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
    @unpack dataset, buffer, inds = reader
    nctimes = nomissing(dataset["time"][:])

    # TODO allow configurable variable names
    precipitation = get_at!(buffer, dataset["P"], nctimes, clock.time)
    vertical.precipitation .= buffer[inds]
    temperature = get_at!(buffer, dataset["TEMP"], nctimes, clock.time)
    vertical.temperature .= buffer[inds]
    potevap = get_at!(buffer, dataset["PET"], nctimes, clock.time)
    vertical.potevap .= buffer[inds]
    return model
end

"prepare an output dataset"
function setup_netcdf(output_path, nclon, nclat, parameters, calendar, time_units)
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
    defVar(
        ds,
        "time",
        Float64,
        ("time",),
        attrib = [
            "units" => time_units,
            "calendar" => calendar,
        ],
    )
    for parameter in parameters
        defVar(
            ds,
            parameter,
            Float32,
            ("lon", "lat", "time"),
            attrib = ["_FillValue" => Float32(NaN)],
        )
    end
    return ds
end

function grow_netcdf!(ds, var::AbstractString, time, A::AbstractArray)
    # index in the time dimension we want to add
    i = length(ds["time"]) + 1
    ds["time"][i] = time
    ds[var][:, :, i] = A
    return ds
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
    buffer::Matrix{T}
    inds::Vector{CartesianIndex{2}}
end

struct NCWriter
    dataset::NCDataset
    parameters::Vector{String}
end

function prepare_reader(path, varname, inds)
    dataset = NCDataset(path)
    var = dataset[varname].var

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
    return NCReader(dataset, buffer, inds)
end

function prepare_writer(config, reader, output_path)
    # TODO remove random string from the filename
    # this makes it easier to develop for now, since we don't run into issues with open files
    base, ext = splitext(output_path)
    randomized_path = string(base, '_', randstring('a':'z', 4), ext)

    nclon = Float64.(nomissing(reader.dataset["lon"][:]))
    nclat = Float64.(nomissing(reader.dataset["lat"][:]))

    output_parameters = config.output.parameters
    calendar = get(config.input, "calendar", "proleptic_gregorian")
    time_units = get(config.input, "time_units", CFTime.DEFAULT_TIME_UNITS)
    ds = Wflow.setup_netcdf(randomized_path, nclon, nclat, output_parameters, calendar, time_units)
    return NCWriter(ds, output_parameters)
end
