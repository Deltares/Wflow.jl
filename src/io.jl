
"Extract a NetCDF variable at a given time"
function get_at(var::NCDatasets.CFVariable, times::AbstractVector{<:TimeType}, t::TimeType)
    dim = findfirst(==("time"), NCDatasets.dimnames(var))
    i = findfirst(>=(t), times)
    i === nothing && throw(DomainError("time $t after dataset end $(last(times))"))
    return selectdim(var, dim, i)
end

"Get dynamic NetCDF input for the given time"
function update_forcing!(model)
    @unpack vertical, clock, reader, inds = model

    nctimes = nomissing(reader["time"][:])
    # TODO allow configurable variable names
    # TODO avoid allocations
    precipitation = get_at(reader["P"], nctimes, clock.time)
    temperature = get_at(reader["TEMP"], nctimes, clock.time)
    potevap = get_at(reader["PET"], nctimes, clock.time)

    # do a mapping from 2d to 1d
    # TODO permute depending on dim ordering like in initialize_sbm_model
    vertical.precipitation .= nomissing(permutedims(precipitation)[inds])
    vertical.temperature .= nomissing(permutedims(temperature)[inds])
    vertical.potevap .= nomissing(permutedims(potevap)[inds])
    return model
end

"prepare an output dataset"
function setup_netcdf(output_path, nclon, nclat)
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
            "units" => CFTime.DEFAULT_TIME_UNITS,
            "calendar" => "proleptic_gregorian",
        ],
    )
    defVar(
        ds,
        "q",
        Float32,
        ("lon", "lat", "time"),
        attrib = ["_FillValue" => Float32(NaN)],
    )
    return ds
end

function grow_netcdf!(ds, var::AbstractString, time, A::AbstractArray)
    # index in the time dimension we want to add
    i = length(ds["time"]) + 1
    ds["time"][i] = time
    ds[var][:, :, i] = A
    return ds
end
