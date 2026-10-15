"""
Minimal, type-stable wrapper around the netCDF-C library.

This replaces the NCDatasets.jl dependency, whose high level API returns types that depend
on the file contents (variable element type, attributes, time decoding), and is therefore
not compatible with `juliac --trim`. All reads here convert to `Float64` (the netCDF-C
library does the type conversion), missing values are represented as `missing` in an
`Array{Union{Missing, Float64}, N}`, with `N` known at compile time.

Only the subset of functionality needed by Wflow is implemented.
"""
module NC

using NetCDF_jll: libnetcdf
using Dates: Dates, DateTime, Millisecond

const NC_NOWRITE = Cint(0x0000)
const NC_WRITE = Cint(0x0001)
const NC_CLOBBER = Cint(0x0000)
const NC_NETCDF4 = Cint(0x1000)
const NC_GLOBAL = Cint(-1)
const NC_UNLIMITED = Csize_t(0)
const NC_ENOTATT = Cint(-43)
const NC_EBADDIM = Cint(-46)
const NC_ENOTVAR = Cint(-49)
const NC_MAX_NAME = 256

const NC_CHAR = Cint(2)
const NC_INT = Cint(4)
const NC_FLOAT = Cint(5)
const NC_DOUBLE = Cint(6)
const NC_INT64 = Cint(10)
const NC_STRING = Cint(12)

nc_type(::Type{Float32}) = NC_FLOAT
nc_type(::Type{Float64}) = NC_DOUBLE
nc_type(::Type{Int32}) = NC_INT
nc_type(::Type{Int64}) = NC_INT64
nc_type(::Type{String}) = NC_STRING

struct NetCDFError <: Exception
    code::Cint
    msg::String
end

Base.showerror(io::IO, err::NetCDFError) =
    print(io, "NetCDF error code ", err.code, ": ", err.msg)

function check(code::Cint, context::String = "")::Nothing
    code == 0 && return nothing
    reason = unsafe_string(ccall((:nc_strerror, libnetcdf), Cstring, (Cint,), code))
    msg = isempty(context) ? reason : string(reason, " (", context, ")")
    throw(NetCDFError(code, msg))
end

"An open netCDF file."
mutable struct Dataset
    ncid::Cint
    path::String
    isopen::Bool
end

"Open an existing netCDF file for reading (`write = false`) or appending."
function Dataset(path::AbstractString; write::Bool = false)
    ncid = Ref{Cint}(0)
    mode = write ? NC_WRITE : NC_NOWRITE
    check(
        ccall((:nc_open, libnetcdf), Cint, (Cstring, Cint, Ptr{Cint}), path, mode, ncid),
        path,
    )
    return Dataset(ncid[], String(path), true)
end

"Open a netCDF file, apply `f` to it and close it afterwards."
function Dataset(f::F, path::AbstractString) where {F}
    ds = Dataset(path)
    try
        return f(ds)
    finally
        close(ds)
    end
end

"Create a new netCDF4 file, overwriting an existing file at `path`."
function create(path::AbstractString)::Dataset
    ncid = Ref{Cint}(0)
    mode = NC_CLOBBER | NC_NETCDF4
    check(
        ccall((:nc_create, libnetcdf), Cint, (Cstring, Cint, Ptr{Cint}), path, mode, ncid),
        path,
    )
    return Dataset(ncid[], String(path), true)
end

function Base.close(ds::Dataset)::Nothing
    if ds.isopen
        check(ccall((:nc_close, libnetcdf), Cint, (Cint,), ds.ncid), ds.path)
        ds.isopen = false
    end
    return nothing
end

Base.isopen(ds::Dataset) = ds.isopen

function sync(ds::Dataset)::Nothing
    check(ccall((:nc_sync, libnetcdf), Cint, (Cint,), ds.ncid), ds.path)
    return nothing
end

###
### Inquiry
###

"Return the variable id, or -1 if the variable does not exist."
function varid(ds::Dataset, name::AbstractString)::Cint
    id = Ref{Cint}(0)
    code = ccall(
        (:nc_inq_varid, libnetcdf),
        Cint,
        (Cint, Cstring, Ptr{Cint}),
        ds.ncid,
        name,
        id,
    )
    code == NC_ENOTVAR && return Cint(-1)
    check(code, name)
    return id[]
end

"Return the dimension id, or -1 if the dimension does not exist."
function dimid(ds::Dataset, name::AbstractString)::Cint
    id = Ref{Cint}(0)
    code = ccall(
        (:nc_inq_dimid, libnetcdf),
        Cint,
        (Cint, Cstring, Ptr{Cint}),
        ds.ncid,
        name,
        id,
    )
    code == NC_EBADDIM && return Cint(-1)
    check(code, name)
    return id[]
end

hasvar(ds::Dataset, name::AbstractString)::Bool = varid(ds, name) >= 0
hasdim(ds::Dataset, name::AbstractString)::Bool = dimid(ds, name) >= 0

function checked_varid(ds::Dataset, name::AbstractString)::Cint
    id = varid(ds, name)
    id < 0 && error("Variable `$name` not found in netCDF file `$(ds.path)`.")
    return id
end

function dimlength(ds::Dataset, id::Cint)::Int
    len = Ref{Csize_t}(0)
    check(
        ccall(
            (:nc_inq_dimlen, libnetcdf),
            Cint,
            (Cint, Cint, Ptr{Csize_t}),
            ds.ncid,
            id,
            len,
        ),
    )
    return Int(len[])
end

function dimlength(ds::Dataset, name::AbstractString)::Int
    id = dimid(ds, name)
    id < 0 && error("Dimension `$name` not found in netCDF file `$(ds.path)`.")
    return dimlength(ds, id)
end

function dimname(ds::Dataset, id::Cint)::String
    buffer = zeros(UInt8, NC_MAX_NAME + 1)
    check(
        ccall(
            (:nc_inq_dimname, libnetcdf),
            Cint,
            (Cint, Cint, Ptr{UInt8}),
            ds.ncid,
            id,
            buffer,
        ),
    )
    return unsafe_string(pointer(buffer))
end

"Names of all dimensions in the file."
function dimnames(ds::Dataset)::Vector{String}
    ndims = Ref{Cint}(0)
    check(ccall((:nc_inq_ndims, libnetcdf), Cint, (Cint, Ptr{Cint}), ds.ncid, ndims))
    ids = zeros(Cint, ndims[])
    check(
        ccall(
            (:nc_inq_dimids, libnetcdf),
            Cint,
            (Cint, Ptr{Cint}, Ptr{Cint}, Cint),
            ds.ncid,
            C_NULL,
            ids,
            Cint(0),
        ),
    )
    return String[dimname(ds, id) for id in ids]
end

"Dimension ids of a variable, in Julia (column major) order."
function var_dimids(ds::Dataset, id::Cint)::Vector{Cint}
    ndims = Ref{Cint}(0)
    check(
        ccall(
            (:nc_inq_varndims, libnetcdf),
            Cint,
            (Cint, Cint, Ptr{Cint}),
            ds.ncid,
            id,
            ndims,
        ),
    )
    ids = zeros(Cint, ndims[])
    check(
        ccall(
            (:nc_inq_vardimid, libnetcdf),
            Cint,
            (Cint, Cint, Ptr{Cint}),
            ds.ncid,
            id,
            ids,
        ),
    )
    return reverse!(ids)
end

"Dimension names of a variable, in Julia (column major) order."
var_dimnames(ds::Dataset, name::AbstractString)::Vector{String} =
    String[dimname(ds, id) for id in var_dimids(ds, checked_varid(ds, name))]

"Size of a variable, in Julia (column major) order."
var_size(ds::Dataset, name::AbstractString)::Vector{Int} =
    Int[dimlength(ds, id) for id in var_dimids(ds, checked_varid(ds, name))]

###
### Attributes
###

"Return the attribute length, or -1 if the attribute does not exist."
function attlength(ds::Dataset, id::Cint, name::AbstractString)::Int
    len = Ref{Csize_t}(0)
    code = ccall(
        (:nc_inq_attlen, libnetcdf),
        Cint,
        (Cint, Cint, Cstring, Ptr{Csize_t}),
        ds.ncid,
        id,
        name,
        len,
    )
    code == NC_ENOTATT && return -1
    check(code, name)
    return Int(len[])
end

function atttype(ds::Dataset, id::Cint, name::AbstractString)::Cint
    type = Ref{Cint}(0)
    check(
        ccall(
            (:nc_inq_atttype, libnetcdf),
            Cint,
            (Cint, Cint, Cstring, Ptr{Cint}),
            ds.ncid,
            id,
            name,
            type,
        ),
        name,
    )
    return type[]
end

"Read a numeric attribute as Float64, returns `nothing` if it does not exist."
function attribute_float(ds::Dataset, id::Cint, name::AbstractString)::Union{Nothing, Float64}
    len = attlength(ds, id, name)
    len < 0 && return nothing
    len == 1 || error("Attribute `$name` should be scalar in `$(ds.path)`.")
    value = Ref{Float64}(0.0)
    check(
        ccall(
            (:nc_get_att_double, libnetcdf),
            Cint,
            (Cint, Cint, Cstring, Ptr{Float64}),
            ds.ncid,
            id,
            name,
            value,
        ),
        name,
    )
    return value[]
end

"Read a text attribute, returns `nothing` if it does not exist."
function attribute_string(
        ds::Dataset,
        id::Cint,
        name::AbstractString,
    )::Union{Nothing, String}
    len = attlength(ds, id, name)
    len < 0 && return nothing
    type = atttype(ds, id, name)
    if type == NC_CHAR
        buffer = zeros(UInt8, len)
        check(
            ccall(
                (:nc_get_att_text, libnetcdf),
                Cint,
                (Cint, Cint, Cstring, Ptr{UInt8}),
                ds.ncid,
                id,
                name,
                buffer,
            ),
            name,
        )
        # strip trailing null characters
        last_index = something(findlast(!iszero, buffer), 0)
        return String(buffer[1:last_index])
    elseif type == NC_STRING
        len == 1 || error("Attribute `$name` should be a single string in `$(ds.path)`.")
        ptr = Ref{Ptr{UInt8}}(C_NULL)
        check(
            ccall(
                (:nc_get_att_string, libnetcdf),
                Cint,
                (Cint, Cint, Cstring, Ptr{Ptr{UInt8}}),
                ds.ncid,
                id,
                name,
                ptr,
            ),
            name,
        )
        str = unsafe_string(ptr[])
        ccall((:nc_free_string, libnetcdf), Cint, (Csize_t, Ptr{Ptr{UInt8}}), 1, ptr)
        return str
    else
        error("Attribute `$name` is not a string in `$(ds.path)`.")
    end
end

attribute_string(ds::Dataset, varname::AbstractString, name::AbstractString) =
    attribute_string(ds, checked_varid(ds, varname), name)

function put_attribute(
        ds::Dataset,
        id::Cint,
        name::AbstractString,
        value::AbstractString,
    )::Nothing
    str = String(value)
    check(
        ccall(
            (:nc_put_att_text, libnetcdf),
            Cint,
            (Cint, Cint, Cstring, Csize_t, Ptr{UInt8}),
            ds.ncid,
            id,
            name,
            sizeof(str),
            str,
        ),
        name,
    )
    return nothing
end

function put_attribute(ds::Dataset, id::Cint, name::AbstractString, value::Float32)::Nothing
    ref = Ref(value)
    check(
        ccall(
            (:nc_put_att_float, libnetcdf),
            Cint,
            (Cint, Cint, Cstring, Cint, Csize_t, Ptr{Float32}),
            ds.ncid,
            id,
            name,
            NC_FLOAT,
            1,
            ref,
        ),
        name,
    )
    return nothing
end

function put_attribute(ds::Dataset, id::Cint, name::AbstractString, value::Float64)::Nothing
    ref = Ref(value)
    check(
        ccall(
            (:nc_put_att_double, libnetcdf),
            Cint,
            (Cint, Cint, Cstring, Cint, Csize_t, Ptr{Float64}),
            ds.ncid,
            id,
            name,
            NC_DOUBLE,
            1,
            ref,
        ),
        name,
    )
    return nothing
end

put_global_attribute(ds::Dataset, name::AbstractString, value::AbstractString) =
    put_attribute(ds, NC_GLOBAL, name, value)

###
### Reading
###

"""
    read_hyperslab(ds, name, start, count)::Vector{Union{Missing, Float64}}

Read a hyperslab of a variable as a flat vector, in Julia (column major) order. `start`
(1-based) and `count` are given in Julia dimension order. Values equal to `_FillValue` or
`missing_value` become `missing`, and `scale_factor` and `add_offset` are applied.
"""
function read_hyperslab(
        ds::Dataset,
        name::AbstractString,
        start::Vector{Int},
        count::Vector{Int},
    )::Vector{Union{Missing, Float64}}
    id = checked_varid(ds, name)
    @assert length(start) == length(count)
    raw = Vector{Float64}(undef, prod(count; init = 1))
    # netCDF-C uses row major order and 0-based indices
    c_start = Csize_t[s - 1 for s in Iterators.reverse(start)]
    c_count = Csize_t[c for c in Iterators.reverse(count)]
    check(
        ccall(
            (:nc_get_vara_double, libnetcdf),
            Cint,
            (Cint, Cint, Ptr{Csize_t}, Ptr{Csize_t}, Ptr{Float64}),
            ds.ncid,
            id,
            c_start,
            c_count,
            raw,
        ),
        name,
    )
    fill_value = attribute_float(ds, id, "_FillValue")
    missing_value = attribute_float(ds, id, "missing_value")
    scale_factor = attribute_float(ds, id, "scale_factor")
    add_offset = attribute_float(ds, id, "add_offset")
    scale = isnothing(scale_factor) ? 1.0 : scale_factor
    offset = isnothing(add_offset) ? 0.0 : add_offset

    data = Vector{Union{Missing, Float64}}(undef, length(raw))
    for (i, x) in enumerate(raw)
        is_missing =
            (!isnothing(fill_value) && isequal(x, fill_value)) ||
            (!isnothing(missing_value) && isequal(x, missing_value))
        data[i] = is_missing ? missing : x * scale + offset
    end
    return data
end

"Read an entire variable as a flat Vector{Union{Missing, Float64}}."
function read_vector(ds::Dataset, name::AbstractString)::Vector{Union{Missing, Float64}}
    count = var_size(ds, name)
    start = ones(Int, length(count))
    return read_hyperslab(ds, name, start, count)
end

"Read an entire 1D variable, erroring on missing values."
function read_vector_nomissing(ds::Dataset, name::AbstractString)::Vector{Float64}
    data = read_vector(ds, name)
    any(ismissing, data) && error("Missing values found in `$name` in `$(ds.path)`.")
    return Float64[x for x in data]
end

###
### CF time
###

const CF_TIME_FACTORS = (
    ("millisecond", 1),
    ("msec", 1),
    ("second", 1000),
    ("sec", 1000),
    ("minute", 60_000),
    ("min", 60_000),
    ("hour", 3_600_000),
    ("hr", 3_600_000),
    ("day", 86_400_000),
)

"""
    parse_time_units(units::AbstractString)::Tuple{Int, DateTime}

Parse CF time units like "days since 2000-01-01 00:00:00", returning the number of
milliseconds per unit and the time origin.
"""
function parse_time_units(units::AbstractString)::Tuple{Int, DateTime}
    parts = split(strip(units), " since "; limit = 2)
    length(parts) == 2 || error("Unsupported time units `$units`.")
    unit_name = lowercase(strip(parts[1]))
    unit_name = endswith(unit_name, 's') ? unit_name[1:(end - 1)] : unit_name
    factor = 0
    for (name, ms) in CF_TIME_FACTORS
        if unit_name == name
            factor = ms
            break
        end
    end
    factor == 0 && error("Unsupported time unit `$(parts[1])` in `$units`.")
    return factor, parse_time_origin(parts[2])
end

function parse_time_origin(str::AbstractString)::DateTime
    str = String(strip(str))
    # drop a trailing time zone indicator, only UTC is supported
    for suffix in ("UTC", "Z")
        if endswith(str, suffix)
            str = String(rstrip(str[1:(end - length(suffix))]))
        end
    end
    str = replace(str, 'T' => ' ')
    date_time = split(str, ' '; keepempty = false)
    ymd = split(date_time[1], '-')
    length(ymd) == 3 || error("Unsupported time origin `$str`.")
    year = parse(Int, ymd[1])
    month = parse(Int, ymd[2])
    day = parse(Int, ymd[3])
    hour = minute = 0
    millisecond = 0
    if length(date_time) > 1
        hms = split(date_time[2], ':')
        hour = parse(Int, hms[1])
        minute = length(hms) > 1 ? parse(Int, hms[2]) : 0
        if length(hms) > 2
            millisecond = round(Int, parse(Float64, hms[3]) * 1000)
        end
    end
    return DateTime(year, month, day, hour, minute) + Millisecond(millisecond)
end

const STANDARD_CALENDARS = ("standard", "gregorian", "proleptic_gregorian")

"Check that the calendar is supported, only standard calendars are."
function check_calendar(calendar::AbstractString)::Nothing
    if lowercase(calendar) ∉ STANDARD_CALENDARS
        error(
            "Calendar `$calendar` is not supported, only $(join(STANDARD_CALENDARS, ", ")).",
        )
    end
    return nothing
end

"Decode CF time values (`Float64`) to DateTime, given the units string."
function decode_time(values::Vector{Float64}, units::AbstractString)::Vector{DateTime}
    factor, origin = parse_time_units(units)
    return DateTime[origin + Millisecond(round(Int64, v * factor)) for v in values]
end

"Encode a DateTime as a CF time value given the units string."
function encode_time(time::DateTime, units::AbstractString)::Float64
    factor, origin = parse_time_units(units)
    return Dates.value(time - origin) / factor
end

"Read a time coordinate variable and decode it to DateTime."
function read_time(ds::Dataset, name::AbstractString)::Vector{DateTime}
    id = checked_varid(ds, name)
    units = attribute_string(ds, id, "units")
    isnothing(units) && error("Time variable `$name` in `$(ds.path)` has no units.")
    calendar = attribute_string(ds, id, "calendar")
    !isnothing(calendar) && check_calendar(calendar)
    return decode_time(read_vector_nomissing(ds, name), units)
end

"Whether a variable has CF time units, of the form `<unit> since <origin>`."
function is_time_variable(ds::Dataset, name::AbstractString)::Bool
    units = attribute_string(ds, checked_varid(ds, name), "units")
    return !isnothing(units) && occursin(" since ", units)
end

###
### Writing
###

"Define a dimension, `len = 0` gives an unlimited dimension."
function def_dim(ds::Dataset, name::AbstractString, len::Integer)::Cint
    id = Ref{Cint}(0)
    check(
        ccall(
            (:nc_def_dim, libnetcdf),
            Cint,
            (Cint, Cstring, Csize_t, Ptr{Cint}),
            ds.ncid,
            name,
            Csize_t(len),
            id,
        ),
        name,
    )
    return id[]
end

"""
    def_var(ds, name, T, dimnames; deflatelevel = 0)::Cint

Define a variable with element type `T` and dimensions `dimnames` in Julia (column major)
order. Dimensions that do not exist yet are created with the length `dimlengths` if given.
"""
function def_var(
        ds::Dataset,
        name::AbstractString,
        ::Type{T},
        dimnames::Vector{String};
        deflatelevel::Int = 0,
    )::Cint where {T}
    dimids = Cint[]
    for dim in Iterators.reverse(dimnames)
        id = dimid(ds, dim)
        id < 0 && error("Dimension `$dim` not defined in `$(ds.path)`.")
        push!(dimids, id)
    end
    id = Ref{Cint}(0)
    check(
        ccall(
            (:nc_def_var, libnetcdf),
            Cint,
            (Cint, Cstring, Cint, Cint, Ptr{Cint}, Ptr{Cint}),
            ds.ncid,
            name,
            nc_type(T),
            length(dimids),
            dimids,
            id,
        ),
        name,
    )
    if deflatelevel > 0
        check(
            ccall(
                (:nc_def_var_deflate, libnetcdf),
                Cint,
                (Cint, Cint, Cint, Cint, Cint),
                ds.ncid,
                id[],
                Cint(1),
                Cint(1),
                Cint(deflatelevel),
            ),
            name,
        )
    end
    return id[]
end

"""
    write_hyperslab(ds, name, data, start, count)

Write `data` to a hyperslab of variable `name`. `start` (1-based) and `count` are given in
Julia (column major) order. The data is converted to the variable type by netCDF-C.
"""
function write_hyperslab(
        ds::Dataset,
        name::AbstractString,
        data::AbstractArray{Float64},
        start::Vector{Int},
        count::Vector{Int},
    )::Nothing
    id = checked_varid(ds, name)
    @assert length(data) == prod(count; init = 1)
    buffer = data isa Array{Float64} ? data : collect(data)
    c_start = Csize_t[s - 1 for s in Iterators.reverse(start)]
    c_count = Csize_t[c for c in Iterators.reverse(count)]
    check(
        ccall(
            (:nc_put_vara_double, libnetcdf),
            Cint,
            (Cint, Cint, Ptr{Csize_t}, Ptr{Csize_t}, Ptr{Float64}),
            ds.ncid,
            id,
            c_start,
            c_count,
            buffer,
        ),
        name,
    )
    return nothing
end

"Write a vector of strings to a 1D NC_STRING variable."
function write_strings(ds::Dataset, name::AbstractString, data::Vector{String})::Nothing
    id = checked_varid(ds, name)
    pointers = Ptr{UInt8}[pointer(s) for s in data]
    GC.@preserve data begin
        check(
            ccall(
                (:nc_put_vara_string, libnetcdf),
                Cint,
                (Cint, Cint, Ptr{Csize_t}, Ptr{Csize_t}, Ptr{Ptr{UInt8}}),
                ds.ncid,
                id,
                Csize_t[0],
                Csize_t[length(data)],
                pointers,
            ),
            name,
        )
    end
    return nothing
end

end # module NC
