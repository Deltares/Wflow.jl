
function lattometres(lat::Float64)
    m1 = 111132.92     # latitude calculation term 1
    m2 = -559.82       # latitude calculation term 2
    m3 = 1.175         # latitude calculation term 3
    m4 = -0.0023       # latitude calculation term 4
    p1 = 111412.84     # longitude calculation term 1
    p2 = -93.5         # longitude calculation term 2
    p3 = 0.118         # longitude calculation term 3

    # Calculate the length of a degree of latitude and longitude in meters
    latlen = m1 + (m2 * cosd(2.0 * lat)) + (m3 * cosd(4.0 * lat)) + (m4 * cosd(6.0 * lat))
    longlen = (p1 * cosd(lat)) + (p2 * cosd(3.0 * lat)) + (p3 * cosd(5.0 * lat))

    return longlen, latlen
end

"""
    ncread(nc, var; <keyword arguments>)

Read parameter `var` from NetCDF file `nc`. Supports various keyword arguments to get
selections of data in desired types, with or without missing values.

# Arguments
- `sel=nothing`: a selection of indices, such as a `Vector{CartesianIndex}` of active cells,
        to return from the NetCDF. By default all cells are returned.
- `defaults=nothing`: a dictionary in which default values are looked up if `var` is not
        in `nc`. By default it gives an error in this case.
- `type=nothing`: type to convert data to after reading. By default no conversion is done.
- `allow_missing=false`: Missing values within `sel` is not allowed by default. Set to
        `true` to allow missing values.
- `fill=nothing`: Missing values are replaced by this fill value if `allow_missing` is `false`.
"""
function ncread(
    nc,
    var;
    sel = nothing,
    defaults = nothing,
    type = nothing,
    allow_missing = false,
    fill = nothing,
)

    if !haskey(nc, var) && !isnothing(defaults)
        # TODO move away from this strategy for defaults
        @warn(string(var, " not found, set to default value ", defaults[var]))
        return Fill(defaults[var], length(sel))
    end

    # Read the entire variable into memory, applying scale, offset and
    # set fill_values to missing.
    A = nc[var][:]

    # Take out only the active cells
    if !isnothing(sel)
        A = A[sel]
    end

    if !allow_missing
        if isnothing(fill)
            # errors if missings are found
            A = nomissing(A)
        else
            # replaces missings with a fill value
            A = nomissing(A, fill)
        end
    end

    # Convert to desired type if needed
    if !isnothing(type)
        if eltype(A) != type
            A = map(type, A)
        end
    end

    return A
end

"""
    set_layerthickness(d::Float64, sl::SVector)

Calculate actual soil thickness of layers based on a reference depth (e.g. soil depth or water table depth) `d`,
a SVector `sl` with cumulative soil depth starting at soil surface (0), and a SVector `tl` with actual thickness
per soil layer.
"""
function set_layerthickness(d::Float64, sl::SVector, tl::SVector)

    act_d = tl .* mv
    for i = 1:length(act_d)
        if d > sl[i+1]
            act_d = setindex(act_d, tl[i], i)
        elseif d - sl[i] > 0.0
            act_d = setindex(act_d, d - sl[i], i)
        end
    end

    nlayers = length(act_d) - sum(isnan.(act_d))
    return act_d, nlayers
end

"""
    detdrainlength(ldd, xl, yl)

Determines the drainaige length for a non square grid. Input `ldd` (drainage network), `xl` (length of cells in x direction),
`yl` (length of cells in y direction). Output is drainage length.
"""
function detdrainlength(ldd, xl, yl)
    # take into account non-square cells
    # if ldd is 8 or 2 use ylength
    # if ldd is 4 or 6 use xlength
    if ldd == 2 || ldd == 8
        yl
    elseif ldd == 4 || ldd == 6
        xl
    else
        hypot(xl, yl)
    end
end

"""
    detdrainwidth(ldd, xl, yl)

Determines the drainaige width for a non square grid. Input `ldd` (drainage network), `xl` (length of cells in x direction),
`yl` (length of cells in y direction). Output is drainage width.
"""
function detdrainwidth(ldd, xl, yl)
    # take into account non-square cells
    # if ldd is 8 or 2 use xlength
    # if ldd is 4 or 6 use ylength
    slantwidth = (xl + yl) * 0.5
    if ldd == 2 || ldd == 8
        xl
    elseif ldd == 4 || ldd == 6
        yl
    else
        slantwidth
    end
end

# 2.5x faster power method
"Faster method for exponentiation"
pow(x, y) = exp(y * log(x))

"Return the sum of the array `A` at indices `inds`"
function sum_at(A, inds)
    v = zero(eltype(A))
    for i in inds
        v += A[i]
    end
    return v
end

"Return the sum of the array `A`, if 'A' is empty return 0.0"
function sum_empty(A)
    if length(A) > 0
        sum(A)
    else
        0.0
    end
end

# https://juliaarrays.github.io/StaticArrays.jl/latest/pages/api/#Arrays-of-static-arrays-1
function svectorscopy(x::Matrix{T}, ::Val{N}) where {T,N}
    size(x, 1) == N || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    copy(reinterpret(SVector{N,T}, vec(x)))
end

"""
    fraction_runoff_toriver(dag, index_river, slope, n)

Determine ratio `frac` between slope river cell `index_river` and slope of each upstream neighbor (based on directed acyclic graph
`dag`).
"""
function fraction_runoff_toriver(dag, index_river, slope, n)
    frac = zeros(n)
    for i in index_river
        nbs = inneighbors(dag, i)
        frac[nbs] = slope[i] ./ (slope[i] .+ slope[nbs])
    end
    return frac
end
