
"Map from PCRaster LDD value to a CartesianIndex"
const pcrdir = [
    CartesianIndex(1, -1),  # 1
    CartesianIndex(1, 0),  # 2
    CartesianIndex(1, 1),  # 3
    CartesianIndex(0, -1),  # 4
    CartesianIndex(0, 0),  # 5
    CartesianIndex(0, 1),  # 6
    CartesianIndex(-1, -1),  # 7
    CartesianIndex(-1, 0),  # 8
    CartesianIndex(-1, 1),  # 9
]

const mv = NaN

# timestep that the parameter units are defined in
const basetimestep = Second(Day(1))

# in case the input data is permuted the other way
permute_indices(inds) = [CartesianIndex(i[2], i[1]) for i in inds]

"Set at indices pit values (default = 5) in a gridded local drainage direction vector"
function set_pit_ldd(pits_2d, ldd, indices; pit = 5)
    pits = pits_2d[indices]
    index = filter(i -> isequal(pits[i], true), 1:length(indices))
    ldd[index] .= pit
    return ldd
end

"Filter upstream neighbors of graph based on logical vector"
function filter_upsteam_nodes(graph, vec_logical::Vector{Bool})
    upstream_nodes = Vector{Int}[]
    for v in topological_sort_by_dfs(graph)
        ups_nodes = inneighbors(graph, v)
        push!(upstream_nodes, filter(i -> !vec_logical[i], ups_nodes))
    end
    return upstream_nodes
end

"""
    active_indices(subcatch_2d, nodata)

Takes a 2D array of the subcatchments. And derive forward and reverse indices.

1: Get a list of `CartesianIndex{2}`` that are active, based on a nodata value.
These map from the 1D internal domain to the 2D external domain.

2: Make a reverse index, a `Matrix{Int}``, which maps from the 2D external domain to
the 1D internal domain, providing an Int which can be used as a linear index. Values of 0
represent inactive cells.
"""
function active_indices(subcatch_2d::AbstractMatrix, nodata)
    A = subcatch_2d
    all_inds = CartesianIndices(size(A))
    indices = filter(i -> !isequal(A[i], nodata), all_inds)

    reverse_indices = zeros(Int, size(A))
    for (i, I) in enumerate(indices)
        reverse_indices[I] = i
    end

    return indices, reverse_indices
end

function active_indices(network, key::Tuple)
    if :reservoir in key
        return network.reservoir.indices_outlet
    elseif :lake in key
        return network.lake.indices_outlet
    elseif :river in key
        return network.river.indices
    elseif :drain in key
        return network.drain.indices
    else
        return network.land.indices
    end
end

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
    set_states(instate_path, model, state_ncnames; <keyword arguments>)

Read states contained in `Dict` `state_ncnames` from NetCDF file located in `instate_path`,
and set states in `model` object. Active cells are selected with the corresponding network's
(`Vector{CartesianIndex}`) from the NetCDF file. 

# Arguments
- `type = nothing`: type to convert data to after reading. By default no conversion is done.
"""
function set_states(instate_path, model, state_ncnames; type = nothing)
    @unpack network = model
    # states in NetCDF include dim time (one value) at index 3 or 4, 3 or 4 dims are allowed
    NCDataset(instate_path) do ds
        for (state, ncname) in state_ncnames
            sel = active_indices(network, state)
            dims = length(dimnames(ds[ncname]))
            # 4 dims, for example (x,y,layer,time) where dim layer is an SVector for soil layers
            if dims == 4
                A = permutedims(ds[ncname][sel, :, 1])
                # note that this array is allowed to have missings, since not every vertical
                # column is `maxlayers` layers deep
                A = replace!(A, missing => NaN)
                # Convert to desired type if needed
                if !isnothing(type)
                    if eltype(A) != type
                        A = map(type, A)
                    end
                end
                # set state in model object
                param(model, state) .= svectorscopy(A, Val{size(A)[1]}())
                # 3 dims (x,y,time)
            elseif dims == 3
                A = ds[ncname][sel, 1]
                A = nomissing(A)
                # Convert to desired type if needed
                if !isnothing(type)
                    if eltype(A) != type
                        A = map(type, A)
                    end
                end
                # set state in model object
                param(model, state) .= A
            else
                error(
                    "Number of state dims should be 3 or 4, number of dims = ",
                    string(dims),
                )
            end
        end
    end
end

"""
    ncread(nc, var; <keyword arguments>)

Read parameter `var` from NetCDF file `nc`. Supports various keyword arguments to get
selections of data in desired types, with or without missing values.

# Arguments
- `key=identity`: `key` is a function which maps `var` to the name of the parameter in the
        NetCDF file. The default `identity` keeps `var` as is.
- `sel=nothing`: a selection of indices, such as a `Vector{CartesianIndex}` of active cells,
        to return from the NetCDF. By default all cells are returned.
- `defaults=nothing`: a dictionary in which default values are looked up if `var` is not
        in `nc`. By default it gives an error in this case.
- `type=nothing`: type to convert data to after reading. By default no conversion is done.
- `allow_missing=false`: Missing values within `sel` is not allowed by default. Set to
        `true` to allow missing values.
- `fill=nothing`: Missing values are replaced by this fill value if `allow_missing` is `false`.
- `dimname` : name of third dimension of parameter `var`. By default no third dimension is expected.
"""
function ncread(
    nc,
    var;
    sel = nothing,
    defaults = nothing,
    type = nothing,
    allow_missing = false,
    fill = nothing,
    dimname = nothing,
)

    if isnothing(var)
        @assert !isnothing(defaults)
        if isnothing(dimname)
            return Base.fill(defaults, length(sel))
        else
            return Base.fill(defaults, (nc.dim[dimname], length(sel)))
        end
    end

    # Read the entire variable into memory, applying scale, offset and
    # set fill_values to missing.
    A = nc[var][:]

    # Take out only the active cells
    if !isnothing(sel)
        if isnothing(dimname)
            A = A[sel]
        else
            dim = findfirst(==(dimname), dimnames(nc[var]))
            if dim == 3
                A = permutedims(A[sel, :])
            elseif dim == 1
                A = A[:, sel]
            end
        end
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
            A = convert(Array{type}, A)
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

"Return the sum of the function `f` at indices `inds`"
function sum_at(f::Function, inds, T)
    v = zero(T)
    for i in inds
        v += f(i)
    end
    return v
end

# https://juliaarrays.github.io/StaticArrays.jl/latest/pages/api/#Arrays-of-static-arrays-1
function svectorscopy(x::Matrix{T}, ::Val{N}) where {T,N}
    size(x, 1) == N || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    copy(reinterpret(SVector{N,T}, vec(x)))
end

"""
    fraction_runoff_toriver(graph, ldd, index_river, slope, n)

Determine ratio `frac` between `slope` river cell `index_river` and `slope` of each upstream neighbor (based on directed acyclic graph
`graph`).
"""
function fraction_runoff_toriver(graph, ldd, index_river, slope, n)
    frac = zeros(n)
    for i in index_river
        nbs = inneighbors(graph, i)
        for j in nbs
            if ldd[j] != ldd[i]
                frac[j] = slope[j] / (slope[i] + slope[j])
            end
        end
    end
    return frac
end

"""
    equal_size_vectors(x)

Used in the structs of arrays to ensure all vectors are of equal length.

`equal_size_vectors(([1,2], [1,2,3]))` would throw an ArgumentError.
`equal_size_vectors(([4,5], [4,5]))` would pass.
`equal_size_vectors((1, [4,5], [4,5]))` would also pass, since `1` is not an AbstractVector.
"""
function equal_size_vectors(x)
    # all vectors in this struct should be the same size
    inds_vec = findall(arg -> isa(arg, AbstractVector), x)
    n = length(x[inds_vec[1]])
    x_vec = x[inds_vec]

    for arr in x_vec
        if length(arr) != n
            throw(ArgumentError("Not all vectors are of equal length"))
        end
    end
    return x
end

"""
    tosecond(x::Period)

Convert a Period into a Float64, which represents the number of seconds. Will fail if this
is not well defined, such as for Month.

# Examples
```julia-repl
julia> tosecond(Day(1))
86400.0
```
"""
tosecond(x::Hour) = Float64(Dates.value(Second(x)))
tosecond(x::Minute) = Float64(Dates.value(Second(x)))
tosecond(x::T) where {T<:DatePeriod} = Float64(Dates.value(Second(x)))
tosecond(x::T) where {T<:TimePeriod} = x / convert(T, Second(1))
