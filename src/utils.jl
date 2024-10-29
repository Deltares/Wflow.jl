
"Map from PCRaster LDD value to a CartesianIndex"
const pcr_dir = [
    CartesianIndex(-1, -1),  # 1
    CartesianIndex(0, -1),  # 2
    CartesianIndex(1, -1),  # 3
    CartesianIndex(-1, 0),  # 4
    CartesianIndex(0, 0),  # 5
    CartesianIndex(1, 0),  # 6
    CartesianIndex(-1, 1),  # 7
    CartesianIndex(0, 1),  # 8
    CartesianIndex(1, 1),  # 9
]

const mv = Float(NaN)

# timestep that the parameter units are defined in
const basetimestep = Second(Day(1))

"""
    scurve(x, a, b, c)

Sigmoid "S"-shaped curve.

# Arguments
- `x::Real`: input
- `a::Real`: determines the centre level
- `b::Real`: determines the amplitude of the curve
- `c::Real`: determines the steepness or "stepwiseness" of the curve.
             The higher c the sharper the function. A negative c reverses the function.
"""
function scurve(x, a, b, c)
    s = one(x) / (b + exp(-c * (x - a)))
    return s
end

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

function active_indices(network::NamedTuple, key::Tuple)
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

function lattometres(lat::Real)
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

function cell_lengths(y::AbstractVector, cellength::Real, sizeinmetres::Bool)
    n = length(y)
    xl = fill(mv, n)
    yl = fill(mv, n)
    if sizeinmetres
        xl .= cellength
        yl .= cellength
    else
        for i in 1:n
            longlen, latlen = lattometres(y[i])
            xl[i] = longlen * cellength
            yl[i] = latlen * cellength
        end
    end
    return xl, yl
end

function river_fraction(
    river::AbstractVector,
    riverlength::AbstractVector,
    riverwidth::AbstractVector,
    xl::AbstractVector,
    yl::AbstractVector,
)
    n = length(river)
    riverfrac = fill(mv, n)
    for i in 1:n
        riverfrac[i] = if river[i]
            min((riverlength[i] * riverwidth[i]) / (xl[i] * yl[i]), 1.0)
        else
            0.0
        end
    end
    return riverfrac
end

"""
    set_states!(instate_path, model, state_ncnames; <keyword arguments>)

Read states contained in `Dict` `state_ncnames` from netCDF file located in `instate_path`,
and set states in `model` object. Active cells are selected with the corresponding network's
(`Vector{CartesianIndex}`) from the netCDF file.

# Arguments
- `type = nothing`: type to convert data to after reading. By default no conversion is done.
"""
function set_states!(instate_path, model; type = nothing, dimname = nothing)
    (; network, config) = model

    # Check if required states are covered
    state_ncnames = check_states(config)

    # states in netCDF include dim time (one value) at index 3 or 4, 3 or 4 dims are allowed
    NCDataset(instate_path) do ds
        for (state, ncname) in state_ncnames
            @info "Setting initial state from netCDF." ncpath = instate_path ncvarname =
                ncname state
            sel = active_indices(network, state)
            n = length(sel)
            dims = length(dimnames(ds[ncname]))
            # 4 dims, for example (x,y,layer,time) where dim layer is an SVector for soil layers
            if dims == 4
                if dimname == :layer
                    dimensions = (x = :, y = :, layer = :, time = 1)
                else
                    error("Unrecognized dimension name $dimname")
                end
                A = read_standardized(ds, ncname, dimensions)
                A = permutedims(A[sel, :])
                # note that this array is allowed to have missings, since not every vertical
                # column is `maxlayers` layers deep
                if dimname == :layer
                    A = replace!(A, missing => NaN)
                end
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
                A = read_standardized(ds, ncname, (x = :, y = :, time = 1))
                A = A[sel]
                A = nomissing(A)
                # Convert to desired type if needed
                if !isnothing(type)
                    if eltype(A) != type
                        A = map(type, A)
                    end
                end
                # set state in model object, only set active cells ([1:n]) (ignore boundary conditions/ghost points)
                param(model, state)[1:n] .= A
            else
                error(
                    "Number of state dims should be 3 or 4, number of dims = ",
                    string(dims),
                )
            end
        end
    end
    return nothing
end

"""
    ncread(nc, config::Config, parameter::AbstractString; <keyword arguments>)

Read a netCDF variable `var` from file `nc`, based on `config` (parsed TOML file) and the
model `parameter` specified in the TOML configuration file. Supports various keyword
arguments to get selections of data in desired types, with or without missing values.

# Arguments
- `alias`=nothing` : An `alias` for the TOML key, by default an `alias` is not expected.
- `optional=true` : By default specifying a model `parameter` in the TOML file is optional.
        Set to false if the model `parameter` is required.
- `sel=nothing`: A selection of indices, such as a `Vector{CartesianIndex}` of active cells,
        to return from the netCDF. By default all cells are returned.
- `defaults=nothing`: A default value if `var` is not in `nc`. By default it gives an error
    in this case.
- `type=nothing`: Type to convert data to after reading. By default no conversion is done.
- `allow_missing=false`: Missing values within `sel` is not allowed by default. Set to
        `true` to allow missing values.
- `fill=nothing`: Missing values are replaced by this fill value if `allow_missing` is
  `false`.
- `dimname` : Name of third dimension of parameter `var`. By default no third dimension is
  expected.
"""
function ncread(
    nc,
    config::Config,
    parameter::AbstractString;
    alias = nothing,
    optional = true,
    sel = nothing,
    defaults = nothing,
    type = nothing,
    allow_missing = false,
    fill = nothing,
    dimname = nothing,
)
    # get var (netCDF variable or type Config) from TOML file.
    # if var has type Config, input parameters can be changed.
    if isnothing(alias)
        if optional
            var = param(config.input, parameter, nothing)
        else
            var = param(config.input, parameter)
        end
    else
        var = get_alias(config.input, parameter, alias, nothing)
    end

    # dim `time` is also included in `dim_sel`: this allows for cyclic parameters (read
    # first timestep), that is later updated with the `update_cyclic!` function.
    if isnothing(dimname)
        dim_sel = (x = :, y = :, time = 1)
    elseif dimname == :layer
        dim_sel = (x = :, y = :, layer = :, time = 1)
    elseif dimname == :flood_depth
        dim_sel = (x = :, y = :, flood_depth = :, time = 1)
    else
        error("Unrecognized dimension name $dimname")
    end

    if isnothing(var)
        @info "Set `$parameter` using default value `$defaults`."
        @assert !isnothing(defaults)
        if !isnothing(type)
            defaults = convert(type, defaults)
        end
        if isnothing(dimname)
            return Base.fill(defaults, length(sel))
        else
            return Base.fill(defaults, (nc.dim[String(dimname)], length(sel)))
        end
    end

    # If var has type Config, input parameters can be changed (through scale, offset and
    # input netCDF var) or set to a uniform value (providing a value). Otherwise, input
    # NetCDF var is read directly.
    var, mod = ncvar_name_modifier(var; config = config)

    if !isnothing(mod.value)
        @info "Set `$parameter` using default value `$(mod.value)`."
        if isnothing(dimname)
            return Base.fill(mod.value, length(sel))
            # set to one uniform value
        elseif length(mod.value) == 1
            return Base.fill(mod.value, (nc.dim[String(dimname)], length(sel)))
            # set to vector of values (should be equal to size dimname)
        elseif length(mod.value) > 1
            @assert length(mod.value) == nc.dim[String(dimname)]
            return repeat(mod.value, 1, length(sel))
        end
    else
        @info "Set `$parameter` using netCDF variable `$var`."
        A = read_standardized(nc, var, dim_sel)
        if !isnothing(mod.index)
            # the modifier index is only set in combination with scale and offset for SVectors,
            # provided through the TOML file.
            if length(mod.index) > 1
                # if index, scale and offset is provided in the TOML as a list.
                for i in 1:length(mod.index)
                    A[:, :, mod.index[i]] =
                        A[:, :, mod.index[i]] .* mod.scale[i] .+ mod.offset[i]
                end
            else
                A[:, :, mod.index] = A[:, :, mod.index] .* mod.scale .+ mod.offset
            end
        elseif mod.scale != 1.0 || mod.offset != 0.0
            A = A .* mod.scale .+ mod.offset
        end
    end

    # Take out only the active cells
    if !isnothing(sel)
        if isnothing(dimname)
            A = A[sel]
        else
            A = permutedims(A[sel, :])
        end
    end

    if !allow_missing
        if isnothing(fill)
            # errors if missings are found
            A = nomissing(A)
            if any(isnan, A)
                error("NaN not allowed in $var")
            end
        else
            # replaces missings with a fill value
            A = nomissing(A, fill)
            # replace also NaN values with the fill value
            replace!(x -> isnan(x) ? fill : x, A)
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
    set_layerthickness(reference_depth::Real, cum_depth::SVector, thickness::SVector)

Calculate actual soil thickness of layers based on a reference depth (e.g. soil depth or
water table depth) `reference_depth`, a SVector `cum_depth` with cumulative soil depth starting
at soil surface (0), and a SVector `thickness` with thickness per soil layer.
"""
function set_layerthickness(reference_depth::Real, cum_depth::SVector, thickness::SVector)
    thicknesslayers = thickness .* mv
    for i in 1:length(thicknesslayers)
        if reference_depth > cum_depth[i + 1]
            thicknesslayers = setindex(thicknesslayers, thickness[i], i)
        elseif reference_depth - cum_depth[i] > 0.0
            thicknesslayers = setindex(thicknesslayers, reference_depth - cum_depth[i], i)
        end
    end
    return thicknesslayers
end

function number_of_active_layers(thickness::SVector)
    nlayers = length(thickness) - sum(isnan.(thickness))
    return nlayers
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

Determines the drainaige width for a non square grid. Input `ldd` (drainage network), `xl`
(length of cells in x direction), `yl` (length of cells in y direction). Output is drainage
width.
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

"""
    det_surfacewidth(ldd, xl, yl)

Determines the surface flow width. Input `dw` (drainage width), `riverwidth` and `river`
(boolean). Output is surface flow width `sw`.
"""
function det_surfacewidth(dw, riverwidth, river)
    sw = river ? max(dw - riverwidth, 0.0) : dw
    return sw
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
function svectorscopy(x::Matrix{T}, ::Val{N}) where {T, N}
    size(x, 1) == N || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    return copy(reinterpret(SVector{N, T}, vec(x)))
end

"""
    fraction_runoff_toriver(graph, ldd, index_river, slope, n)

Determine ratio `frac` between `slope` river cell `index_river` and `slope` of each upstream
neighbor (based on directed acyclic graph `graph`).
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
tosecond(x::T) where {T <: DatePeriod} = Float64(Dates.value(Second(x)))
tosecond(x::T) where {T <: TimePeriod} = x / convert(T, Second(1))

"""
    adjacent_nodes_at_link(graph)

Return the source node `src` and destination node `dst` of each link of a directed `graph`.
"""
function adjacent_nodes_at_link(graph)
    links = collect(edges(graph))
    return (src = src.(links), dst = dst.(links))
end

"""
    adjacent_links_at_node(graph, nodes_at_link)

Return the source link `src` and destination link `dst` of each node of a directed `graph`.
"""
function adjacent_links_at_node(graph, nodes_at_link)
    nodes = vertices(graph)
    src_link = Vector{Int}[]
    dst_link = copy(src_link)
    for i in 1:nv(graph)
        push!(src_link, findall(isequal(nodes[i]), nodes_at_link.dst))
        push!(dst_link, findall(isequal(nodes[i]), nodes_at_link.src))
    end
    return (src = src_link, dst = dst_link)
end

"Add `vertex` and `edge` to `pits` of a directed `graph`"
function add_vertex_edge_graph!(graph, pits)
    n = nv(graph)
    for (i, v) in enumerate(pits)
        add_vertex!(graph)
        add_edge!(graph, v, n + i)
    end
    return nothing
end

"""
    set_effective_flowwidth!(we_x, we_y, indices, graph_riv, riverwidth, ldd_riv, inds_rev_riv)

For river cells (D8 flow direction) in a staggered grid the effective flow width at cell
edges (floodplain) `we_x` in the x-direction and `we_y` in the y-direction is corrected by
subtracting the river width `riverwidth` from the cell edges. For diagonal directions, the
`riverwidth` is split between the two adjacent cell edges. A cell edge at linear index `idx`
is defined as the edge between node `idx` and the adjacent node (+ CartesianIndex(1, 0)) for
x and (+ CartesianIndex(0, 1)) for y. For cells that contain a `waterbody` (reservoir or
lake), the effective flow width is set to zero.
"""
function set_effective_flowwidth!(
    we_x,
    we_y,
    indices,
    graph_riv,
    riverwidth,
    ldd_riv,
    waterbody,
    inds_rev_riv,
)
    toposort = topological_sort_by_dfs(graph_riv)
    n = length(we_x)
    for v in toposort
        dst = outneighbors(graph_riv, v)
        isempty(dst) && continue
        w = min(riverwidth[v], riverwidth[only(dst)])
        dir = pcr_dir[ldd_riv[v]]
        idx = inds_rev_riv[v]
        # loop over river D8 directions
        if dir == CartesianIndex(1, 1)
            we_x[idx] = waterbody[v] ? 0.0 : max(we_x[idx] - 0.5 * w, 0.0)
            we_y[idx] = waterbody[v] ? 0.0 : max(we_y[idx] - 0.5 * w, 0.0)
        elseif dir == CartesianIndex(-1, -1)
            if indices.xd[idx] <= n
                we_y[indices.xd[idx]] =
                    waterbody[v] ? 0.0 : max(we_y[indices.xd[idx]] - 0.5 * w, 0.0)
            end
            if indices.yd[idx] <= n
                we_x[indices.yd[idx]] =
                    waterbody[v] ? 0.0 : max(we_x[indices.yd[idx]] - 0.5 * w, 0.0)
            end
        elseif dir == CartesianIndex(1, 0)
            we_y[idx] = waterbody[v] ? 0.0 : max(we_y[idx] - w, 0.0)
        elseif dir == CartesianIndex(0, 1)
            we_x[idx] = waterbody[v] ? 0.0 : max(we_x[idx] - w, 0.0)
        elseif dir == CartesianIndex(-1, 0)
            if indices.xd[idx] <= n
                we_y[indices.xd[idx]] =
                    waterbody[v] ? 0.0 : max(we_y[indices.xd[idx]] - w, 0.0)
            end
        elseif dir == CartesianIndex(0, -1)
            if indices.yd[idx] <= n
                we_x[indices.yd[idx]] =
                    waterbody[v] ? 0.0 : max(we_x[indices.yd[idx]] - w, 0.0)
            end
        elseif dir == CartesianIndex(1, -1)
            we_y[idx] = max(we_y[idx] - 0.5 * w, 0.0)
            if indices.yd[idx] <= n
                we_x[indices.yd[idx]] =
                    waterbody[v] ? 0.0 : max(we_x[indices.yd[idx]] - 0.5 * w, 0.0)
            end
        elseif dir == CartesianIndex(-1, 1)
            if indices.xd[idx] <= n
                we_y[indices.xd[idx]] =
                    waterbody[v] ? 0.0 : max(we_y[indices.xd[idx]] - 0.5 * w, 0.0)
            end
            we_x[idx] = waterbody[v] ? 0.0 : max(we_x[idx] - 0.5 * w, 0.0)
        end
    end
    return nothing
end

"Return julian day of year (leap days are not counted)"
function julian_day(time)
    # for all years February 28 is day 59 and March 1 is day 60.
    day = dayofyear(time) - (isleapyear(time) && dayofyear(time) > 60)
    return day
end

"Partition indices with at least size `basesize`"
function _partition(xs::Integer, basesize::Integer)
    n = Int(max(1, xs ÷ basesize))
    return (Int(1 + ((i - 1) * xs) ÷ n):Int((i * xs) ÷ n) for i in 1:n)
end

"""
    threaded_foreach(f, x::AbstractArray; basesize::Integer)

Run function `f` in parallel by spawning tasks (nthreads <= 8), each task iterates over a
chunk of size `basesize`. For nthreads > 8 run function `f` in parallel with
`Polyester@batch` with `minbatch` equal to `basesize`.
"""
function threaded_foreach(f, x::AbstractArray; basesize::Integer)
    if Threads.nthreads() <= 8
        len = length(x)
        partitions = _partition(len, basesize)
        if length(partitions) > 1 && Threads.nthreads() > 1
            @sync for p in partitions
                Threads.@spawn begin
                    for i in eachindex(p)
                        f(@inbounds p[i])
                    end
                end
            end
        else
            for i in eachindex(x)
                f(@inbounds x[i])
            end
        end
    else
        @batch per = thread minbatch = basesize for i in eachindex(x)
            f(@inbounds x[i])
        end
    end
    return nothing
end

"""
    hydraulic_conductivity_at_depth(p::KvExponential, kvfrac, z, i, n)
    hydraulic_conductivity_at_depth(p::KvExponentialConstant, kvfrac, z, i, n)
    hydraulic_conductivity_at_depth(p::KvLayered, kvfrac, z, i, n)
    hydraulic_conductivity_at_depth(p::KvLayeredExponential, kvfrac, z, i, n)

Return vertical hydraulic conductivity `kv_z` at depth `z` for index `i` using muliplication
factor `kv_frac` at soil layer `n` and vertical hydraulic conductivity profile `p`.
"""
function hydraulic_conductivity_at_depth(p::KvExponential, kvfrac, z, i, n)
    kv_z = kvfrac[i][n] * p.kv_0[i] * exp(-p.f[i] * z)
    return kv_z
end

function hydraulic_conductivity_at_depth(p::KvExponentialConstant, kvfrac, z, i, n)
    (; kv_0, f) = p.exponential
    if z < p.z_exp[i]
        kv_z = kvfrac[i][n] * kv_0[i] * exp(-f[i] * z)
    else
        kv_z = kvfrac[i][n] * kv_0[i] * exp(-f[i] * p.z_exp[i])
    end
    return kv_z
end

function hydraulic_conductivity_at_depth(p::KvLayered, kvfrac, z, i, n)
    kv_z = kvfrac[i][n] * p.kv[i][n]
    return kv_z
end

function hydraulic_conductivity_at_depth(p::KvLayeredExponential, kvfrac, z, i, n)
    return if z < p.z_layered[i]
        kvfrac[i][n] * p.kv[i][n]
    else
        n = p.nlayers_kv[i]
        kvfrac[i][n] * p.kv[i][n] * exp(-p.f[i] * (z - p.z_layered[i]))
    end
end

"""
    kh_layered_profile!(soil::SbmSoilModel, subsurface::LateralSSF, kv_profile::KvLayered, dt)
    kh_layered_profile!(soil::SbmSoilModel, subsurface::LateralSSF, kv_profile::KvLayeredExponential, dt)

Compute equivalent horizontal hydraulic conductivity `kh` [m d⁻¹] using vertical hydraulic
conductivity profile `kv_profile`.
"""
function kh_layered_profile!(
    soil::SbmSoilModel,
    subsurface::LateralSSF,
    kv_profile::KvLayered,
    dt,
)
    (; nlayers, sumlayers, act_thickl, soilthickness) = soil.parameters
    (; n_unsatlayers, zi) = soil.variables
    (; kh) = subsurface.kh_profile
    (; khfrac) = subsurface

    t_factor = (tosecond(basetimestep) / dt)
    for i in eachindex(kh)
        m = nlayers[i]

        if soilthickness[i] > zi[i]
            transmissivity = 0.0
            _sumlayers = @view sumlayers[i][2:end]
            n = max(n_unsatlayers[i], 1)
            transmissivity += (_sumlayers[n] - zi[i]) * kv_profile.kv[i][n]
            n += 1
            while n <= m
                transmissivity += act_thickl[i][n] * kv_profile.kv[i][n]
                n += 1
            end
            # convert units for kh [m d⁻¹] computation (transmissivity [mm² Δt⁻¹], soilthickness
            # [mm] and zi [mm])
            kh[i] =
                0.001 * (transmissivity / (soilthickness[i] - zi[i])) * t_factor * khfrac[i]
        else
            kh[i] = 0.001 * kv_profile.kv[i][m] * t_factor * khfrac[i]
        end
    end
    return nothing
end

function kh_layered_profile!(
    soil::SbmSoilModel,
    subsurface::LateralSSF,
    kv_profile::KvLayeredExponential,
    dt,
)
    (; nlayers, sumlayers, act_thickl, soilthickness) = soil.parameters
    (; nlayers_kv, z_layered, kv, f) = kv_profile
    (; n_unsatlayers, zi) = soil.variables
    (; kh) = subsurface.kh_profile
    (; khfrac) = subsurface
    t_factor = (tosecond(basetimestep) / dt)

    for i in eachindex(kh)
        m = nlayers[i]

        if soilthickness[i] > zi[i]
            transmissivity = 0.0
            n = max(n_unsatlayers[i], 1)
            if zi[i] >= z_layered[i]
                zt = soilthickness[i] - z_layered[i]
                j = nlayers_kv[i]
                transmissivity +=
                    kv[i][j] / f[i] *
                    (exp(-f[i] * (zi[i] - z_layered[i])) - exp(-f[i] * zt))
                n = m
            else
                _sumlayers = @view sumlayers[i][2:end]
                transmissivity += (_sumlayers[n] - zi[i]) * kv[i][n]
            end
            n += 1
            while n <= m
                if n > nlayers_kv[i]
                    zt = soilthickness[i] - z_layered[i]
                    j = nlayers_kv[i]
                    transmissivity += kv[i][j] / f[i] * (1.0 - exp(-f[i] * zt))
                    n = m
                else
                    transmissivity += act_thickl[i][n] * kv[i][n]
                end
                n += 1
            end
            # convert units for kh [m d⁻¹] computation (transmissivity [mm² Δt⁻¹], soilthickness
            # [mm] and zi [mm])
            kh[i] =
                0.001 * (transmissivity / (soilthickness[i] - zi[i])) * t_factor * khfrac[i]
        else
            if zi[i] >= z_layered[i]
                j = nlayers_kv[i]
                kh[i] =
                    0.001 *
                    kv[i][j] *
                    exp(-f[i] * (zi[i] - z_layered[i])) *
                    khfrac[i] *
                    t_factor
            else
                kh[i] = 0.001 * kv[i][m] * t_factor * khfrac[i]
            end
        end
    end
    return nothing
end

kh_layered_profile!(
    soil::SbmSoilModel,
    subsurface::LateralSSF,
    kv_profile::Union{KvExponential, KvExponentialConstant},
    dt,
) = nothing

"""
    initialize_lateralssf!(model::LateralSSF, kh_profile::KhExponential)
    initialize_lateralssf!(model::LateralSSF, kh_profile::KhExponentialConstant)

Initialize lateral subsurface variables `ssf` and `ssfmax` using horizontal hydraulic
conductivity profile `kh_profile`.
"""
function initialize_lateralssf!(model::LateralSSF, kh_profile::KhExponential)
    (; kh_0, f) = kh_profile
    (; ssf, ssfmax, zi, slope, soilthickness, dw) = model

    @. ssfmax = ((kh_0 * slope) / f) * (1.0 - exp(-f * soilthickness))
    @. ssf = ((kh_0 * slope) / f) * (exp(-f * zi) - exp(-f * soilthickness)) * dw
    return nothing
end

function initialize_lateralssf!(model::LateralSSF, kh_profile::KhExponentialConstant)
    (; kh_0, f) = kh_profile.exponential
    (; z_exp) = kh_profile
    (; ssf, ssfmax, zi, slope, soilthickness, dw) = model
    ssf_constant = @. kh_0 * exp(-f * z_exp) * slope * (soilthickness - z_exp)
    for i in eachindex(ssf)
        ssfmax[i] =
            ((kh_0[i] * slope[i]) / f[i]) * (1.0 - exp(-f[i] * z_exp[i])) + ssf_constant[i]
        if zi[i] < z_exp[i]
            ssf[i] =
                (
                    ((kh_0[i] * slope[i]) / f[i]) *
                    (exp(-f[i] * zi[i]) - exp(-f[i] * z_exp[i])) + ssf_constant[i]
                ) * dw[i]
        else
            ssf[i] =
                kh_0[i] * exp(-f[i] * zi[i]) * slope[i] * (soilthickness[i] - zi[i]) * dw[i]
        end
    end
    return nothing
end

"""
    initialize_lateralssf!(subsurface::LateralSSF, soil::SbmSoilModel, kv_profile::KvLayered, dt)
    initialize_lateralssf!(subsurface::LateralSSF, soil::SbmSoilModel, kv_profile::KvLayeredExponential, dt)

Initialize lateral subsurface variables `ssf` and `ssfmax` using  vertical hydraulic
conductivity profile `kv_profile`.
"""
function initialize_lateralssf!(
    subsurface::LateralSSF,
    soil::SbmSoilModel,
    kv_profile::KvLayered,
    dt,
)
    (; kh) = subsurface.kh_profile
    (; nlayers, act_thickl) = soil.parameters
    (; ssf, ssfmax, zi, khfrac, soilthickness, slope, dw) = subsurface
    kh_layered_profile!(soil, subsurface, kv_profile, dt)
    for i in eachindex(ssf)
        ssf[i] = kh[i] * (soilthickness[i] - zi[i]) * slope[i] * dw[i]
        kh_max = 0.0
        for j in 1:nlayers[i]
            kh_max += kv_profile.kv[i][j] * act_thickl[i][j]
        end
        kh_max = kh_max * khfrac[i] * 0.001 * 0.001
        ssfmax[i] = kh_max * slope[i]
    end
    return nothing
end

function initialize_lateralssf!(
    subsurface::LateralSSF,
    soil::SbmSoilModel,
    kv_profile::KvLayeredExponential,
    dt,
)
    (; ssf, ssfmax, zi, khfrac, soilthickness, slope, dw) = subsurface
    (; nlayers, act_thickl) = soil.parameters
    (; kh) = subsurface.kh_profile
    (; kv, f, nlayers_kv, z_layered) = kv_profile

    kh_layered_profile!(soil, subsurface, kv_profile, dt)
    for i in eachindex(ssf)
        ssf[i] = kh[i] * (soilthickness[i] - zi[i]) * slope[i] * dw[i]
        kh_max = 0.0
        for j in 1:nlayers[i]
            if j <= nlayers_kv[i]
                kh_max += kv[i][j] * act_thickl[i][j]
            else
                zt = soil.parameters.soilthickness[i] - z_layered[i]
                k = max(j - 1, 1)
                kh_max += kv[i][k] / f[i] * (1.0 - exp(-f[i] * zt))
                break
            end
        end
        kh_max = kh_max * khfrac[i] * 0.001 * 0.001
        ssfmax[i] = kh_max * slope[i]
    end
    return nothing
end

"""
    bounded_divide(x, y; max = 1.0, default = 0.0)

Return the division of `x` by `y`, bounded by a maximum value `max`, when `y` > 0.0.
Otherwise return a `default` value.
"""
function bounded_divide(x, y; max = 1.0, default = 0.0)
    z = y > 0.0 ? min(x / y, max) : default
    return z
end
