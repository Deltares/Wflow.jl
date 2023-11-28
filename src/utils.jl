
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
        for i = 1:n
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
    for i = 1:n
        riverfrac[i] = if river[i]
            min((riverlength[i] * riverwidth[i]) / (xl[i] * yl[i]), 1.0)
        else
            0.0
        end
    end
    return riverfrac
end

"""
    set_states(instate_path, model, state_ncnames; <keyword arguments>)

Read states contained in `Dict` `state_ncnames` from NetCDF file located in `instate_path`,
and set states in `model` object. Active cells are selected with the corresponding network's
(`Vector{CartesianIndex}`) from the NetCDF file.

# Arguments
- `type = nothing`: type to convert data to after reading. By default no conversion is done.
"""
function set_states(instate_path, model, state_ncnames; type = nothing, dimname = nothing)
    @unpack network = model
    # states in NetCDF include dim time (one value) at index 3 or 4, 3 or 4 dims are allowed
    NCDataset(instate_path) do ds
        for (state, ncname) in state_ncnames
            @info "Setting initial state from NetCDF." ncpath = instate_path ncvarname =
                ncname state
            sel = active_indices(network, state)
            n = length(sel)
            dims = length(dimnames(ds[ncname]))
            # 4 dims, for example (x,y,layer,time) where dim layer is an SVector for soil layers
            if dims == 4
                if dimname == :layer
                    dimensions = (x = :, y = :, layer = :, time = 1)
                elseif dimname == :classes
                    dimensions = (x = :, y = :, classes = :, time = 1)
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
end

"""
    ncread(nc, config::Config, parameter::AbstractString; <keyword arguments>)

Read a NetCDF variable `var` from file `nc`, based on `config` (parsed TOML file) and the
model `parameter` specified in the TOML configuration file. Supports various keyword
arguments to get selections of data in desired types, with or without missing values.

# Arguments
- `alias`=nothing` : An `alias` for the TOML key, by default an `alias` is not expected.
- `optional=true` : By default specifying a model `parameter` in the TOML file is optional.
        Set to false if the model `parameter` is required.
- `sel=nothing`: A selection of indices, such as a `Vector{CartesianIndex}` of active cells,
        to return from the NetCDF. By default all cells are returned.
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
    # get var (NetCDF variable or type Config) from TOML file.
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
    elseif dimname == :classes
        dim_sel = (x = :, y = :, classes = :, time = 1)
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
    # input NetCDF var) or set to a uniform value (providing a value). Otherwise, input
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
        @info "Set `$parameter` using NetCDF variable `$var`."
        A = read_standardized(nc, var, dim_sel)
        if !isnothing(mod.index)
            # the modifier index is only set in combination with scale and offset for SVectors,
            # provided through the TOML file.
            if length(mod.index) > 1
                # if index, scale and offset is provided in the TOML as a list.
                for i = 1:length(mod.index)
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
    set_layerthickness(d::Real, sl::SVector)

Calculate actual soil thickness of layers based on a reference depth (e.g. soil depth or water table depth) `d`,
a SVector `sl` with cumulative soil depth starting at soil surface (0), and a SVector `tl` with actual thickness
per soil layer.
"""
function set_layerthickness(d::Real, sl::SVector, tl::SVector)

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
function svectorscopy(x::Matrix{T}, ::Val{N}) where {T,N}
    size(x, 1) == N || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    copy(reinterpret(SVector{N,T}, vec(x)))
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
tosecond(x::T) where {T<:DatePeriod} = Float64(Dates.value(Second(x)))
tosecond(x::T) where {T<:TimePeriod} = x / convert(T, Second(1))

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
    for i = 1:nv(graph)
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
    return (Int(1 + ((i - 1) * xs) ÷ n):Int((i * xs) ÷ n) for i = 1:n)
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
    hydraulic_conductivity_at_depth(sbm::SBM, z, i, ksat_profile)

Return vertical hydraulic conductivity `kv_z` for soil layer `n` at depth `z` for vertical
concept `SBM` (at index `i`) based on hydraulic conductivity profile `ksat_profile`.
"""
function hydraulic_conductivity_at_depth(sbm::SBM, z, i, n, ksat_profile)
    if ksat_profile == "exponential"
        kv_z = sbm.kvfrac[i][n] * sbm.kv₀[i] * exp(-sbm.f[i] * z)
    elseif ksat_profile == "exponential_constant"
        if sbm.zi[i] < sbm.z_exp[i]
            kv_z = sbm.kvfrac[i][n] * sbm.kv₀[i] * exp(-sbm.f[i] * z)
        else
            kv_z = sbm.kvfrac[i][n] * sbm.kv₀[i] * exp(-sbm.f[i] * sbm.z_exp[i])
        end
    elseif ksat_profile == "layered"
        kv_z = sbm.kvfrac[i][n] * sbm.kv[i][n]
    elseif ksatprofile == "layered_exponential"
        if sbm.zi[i] < sbm.z_exp[i]
            kv_z = sbm.kvfrac[i][n] * sbm.kv[i][n]
        else
            n = sbm.nlayers_kv[i]
            kv_z = sbm.kvfrac[i][n] * sbm.kv[i][n] * exp(-sbm.f[i] * (z - sbm.z_exp[i]))
        end
    end
    return kv_z
end
