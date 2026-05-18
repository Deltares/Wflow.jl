
"Map from PCRaster LDD value to a CartesianIndex"
const PCR_DIR = [
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

# timestep that the parameter units are defined in
const BASETIMESTEP = Second(Day(1))

"""
    scurve(x, a, b, c)

Sigmoid "S"-shaped curve.

# Arguments
- `x::Real`: input
- `a::Real`: determines the center level
- `b::Real`: determines the amplitude of the curve (range: (0, b⁻¹))
- `c::Real`: determines the steepness or "stepwiseness" of the curve.
             The higher c the sharper the function. A negative c reverses the function.
"""
function scurve(x::Real, a::Real, b::Real, c::Real)::Real
    s = inv(b + exp(-c * (x - a)))
    return s
end

function to_enumx(T, i::Int)
    options = instances(T)
    n_options = length(options)
    if 1 ≤ i ≤ n_options
        return options[i]
    else
        options_repr = repr(MIME("text/plain"), T)
        throw(
            error(
                "Cannot convert $i to $T, there are only $n_options options:\n$options_repr.",
            ),
        )
    end
end

"Set at indices pit values (default = 5) in a gridded local drainage direction vector"
function set_pit_ldd(
    pits_2d::AbstractMatrix{Bool},
    ldd::Vector{UInt8},
    indices::Vector{CartesianIndex{2}};
    pit::Integer = 5,
)::Vector{UInt8}
    pits = pits_2d[indices]
    index = filter(i -> isequal(pits[i], true), 1:length(indices))
    ldd[index] .= UInt8(pit)
    return ldd
end

"Filter upstream neighbors of graph based on logical vector"
function filter_upstream_nodes(
    graph::SimpleDiGraph{Int},
    vec_logical::Vector{Bool},
)::Vector{Vector{Int}}
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
function active_indices(
    subcatch_2d::AbstractMatrix,
    nodata,
)::Tuple{Vector{CartesianIndex{2}}, Matrix{Int}}
    A = subcatch_2d
    all_inds = CartesianIndices(size(A))
    indices = filter(i -> !isequal(A[i], nodata), all_inds)

    reverse_indices = zeros(Int, size(A))
    for (linear_idx, cartesian_idx) in enumerate(indices)
        reverse_indices[cartesian_idx] = linear_idx
    end

    return indices, reverse_indices
end

function active_indices(domain::Domain, key::AbstractString)::Vector{CartesianIndex{2}}
    if occursin("reservoir", key)
        return domain.reservoir.network.outlet_indices_2d
    elseif occursin("river", key) || occursin("floodplain", key)
        return domain.river.network.river_indices_2d
    elseif occursin("drain", key)
        return domain.drain.network.drain_indices_2d
    else
        return domain.land.network.land_indices_2d
    end
end

function lattometres(lat::Real)::Tuple{Float64, Float64}
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

function cell_lengths(
    y::AbstractVector{<:Real},
    celllength::Real,
    cell_length_in_meter::Bool,
)::Tuple{Vector{Float64}, Vector{Float64}}
    n_cells = length(y)
    xl = fill(MISSING_VALUE, n_cells)
    yl = fill(MISSING_VALUE, n_cells)
    if cell_length_in_meter
        xl .= celllength
        yl .= celllength
    else
        for cell_idx in 1:n_cells
            longlen, latlen = lattometres(y[cell_idx])
            xl[cell_idx] = longlen * celllength
            yl[cell_idx] = latlen * celllength
        end
    end
    return xl, yl
end

"""
    set_states!(instate_path, model, state_ncnames; <keyword arguments>)

Read states contained in `Dict` `state_ncnames` from netCDF file located in `instate_path`,
and set states in `model` object. Active cells are selected with the corresponding network's
(`Vector{CartesianIndex}`) from the netCDF file.

# Arguments
- `type = nothing`: type to convert data to after reading. By default no conversion is done.
"""
function set_states!(
    instate_path::AbstractString,
    model;
    type = nothing,
    dimname = nothing,
)::Nothing
    (; domain, land, config) = model

    # Check if required states are covered
    state_ncnames = check_states(config)

    # states in netCDF include dim time (one value) at index 3 or 4, 3 or 4 dims are allowed
    NCDataset(instate_path) do ds
        for (state, ncname) in state_ncnames
            @info "Setting initial state from netCDF." ncpath = instate_path ncvarname =
                ncname state
            sel = active_indices(domain, state)
            n_active_cells = length(sel)
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
                # note that this array is allowed to have missing, since not every land
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
                get_field_in_model(model, state)[1] .= svectorscopy(A, Val{size(A)[1]}())
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
                # set state in model object, only set active cells ([1:n_active_cells]) (ignore boundary conditions/ghost points)
                lens = get_metadata(state, typeof(land), Routing; model).lens
                lens(model)[1:n_active_cells] .= A
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

function get_var(config::Config, parameter::AbstractString; optional = true)
    if hasfield(InputSection, Symbol(parameter))
        var = getfield(config.input, Symbol(parameter))
    elseif haskey(config.input.location_maps, parameter)
        var = config.input.location_maps[parameter]
    elseif haskey(config.input.static, parameter)
        var = config.input.static[parameter]
    elseif haskey(config.input.cyclic, parameter)
        var = config.input.cyclic[parameter]
    elseif optional
        var = nothing
    else
        error(
            "Required input model parameter with standard name '$parameter' not set in TOML file",
        )
    end
    return var
end

"""
Apply the affine transform in `var` to the incoming array `A` in place element-wise.
The affine transform consists of a scaling by `scale` and a translation by `offset`.
These operations are only applied when non-trivial.
"""
function apply_affine_transform!(A::AbstractArray, var::InputEntry)
    (; do_scaling, scale_scalar, scale, do_offsetting, offset_scalar, offset) = var
    if do_scaling
        if scale_scalar
            A .*= only(scale)
        else
            A .*= scale
        end
    end
    if do_offsetting
        if offset_scalar
            A .+= only(offset)
        else
            A .+= offset
        end
    end
    return A
end

"""
    ncread(nc, config::Config, parameter::AbstractString, model_type; sel = nothing)

Read a netCDF variable `var` from file `nc`, based on `config` (parsed TOML file) and the
model `parameter` (standard name) specified in the TOML configuration file. Supports various
keyword arguments to get selections of data in desired types, with or without missing
values.

# Arguments
- `model_type`: The model type (e.g., LandHydrologySBM, SoilLoss, Domain, Routing) used to
        determine the appropriate standard name mapping.
- `sel=nothing`: A selection of indices, such as a `Vector{CartesianIndex}` of active cells,
        to return from the netCDF. By default all cells are returned.
- `logging=true`: Generate a logging message when reading a netCDF variable.
"""
function ncread(
    nc,
    config::Config,
    parameter::AbstractString,
    model_type;
    sel = nothing,
    logging = true,
    metadata = get_metadata(parameter, model_type),
)
    (; default, fill, type, allow_missing, dimname) = metadata
    var = get_var(config, parameter; optional = !isnothing(default))

    # for optional parameters default values are used.
    if isnothing(var)
        @info "Set `$parameter` using default value `$default`."
        @assert !isnothing(default) "Default value required but not available for $parameter (if you see this as a user please open an issue)."
        if isnothing(dimname)
            return Base.fill(default, length(sel))
        else
            return Base.fill(default, (nc.dim[String(dimname)], length(sel)))
        end
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

    var = convert(InputEntry, var)
    (; value, layer, scale, offset) = var
    variable_info(var)

    if !isnothing(value)
        @info "Set `$parameter` using uniform value `$value` from TOML file."
        if isnothing(dimname)
            # set to one uniform value
            return Base.fill(only(value), length(sel))
        elseif length(value) == 1
            # set to one uniform value (parameter with third dimension of size 1)
            return Base.fill(only(value), (nc.dim[String(dimname)], length(sel)))
        elseif length(value) > 1
            # set to multiple uniform values (parameter with third dimension of size > 1)
            @assert length(value) == nc.dim[String(dimname)]
            return repeat(value, 1, length(sel))
        end
    else
        if logging
            @info "Set `$parameter` using netCDF variable `$var`."
        end
        A = read_standardized(nc, variable_name(var), dim_sel)
        if !isnothing(layer)
            # the modifier index is only set in combination with scale and offset for SVectors,
            # provided through the TOML file.
            # if index, scale and offset is provided in the TOML as a list.
            for layer_idx in eachindex(layer)
                A[:, :, layer[layer_idx]] =
                    A[:, :, layer[layer_idx]] .* scale[layer_idx] .+ offset[layer_idx]
            end
        else
            apply_affine_transform!(A, var)
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

    if allow_missing
        # Convert to desired type if needed
        A = map(x -> ismissing(x) ? x : type(x), A)
    else
        if isnothing(fill)
            # errors if missing are found
            A = nomissing(A)
            if any(isnan, A)
                error("NaN not allowed in $var")
            end
        else
            # replaces missing with a fill value
            A = nomissing(A, fill)
            # replace also NaN values with the fill value
            replace!(x -> isnan(x) ? fill : x, A)
        end

        # Convert to desired type if needed
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
function set_layerthickness(
    reference_depth::Real,
    cum_depth::SVector,
    thickness::SVector,
)::SVector
    n_cells = length(thickness)
    thicknesslayers = thickness .* MISSING_VALUE
    for cell_idx in 1:n_cells
        if reference_depth > cum_depth[cell_idx + 1]
            thicknesslayers = setindex(thicknesslayers, thickness[cell_idx], cell_idx)
        elseif reference_depth - cum_depth[cell_idx] > 0.0
            thicknesslayers =
                setindex(thicknesslayers, reference_depth - cum_depth[cell_idx], cell_idx)
        end
    end
    return thicknesslayers
end

function number_of_active_layers(thickness::SVector)::Int
    nlayers = length(thickness) - sum(isnan.(thickness))
    return nlayers
end

"""
    get_flow_length(ldd, x_length, y_length)

Return the flow length for a non square grid. Input `ldd` (drainage network), `x_length`
(length of cells in x direction), `y_length` (length of cells in y direction). Output is
flow length.
"""
function get_flow_length(ldd::UInt8, x_length::Real, y_length::Real)::Real
    # take into account non-square cells
    # if ldd is 8 or 2 use y_length
    # if ldd is 4 or 6 use x_length
    if ldd == 2 || ldd == 8
        y_length
    elseif ldd == 4 || ldd == 6
        x_length
    else
        hypot(x_length, y_length)
    end
end

"""
    get_flow_width(ldd, x_length, y_length)

Return the flow width for a non square grid. Input `ldd` (drainage network), `x_length`
(length of cells in x direction), `y_length` (length of cells in y direction). Output is
flow width.
"""
function get_flow_width(ldd::UInt8, x_length::Real, y_length::Real)::Real
    # take into account non-square cells
    # if ldd is 8 or 2 use x_length
    # if ldd is 4 or 6 use y_length
    if ldd == 2 || ldd == 8
        x_length
    elseif ldd == 4 || ldd == 6
        y_length
    else
        (x_length * y_length) / hypot(x_length, y_length)
    end
end

"""
    get_surface_width(flow_width, flow_length, land_area, river_location)

Return the surface flow width. Input `flow_width` (flow width), `flow_length` (flow length),
`land_area` (area covered by land (excluding river coverage)) and `river_location` (river
cell, boolean). Output is surface flow width `surface_width`.
"""
function get_surface_width(
    flow_width::Real,
    flow_length::Real,
    land_area::Real,
    river_location::Bool,
)::Real
    surface_width = river_location ? land_area / flow_length : flow_width
    return surface_width
end

# 2.5x faster power method
"Faster method for exponentiation"
pow(x::Real, y::Real)::Real = exp(y * log(x))

function sum_at(A::AbstractVector{T}, inds::AbstractVector{Int})::T where {T}
    mapreduce(i -> A[i], +, inds; init = zero(T))
end

sum_at(f::Function, inds::AbstractVector{Int}; T::Type{<:Number} = Float64) =
    mapreduce(f, +, inds; init = zero(T))

# https://juliaarrays.github.io/StaticArrays.jl/latest/pages/api/#Arrays-of-static-arrays-1
function svectorscopy(x::Matrix{T}, ::Val{N})::Vector{SVector{N, T}} where {T, N}
    size(x, 1) == N || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    return copy(reinterpret(SVector{N, T}, vec(x)))
end

"""
    get_flow_fraction_to_river(graph, ldd, inds_river, slope)

Return flow `fraction` to a river cell (at `neighbor_idx`) based on the ratio of the land surface
`slope` at `neighbor_idx` to the sum of the land surface `slope` at `neighbor_idx` and at river cell
`river_cell_idx`.
"""
function get_flow_fraction_to_river(
    graph::SimpleDiGraph{Int},
    ldd::Vector{UInt8},
    inds_river::Vector{Int},
    slope::Vector{<:Real},
)::Vector{Float64}
    n_cells = length(slope)
    fraction = zeros(n_cells)
    for river_cell_idx in inds_river
        nbs = inneighbors(graph, river_cell_idx)
        for neighbor_idx in nbs
            if ldd[neighbor_idx] != ldd[river_cell_idx]
                fraction[neighbor_idx] =
                    slope[neighbor_idx] / (slope[river_cell_idx] + slope[neighbor_idx])
            end
        end
    end
    return fraction
end

"""
    equal_size_vectors(x)

Used in the structs of arrays to ensure all vectors are of equal length.

`equal_size_vectors(([1,2], [1,2,3]))` would throw an ArgumentError.
`equal_size_vectors(([4,5], [4,5]))` would pass.
`equal_size_vectors((1, [4,5], [4,5]))` would also pass, since `1` is not an AbstractVector.
"""
function equal_size_vectors(x::Tuple)
    # all vectors in this struct should be the same size
    inds_vec = findall(arg -> isa(arg, AbstractVector), x)
    expected_length = length(x[inds_vec[1]])
    x_vec = x[inds_vec]

    for arr in x_vec
        if length(arr) != expected_length
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
    adjacent_nodes_at_edge(graph)

Return the source node `src` and destination node `dst` of each edge of a directed `graph`.
"""
function adjacent_nodes_at_edge(
    graph::SimpleDiGraph{Int},
)::NamedTuple{(:src, :dst), Tuple{Vector{Int}, Vector{Int}}}
    _edges = collect(edges(graph))
    return (src = src.(_edges), dst = dst.(_edges))
end

"""
    adjacent_edges_at_node(graph, nodes_at_edge)

Return the source edge `src` and destination edge `dst` of each node of a directed `graph`.
"""
function adjacent_edges_at_node(
    graph::SimpleDiGraph{Int},
    nodes_at_edge,
)::NamedTuple{(:src, :dst), Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}}
    nodes = vertices(graph)
    src_edge = Vector{Int}[]
    dst_edge = copy(src_edge)
    for node_idx in 1:nv(graph)
        push!(src_edge, findall(isequal(nodes[node_idx]), nodes_at_edge.dst))
        push!(dst_edge, findall(isequal(nodes[node_idx]), nodes_at_edge.src))
    end
    return (src = src_edge, dst = dst_edge)
end

"Add `vertex` and `edge` to `pits` of a directed `graph`"
function add_vertex_edge_graph!(graph::SimpleDiGraph{Int}, pits::Vector{Int})::Nothing
    n_nodes = nv(graph)
    for (pit_idx, pit_node) in enumerate(pits)
        add_vertex!(graph)
        add_edge!(graph, pit_node, n_nodes + pit_idx)
    end
    return nothing
end

"""
    set_effective_flowwidth!(we_x::Vector{Float64}, we_y::Vector{Float64}, domain::Domain)

For river cells (D8 flow direction) in a staggered grid the effective flow width at cell
edges (floodplain) `we_x` in the x-direction and `we_y` in the y-direction is corrected by
subtracting the river width `flow_width` from the cell edges. For diagonal directions, the
`flow_width` is split between the two adjacent cell edges. A cell edge at linear index `idx`
is defined as the edge between node `idx` and the adjacent node (+ CartesianIndex(1, 0)) for
x and (+ CartesianIndex(0, 1)) for y. For cells that contain a `reservoir_outlet`
(reservoir), the effective flow width is set to zero.
"""
function set_effective_flowwidth!(
    we_x::Vector{Float64},
    we_y::Vector{Float64},
    domain::Domain,
)::Nothing
    (; local_drain_direction, river_indices_2d) = domain.river.network
    (; edge_indices, reverse_indices) = domain.land.network
    (; flow_width, reservoir_outlet) = domain.river.parameters
    reverse_indices = reverse_indices[river_indices_2d]

    graph = flowgraph(local_drain_direction, river_indices_2d, PCR_DIR)
    toposort = topological_sort_by_dfs(graph)
    n_cells = length(we_x)
    for cell_idx_2d in toposort
        dst = outneighbors(graph, cell_idx_2d)
        isempty(dst) && continue
        w = min(flow_width[cell_idx_2d], flow_width[only(dst)])
        dir = PCR_DIR[local_drain_direction[cell_idx_2d]]
        cell_idx = reverse_indices[cell_idx_2d]
        # loop over river D8 directions
        if dir == CartesianIndex(1, 1)
            we_x[cell_idx] =
                reservoir_outlet[cell_idx_2d] ? 0.0 : max(we_x[cell_idx] - 0.5 * w, 0.0)
            we_y[cell_idx] =
                reservoir_outlet[cell_idx_2d] ? 0.0 : max(we_y[cell_idx] - 0.5 * w, 0.0)
        elseif dir == CartesianIndex(-1, -1)
            if edge_indices.xd[cell_idx] <= n_cells
                we_y[edge_indices.xd[cell_idx]] =
                    reservoir_outlet[cell_idx_2d] ? 0.0 :
                    max(we_y[edge_indices.xd[cell_idx]] - 0.5 * w, 0.0)
            end
            if edge_indices.yd[cell_idx] <= n_cells
                we_x[edge_indices.yd[cell_idx]] =
                    reservoir_outlet[cell_idx_2d] ? 0.0 :
                    max(we_x[edge_indices.yd[cell_idx]] - 0.5 * w, 0.0)
            end
        elseif dir == CartesianIndex(1, 0)
            we_y[cell_idx] =
                reservoir_outlet[cell_idx_2d] ? 0.0 : max(we_y[cell_idx] - w, 0.0)
        elseif dir == CartesianIndex(0, 1)
            we_x[cell_idx] =
                reservoir_outlet[cell_idx_2d] ? 0.0 : max(we_x[cell_idx] - w, 0.0)
        elseif dir == CartesianIndex(-1, 0)
            if edge_indices.xd[cell_idx] <= n_cells
                we_y[edge_indices.xd[cell_idx]] =
                    reservoir_outlet[cell_idx_2d] ? 0.0 :
                    max(we_y[edge_indices.xd[cell_idx]] - w, 0.0)
            end
        elseif dir == CartesianIndex(0, -1)
            if edge_indices.yd[cell_idx] <= n_cells
                we_x[edge_indices.yd[cell_idx]] =
                    reservoir_outlet[cell_idx_2d] ? 0.0 :
                    max(we_x[edge_indices.yd[cell_idx]] - w, 0.0)
            end
        elseif dir == CartesianIndex(1, -1)
            we_y[cell_idx] = max(we_y[cell_idx] - 0.5 * w, 0.0)
            if edge_indices.yd[cell_idx] <= n_cells
                we_x[edge_indices.yd[cell_idx]] =
                    reservoir_outlet[cell_idx_2d] ? 0.0 :
                    max(we_x[edge_indices.yd[cell_idx]] - 0.5 * w, 0.0)
            end
        elseif dir == CartesianIndex(-1, 1)
            if edge_indices.xd[cell_idx] <= n_cells
                we_y[edge_indices.xd[cell_idx]] =
                    reservoir_outlet[cell_idx_2d] ? 0.0 :
                    max(we_y[edge_indices.xd[cell_idx]] - 0.5 * w, 0.0)
            end
            we_x[cell_idx] =
                reservoir_outlet[cell_idx_2d] ? 0.0 : max(we_x[cell_idx] - 0.5 * w, 0.0)
        end
    end
    return nothing
end

"Return julian day of year (leap days are not counted)"
function julian_day(time::TimeType)::Int
    # for all years February 28 is day 59 and March 1 is day 60.
    day = dayofyear(time) - (isleapyear(time) && dayofyear(time) > 60)
    return day
end

"Partition indices with at least size `basesize`"
function _partition(xs::Integer, basesize::Integer)
    n_partitions = Int(max(1, xs ÷ basesize))
    return (
        Int(1 + ((partition_idx - 1) * xs) ÷ n_partitions):Int(
            (partition_idx * xs) ÷ n_partitions,
        ) for partition_idx in 1:n_partitions
    )
end

"""
    threaded_foreach(f, x::AbstractArray; basesize::Integer)

Run function `f` in parallel by spawning tasks (nthreads <= 8), each task iterates over a
chunk of size `basesize`. For nthreads > 8 run function `f` in parallel with
`Polyester@batch` with `minbatch` equal to `basesize`.
"""
function threaded_foreach(f::Function, x::AbstractArray; basesize::Integer)::Nothing
    if Threads.nthreads() <= 8
        len = length(x)
        partitions = _partition(len, basesize)
        if length(partitions) > 1 && Threads.nthreads() > 1
            @sync for p in partitions
                Threads.@spawn begin
                    for element_idx in eachindex(p)
                        f(@inbounds p[element_idx])
                    end
                end
            end
        else
            for element_idx in eachindex(x)
                f(@inbounds x[element_idx])
            end
        end
    else
        @batch per = thread minbatch = basesize for element_idx in eachindex(x)
            f(@inbounds x[element_idx])
        end
    end
    return nothing
end

"""
    hydraulic_conductivity_at_depth(p::KvExponential, kvfrac, z, cell_idx, soil_layer_idx)
    hydraulic_conductivity_at_depth(p::KvExponentialConstant, kvfrac, z, cell_idx, soil_layer_idx)
    hydraulic_conductivity_at_depth(p::KvLayered, kvfrac, z, cell_idx, soil_layer_idx)
    hydraulic_conductivity_at_depth(p::KvLayeredExponential, kvfrac, z, cell_idx, soil_layer_idx)

Return vertical hydraulic conductivity `kv_z` at depth `z` at `cell_idx` using multiplication
factor `kv_frac` at `soil_layer_idx` and vertical hydraulic conductivity profile `p`.
"""
function hydraulic_conductivity_at_depth(
    p::KvExponential,
    kvfrac,
    z,
    cell_idx,
    soil_layer_idx,
)
    kv_z = kvfrac[cell_idx][soil_layer_idx] * p.kv_0[cell_idx] * exp(-p.f[cell_idx] * z)
    return kv_z
end

function hydraulic_conductivity_at_depth(
    p::KvExponentialConstant,
    kvfrac,
    z,
    cell_idx,
    soil_layer_idx,
)
    (; kv_0, f) = p.exponential
    if z < p.z_exp[cell_idx]
        kv_z = kvfrac[cell_idx][soil_layer_idx] * kv_0[cell_idx] * exp(-f[cell_idx] * z)
    else
        kv_z =
            kvfrac[cell_idx][soil_layer_idx] *
            kv_0[cell_idx] *
            exp(-f[cell_idx] * p.z_exp[cell_idx])
    end
    return kv_z
end

function hydraulic_conductivity_at_depth(p::KvLayered, kvfrac, z, cell_idx, soil_layer_idx)
    kv_z = kvfrac[cell_idx][soil_layer_idx] * p.kv[cell_idx][soil_layer_idx]
    return kv_z
end

function hydraulic_conductivity_at_depth(
    p::KvLayeredExponential,
    kvfrac,
    z,
    cell_idx,
    soil_layer_idx,
)
    return if z < p.z_layered[cell_idx]
        kvfrac[cell_idx][soil_layer_idx] * p.kv[cell_idx][soil_layer_idx]
    else
        soil_layer_idx = p.nlayers_kv[cell_idx]
        kvfrac[cell_idx][soil_layer_idx] *
        p.kv[cell_idx][soil_layer_idx] *
        exp(-p.f[cell_idx] * (z - p.z_layered[cell_idx]))
    end
end

"""
    kh_layered_profile!(soil_model::SbmSoilModel, subsurface_flow_model::LateralSSFModel, kv_profile::KvLayered, dt)
    kh_layered_profile!(soil_model::SbmSoilModel, subsurface_flow_model::LateralSSFModel, kv_profile::KvLayeredExponential, dt)

Compute equivalent horizontal hydraulic conductivity `kh` [m d⁻¹] using vertical hydraulic
conductivity profile `kv_profile`.
"""
function kh_layered_profile!(
    soil_model::SbmSoilModel,
    subsurface_flow_model::LateralSSFModel,
    kv_profile::KvLayered,
    dt,
)
    (; n_cells) = soil_model
    (; nlayers, sumlayers, act_thickl, soilthickness) = soil_model.parameters
    (; n_unsatlayers, zi) = soil_model.variables
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; khfrac) = subsurface_flow_model.parameters

    t_factor = (tosecond(BASETIMESTEP) / dt)
    for cell_idx in 1:n_cells
        n_soil_layers = nlayers[cell_idx]

        if soilthickness[cell_idx] > zi[cell_idx]
            transmissivity = 0.0
            _sumlayers = @view sumlayers[cell_idx][2:end]
            soil_layer_idx = max(n_unsatlayers[cell_idx], 1)
            transmissivity +=
                (_sumlayers[soil_layer_idx] - zi[cell_idx]) *
                kv_profile.kv[cell_idx][soil_layer_idx]
            soil_layer_idx += 1
            while soil_layer_idx <= n_soil_layers
                transmissivity +=
                    act_thickl[cell_idx][soil_layer_idx] *
                    kv_profile.kv[cell_idx][soil_layer_idx]
                soil_layer_idx += 1
            end
            # convert units for kh [m d⁻¹] computation (transmissivity [mm² Δt⁻¹], soilthickness
            # [mm] and zi [mm])
            kh[cell_idx] =
                0.001 *
                (transmissivity / (soilthickness[cell_idx] - zi[cell_idx])) *
                t_factor *
                khfrac[cell_idx]
        else
            kh[cell_idx] =
                0.001 * kv_profile.kv[cell_idx][n_soil_layers] * t_factor * khfrac[cell_idx]
        end
    end
    return nothing
end

function kh_layered_profile!(
    soil_model::SbmSoilModel,
    subsurface_flow_model::LateralSSFModel,
    kv_profile::KvLayeredExponential,
    dt,
)
    (; n_cells) = soil_model
    (; nlayers, sumlayers, act_thickl, soilthickness) = soil_model.parameters
    (; nlayers_kv, z_layered, kv, f) = kv_profile
    (; n_unsatlayers, zi) = soil_model.variables
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; khfrac) = subsurface_flow_model.parameters
    t_factor = (tosecond(BASETIMESTEP) / dt)

    for cell_idx in 1:n_cells
        n_soil_layers = nlayers[cell_idx]

        if soilthickness[cell_idx] > zi[cell_idx]
            transmissivity = 0.0
            soil_layer_idx = max(n_unsatlayers[cell_idx], 1)
            if zi[cell_idx] >= z_layered[cell_idx]
                zt = soilthickness[cell_idx] - z_layered[cell_idx]
                soil_layer_idx = nlayers_kv[cell_idx]
                transmissivity +=
                    kv[cell_idx][soil_layer_idx] / f[cell_idx] * (
                        exp(-f[cell_idx] * (zi[cell_idx] - z_layered[cell_idx])) -
                        exp(-f[cell_idx] * zt)
                    )
                soil_layer_idx = n_soil_layers
            else
                _sumlayers = @view sumlayers[cell_idx][2:end]
                transmissivity +=
                    (_sumlayers[soil_layer_idx] - zi[cell_idx]) *
                    kv[cell_idx][soil_layer_idx]
            end
            soil_layer_idx += 1
            while soil_layer_idx <= n_soil_layers
                if soil_layer_idx > nlayers_kv[cell_idx]
                    zt = soilthickness[cell_idx] - z_layered[cell_idx]
                    j = nlayers_kv[cell_idx]
                    transmissivity +=
                        kv[cell_idx][j] / f[cell_idx] * (1.0 - exp(-f[cell_idx] * zt))
                    soil_layer_idx = n_soil_layers
                else
                    transmissivity +=
                        act_thickl[cell_idx][soil_layer_idx] * kv[cell_idx][soil_layer_idx]
                end
                soil_layer_idx += 1
            end
            # convert units for kh [m d⁻¹] computation (transmissivity [mm² Δt⁻¹], soilthickness
            # [mm] and zi [mm])
            kh[cell_idx] =
                0.001 *
                (transmissivity / (soilthickness[cell_idx] - zi[cell_idx])) *
                t_factor *
                khfrac[cell_idx]
        else
            if zi[cell_idx] >= z_layered[cell_idx]
                j = nlayers_kv[cell_idx]
                kh[cell_idx] =
                    0.001 *
                    kv[cell_idx][j] *
                    exp(-f[cell_idx] * (zi[cell_idx] - z_layered[cell_idx])) *
                    khfrac[cell_idx] *
                    t_factor
            else
                kh[cell_idx] =
                    0.001 * kv[cell_idx][n_soil_layers] * t_factor * khfrac[cell_idx]
            end
        end
    end
    return nothing
end

kh_layered_profile!(
    soil_model::SbmSoilModel,
    subsurface_flow_model::LateralSSFModel,
    kv_profile::Union{KvExponential, KvExponentialConstant},
    dt,
) = nothing

"""
    initialize_lateral_ssf_model!(subsurface_flow_model::LateralSSFModel, parameters::LandParameters, kh_profile::KhExponential)
    initialize_lateral_ssf_model!(subsurface_flow_model::LateralSSFModel, parameters::LandParameters, kh_profile::KhExponentialConstant)

Initialize lateral subsurface variables `q` and `q_max` using horizontal hydraulic
conductivity profile `kh_profile`.
"""
function initialize_lateral_ssf_model!(
    subsurface_flow_model::LateralSSFModel,
    parameters::LandParameters,
    kh_profile::KhExponential,
)
    (; kh_0, f) = kh_profile
    (; q, q_max, zi) = subsurface_flow_model.variables
    (; soilthickness) = subsurface_flow_model.parameters
    (; slope, flow_width) = parameters

    @. q_max = ((kh_0 * slope) / f) * (1.0 - exp(-f * soilthickness))
    @. q = ((kh_0 * slope) / f) * (exp(-f * zi) - exp(-f * soilthickness)) * flow_width
    return nothing
end

function initialize_lateral_ssf_model!(
    subsurface_flow_model::LateralSSFModel,
    parameters::LandParameters,
    kh_profile::KhExponentialConstant,
)
    (; kh_0, f) = kh_profile.exponential
    (; z_exp) = kh_profile
    (; q, q_max, zi, n_cells) = subsurface_flow_model.variables
    (; soilthickness) = subsurface_flow_model.parameters
    (; slope, flow_width) = parameters

    q_constant = @. kh_0 * exp(-f * z_exp) * slope * (soilthickness - z_exp)
    for cell_idx in 1:n_cells
        q_max[cell_idx] =
            ((kh_0[cell_idx] * slope[cell_idx]) / f[cell_idx]) *
            (1.0 - exp(-f[cell_idx] * z_exp[cell_idx])) + q_constant[cell_idx]
        if zi[cell_idx] < z_exp[cell_idx]
            q[cell_idx] =
                (
                    ((kh_0[cell_idx] * slope[cell_idx]) / f[cell_idx]) * (
                        exp(-f[cell_idx] * zi[cell_idx]) -
                        exp(-f[cell_idx] * z_exp[cell_idx])
                    ) + q_constant[cell_idx]
                ) * flow_width[cell_idx]
        else
            q[cell_idx] =
                kh_0[cell_idx] *
                exp(-f[cell_idx] * zi[cell_idx]) *
                slope[cell_idx] *
                (soilthickness[cell_idx] - zi[cell_idx]) *
                flow_width[cell_idx]
        end
    end
    return nothing
end

"""
    initialize_lateral_ssf_model!(subsurface_flow_model::LateralSSFModel, soil_model::SbmSoilModel, parameters::LandParameters, kv_profile::KvLayered, dt)
    initialize_lateral_ssf_model!(subsurface_flow_model::LateralSSFModel, soil_model::SbmSoilModel, parameters::LandParameters, kv_profile::KvLayeredExponential, dt)

Initialize lateral subsurface variables `q` and `q_max` using  vertical hydraulic
conductivity profile `kv_profile`.
"""
function initialize_lateral_ssf_model!(
    subsurface_flow_model::LateralSSFModel,
    soil_model::SbmSoilModel,
    parameters::LandParameters,
    kv_profile::KvLayered,
    dt,
)
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; nlayers, act_thickl) = soil_model.parameters
    (; q, q_max, zi, n_cells) = subsurface_flow_model.variables
    (; khfrac, soilthickness) = subsurface_flow_model.parameters
    (; slope, flow_width) = parameters

    kh_layered_profile!(soil_model, subsurface_flow_model, kv_profile, dt)
    for cell_idx in 1:n_cells
        q[cell_idx] =
            kh[cell_idx] *
            (soilthickness[cell_idx] - zi[cell_idx]) *
            slope[cell_idx] *
            flow_width[cell_idx]
        kh_max = 0.0
        for soil_layer_idx in 1:nlayers[cell_idx]
            kh_max +=
                kv_profile.kv[cell_idx][soil_layer_idx] *
                act_thickl[cell_idx][soil_layer_idx]
        end
        kh_max *= khfrac[cell_idx] * 0.001 * 0.001
        q_max[cell_idx] = kh_max * slope[cell_idx]
    end
    return nothing
end

function initialize_lateral_ssf_model!(
    subsurface_flow_model::LateralSSFModel,
    soil_model::SbmSoilModel,
    parameters::LandParameters,
    kv_profile::KvLayeredExponential,
    dt,
)
    (; q, q_max, zi, n_cells) = subsurface_flow_model.variables
    (; khfrac, soilthickness) = subsurface_flow_model.parameters
    (; slope, flow_width) = parameters
    (; nlayers, act_thickl) = soil_model.parameters
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; kv, f, nlayers_kv, z_layered) = kv_profile

    kh_layered_profile!(soil_model, subsurface_flow_model, kv_profile, dt)
    for cell_idx in 1:n_cells
        q[cell_idx] =
            kh[cell_idx] *
            (soilthickness[cell_idx] - zi[cell_idx]) *
            slope[cell_idx] *
            flow_width[cell_idx]
        kh_max = 0.0
        for soil_layer_idx in 1:nlayers[cell_idx]
            if soil_layer_idx <= nlayers_kv[cell_idx]
                kh_max +=
                    kv[cell_idx][soil_layer_idx] * act_thickl[cell_idx][soil_layer_idx]
            else
                zt = soil_model.parameters.soilthickness[cell_idx] - z_layered[cell_idx]
                kh_max +=
                    kv[cell_idx][max(soil_layer_idx - 1, 1)] / f[cell_idx] *
                    (1.0 - exp(-f[cell_idx] * zt))
                break
            end
        end
        kh_max = kh_max * khfrac[cell_idx] * 0.001 * 0.001
        q_max[cell_idx] = kh_max * slope[cell_idx]
    end
    return nothing
end

"""
    bounded_divide(x, y; max = 1.0, default = 0.0)

Return the division of `x` by `y`, bounded by a maximum value `max`, when `y` > 0.0.
Otherwise return a `default` value.
"""
function bounded_divide(x::Real, y::Real; max::Real = 1.0, default::Real = 0.0)::Real
    z = y > 0.0 ? min(x / y, max) : default
    return z
end

"""
The sine of the slope in radians;
sin(arctan(x)) = x / √(1 + x²)
"""
sin_slope(slope) = slope / sqrt(1 + slope^2)

"""
Return water table change `dh` and exfiltration rate `exfilt`. For a falling water table
`dh` is based on subsurface net flux `net_flux` and specific yield `specific_yield`. For a
rising water table `dh` is based on `net_flux` and the unsaturated store capacity (per soil
layer). For a rising water table a dynamic specific yield is computed.
"""
function water_table_change(
    soil_model::SbmSoilModel,
    net_flux::Float64,
    specific_yield::Float64,
    cell_idx::Int,
)
    (; n_unsatlayers, ustorelayerthickness, ustorelayerdepth) = soil_model.variables
    (; theta_s, theta_r) = soil_model.parameters

    # effective porosity (difference between saturated and residual water content)
    theta_e = theta_s[cell_idx] - theta_r[cell_idx]

    if net_flux <= 0.0
        dh = net_flux / specific_yield
    else
        dh = 0.0
        f_conv = 0.001 # convert units from [mm] to [m]
        for soil_layer_idx in n_unsatlayers[cell_idx]:-1:1
            capacity = max(
                f_conv * (
                    ustorelayerthickness[cell_idx][soil_layer_idx] * theta_e -
                    ustorelayerdepth[cell_idx][soil_layer_idx]
                ),
                0.0,
            )
            flux_layer = min(net_flux, capacity)
            if capacity <= net_flux
                # if unsaturated layer is fully saturated dh equals layer thickness
                dh += f_conv * ustorelayerthickness[cell_idx][soil_layer_idx]
            else
                sy =
                    theta_e - (
                        ustorelayerdepth[cell_idx][soil_layer_idx] /
                        ustorelayerthickness[cell_idx][soil_layer_idx]
                    )
                dh += flux_layer / sy
            end
            net_flux -= flux_layer
            net_flux == 0.0 && break
        end
    end
    exfilt = max(net_flux, 0.0)
    return dh, exfilt
end

"Set lower bound for drainable porosity"
function lower_bound_drainable_porosity(theta_s, theta_fc; lower_bound = 0.02)
    return max(theta_s - theta_fc, lower_bound)
end
