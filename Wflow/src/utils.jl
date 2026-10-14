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
    for (i, I) in enumerate(indices)
        reverse_indices[I] = i
    end

    return indices, reverse_indices
end

"Get active indices of `key` (standard name or model path) by prefix string matching."
function active_indices(domain::Domain, key::AbstractString)::Vector{CartesianIndex{2}}
    if startswith(key, domain_parameter_map["reservoir"])
        return domain.reservoir.network.indices_outlet
    elseif startswith(key, domain_parameter_map["river"])
        return domain.river.network.indices
    elseif startswith(key, domain_parameter_map["drain"])
        return domain.drain.network.indices
    else
        return domain.land.network.indices
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
    n = length(y)
    xl = fill(MISSING_VALUE, n)
    yl = fill(MISSING_VALUE, n)
    if cell_length_in_meter
        xl .= celllength
        yl .= celllength
    else
        for i in 1:n
            longlen, latlen = lattometres(y[i])
            xl[i] = longlen * celllength
            yl[i] = latlen * celllength
        end
    end
    return xl, yl
end

"Convert an array with missing values to an array without, erroring if any are missing."
function nomissing(A::AbstractArray{Union{Missing, T}, N})::Array{T, N} where {T, N}
    any(ismissing, A) && error("Arrays with missing values are not supported.")
    return convert(Array{T, N}, A)
end

nomissing(A::AbstractArray{T}) where {T <: Number} = A

"Replace missing values in an array with `fill`."
nomissing(A::AbstractArray{Union{Missing, T}, N}, fill::T) where {T, N} =
    Array{T, N}(coalesce.(A, fill))

"""
    set_states!(instate_path, model, state_ncnames; <keyword arguments>)

Read states contained in `Dict` `state_ncnames` from netCDF file located in `instate_path`,
and set states in `model` object. Active cells are selected with the corresponding network's
(`Vector{CartesianIndex}`) from the netCDF file.

# Arguments
- `type = nothing`: type to convert data to after reading. By default no conversion is done.
"""
function set_states!(instate_path::AbstractString, model; dimname = nothing)::Nothing
    (; domain, config, clock) = model
    dt_val = tosecond(clock.dt)

    # Check if required states are covered
    state_ncnames = check_states(config)

    # states in netCDF include dim time (one value) at index 3 or 4, 3 or 4 dims are allowed
    ds = NC.Dataset(instate_path)
    try
        for (state, ncname) in state_ncnames
            sel = active_indices(domain, state)
            dims = length(NC.var_dimnames(ds, ncname))
            with_standard_name(model, state) do field, metadata
                (; unit) = metadata
                @info "Setting initial state from netCDF." ncpath = instate_path ncvarname =
                    ncname state unit
                set_state!(field, ds, ncname, sel, dims, dimname, unit, dt_val)
            end
        end
    finally
        close(ds)
    end
    return nothing
end

"Set a layered state (x, y, layer, time) in a vector of SVectors."
function set_state!(
        field::Vector{SVector{N, Float64}},
        ds::NC.Dataset,
        ncname::String,
        sel::Vector{CartesianIndex{2}},
        dims::Int,
        dimname,
        unit::Unit,
        dt_val::Float64,
    )::Nothing where {N}
    # 4 dims, for example (x,y,layer,time) where dim layer is an SVector for soil layers
    dims == 4 || error("Number of dims for layered state $ncname should be 4, got $dims")
    dimname == :layer || error("Unrecognized dimension name $dimname")
    A = read_standardized(ds, ncname, (x = :, y = :, layer = :, time = 1))
    # note that this array is allowed to have missing, since not every land
    # column is `maximum_number_of_layers` layers deep
    B = Matrix{Float64}(coalesce.(permutedims(A[sel, :]), NaN))
    to_SI!(B, unit; dt_val)
    # set state in model object
    field .= svectorscopy(B, Val(N))
    return nothing
end

"Set a state (x, y, time) in a vector, only for the active cells."
function set_state!(
        field::AbstractVector{Float64},
        ds::NC.Dataset,
        ncname::String,
        sel::Vector{CartesianIndex{2}},
        dims::Int,
        dimname,
        unit::Unit,
        dt_val::Float64,
    )::Nothing
    dims == 3 || error("Number of dims for state $ncname should be 3, got $dims")
    A = read_standardized(ds, ncname, (x = :, y = :, time = 1))
    B = nomissing(A[sel])
    to_SI!(B, unit; dt_val)
    # set state in model object, only set active cells ([1:n]) (ignore boundary conditions/ghost points)
    n = length(sel)
    field[1:n] .= B
    return nothing
end

set_state!(field, ds, ncname, sel, dims, dimname, unit, dt_val) =
    error("Unsupported state type $(typeof(field)) for $ncname")

function get_var(config::Config, parameter::AbstractString; optional = true)
    if hasfield(InputSection, Symbol(parameter))
        var = getfield(config.input, Symbol(parameter))
    elseif haskey(config.input._location_maps, parameter)
        var = config.input._location_maps[parameter]
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
    (; _do_scaling, _scale_scalar, scale, _do_offsetting, _offset_scalar, offset) = var
    if _do_scaling
        if _scale_scalar
            A .*= only(scale)
        else
            A .*= scale
        end
    end
    if _do_offsetting
        if _offset_scalar
            A .+= only(offset)
        else
            A .+= offset
        end
    end
    return A
end

"""
    ncread(nc, config::Config, parameter::AbstractString, model_type; <keyword arguments>)

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
- `metadata`: The metadata of the read parameter or state, obtained from the parameter name by default.
"""
Base.@constprop :aggressive function ncread(
        nc::NC.Dataset,
        config::Config,
        parameter::AbstractString,
        model_type::Type{M};
        sel = nothing,
        logging::Bool = true,
        metadata::Union{Nothing, ParameterMetadata} = nothing,
    ) where {M}
    # The parameter name is a literal at nearly all call sites, so with aggressive constant
    # propagation and a foldable metadata lookup, the metadata is a compile time constant.
    # This makes the return type concrete without annotations at the call sites.
    md = isnothing(metadata) ? get_metadata_foldable(parameter, model_type) : metadata
    return ncread_dispatch(md, nc, config, parameter, sel, logging)
end

"""
Call `ncread_typed` with the element type, extra dimension and missing value settings of
the metadata as compile time constants. For constant metadata this folds to a single call,
otherwise it is a static branch over all supported combinations.
"""
Base.@constprop :aggressive function ncread_dispatch(
        @nospecialize(md::ParameterMetadata),
        nc::NC.Dataset,
        config::Config,
        parameter::AbstractString,
        sel,
        logging::Bool,
    )
    type = md.type
    no_extra_dim = isnothing(md.dimname)
    allow_missing = md.allow_missing
    return if type === Int
        ncread_dispatch(Int, no_extra_dim, allow_missing, nc, config, parameter, md, sel, logging)
    elseif type === Bool
        ncread_dispatch(Bool, no_extra_dim, allow_missing, nc, config, parameter, md, sel, logging)
    else
        ncread_dispatch(
            Float64,
            no_extra_dim,
            allow_missing,
            nc,
            config,
            parameter,
            md,
            sel,
            logging,
        )
    end
end

Base.@constprop :aggressive function ncread_dispatch(
        ::Type{T},
        no_extra_dim::Bool,
        allow_missing::Bool,
        nc::NC.Dataset,
        config::Config,
        parameter::AbstractString,
        @nospecialize(md::ParameterMetadata),
        sel,
        logging::Bool,
    ) where {T}
    return if no_extra_dim
        if allow_missing
            ncread_typed(T, Val(true), Val(true), nc, config, parameter, md, sel, logging)
        else
            ncread_typed(T, Val(true), Val(false), nc, config, parameter, md, sel, logging)
        end
    else
        if allow_missing
            ncread_typed(T, Val(false), Val(true), nc, config, parameter, md, sel, logging)
        else
            ncread_typed(T, Val(false), Val(false), nc, config, parameter, md, sel, logging)
        end
    end
end

"""
    get_metadata_foldable(name, model_type)

Metadata lookup that may be evaluated at compile time for a constant `name`. The standard
name maps are constant after module initialization, so the lookup is consistent.
"""
Base.@assume_effects :foldable get_metadata_foldable(
    name::AbstractString,
    model_type::Type{M},
) where {M} = something(get_metadata(name, model_type))::ParameterMetadata

"Convert a parameter value to type `T`, used for defaults and fill values."
convert_parameter(::Type{T}, x::Union{Float64, Int, Bool}) where {T} = T(x)

"Uniform value array, for a parameter without (`V = true`) or with an extra dimension."
uniform_array(value::T, ::Val{true}, nc, dimname, sel) where {T} = Base.fill(value, length(sel))
uniform_array(value::T, ::Val{false}, nc, dimname, sel) where {T} =
    Base.fill(value, (NC.dimlength(nc, String(dimname)), length(sel)))

"Convert an array to allow missing values if `M = true`."
missing_eltype(A::Array{T, N}, ::Val{true}) where {T, N} = Array{Union{Missing, T}, N}(A)
missing_eltype(A::Array, ::Val{false}) = A

function ncread_typed(
        ::Type{T},
        no_extra_dim::Val{V},
        allow_missing::Val{M},
        nc::NC.Dataset,
        config::Config,
        parameter::AbstractString,
        @nospecialize(metadata::ParameterMetadata),
        sel,
        logging::Bool,
    ) where {T, V, M}
    (; dimname, unit) = metadata
    default = metadata.default::Union{Nothing, Float64, Int, Bool}
    fill = metadata.fill::Union{Nothing, Float64, Int, Bool}
    var = get_var(config, parameter; optional = !isnothing(default))
    # the timestep is set from the forcing if not given, before reading parameters
    dt_val = something(config.time.timestepsecs)

    # for optional parameters default values are used.
    if isnothing(var)
        @info "Set `$parameter [$unit]` using default value `$default $unit`."
        @assert !isnothing(default) "Default value required but not available for $parameter (if you see this as a user please open an issue)."
        default_SI = to_SI(convert_parameter(T, default), unit; dt_val)
        A = uniform_array(default_SI, no_extra_dim, nc, dimname, sel)
        return missing_eltype(A, allow_missing)
    end

    if var isa Float64 || var isa Int
        var = InputEntry(; value = [Float64(var)])
    elseif var isa String
        var = InputEntry(; external_name = var)
    else
        var = var::InputEntry
    end

    (; value, layer, scale, offset) = var
    variable_info(var)

    if !isnothing(value)
        @info "Set `$parameter [$unit]` using uniform value `$(join(value, ", ")) $unit` from TOML file."
        values_T = T[convert_parameter(T, x) for x in value]
        A = if V || length(values_T) == 1
            # set to one uniform value (also for a parameter with third dimension of size 1)
            uniform_array(only(values_T), no_extra_dim, nc, dimname, sel)
        else
            # set to multiple uniform values (parameter with third dimension of size > 1)
            @assert length(values_T) == NC.dimlength(nc, String(dimname))
            repeat(values_T, 1, length(sel))
        end
        return missing_eltype(to_SI!(A, unit; dt_val), allow_missing)
    end

    if logging
        @info "Set `$parameter [$unit]` using netCDF variable `$var`."
    end
    ncname = something(variable_name(var))
    return if V
        # dim `time` is also included in `dim_sel`: this allows for cyclic parameters (read
        # first timestep), that is later updated with the `update_cyclic!` function.
        A = read_standardized(nc, ncname, (x = :, y = :, time = 1))
        apply_affine_transform!(A, var)
        # Take out only the active cells
        finalize_parameter(T, select_cells(A, sel), var, fill, allow_missing, unit, dt_val)
    else
        A = if dimname == :layer
            read_standardized(nc, ncname, (x = :, y = :, layer = :, time = 1))
        elseif dimname == :flood_depth
            read_standardized(nc, ncname, (x = :, y = :, flood_depth = :, time = 1))
        else
            error("Unrecognized dimension name $dimname")
        end
        if !isnothing(layer)
            # the modifier index is only set in combination with scale and offset for SVectors,
            # provided through the TOML file.
            # if index, scale and offset is provided in the TOML as a list.
            for i in eachindex(layer)
                A[:, :, layer[i]] = A[:, :, layer[i]] .* scale[i] .+ offset[i]
            end
        else
            apply_affine_transform!(A, var)
        end
        # Take out only the active cells
        finalize_parameter(
            T,
            permutedims(A[sel, :]),
            var,
            fill,
            allow_missing,
            unit,
            dt_val,
        )
    end
end

"Select active cells from a 2D input array, or return the full array for `sel = nothing`."
select_cells(A::Matrix, sel::AbstractVector{CartesianIndex{2}}) = A[sel]
select_cells(A::Matrix, ::Nothing) = A

"""
Handle missing values, convert to element type `T` and convert to SI units, for a parameter
read from netCDF.
"""
function finalize_parameter(
        ::Type{T},
        A::Array{Union{Missing, Float64}, N},
        var::InputEntry,
        fill::Union{Nothing, Float64, Int, Bool},
        ::Val{true},
        unit::Unit,
        dt_val,
    ) where {T, N}
    # Convert to desired type if needed
    B = Array{Union{Missing, T}, N}(undef, size(A))
    for i in eachindex(A, B)
        x = A[i]
        B[i] = ismissing(x) ? missing : convert_parameter(T, x)
    end
    return to_SI!(B, unit; dt_val)
end

function finalize_parameter(
        ::Type{T},
        A::Array{Union{Missing, Float64}, N},
        var::InputEntry,
        fill::Union{Nothing, Float64, Int, Bool},
        ::Val{false},
        unit::Unit,
        dt_val,
    ) where {T, N}
    if isnothing(fill)
        # errors if missing are found
        C = nomissing(A)
        if any(isnan, C)
            error("NaN not allowed in $var")
        end
    else
        # replaces missing with a fill value
        fill_value = Float64(fill)
        C = nomissing(A, fill_value)
        # replace also NaN values with the fill value
        replace!(x -> isnan(x) ? fill_value : x, C)
    end
    D = T === Float64 ? C : Array{T, N}(convert_parameter.(T, C))
    return to_SI!(D, unit; dt_val)
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
        thickness::SVector{N, Float64},
    )::SVector{N, Float64} where {N}
    thicknesslayers = thickness .* MISSING_VALUE
    for i in 1:length(thicknesslayers)
        if reference_depth > cum_depth[i + 1]
            thicknesslayers = setindex(thicknesslayers, thickness[i], i)
        elseif reference_depth - cum_depth[i] > 0.0
            thicknesslayers = setindex(thicknesslayers, reference_depth - cum_depth[i], i)
        end
    end
    return thicknesslayers
end

function number_of_active_layers(thickness::SVector)::Int
    number_of_layers = length(thickness) - sum(isnan.(thickness))
    return number_of_layers
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
        return y_length
    elseif ldd == 4 || ldd == 6
        return x_length
    else
        return hypot(x_length, y_length)
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
        return x_length
    elseif ldd == 4 || ldd == 6
        return y_length
    else
        return (x_length * y_length) / hypot(x_length, y_length)
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
    return mapreduce(i -> A[i], +, inds; init = zero(T))
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

Return flow `fraction` to a river cell (at index `j`) based on the ratio of the land surface
`slope` at index `j` to the sum of the land surface `slope` at index `j` and at river cell
index `i`.
"""
function get_flow_fraction_to_river(
        graph::SimpleDiGraph{Int},
        ldd::Vector{UInt8},
        inds_river::Vector{Int},
        slope::Vector{<:Real},
    )::Vector{Float64}
    n = length(slope)
    fraction = zeros(n)
    for i in inds_river
        nbs = inneighbors(graph, i)
        for j in nbs
            if ldd[j] != ldd[i]
                fraction[j] = slope[j] / (slope[i] + slope[j])
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
    for i in 1:nv(graph)
        push!(src_edge, findall(isequal(nodes[i]), nodes_at_edge.dst))
        push!(dst_edge, findall(isequal(nodes[i]), nodes_at_edge.src))
    end
    return (src = src_edge, dst = dst_edge)
end

"Add `vertex` and `edge` to `pits` of a directed `graph`"
function add_vertex_edge_graph!(graph::SimpleDiGraph{Int}, pits::Vector{Int})::Nothing
    n = nv(graph)
    for (i, v) in enumerate(pits)
        add_vertex!(graph)
        add_edge!(graph, v, n + i)
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
    (; local_drain_direction, indices) = domain.river.network
    (; edge_indices, reverse_indices) = domain.land.network
    (; flow_width, reservoir_outlet) = domain.river.parameters
    reverse_indices = reverse_indices[indices]

    graph, local_drain_direction = flowgraph(local_drain_direction, indices, PCR_DIR)
    toposort = topological_sort_by_dfs(graph)
    n = length(we_x)
    for v in toposort
        dst = outneighbors(graph, v)
        isempty(dst) && continue
        w = min(flow_width[v], flow_width[only(dst)])
        dir = PCR_DIR[local_drain_direction[v]]
        idx = reverse_indices[v]
        # loop over river D8 directions
        if dir == CartesianIndex(1, 1)
            we_x[idx] = reservoir_outlet[v] ? 0.0 : max(we_x[idx] - 0.5 * w, 0.0)
            we_y[idx] = reservoir_outlet[v] ? 0.0 : max(we_y[idx] - 0.5 * w, 0.0)
        elseif dir == CartesianIndex(-1, -1)
            if edge_indices.ind_x_down[idx] <= n
                we_y[edge_indices.ind_x_down[idx]] =
                    reservoir_outlet[v] ? 0.0 :
                    max(we_y[edge_indices.ind_x_down[idx]] - 0.5 * w, 0.0)
            end
            if edge_indices.ind_y_down[idx] <= n
                we_x[edge_indices.ind_y_down[idx]] =
                    reservoir_outlet[v] ? 0.0 :
                    max(we_x[edge_indices.ind_y_down[idx]] - 0.5 * w, 0.0)
            end
        elseif dir == CartesianIndex(1, 0)
            we_y[idx] = reservoir_outlet[v] ? 0.0 : max(we_y[idx] - w, 0.0)
        elseif dir == CartesianIndex(0, 1)
            we_x[idx] = reservoir_outlet[v] ? 0.0 : max(we_x[idx] - w, 0.0)
        elseif dir == CartesianIndex(-1, 0)
            if edge_indices.ind_x_down[idx] <= n
                we_y[edge_indices.ind_x_down[idx]] =
                    reservoir_outlet[v] ? 0.0 :
                    max(we_y[edge_indices.ind_x_down[idx]] - w, 0.0)
            end
        elseif dir == CartesianIndex(0, -1)
            if edge_indices.ind_y_down[idx] <= n
                we_x[edge_indices.ind_y_down[idx]] =
                    reservoir_outlet[v] ? 0.0 :
                    max(we_x[edge_indices.ind_y_down[idx]] - w, 0.0)
            end
        elseif dir == CartesianIndex(1, -1)
            we_y[idx] = max(we_y[idx] - 0.5 * w, 0.0)
            if edge_indices.ind_y_down[idx] <= n
                we_x[edge_indices.ind_y_down[idx]] =
                    reservoir_outlet[v] ? 0.0 :
                    max(we_x[edge_indices.ind_y_down[idx]] - 0.5 * w, 0.0)
            end
        elseif dir == CartesianIndex(-1, 1)
            if edge_indices.ind_x_down[idx] <= n
                we_y[edge_indices.ind_x_down[idx]] =
                    reservoir_outlet[v] ? 0.0 :
                    max(we_y[edge_indices.ind_x_down[idx]] - 0.5 * w, 0.0)
            end
            we_x[idx] = reservoir_outlet[v] ? 0.0 : max(we_x[idx] - 0.5 * w, 0.0)
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
    n = Int(max(1, xs ÷ basesize))
    return (Int(1 + ((i - 1) * xs) ÷ n):Int((i * xs) ÷ n) for i in 1:n)
end

"""
    threaded_foreach(f, x::AbstractArray; basesize::Integer)

Run function `f` in parallel by spawning tasks (nthreads <= 8), each task iterates over a
chunk of size `basesize`. For nthreads > 8 run function `f` in parallel with
`Polyester@batch` with `minbatch` equal to `basesize`.
"""
function threaded_foreach(f::F, x::AbstractArray; basesize::Integer)::Nothing where {F}
    if Threads.nthreads() <= 8
        len = length(x)
        partitions = _partition(len, basesize)
        if length(partitions) > 1 && Threads.nthreads() > 1
            # collect the tasks in a typed vector rather than using `@sync`, which waits on
            # a `Channel{Any}` and therefore needs dynamic dispatch
            tasks = Task[]
            for p in partitions
                task = Threads.@spawn begin
                    for i in p
                        f(@inbounds x[i])
                    end
                end
                push!(tasks, task)
            end
            foreach(wait, tasks)
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
    hydraulic_conductivity_at_depth(p::KvExponential, vertical_hydraulic_conductivity_factor, z, i, n)
    hydraulic_conductivity_at_depth(p::KvExponentialConstant, vertical_hydraulic_conductivity_factor, z, i, n)
    hydraulic_conductivity_at_depth(p::KvLayered, vertical_hydraulic_conductivity_factor, z, i, n)
    hydraulic_conductivity_at_depth(p::KvLayeredExponential, vertical_hydraulic_conductivity_factor, z, i, n)

Return vertical hydraulic conductivity `kv_z` at depth `z` for index `i` using multiplication
factor `kv_frac` at soil layer `n` and vertical hydraulic conductivity profile `p`.
"""
function hydraulic_conductivity_at_depth(
        p::KvExponential,
        vertical_hydraulic_conductivity_factor,
        z,
        i,
        n,
    )
    kv_z =
        vertical_hydraulic_conductivity_factor[i][n] *
        p.kv_0[i] *
        exp(-p.hydraulic_conductivity_scale_parameter[i] * z)
    return kv_z
end

function hydraulic_conductivity_at_depth(
        p::KvExponentialConstant,
        vertical_hydraulic_conductivity_factor,
        z,
        i,
        n,
    )
    (; kv_0, hydraulic_conductivity_scale_parameter) = p.exponential
    if z < p.z_exp[i]
        kv_z =
            vertical_hydraulic_conductivity_factor[i][n] *
            kv_0[i] *
            exp(-hydraulic_conductivity_scale_parameter[i] * z)
    else
        kv_z =
            vertical_hydraulic_conductivity_factor[i][n] *
            kv_0[i] *
            exp(-hydraulic_conductivity_scale_parameter[i] * p.z_exp[i])
    end
    return kv_z
end

function hydraulic_conductivity_at_depth(
        p::KvLayered,
        vertical_hydraulic_conductivity_factor,
        z,
        i,
        n,
    )
    kv_z = vertical_hydraulic_conductivity_factor[i][n] * p.kv[i][n]
    return kv_z
end

function hydraulic_conductivity_at_depth(
        p::KvLayeredExponential,
        vertical_hydraulic_conductivity_factor,
        z,
        i,
        n,
    )
    kv_z = if z < p.z_layered[i]
        vertical_hydraulic_conductivity_factor[i][n] * p.kv[i][n]
    else
        n = p.nlayers_kv[i]
        vertical_hydraulic_conductivity_factor[i][n] *
            p.kv[i][n] *
            exp(-p.hydraulic_conductivity_scale_parameter[i] * (z - p.z_layered[i]))
    end
    return kv_z
end

"""
    kh_layered_profile!(soil_model::SbmSoilModel, subsurface_flow_model::LateralSSFModel, kv_profile::KvLayered)
    kh_layered_profile!(soil_model::SbmSoilModel, subsurface_flow_model::LateralSSFModel, kv_profile::KvLayeredExponential)

Compute equivalent horizontal hydraulic conductivity `kh` [m s⁻¹] using vertical hydraulic
conductivity profile `kv_profile`.
"""
function kh_layered_profile!(
        soil_model::SbmSoilModel,
        subsurface_flow_model::LateralSSFModel,
        kv_profile::KvLayered,
    )
    (; number_of_layers, cumulative_layer_depth, actual_layer_thickness, soil_thickness) =
        soil_model.parameters
    (; n_unsatlayers, water_table_depth) = soil_model.variables
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; horizontal_to_vertical_hydraulic_conductivity_ratio) =
        subsurface_flow_model.parameters

    for i in eachindex(kh)
        m = number_of_layers[i]

        if soil_thickness[i] > water_table_depth[i]
            transmissivity = 0.0
            _sumlayers = @view cumulative_layer_depth[i][2:end]
            n = max(n_unsatlayers[i], 1)
            transmissivity += (_sumlayers[n] - water_table_depth[i]) * kv_profile.kv[i][n]
            n += 1
            while n <= m
                transmissivity += actual_layer_thickness[i][n] * kv_profile.kv[i][n]
                n += 1
            end
            kh[i] =
                (transmissivity / (soil_thickness[i] - water_table_depth[i])) *
                horizontal_to_vertical_hydraulic_conductivity_ratio[i]
        else
            kh[i] =
                kv_profile.kv[i][m] * horizontal_to_vertical_hydraulic_conductivity_ratio[i]
        end
    end
    return nothing
end

function kh_layered_profile!(
        soil_model::SbmSoilModel,
        subsurface_flow_model::LateralSSFModel,
        kv_profile::KvLayeredExponential,
    )
    (; number_of_layers, cumulative_layer_depth, actual_layer_thickness, soil_thickness) =
        soil_model.parameters
    (; nlayers_kv, z_layered, kv, hydraulic_conductivity_scale_parameter) = kv_profile
    (; n_unsatlayers, water_table_depth) = soil_model.variables
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; horizontal_to_vertical_hydraulic_conductivity_ratio) =
        subsurface_flow_model.parameters

    for i in eachindex(kh)
        m = number_of_layers[i]

        if soil_thickness[i] > water_table_depth[i]
            transmissivity = 0.0
            n = max(n_unsatlayers[i], 1)
            if water_table_depth[i] >= z_layered[i]
                zt = soil_thickness[i] - z_layered[i]
                j = nlayers_kv[i]
                transmissivity +=
                    kv[i][j] / hydraulic_conductivity_scale_parameter[i] * (
                    exp(
                        -hydraulic_conductivity_scale_parameter[i] *
                            (water_table_depth[i] - z_layered[i]),
                    ) - exp(-hydraulic_conductivity_scale_parameter[i] * zt)
                )
                n = m
            else
                _sumlayers = @view cumulative_layer_depth[i][2:end]
                transmissivity += (_sumlayers[n] - water_table_depth[i]) * kv[i][n]
            end
            n += 1
            while n <= m
                if n > nlayers_kv[i]
                    zt = soil_thickness[i] - z_layered[i]
                    j = nlayers_kv[i]
                    transmissivity +=
                        kv[i][j] / hydraulic_conductivity_scale_parameter[i] *
                        (1.0 - exp(-hydraulic_conductivity_scale_parameter[i] * zt))
                    n = m
                else
                    transmissivity += actual_layer_thickness[i][n] * kv[i][n]
                end
                n += 1
            end
            kh[i] =
                (transmissivity / (soil_thickness[i] - water_table_depth[i])) *
                horizontal_to_vertical_hydraulic_conductivity_ratio[i]
        else
            if water_table_depth[i] >= z_layered[i]
                j = nlayers_kv[i]
                kh[i] =
                    kv[i][j] *
                    exp(
                    -hydraulic_conductivity_scale_parameter[i] *
                        (water_table_depth[i] - z_layered[i]),
                ) *
                    horizontal_to_vertical_hydraulic_conductivity_ratio[i]
            else
                kh[i] = kv[i][m] * horizontal_to_vertical_hydraulic_conductivity_ratio[i]
            end
        end
    end
    return nothing
end

kh_layered_profile!(
    soil_model::SbmSoilModel,
    subsurface_flow_model::LateralSSFModel,
    kv_profile::Union{KvExponential, KvExponentialConstant},
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
    (; kh_0, hydraulic_conductivity_scale_parameter) = kh_profile
    (; q, q_max, water_table_depth) = subsurface_flow_model.variables
    (; soil_thickness, slope) = subsurface_flow_model.parameters
    (; flow_width) = parameters

    @. q_max =
        ((kh_0 * slope) / hydraulic_conductivity_scale_parameter) *
        (1.0 - exp(-hydraulic_conductivity_scale_parameter * soil_thickness))
    @. q =
        ((kh_0 * slope) / hydraulic_conductivity_scale_parameter) *
        (
        exp(-hydraulic_conductivity_scale_parameter * water_table_depth) -
            exp(-hydraulic_conductivity_scale_parameter * soil_thickness)
    ) *
        flow_width
    return nothing
end

function initialize_lateral_ssf_model!(
        subsurface_flow_model::LateralSSFModel,
        parameters::LandParameters,
        kh_profile::KhExponentialConstant,
    )
    (; kh_0, hydraulic_conductivity_scale_parameter) = kh_profile.exponential
    (; z_exp) = kh_profile
    (; q, q_max, water_table_depth) = subsurface_flow_model.variables
    (; soil_thickness, slope) = subsurface_flow_model.parameters
    (; flow_width) = parameters

    q_constant = @. kh_0 *
        exp(-hydraulic_conductivity_scale_parameter * z_exp) *
        slope *
        (soil_thickness - z_exp)
    for i in eachindex(q)
        q_max[i] =
            ((kh_0[i] * slope[i]) / hydraulic_conductivity_scale_parameter[i]) *
            (1.0 - exp(-hydraulic_conductivity_scale_parameter[i] * z_exp[i])) +
            q_constant[i]
        if water_table_depth[i] < z_exp[i]
            q[i] =
                (
                ((kh_0[i] * slope[i]) / hydraulic_conductivity_scale_parameter[i]) * (
                    exp(
                        -hydraulic_conductivity_scale_parameter[i] *
                            water_table_depth[i],
                    ) - exp(-hydraulic_conductivity_scale_parameter[i] * z_exp[i])
                ) + q_constant[i]
            ) * flow_width[i]
        else
            q[i] =
                kh_0[i] *
                exp(-hydraulic_conductivity_scale_parameter[i] * water_table_depth[i]) *
                slope[i] *
                (soil_thickness[i] - water_table_depth[i]) *
                flow_width[i]
        end
    end
    return nothing
end

"""
    initialize_lateral_ssf_model!(subsurface_flow_model::LateralSSFModel, soil_model::SbmSoilModel, parameters::LandParameters, kv_profile::KvLayered)
    initialize_lateral_ssf_model!(subsurface_flow_model::LateralSSFModel, soil_model::SbmSoilModel, parameters::LandParameters, kv_profile::KvLayeredExponential)

Initialize lateral subsurface variables `q` and `q_max` using  vertical hydraulic
conductivity profile `kv_profile`.
"""
function initialize_lateral_ssf_model!(
        subsurface_flow_model::LateralSSFModel,
        soil_model::SbmSoilModel,
        parameters::LandParameters,
        kv_profile::KvLayered,
    )
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; number_of_layers, actual_layer_thickness) = soil_model.parameters
    (; q, q_max, water_table_depth) = subsurface_flow_model.variables
    (; horizontal_to_vertical_hydraulic_conductivity_ratio, soil_thickness, slope) =
        subsurface_flow_model.parameters
    (; flow_width) = parameters

    kh_layered_profile!(soil_model, subsurface_flow_model, kv_profile)
    for i in eachindex(q)
        q[i] = kh[i] * (soil_thickness[i] - water_table_depth[i]) * slope[i] * flow_width[i]
        kh_max = 0.0
        for j in 1:number_of_layers[i]
            kh_max += kv_profile.kv[i][j] * actual_layer_thickness[i][j]
        end
        kh_max *= horizontal_to_vertical_hydraulic_conductivity_ratio[i]
        q_max[i] = kh_max * slope[i]
    end
    return nothing
end

function initialize_lateral_ssf_model!(
        subsurface_flow_model::LateralSSFModel,
        soil_model::SbmSoilModel,
        parameters::LandParameters,
        kv_profile::KvLayeredExponential,
    )
    (; q, q_max, water_table_depth) = subsurface_flow_model.variables
    (; horizontal_to_vertical_hydraulic_conductivity_ratio, soil_thickness, slope) =
        subsurface_flow_model.parameters
    (; flow_width) = parameters
    (; number_of_layers, actual_layer_thickness) = soil_model.parameters
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; kv, hydraulic_conductivity_scale_parameter, nlayers_kv, z_layered) = kv_profile

    kh_layered_profile!(soil_model, subsurface_flow_model, kv_profile)
    for i in eachindex(q)
        q[i] = kh[i] * (soil_thickness[i] - water_table_depth[i]) * slope[i] * flow_width[i]
        kh_max = 0.0
        for j in 1:number_of_layers[i]
            if j <= nlayers_kv[i]
                kh_max += kv[i][j] * actual_layer_thickness[i][j]
            else
                zt = soil_model.parameters.soil_thickness[i] - z_layered[i]
                k = max(j - 1, 1)
                kh_max +=
                    kv[i][k] / hydraulic_conductivity_scale_parameter[i] *
                    (1.0 - exp(-hydraulic_conductivity_scale_parameter[i] * zt))
                break
            end
        end
        kh_max = kh_max * horizontal_to_vertical_hydraulic_conductivity_ratio[i]
        q_max[i] = kh_max * slope[i]
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
    bounded_power(base, power)

Computes min(base^power, 1) without computing the power
if the result is known to be larger than 1.
Assumes base, power > 0
"""
function bounded_power(base::T, power) where {T}
    if base > 1
        return one(T)
    else
        return pow(base, power)
    end
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
        i::Int,
        dt::Float64,
    )
    (; n_unsatlayers, unsaturated_layer_thickness, unsaturated_layer_depth) =
        soil_model.variables
    (; theta_s, theta_r) = soil_model.parameters

    # effective porosity (difference between saturated and residual water content)
    theta_e = theta_s[i] - theta_r[i]

    if net_flux <= 0.0
        dh = net_flux * dt / specific_yield
    else
        dh = 0.0
        for k in n_unsatlayers[i]:-1:1
            capacity =
                max(
                unsaturated_layer_thickness[i][k] * theta_e -
                    unsaturated_layer_depth[i][k],
                0.0,
            ) / dt
            flux_layer = min(net_flux, capacity)
            if capacity <= net_flux
                # if unsaturated layer is fully saturated dh equals layer thickness
                dh += unsaturated_layer_thickness[i][k]
            else
                sy =
                    theta_e -
                    (unsaturated_layer_depth[i][k] / unsaturated_layer_thickness[i][k])
                dh += flux_layer * dt / sy
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
