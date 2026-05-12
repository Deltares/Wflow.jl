"""
    DataLookup

Dict-based lookup for model parameters and variables, keyed by standard name strings.
Stores references to the underlying arrays so that downstream code can access model
data by standard name without knowing the struct layout.

Separate dictionaries are kept per element type for type stability.
"""
struct DataLookup{N}
    data_float::Dict{String, Vector{Float64}}
    data_int::Dict{String, Vector{Int}}
    data_bool::Dict{String, Vector{Bool}}
    data_svector::Dict{String, Vector{SVector{N, Float64}}}
end

function DataLookup{N}() where {N}
    DataLookup{N}(
        Dict{String, Vector{Float64}}(),
        Dict{String, Vector{Int}}(),
        Dict{String, Vector{Bool}}(),
        Dict{String, Vector{SVector{N, Float64}}}(),
    )
end

"Convenience constructor without specifying N (defaults to N=1)."
DataLookup() = DataLookup{1}()

"""
    DataLookup(config::Config)

Construct a `DataLookup{N}` where `N` is derived from the number of soil layers
configured in `config.model.soil_layer__thickness` (N = length + 1).
"""
function DataLookup(config::Config)
    N = length(config.model.soil_layer__thickness) + 1
    DataLookup{N}()
end

function Base.setindex!(
    data_lookup::DataLookup,
    value::Vector{Float64},
    key::AbstractString,
)
    data_lookup.data_float[key] = value
end
function Base.setindex!(data_lookup::DataLookup, value::Vector{Int}, key::AbstractString)
    data_lookup.data_int[key] = value
end
function Base.setindex!(data_lookup::DataLookup, value::Vector{Bool}, key::AbstractString)
    data_lookup.data_bool[key] = value
end
function Base.setindex!(
    data_lookup::DataLookup{N},
    value::Vector{<:SVector{N, Float64}},
    key::AbstractString,
) where {N}
    data_lookup.data_svector[key] = value
end
function Base.setindex!(
    data_lookup::DataLookup,
    value::Vector{<:EnumX.Enum},
    key::AbstractString,
)
    data_lookup.data_int[key] = Int.(value)
end
function Base.setindex!(data_lookup::DataLookup, value::Vector, key::AbstractString)
    # fallback: convert to Float64 if possible, otherwise error
    data_lookup.data_float[key] = convert(Vector{Float64}, value)
end

function Base.show(io::IO, data_lookup::DataLookup{N}) where {N}
    nf = length(data_lookup.data_float)
    ni = length(data_lookup.data_int)
    nb = length(data_lookup.data_bool)
    ns = length(data_lookup.data_svector)
    total = nf + ni + nb + ns
    println(io, "DataLookup{$N}($total entries)")
    for (label, dict) in (
        ("Float64", data_lookup.data_float),
        ("Int", data_lookup.data_int),
        ("Bool", data_lookup.data_bool),
        ("SVector{$N,Float64}", data_lookup.data_svector),
    )
        isempty(dict) && continue
        println(io, "  $label ($(length(dict))):")
        for key in sort!(collect(keys(dict)))
            println(io, "    ", key)
        end
    end
end

function Base.getindex(data_lookup::DataLookup, key::AbstractString)
    haskey(data_lookup.data_float, key) && return data_lookup.data_float[key]
    haskey(data_lookup.data_int, key) && return data_lookup.data_int[key]
    haskey(data_lookup.data_bool, key) && return data_lookup.data_bool[key]
    haskey(data_lookup.data_svector, key) && return data_lookup.data_svector[key]
    throw(KeyError(key))
end

function Base.haskey(data_lookup::DataLookup, key::AbstractString)
    return haskey(data_lookup.data_float, key) ||
           haskey(data_lookup.data_int, key) ||
           haskey(data_lookup.data_bool, key) ||
           haskey(data_lookup.data_svector, key)
end

"""
    @with_data_lookup struct Name
        "standard_name"
        field::Type = default
        ...
    end

Wraps `@kwdef` to create a keyword-argument constructor with default values.
Additionally generates a method `Name(data_lookup::DataLookup; kwargs...)` that
constructs the struct via `@kwdef`, then registers every field preceded by a
standard-name string literal into `data_lookup`.

Fields **without** a preceding string annotation behave exactly as with `@kwdef`.

# Example

```julia
@with_data_lookup struct GlacierVariables
    n::Int
    "glacier_ice__leq_depth"
    glacier_store::Vector{Float64}
    "glacier_ice__melt_volume_flux"
    glacier_melt::Vector{Float64} = fill(MISSING_VALUE, n)
end
```

Calling `GlacierVariables(; n, glacier_store)` works identically to `@kwdef`.
Calling `GlacierVariables(data_lookup; n, glacier_store)` additionally runs:
```julia
data_lookup["glacier_ice__leq_depth"] = obj.glacier_store
data_lookup["glacier_ice__melt_volume_flux"] = obj.glacier_melt
```
"""
macro with_data_lookup(struct_expr)
    return esc(_build_data_lookup_struct(struct_expr))
end

function _build_data_lookup_struct(raw_expr::Expr)
    raw_expr.head === :struct ||
        error("@with_data_lookup must be applied to a struct definition")

    # ── Collect string annotations and strip them from the body ─────────
    body = raw_expr.args[3]::Expr
    lookup_pairs = Pair{String, Symbol}[]     # "standard_name" => field_symbol
    pending_annotation::Union{String, Nothing} = nothing
    cleaned_items = []

    for item in body.args
        if item isa String
            pending_annotation = item
            continue
        end
        push!(cleaned_items, item)
        if !isnothing(pending_annotation)
            fname = _field_name(item)
            if !isnothing(fname)
                push!(lookup_pairs, pending_annotation => fname)
            end
            pending_annotation = nothing
        end
    end

    # Rebuild body without string annotations
    raw_expr.args[3] = Expr(:block, cleaned_items...)

    struct_name = _extract_struct_name(raw_expr.args[2])

    # ── Let @kwdef handle the struct definition + keyword constructor ───
    kwdef_result = :(@kwdef $raw_expr)

    if isempty(lookup_pairs)
        return kwdef_result
    end

    # ── Generate wrapper: Name(data_lookup::DataLookup; kwargs...) ──────
    # Forwards kwargs to the @kwdef constructor, then registers fields.
    # Always use the bare struct name (no type parameters) so that the
    # @kwdef-generated keyword constructor can infer type parameters from
    # the field values.  Using Name{N,M,...}(data_lookup; ...) would fail
    # because Julia cannot infer the type parameters from `data_lookup`.
    registrations = [:(data_lookup[$sn] = obj.$fs) for (sn, fs) in lookup_pairs]

    fn_body = quote
        obj = $struct_name(; kwargs...)
        $(registrations...)
        return obj
    end

    fn_call = :($struct_name(data_lookup::DataLookup; kwargs...))
    wrapper = Expr(:function, fn_call, fn_body)

    return quote
        $kwdef_result
        $wrapper
    end
end

# ── helpers ────────────────────────────────────────────────────────────────

"Extract the bare struct name (a `Symbol`) from the header expression."
function _extract_struct_name(header)
    if header isa Symbol
        return header
    elseif header isa Expr
        if header.head === :curly
            return header.args[1]::Symbol
        elseif header.head === :<:
            return _extract_struct_name(header.args[1])
        end
    end
    error("@with_data_lookup: cannot parse struct name from: $header")
end

"Extract the field name `Symbol` from a body expression, or `nothing`."
function _field_name(ex)
    ex isa Expr || return nothing
    if ex.head === :(::) && length(ex.args) == 2 && ex.args[1] isa Symbol
        return ex.args[1]
    elseif ex.head === :(=) && ex.args[1] isa Expr
        return _field_name(ex.args[1])
    end
    return nothing
end
