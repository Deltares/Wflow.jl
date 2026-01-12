# wrapper methods for standard name mapping
standard_name_map(model::T) where {T} = standard_name_map(T)
standard_name_map(::Type{<:LandHydrologySBM}) = sbm_standard_name_map
standard_name_map(::Type{<:SoilLoss}) = sediment_standard_name_map
standard_name_map(::Type{<:Domain}) = domain_standard_name_map
standard_name_map(::Type{<:Routing}) = routing_standard_name_map

const PARAMETER_TYPES = Union{Float64, Int, Bool, Nothing}

"""
Metadata associated with parameters and variables.

# Arguments
- `lens`: The path in the model data structure to the parameter/variable if it exists
- `unit`: The unit of the parameter/variable in the Wflow input
- `default`: The default (initial) value of the parameter/variable if it exists
- `fill`: Missing input values are replaced by this value if allow_missing == false
- `type`: The output type of the data. Assumed to be `Float64` if it is not provided and cannot be derived
    from `default` or `fill`
- `description`: The description of the parameter/variable provided in the Wflow docs
- `allow_missing`: Whether the parameter/variable is allowed to have missing entries
- `dimname`: The name of the third dimension of the parameter/variable if it exists
- `tags`: Identifiers to filter parameters/variables for specific tables in the docs
"""
@kwdef struct ParameterMetadata{
    L,
    D <: PARAMETER_TYPES,
    F <: PARAMETER_TYPES,
    T <: PARAMETER_TYPES,
    N <: Union{Symbol, Nothing},
}
    lens::L = nothing
    unit::Unit = EMPTY_UNIT
    default::D = nothing
    fill::F = nothing
    type::Type{T} = nothing
    description::String = ""
    allow_missing::Bool = false
    dimname::N = nothing
    tags::Vector{Symbol} = []
    function ParameterMetadata(
        lens::L,
        unit,
        default::D,
        fill::F,
        type,
        description,
        allow_missing,
        dimname::N,
        flags,
    ) where {L, D, F, N}
        if isnothing(type)
            type = if !isnothing(default)
                D
            elseif !isnothing(fill)
                F
            else
                # Assume the type is Float64 if it is not provided and cannot be derived
                # from the default or fill
                Float64
            end
        end
        return new{L, D, F, type, N}(
            lens,
            unit,
            default,
            fill,
            type,
            description,
            allow_missing,
            dimname,
            flags,
        )
    end
end

get_metadata(name::AbstractString, ::Type{<:Writer}) = ParameterMetadata()
get_metadata(name::AbstractString, model) = get_metadata(name, typeof(model))
get_metadata(name::AbstractString, L::Type) = standard_name_map(L)[name]

function get_metadata(name::AbstractString, types::Vararg{Type})
    for type in types
        metadata = get(standard_name_map(type), name, nothing)
        !isnothing(metadata) && return metadata
    end
    return nothing
end

function get_metadata(name::AbstractString, land::AbstractLandModel)
    # Check whether it is a land variable first
    metadata = get(standard_name_map(land), name, nothing)

    if isnothing(metadata)
        # Then check other variable types
        for (name_map, _standard_name_map) in standard_name_maps
            (name_map ∈ ("sbm", "sediment")) && continue
            metadata = get(_standard_name_map, name, nothing)
            !isnothing(metadata) && break
        end
    end
    return metadata
end

get_metadata(name::AbstractString) =
    get_metadata(name, LandHydrologySBM, Routing, SoilLoss, Domain)

get_lens(name::AbstractString, model::T) where {T} = get_lens(name, T)
get_lens(name::AbstractString, ::Type{T}) where {T} = get_metadata(name, T).lens
