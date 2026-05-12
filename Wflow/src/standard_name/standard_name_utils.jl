# wrapper methods for standard name mapping
get_standard_name_map(model::T) where {T} = get_standard_name_map(T)
get_standard_name_map(::Type{<:LandHydrologySBM}) = sbm_standard_name_map
get_standard_name_map(::Type{<:SoilLossModel}) = sediment_standard_name_map
get_standard_name_map(::Type{<:Domain}) = domain_standard_name_map
get_standard_name_map(::Type{<:Routing}) = routing_standard_name_map

const PARAMETER_TYPES = Union{Float64, Int, Bool, Nothing}

"""
Metadata associated with parameters and variables.

# Arguments
- `unit`: The unit of the parameter/variable in the Wflow input
- `default`: The default (initial) value of the parameter/variable if it exists
- `fill`: Missing input values are replaced by this value if allow_missing == false
- `type`: The output type of the data. Assumed to be `Float64` if it is not provided and cannot be derived
    from `default` or `fill`
- `description`: The description of the parameter/variable provided in the Wflow docs
- `allow_missing`: Whether the parameter/variable is allowed to have missing entries
- `allow_dynamic_input`: Allow updating this parameter from input via cyclic/forcing
- `dimname`: The name of the third dimension of the parameter/variable if it exists
- `tags`: Identifiers to filter parameters/variables for specific tables in the docs
"""
@kwdef struct ParameterMetadata{
    D <: PARAMETER_TYPES,
    F <: PARAMETER_TYPES,
    T <: PARAMETER_TYPES,
    N <: Union{Symbol, Nothing},
}
    unit::Unit = EMPTY_UNIT
    default::D = nothing
    fill::F = nothing
    type::Type{T} = nothing
    description::String = ""
    allow_missing::Bool = false
    allow_dynamic_input::Bool = false
    dimname::N = nothing
    tags::Vector{Symbol} = []
    function ParameterMetadata(
        unit,
        default::D,
        fill::F,
        type,
        description,
        allow_missing,
        allow_dynamic_input,
        dimname::N,
        flags,
    ) where {D, F, N}
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
        return new{D, F, type, N}(
            unit,
            default,
            fill,
            type,
            description,
            allow_missing,
            allow_dynamic_input,
            dimname,
            flags,
        )
    end
end

function get_metadata(
    name::AbstractString,
    types::Vararg{Type};
    kwargs...,
)::Union{ParameterMetadata, Nothing}
    metadata = nothing
    for type in types
        standard_name_map = get_standard_name_map(type)
        metadata_candidate = get(standard_name_map, name, nothing)

        if !isnothing(metadata_candidate)
            !isnothing(metadata) && error(
                "Ambiguity found for obtaining metadata for '$name'; this key is in the standard name map for at least 2 of $types.",
            )
            metadata = metadata_candidate
        end
    end
    return metadata
end

function get_metadata(
    name::AbstractString,
    land::AbstractLandModel;
    kwargs...,
)::Union{ParameterMetadata, Nothing}
    # Check whether it is a land variable first
    metadata = get(get_standard_name_map(land), name, nothing)

    if isnothing(metadata)
        # Then check other variable types
        for (name_map, standard_name_map) in STANDARD_NAME_MAPS
            (name_map ∈ ("sbm", "sediment")) && continue
            metadata = get(standard_name_map, name, nothing)
            !isnothing(metadata) && break
        end
    end
    return metadata
end

get_metadata(name::AbstractString, model) = get_metadata(name, typeof(model))
get_metadata(name::AbstractString, L::Type) = get_standard_name_map(L)[name]

# When no model or model type is specified, search all standard name maps
get_metadata(name::AbstractString) =
    get_metadata(name, map(d -> d[3], Wflow.STANDARD_NAME_MAPS)...)

function get_field_in_model(model, name::AbstractString; check_allow_dynamic_input = false)
    # Scope metadata lookup to the model's land type to avoid ambiguity
    # between standard name maps (e.g. routing vs sediment)
    land = hasproperty(model, :land) ? model.land : nothing
    metadata = if !isnothing(land) && land isa AbstractLandModel
        get_metadata(name, land)
    else
        get_metadata(name)
    end

    if check_allow_dynamic_input && !isnothing(metadata) && !metadata.allow_dynamic_input
        error(
            "Tried to set '$name' dynamically via cyclic/forcing input, which is not allowed.",
        )
    end

    # Try data_lookup first (standard name registered at construction time)
    if hasproperty(model, :data_lookup) && haskey(model.data_lookup, name)
        return model.data_lookup[name], metadata
    end

    # Fall back to param (path-based access into the model struct tree)
    try
        return param(model, name), metadata
    catch
        error("Couldn't obtain a field from this model specified by '$name'.")
    end
end
