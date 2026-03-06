# wrapper methods for standard name mapping
get_standard_name_map(model::T) where {T} = get_standard_name_map(T)
get_standard_name_map(::Type{<:LandHydrologySBM}) = sbm_standard_name_map
get_standard_name_map(::Type{<:SoilLoss}) = sediment_standard_name_map
get_standard_name_map(::Type{<:Domain}) = domain_standard_name_map
get_standard_name_map(::Type{<:Routing}) = routing_standard_name_map

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
- `allow_dynamic_input`: Allow updating this parameter from input via cyclic/forcing
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
    allow_dynamic_input::Bool = false
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
        allow_dynamic_input,
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
            allow_dynamic_input,
            dimname,
            flags,
        )
    end
end

function metadata_from_lens_string(
    lens_string::AbstractString,
    standard_name_map::OrderedDict{String, ParameterMetadata},
)::Union{ParameterMetadata, Nothing}
    for metadata_candidate in values(standard_name_map)
        if string(metadata_candidate.lens)[7:(end - 1)] == lens_string
            return metadata_candidate
        end
    end
    return nothing
end

function get_metadata(
    name::AbstractString,
    types::Vararg{Type};
    model = nothing,
)::Union{ParameterMetadata, Nothing}
    metadata = nothing
    for type in types
        standard_name_map = get_standard_name_map(type)
        # First see whether 'name' is a standard name within the standard name map
        # corresponding to 'type'
        metadata_candidate = get(standard_name_map, name, nothing)
        # If not, see whether 'name' is a path in the model object which matches
        # a lens in the standard name map
        if isnothing(metadata_candidate)
            metadata_candidate = metadata_from_lens_string(name, standard_name_map)
        end

        if !isnothing(metadata_candidate)
            # Metadata was found; if a model was provided check whether
            # the lens matches
            if !isnothing(model) && !isnothing(metadata_candidate.lens)
                found_matching_lens = false
                try
                    metadata_candidate.lens(model)
                    found_matching_lens = true
                catch
                end
                if found_matching_lens
                    !isnothing(metadata) && error(
                        "Ambiguity found for obtaining metadata for '$name'; this key is in the standard nampe map for at least 2 of $types with a fitting lens.",
                    )
                    metadata = metadata_candidate
                end
            else
                # If model or lens was not provided assume that the metadata matches
                !isnothing(metadata) && error(
                    "Ambiguity found for obtaining metadata for '$name'; this key is in the standard name map for at least 2 of $types and there was no model provided to disambiguate.",
                )
                metadata = metadata_candidate
            end
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

get_metadata(name::AbstractString, model) = get_metadata(name, typeof(model); model)
get_metadata(name::AbstractString, L::Type) = get_standard_name_map(L)[name]

# When no model or model type is specified, search all standard name maps
get_metadata(name::AbstractString; kwargs...) =
    get_metadata(name, map(d -> d[3], Wflow.STANDARD_NAME_MAPS)...; kwargs...)

function get_field_in_model(model, name::AbstractString; check_allow_dynamic_input = false)
    metadata = get_metadata(name; model)

    return if !isnothing(metadata)
        # If metadata was found, `str` is a standard name or a path in the model object which matches a lens
        if check_allow_dynamic_input && !metadata.allow_dynamic_input
            error(
                "Tried to set '$name' dynamically via cyclic/forcing input, which is not allowed.",
            )
        end
        metadata.lens(model)
    else
        # If no metadata was found, `str` is either a path in the model object that doesn't match a lens or is invalid
        try
            param(model, name)
        catch
            error("Couldn't obtain a field from this model specified by '$name'.")
        end
    end
end
