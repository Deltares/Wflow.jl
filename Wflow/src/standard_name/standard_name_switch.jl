"""
    lens_property_names(lens)::Union{Nothing, Vector{Symbol}}

Return the property names of a lens consisting only of property access, e.g.
`[:land, :soil, :variables, :q]` for `@optic(_.land.soil.variables.q)`. Returns `nothing`
for other lenses.
"""
function lens_property_names(lens)::Union{Nothing, Vector{Symbol}}
    isnothing(lens) && return nothing
    names = Symbol[]
    # decompose returns the outermost optic first
    for optic in reverse(collect(decompose(lens)))
        optic isa PropertyLens || return nothing
        push!(names, typeof(optic).parameters[1])
    end
    return names
end

"""
    without_lens(metadata::ParameterMetadata)

Return a copy of the metadata without lens. The lens has a unique type per parameter, so
dropping it limits the number of method specializations of functions that take metadata.
"""
without_lens(metadata::ParameterMetadata) = ParameterMetadata(
    nothing,
    metadata.unit,
    metadata.default,
    metadata.fill,
    metadata.type,
    metadata.description,
    metadata.allow_missing,
    metadata.allow_dynamic_input,
    metadata.dimname,
    metadata.tags,
)

"Whether the property chain `names` exists for (concrete) type `T`."
function has_property_chain(T::Type, names::Vector{Symbol})::Bool
    for name in names
        isconcretetype(T) || return false
        hasfield(T, name) || return false
        T = fieldtype(T, name)
    end
    return true
end

"""
    with_standard_name(f, model, name::String)

Call `f(field, metadata)` with the model `field` and `ParameterMetadata` corresponding to
the standard name `name`, or the path in the model, e.g. "land.soil.variables.q". The
lens is removed from the metadata, see `without_lens`. Unlike `get_field_in_model`, this is type stable and avoids
dynamic dispatch: the lens lookup is unrolled at compile time into a chain of string
comparisons, where each branch has a concrete field and metadata type. Only standard names
whose lens exists in the concrete type of `model` are included. This makes it compatible
with `juliac --trim`.
"""
@generated function with_standard_name(f::F, model::M, name::String) where {F, M}
    body = Expr(:block)
    seen = Set{String}()
    for (_, standard_name_map, _) in STANDARD_NAME_MAPS
        for (standard_name, metadata) in standard_name_map
            names = lens_property_names(metadata.lens)
            isnothing(names) && continue
            has_property_chain(M, names) || continue
            getter = foldl(names; init = :model) do ex, property_name
                :(getproperty($ex, $(QuoteNode(property_name))))
            end
            metadata_without_lens = without_lens(metadata)
            # support both the standard name and the path in the model, like
            # "land.soil.variables.q"
            for key in (standard_name, join(names, '.'))
                key in seen && continue
                push!(seen, key)
                push!(body.args, :(name == $key && return f($getter, $metadata_without_lens)))
            end
        end
    end
    push!(body.args, :(throw(ArgumentError(string("Unknown standard name: ", name)))))
    return body
end

"""
    has_standard_name(model, name::String)::Bool

Whether `with_standard_name` supports `name` for this `model`.
"""
@generated function has_standard_name(model::M, name::String)::Bool where {M}
    body = Expr(:block)
    for (_, standard_name_map, _) in STANDARD_NAME_MAPS
        for (standard_name, metadata) in standard_name_map
            names = lens_property_names(metadata.lens)
            isnothing(names) && continue
            has_property_chain(M, names) || continue
            for key in (standard_name, join(names, '.'))
                push!(body.args, :(name == $key && return true))
            end
        end
    end
    push!(body.args, :(return false))
    return body
end

"""
    allow_dynamic_field_access(model)

Whether model fields can be accessed by a path that is not a standard name, e.g.
"routing.river_flow.variables.q". This needs dynamic dispatch, and is therefore not
supported for builds with a fixed build profile, see `SbmModel`.
"""
allow_dynamic_field_access(::AbstractModel{T}) where {T} = allow_dynamic_field_access(T)
allow_dynamic_field_access(modelmap::NamedTuple) =
    haskey(modelmap, :type) ? allow_dynamic_field_access(typeof(modelmap.type)) : true
allow_dynamic_field_access(::Type{SbmModel{P}}) where {P} = P === nothing
allow_dynamic_field_access(::Type{<:AbstractModelType}) = true

"""
    with_model_field(f, model, name::String)

Call `f(field, metadata)` for the model field corresponding to `name`, a standard name or a
path in the model. Standard names and paths matching the lens of a standard name are looked
up type stably with `with_standard_name`. Other paths are looked up dynamically if allowed,
see `allow_dynamic_field_access`, with default metadata.
"""
function with_model_field(f::F, model, name::String) where {F}
    if has_standard_name(model, name)
        return with_standard_name(f, model, name)
    elseif allow_dynamic_field_access(model)
        field = try
            param(model, name)
        catch
            error("Couldn't obtain a field from this model specified by '$name'.")
        end
        @warn "No metadata was found for $name, so the output will be expressed in standard SI units ($(join(STANDARD_UNITS, ", "))) and might fail."
        return f(field, ParameterMetadata())
    else
        error("Couldn't obtain a field from this model specified by '$name'.")
    end
end
