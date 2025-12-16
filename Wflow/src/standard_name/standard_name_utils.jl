# wrapper methods for standard name mapping
standard_name_map(model) = standard_name_map(typeof(model))
standard_name_map(::Type{<:LandHydrologySBM}) = sbm_standard_name_map
standard_name_map(::Type{<:SoilLoss}) = sediment_standard_name_map
standard_name_map(::Type{<:Domain}) = domain_standard_name_map
standard_name_map(::Type{<:Routing}) = routing_standard_name_map

get_lens(name::AbstractString, model) = get_lens(name, typeof(model))
get_lens(name::AbstractString, L::Type) = standard_name_map(L)[name].lens

@kwdef struct VariableMetadata{L, D}
    lens::L = nothing
    unit::Unit = EMPTY_UNIT
    default::D = missing
    description::String
end
