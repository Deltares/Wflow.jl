"Struct to store atmospheric forcing variables"
@get_units @with_kw struct AtmosphericForcing{T}
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{T}
    # Potential reference evapotranspiration [mm Δt⁻¹]
    potential_evaporation::Vector{T}
    # Temperature [ᵒC]
    temperature::Vector{T} | "°C"
end

"Initialize atmospheric forcing"
function AtmosphericForcing(
    n;
    precipitation::Vector{T} = fill(mv, n),
    potential_evaporation::Vector{T} = fill(mv, n),
    temperature::Vector{T} = fill(mv, n),
) where {T}
    return AtmosphericForcing{T}(;
        precipitation = precipitation,
        potential_evaporation = potential_evaporation,
        temperature = temperature,
    )
end