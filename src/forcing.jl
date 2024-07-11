
@get_units @with_kw struct AtmosphericForcing{T}
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{T}
    # Potential reference evapotranspiration [mm Δt⁻¹]
    potential_evaporation::Vector{T}
    # Temperature [ᵒC]
    temperature::Vector{T} | "°C"
end

function initialize_atmospheric_forcing(n)
    atmospheric_forcing = AtmosphericForcing(
        precipitation = fill(mv, n),
        potential_evaporation = fill(mv, n),
        temperature = fill(mv, n),
    )
    return atmospheric_forcing
end
