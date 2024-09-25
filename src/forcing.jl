
@get_units @with_kw struct AtmosphericForcing{T}
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{T}
    # Potential reference evapotranspiration [mm Δt⁻¹]
    potential_evaporation::Vector{T}
    # Temperature [ᵒC]
    temperature::Vector{T} | "°C"
end

function AtmosphericForcing(n)
    atmospheric_forcing = AtmosphericForcing(;
        precipitation = fill(mv, n),
        potential_evaporation = fill(mv, n),
        temperature = fill(mv, n),
    )
    return atmospheric_forcing
end

@get_units @with_kw struct HydrometeoForcing{T}
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{T}
    # Overland flow depth [m]
    waterlevel_land::Vector{T} | "m"
    # Overland flow discharge [m3 s-1]
    q_land::Vector{T} | "m3 s-1"
    # River depth [m]
    waterlevel_river::Vector{T} | "m"
    # River discharge [m3 s-1]
    q_river::Vector{T} | "m3 s-1"
end

function HydrometeoForcing(n)
    hydrometeo_forcing = HydrometeoForcing(;
        precipitation = fill(mv, n),
        waterlevel_land = fill(mv, n),
        q_land = fill(mv, n),
        waterlevel_river = fill(mv, n),
        q_river = fill(mv, n),
    )
    return hydrometeo_forcing
end