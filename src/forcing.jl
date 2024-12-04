"Struct to store atmospheric forcing variables"
@get_units @grid_loc @with_kw struct AtmosphericForcing{T}
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
    return AtmosphericForcing{T}(; precipitation, potential_evaporation, temperature)
end

"Struct to store hydrological forcing variables"
@get_units @grid_loc @with_kw struct HydrologicalForcing{T}
    # Overland flow depth [m]
    waterlevel_land::Vector{T} | "m"
    # Overland flow discharge [m3 s-1]
    q_land::Vector{T} | "m3 s-1"
    # River depth [m]
    waterlevel_river::Vector{T} | "m"
    # River discharge [m3 s-1]
    q_river::Vector{T} | "m3 s-1"
end

"Initialize hydrological forcing"
function HydrologicalForcing(
    n;
    waterlevel_land::Vector{T} = fill(mv, n),
    q_land::Vector{T} = fill(mv, n),
    waterlevel_river::Vector{T} = fill(mv, n),
    q_river::Vector{T} = fill(mv, n),
) where {T}
    return HydrologicalForcing{T}(; waterlevel_land, q_land, waterlevel_river, q_river)
end