"Struct to store atmospheric forcing variables"
@with_kw struct AtmosphericForcing{T}
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{T}
    # Potential reference evapotranspiration [mm Δt⁻¹]
    potential_evaporation::Vector{T}
    # Temperature [ᵒC]
    temperature::Vector{T}
end

"Initialize atmospheric forcing"
function AtmosphericForcing(
    n;
    precipitation::Vector{T} = fill(MISSING_VALUE, n),
    potential_evaporation::Vector{T} = fill(MISSING_VALUE, n),
    temperature::Vector{T} = fill(MISSING_VALUE, n),
) where {T}
    return AtmosphericForcing{T}(; precipitation, potential_evaporation, temperature)
end

"Struct to store hydrological forcing variables"
@with_kw struct HydrologicalForcing{T}
    # Rainfall interception by the vegetation [mm]
    interception::Vector{T}
    # Overland flow depth [m]
    waterlevel_land::Vector{T}
    # Overland flow discharge [m3 s-1]
    q_land::Vector{T}
    # River depth [m]
    waterlevel_river::Vector{T}
    # River discharge [m3 s-1]
    q_river::Vector{T}
end

"Initialize hydrological forcing"
function HydrologicalForcing(
    n;
    interception::Vector{T} = fill(MISSING_VALUE, n),
    waterlevel_land::Vector{T} = fill(MISSING_VALUE, n),
    q_land::Vector{T} = fill(MISSING_VALUE, n),
    waterlevel_river::Vector{T} = fill(MISSING_VALUE, n),
    q_river::Vector{T} = fill(MISSING_VALUE, n),
) where {T}
    return HydrologicalForcing{T}(;
        interception,
        waterlevel_land,
        q_land,
        waterlevel_river,
        q_river,
    )
end