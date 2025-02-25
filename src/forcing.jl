"Struct to store atmospheric forcing variables"
@with_kw struct AtmosphericForcing
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{Float64}
    # Potential reference evapotranspiration [mm Δt⁻¹]
    potential_evaporation::Vector{Float64}
    # Temperature [ᵒC]
    temperature::Vector{Float64}
end

"Initialize atmospheric forcing"
function AtmosphericForcing(
    n::Int;
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n),
    potential_evaporation::Vector{Float64} = fill(MISSING_VALUE, n),
    temperature::Vector{Float64} = fill(MISSING_VALUE, n),
)
    return AtmosphericForcing(; precipitation, potential_evaporation, temperature)
end

"Struct to store hydrological forcing variables"
@with_kw struct HydrologicalForcing
    # Rainfall interception by the vegetation [mm]
    interception::Vector{Float64}
    # Overland flow depth [m]
    waterlevel_land::Vector{Float64}
    # Overland flow discharge [m3 s-1]
    q_land::Vector{Float64}
    # River depth [m]
    waterlevel_river::Vector{Float64}
    # River discharge [m3 s-1]
    q_river::Vector{Float64}
end

"Initialize hydrological forcing"
function HydrologicalForcing(
    n::Int;
    interception::Vector{Float64} = fill(MISSING_VALUE, n),
    waterlevel_land::Vector{Float64} = fill(MISSING_VALUE, n),
    q_land::Vector{Float64} = fill(MISSING_VALUE, n),
    waterlevel_river::Vector{Float64} = fill(MISSING_VALUE, n),
    q_river::Vector{Float64} = fill(MISSING_VALUE, n),
)
    return HydrologicalForcing(;
        interception,
        waterlevel_land,
        q_land,
        waterlevel_river,
        q_river,
    )
end