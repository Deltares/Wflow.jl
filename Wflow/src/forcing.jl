"Struct to store atmospheric forcing variables"
@with_kw struct AtmosphericForcing
    n::Int
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential reference evapotranspiration [mm Δt⁻¹]
    potential_evaporation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Temperature [ᵒC]
    temperature::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store hydrological forcing variables"
@with_kw struct HydrologicalForcing
    n::Int
    # Rainfall interception by the vegetation [mm]
    interception::Vector{Float64} = fill(MISSING_VALUE, n)
    # Overland flow depth [m]
    waterlevel_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # Overland flow discharge [m3 s-1]
    q_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # River depth [m]
    waterlevel_river::Vector{Float64} = fill(MISSING_VALUE, n)
    # River discharge [m3 s-1]
    q_river::Vector{Float64} = fill(MISSING_VALUE, n)
end
