"Struct to store atmospheric forcing variables"
@with_kw struct AtmosphericForcing
    n::Int
    # Precipitation [m s⁻¹]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential reference evapotranspiration [m s⁻¹]
    potential_evaporation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Temperature [K]
    temperature::Vector{Float64} = fill(MISSING_VALUE, n)
    # Downward shortwave radiation [W m-2]
    shortwave_radiation_in::Vector{Float64} = Float64[]
    # Wind speed [m s-1]
    wind_speed::Vector{Float64} = Float64[]
    # Net radiation [W m-2]
    net_radiation::Vector{Float64} = Float64[]
end

"Struct to store hydrological forcing variables"
@with_kw struct HydrologicalForcing
    n::Int
    # Rainfall interception by the vegetation [m s⁻¹]
    interception::Vector{Float64} = fill(MISSING_VALUE, n)
    # Overland flow depth [m]
    waterlevel_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # Overland flow discharge [m³ s⁻¹]
    q_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # River depth [m]
    waterlevel_river::Vector{Float64} = fill(MISSING_VALUE, n)
    # River discharge [m³ s⁻¹]
    q_river::Vector{Float64} = fill(MISSING_VALUE, n)
end
