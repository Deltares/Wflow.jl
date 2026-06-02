"Struct to store atmospheric forcing variables"
@with_data_lookup struct AtmosphericForcing
    n::Int
    # Precipitation [m s⁻¹]
    "atmosphere_water__precipitation_volume_flux"
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential reference evapotranspiration [m s⁻¹]
    "land_surface_water__potential_evaporation_volume_flux"
    potential_evaporation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Temperature [K]
    "atmosphere_air__temperature"
    temperature::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store hydrological forcing variables"
@kwdef struct HydrologicalForcing
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
