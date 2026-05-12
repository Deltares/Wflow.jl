"Struct to store atmospheric forcing variables"
@with_data_lookup struct AtmosphericForcing
    n::Int
    # Precipitation [mm Δt⁻¹]
    "atmosphere_water__precipitation_volume_flux"
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential reference evapotranspiration [mm Δt⁻¹]
    "land_surface_water__potential_evaporation_volume_flux"
    potential_evaporation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Temperature [ᵒC]
    "atmosphere_air__temperature"
    temperature::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store hydrological forcing variables"
@with_data_lookup struct HydrologicalForcing
    n::Int
    # Rainfall interception by the vegetation [mm]
    "vegetation_canopy_water__interception_volume_flux"
    interception::Vector{Float64} = fill(MISSING_VALUE, n)
    # Overland flow depth [m]
    "land_surface_water__depth"
    waterlevel_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # Overland flow discharge [m3 s-1]
    "land_surface_water__volume_flow_rate"
    q_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # River depth [m]
    "river_water__depth"
    waterlevel_river::Vector{Float64} = fill(MISSING_VALUE, n)
    # River discharge [m3 s-1]
    "river_water__volume_flow_rate"
    q_river::Vector{Float64} = fill(MISSING_VALUE, n)
end
