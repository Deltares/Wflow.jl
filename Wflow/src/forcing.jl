"Struct to store atmospheric forcing variables"
@with_kw struct AtmosphericForcing
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{Float64} = Float64[]
    # Potential reference evapotranspiration [mm Δt⁻¹]
    potential_evaporation::Vector{Float64} = Float64[]
    # Temperature [ᵒC]
    temperature::Vector{Float64} = Float64[]
    # Downward shortwave radiation [W m-2]
    shortwave_radiation_in::Vector{Float64} = Float64[]
    # Wind speed [m s-1]
    wind_speed::Vector{Float64} = Float64[]
    # Net radiation [W m-2] (pre-calculated in hydromt_wflow preprocessing)
    net_radiation::Vector{Float64} = Float64[]
end

"Initialize atmospheric forcing"
function AtmosphericForcing(config::Config, n::Int)
    do_land_surface_temperature =
        get(config.model, "land_surface_temperature__flag", false)::Bool

    atmos_forcing = AtmosphericForcing(;
        precipitation = fill(MISSING_VALUE, n),
        potential_evaporation = fill(MISSING_VALUE, n),
        temperature = fill(MISSING_VALUE, n),
    )

    if do_land_surface_temperature
        @reset atmos_forcing.shortwave_radiation_in = fill(MISSING_VALUE, n)
        @reset atmos_forcing.wind_speed = fill(MISSING_VALUE, n)
        @reset atmos_forcing.net_radiation = fill(MISSING_VALUE, n)
    end

    return atmos_forcing
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