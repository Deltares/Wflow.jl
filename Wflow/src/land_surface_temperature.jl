# after Devi Purnamasari et al. 2025 
# https://doi.org/10.5194/hess-29-1483-2025

abstract type AbstractLandSurfaceTemperatureModel end
struct NoLandSurfaceTemperatureModel <: AbstractLandSurfaceTemperatureModel end

"Struct for storing LST model variables"
@with_kw struct LandSurfaceTemperatureVariables
    aerodynamic_resistance::Vector{Float64} = Float64[]      # Aerodynamic resistance (s/m)
    latent_heat_flux::Vector{Float64} = Float64[]      # Latent heat flux (W/m2)
    sensible_heat_flux::Vector{Float64} = Float64[]       # Sensible heat flux (W/m2)
    latent_heat_of_vaporization::Vector{Float64} = Float64[]  # Latent heat of vaporization (J/kg)
    land_surface_temperature::Vector{Float64} = Float64[]     # Land surface temperature (K)
    net_radiation::Vector{Float64} = Float64[]               # Net radiation (W/m2)
    net_shortwave_radiation::Vector{Float64} = Float64[]     # Net shortwave radiation (W/m2)
    net_longwave_radiation::Vector{Float64} = Float64[]      # Net longwave radiation (W/m2)
end

"Initialize Land Surface Temperature model variables"
function LandSurfaceTemperatureVariables(n::Int)
    return LandSurfaceTemperatureVariables(;
        aerodynamic_resistance = fill(MISSING_VALUE, n),
        latent_heat_flux = fill(MISSING_VALUE, n),
        sensible_heat_flux = fill(MISSING_VALUE, n),
        latent_heat_of_vaporization = fill(MISSING_VALUE, n),
        land_surface_temperature = fill(MISSING_VALUE, n),
        net_radiation = fill(MISSING_VALUE, n),
        net_shortwave_radiation = fill(MISSING_VALUE, n),
        net_longwave_radiation = fill(MISSING_VALUE, n),
    )
end

@with_kw struct LandSurfaceTemperatureModel <: AbstractLandSurfaceTemperatureModel
    variables::LandSurfaceTemperatureVariables
end

"""
Initialize LandSurfaceTemperatureModel with parameters from Gridded parameters file
"""
function LandSurfaceTemperatureModel(n::Int)
    vars = LandSurfaceTemperatureVariables(n)
    return LandSurfaceTemperatureModel(; variables = vars)
end

"""
Update LST model for a single timestep using pre-calculated net radiation from forcing
"""
function update_land_surface_temperature(
    land_surface_temperature_model::LandSurfaceTemperatureModel,
    soil_model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    vegetation_parameters::VegetationParameters,
    config::Config,
)
    n = length(land_surface_temperature_model.variables.land_surface_temperature)

    # Get wind measurement height from config (default to 10m if not specified)
    wind_measurement_height = Float64(get(config.input, "wind_altitude", 10.0))

    for i in 1:n
        # Use pre-calculated net radiation from forcing (calculated in hydromt_wflow preprocessing)
        land_surface_temperature_model.variables.latent_heat_of_vaporization[i] =
            compute_latent_heat_of_vaporization(atmospheric_forcing.temperature[i])

        land_surface_temperature_model.variables.latent_heat_flux[i] =
            compute_latent_heat_flux(
                atmospheric_forcing.temperature[i],
                soil_model.variables.actevap[i],
                config,
            )

        # Calculate sensible heat flux
        land_surface_temperature_model.variables.sensible_heat_flux[i] =
            compute_sensible_heat_flux(
                atmospheric_forcing.net_radiation[i],
                land_surface_temperature_model.variables.latent_heat_flux[i],
            )

        # Calculate aerodynamic resistance using wind speed at canopy height
        land_surface_temperature_model.variables.aerodynamic_resistance[i] =
            wind_and_aero_resistance(
                atmospheric_forcing.wind_speed[i],
                wind_measurement_height,
                vegetation_parameters.canopy_height[i],
            )

        # Calculate LST
        land_surface_temperature_model.variables.land_surface_temperature[i] =
            compute_land_surface_temperature(
                land_surface_temperature_model.variables.sensible_heat_flux[i],
                land_surface_temperature_model.variables.aerodynamic_resistance[i],
                atmospheric_forcing.temperature[i],
            )
    end

    return nothing
end

"""
Update LST model for a single timestep
"""
function update!(
    land_surface_temperature_model::LandSurfaceTemperatureModel,
    soil_model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    vegetation_parameters::VegetationParameters,
    config::Config,
)
    update_land_surface_temperature(
        land_surface_temperature_model,
        soil_model,
        atmospheric_forcing,
        vegetation_parameters,
        config,
    )
    return nothing
end

function update!(
    model::NoLandSurfaceTemperatureModel,
    soil_model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    vegetation_parameters::VegetationParameters,
    config::Config,
)
    return nothing
end

# wrapper methods
get_LandSurfaceTemperature(model::NoLandSurfaceTemperatureModel) = 0.0
get_LandSurfaceTemperature(model::AbstractLandSurfaceTemperatureModel) =
    model.variables.land_surface_temperature

""" 'latent heat of vaporization' :: λ=2501−2.375Ta.(A1) """
function compute_latent_heat_of_vaporization(air_temperature::Float64)
    return 2501.0 - 2.375 * air_temperature
end
""" 'latent heat flux' :: LE=λ×ρwater×ET,(3)"""
function compute_latent_heat_flux(
    air_temperature::Float64,
    actual_evapotranspiration::Float64,
    config::Config,
)
    latent_heat_of_vaporization = compute_latent_heat_of_vaporization(air_temperature)
    # Convert actual_evapotranspiration from mm/Δt to m/s
    actual_evapotranspiration_ms =
        (actual_evapotranspiration / 1000.0) / config.time.timestepsecs
    latent_heat_flux =
        latent_heat_of_vaporization * density_water * actual_evapotranspiration_ms
    return latent_heat_flux
end
""" 'sensible heat flux' :: H  ≈ RNet - LE """
function compute_sensible_heat_flux(net_radiation::Float64, latent_heat_flux::Float64)
    # Handle NaN values in net radiation
    if isnan(net_radiation)
        @warn "Net radiation is NaN, setting sensible heat flux to 0"
        return 0.0
    end
    sensible_heat_flux = net_radiation - latent_heat_flux
    return sensible_heat_flux
end
""" 
'aerodynamic resistance' :: ra = (ln(z/z0m) - psi_m) / (k^2 * u) 
no clean way yet to deal with variable canopy height empirically

| Cover type        | Typical d/h   | Reference                                    |
| ----------------- | ------------- | -------------------------------------------- |
| Short grass       | 0.67          | Allen et al. (1998), Brutsaert (1982)        |
| Wheat, shrubs     | 0.65          | Brutsaert (1982); Thom (1975); Shuttleworth  |
| Tall crops (corn) | 0.6           | Monteith & Unsworth (1990)                   |
| Forest            | 0.5 (capped)  | Garratt (1992); Shuttleworth & Gurney (1990) |

Empirical averages; real sites vary

Sensitive to canopy density & LAI

Seasonal and structural variability

Forest cap avoids z - d < 0 issues
"""
#EQN: A10-13 adapted with dh ratio to ensure stability
function wind_and_aero_resistance(
    wind_speed_measured::Float64,
    z_measured::Float64,
    canopy_height::Float64;
    z_target::Float64 = 2.0,
    k::Float64 = 0.41,
)
    # Handle zero or negative wind speed
    if wind_speed_measured <= 0
        @warn "Wind speed is zero or negative: $wind_speed_measured, using minimum value"
        wind_speed_measured = 0.5  # Minimum wind speed for stability
    end
    # Empirical d/h ratios
    if canopy_height < 0.2
        dh_ratio = 0.67  # grass
    elseif canopy_height < 2.0
        dh_ratio = 0.65  # wheat/shrubs
    elseif canopy_height < 5.0
        dh_ratio = 0.6   # corn, tall crops
    else
        dh_ratio = 0.5   # forest, capped
    end

    d = dh_ratio * canopy_height
    # Cap d to 0.9 * z_measured for tall canopies
    if d > 0.9 * z_measured
        d = 0.9 * z_measured
    end

    z0m = max(0.005, 0.123 * canopy_height) # in case the canopy height is too small
    z0h = 0.1 * z0m
    zm = canopy_height + z_target
    zh = zm

    # Only use log-profile if both heights are above d
    if (z_measured > d + 0.1) && (zm > d + 0.1) && (zh > d + 0.1)
        wind_speed_canopy =
            wind_speed_measured * (log((zm - d) / z0m) / log((z_measured - d) / z0m))
        ra = (log((zm - d) / z0m) * log((zh - d) / z0h)) / (k^2 * wind_speed_canopy)

        # Check for invalid results
        if isnan(ra) || isinf(ra) || ra <= 0
            @warn "Invalid aerodynamic resistance from log-profile: $ra, using fallback"
            ra = 208 / max(wind_speed_measured, 0.5)  # FAO56 reference value
        end
        return ra
    end

    # Fallback: use measured wind, empirical ra
    @warn "Aerodynamic resistance: log-profile not valid for canopy_height=$canopy_height, z_measured=$z_measured, d=$d. Using measured wind and empirical ra."
    ra = 208 / max(wind_speed_measured, 0.5)  # FAO56 reference value
    return ra
end

""" 'land surface temperature' :: Ts=(H ra) /(ρacp)+Ta,(4)"""
function compute_land_surface_temperature(
    sensible_heat_flux::Float64,
    aerodynamic_resistance::Float64,
    air_temperature::Float64,
)
    # Handle invalid aerodynamic resistance values
    if isnan(aerodynamic_resistance) ||
       isinf(aerodynamic_resistance) ||
       aerodynamic_resistance <= 0
        @warn "Invalid aerodynamic resistance: $aerodynamic_resistance, using air temperature as land surface temperature"
        return air_temperature
    end

    land_surface_temperature =
        (sensible_heat_flux * aerodynamic_resistance) /
        (density_air * specific_heat_capacity_air) + air_temperature
    return land_surface_temperature
end
