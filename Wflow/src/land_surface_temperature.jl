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
        (actual_evapotranspiration / 1000.0) / Float64(config.time.timestepsecs)
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

Sensitive to canopy density & LAI

Seasonal and structural variability

Forest cap avoids z - d < 0 issues
Alternative aerodynamic conductance calculation from AWRA05
https://www.researchgate.net/publication/233757155_AWRA_Technical_Report_3_Landscape_Model_version_05_Technical_Description
    
"""

function wind_and_aero_resistance(
    wind_speed_measured::Float64,
    z_measured::Float64,
    canopy_height::Float64;
    zm_ref::Float64 = 2.0,
    k::Float64 = 0.41,
)
    # Handle measurement height below canopy
    if z_measured < canopy_height
        z_measured = canopy_height
    end

    # Ensure minimum canopy height using Allen FAO reference
    canopy_height = max(canopy_height, 0.12)

    # Simplified empirical d/h ratios and roughness height adjustments
    if canopy_height < 1
        ref_h = 0.12
        dh_ratio = 2.0 / 3.0
        z0m_ratio = 1.23e-1
        z0h_ratio = 0.1

    elseif canopy_height >= 1
        z0m_ratio = 1.23e-1 * (canopy_height / 2.0) #z0m increases with canopy height
        ref_h = 0.12
        dh_ratio = 2.0 / 3.0
        z0h_ratio = 0.2
    else
        error("Canopy height $canopy_height is not supported")
    end

    d = dh_ratio * ref_h
    z0m = z0m_ratio * ref_h
    z0h = z0h_ratio * z0m  # Canopy sublayer roughness

    # Wind speed conversion to reference height (2m) for consistency with FAO-56
    wind_speed_ref =
        max(wind_speed_measured * (log(zm_ref / z0m) / log(z_measured / z0m)), 0.5)

    if canopy_height < 1
        # Aerodynamic resistance using reference height
        # ra = (log((zm_ref - d) / z0m)) * (log((zm_ref - d) / z0h)) / (k^2 * wind_speed_ref)
        ra = log((zm_ref - d) / z0m) / (k^2 * wind_speed_ref)
    elseif canopy_height >= 1
        #AWRA05
        f_h = log((813 / canopy_height) - 5.45)
        ku = 0.305 / (f_h * (f_h + 2.3))
        ga = ku * wind_speed_ref
        ra = (1 / ga) / 10.0
    end

    # Ensure positive aerodynamic resistance
    if ra <= 0
        error("Aerodynamic resistance is negative: $ra")
    end

    return ra, wind_speed_ref
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
