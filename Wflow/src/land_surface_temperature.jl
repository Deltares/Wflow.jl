# after Devi Purnamasari et al. 2025
# https://doi.org/10.5194/hess-29-1483-2025

abstract type AbstractLandSurfaceTemperatureModel end
struct NoLandSurfaceTemperatureModel <: AbstractLandSurfaceTemperatureModel end

"Struct for storing land surface temperature model variables"
@with_kw struct LandSurfaceTemperatureVariables
    n::Int
    aerodynamic_resistance::Vector{Float64} = fill(MISSING_VALUE, n)        # Aerodynamic resistance (s/m)
    latent_heat_flux::Vector{Float64} = fill(MISSING_VALUE, n)              # Latent heat flux (W/m2)
    sensible_heat_flux::Vector{Float64} = fill(MISSING_VALUE, n)            # Sensible heat flux (W/m2)
    latent_heat_of_vaporization::Vector{Float64} = fill(MISSING_VALUE, n)   # Latent heat of vaporization (J/kg)
    land_surface_temperature::Vector{Float64} = fill(MISSING_VALUE, n)      # Land surface temperature (K)
    net_radiation::Vector{Float64} = fill(MISSING_VALUE, n)                 # Net radiation (W/m2)
    net_shortwave_radiation::Vector{Float64} = fill(MISSING_VALUE, n)       # Net shortwave radiation (W/m2)
    net_longwave_radiation::Vector{Float64} = fill(MISSING_VALUE, n)        # Net longwave radiation (W/m2)
end

@with_kw struct LandSurfaceTemperatureModel <: AbstractLandSurfaceTemperatureModel
    variables::LandSurfaceTemperatureVariables
end

"Initialize land surface temperature model."
function LandSurfaceTemperatureModel(n::Int)
    variables = LandSurfaceTemperatureVariables(; n)
    lst_model = LandSurfaceTemperatureModel(; variables)
    return lst_model
end

"Update land surface temperarure model for a single timestep."
function update_land_surface_temperature!(
    land_surface_temperature_model::LandSurfaceTemperatureModel,
    soil_model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    vegetation_parameters::VegetationParameters,
    wind_measurement_height::Float64,
    dt::Float64,
)
    n = length(land_surface_temperature_model.variables.land_surface_temperature)

    for i in 1:n
        # Use pre-calculated net radiation from forcing
        land_surface_temperature_model.variables.latent_heat_of_vaporization[i] =
            compute_latent_heat_of_vaporization(atmospheric_forcing.temperature[i])

        land_surface_temperature_model.variables.latent_heat_flux[i] =
            compute_latent_heat_flux(
                atmospheric_forcing.temperature[i],
                soil_model.variables.actevap[i],
                dt,
            )

        # Calculate sensible heat flux
        land_surface_temperature_model.variables.sensible_heat_flux[i] =
            compute_sensible_heat_flux(
                atmospheric_forcing.net_radiation[i],
                land_surface_temperature_model.variables.latent_heat_flux[i],
            )

        # Calculate aerodynamic resistance using wind speed at canopy height
        canopy_height = max(vegetation_parameters.canopy_height[i], 0.12)
        land_surface_temperature_model.variables.aerodynamic_resistance[i] =
            wind_and_aero_resistance(
                atmospheric_forcing.wind_speed[i],
                wind_measurement_height,
                canopy_height,
            )

        # Calculate land surface temperature
        land_surface_temperature_model.variables.land_surface_temperature[i] =
            compute_land_surface_temperature(
                land_surface_temperature_model.variables.sensible_heat_flux[i],
                land_surface_temperature_model.variables.aerodynamic_resistance[i],
                atmospheric_forcing.temperature[i],
            )
    end

    return nothing
end

function update_land_surface_temperature!(
    model::NoLandSurfaceTemperatureModel,
    soil_model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    vegetation_parameters::VegetationParameters,
    wind_measurement_height::Float64,
    dt::Float64,
)
    return nothing
end

""" 'latent heat of vaporization' :: λ=2501 - 2.375 Ta (A1) """
function compute_latent_heat_of_vaporization(air_temperature::Float64)
    return (2501.0 - 2.375 * air_temperature) * 1000.0 # J/kg converted fro kj.kg
end

""" 'latent heat flux' :: LE=λ x ρwater x ET (3)"""
function compute_latent_heat_flux(
    air_temperature::Float64,
    actual_evapotranspiration::Float64,
    dt::Float64,
)
    latent_heat_of_vaporization = compute_latent_heat_of_vaporization(air_temperature)
    # Convert actual_evapotranspiration from mm/Δt to m/s
    actual_evapotranspiration_ms = (actual_evapotranspiration / 1000.0) / dt
    latent_heat_flux =
        latent_heat_of_vaporization * WATER_DENSITY * actual_evapotranspiration_ms
    return latent_heat_flux
end

""" 'sensible heat flux' :: H  ≈ RNet - LE - G"""
function compute_sensible_heat_flux(net_radiation::Float64, latent_heat_flux::Float64)
    # Handle NaN values in net radiation
    if isnan(net_radiation)
        return 0.0
    end
    #TODO:run the snow module assimilate soil temperature
    # allowing a better estimate for G, currently G is daytime proportional to (0.1 nighttime, 0.5 daytime)
    G = 0.1 * net_radiation
    sensible_heat_flux = net_radiation - latent_heat_flux - G
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
    zm_ref::Float64 = 2.0, # reference height for wind speed (m)
    k::Float64 = 0.41, # von Kármán constant
)
    # Handle measurement height below canopy
    if z_measured < canopy_height
        z_measured = canopy_height
    end

    # Simplified empirical d/h ratios and roughness height adjustments
    if canopy_height < 1.0
        ref_h = 0.12
        dh_ratio = 2.0 / 3.0
        z0m_ratio = 1.23e-1
        z0h_ratio = 0.1

    elseif canopy_height >= 1.0
        z0m_ratio = 1.23e-1 * (canopy_height / 2.0) #z0m increases with canopy height
        ref_h = 0.33
        dh_ratio = 2.0 / 3.0
        z0h_ratio = 0.2
    end

    # Calculate canopy height and roughness height
    d = dh_ratio * ref_h
    z0m = z0m_ratio * ref_h
    z0h = z0h_ratio * z0m  # Canopy sublayer roughness

    # Wind speed conversion to reference height (2m) for consistency with FAO-56
    wind_speed_ref =
        max(wind_speed_measured * (log(zm_ref / z0m) / log(z_measured / z0m)), 0.5)

    if canopy_height < 1.0
        # Aerodynamic resistance using reference height
        # ra = (log((zm_ref - d) / z0m)) * (log((zm_ref - d) / z0h)) / (k^2 * wind_speed_ref)
        ra = log((zm_ref - d) / z0m) / (k^2 * wind_speed_ref)
    elseif canopy_height >= 1.0
        #AWRA05
        f_h = log((813 / canopy_height) - 5.45)
        ku = 0.305 / (f_h * (f_h + 2.3))
        ga = ku * wind_speed_ref
        ra = 1 / ga
    end

    return max(ra, 10.0)
end

""" 'land surface temperature' :: Ts=(H ra) /(ρacp)+Ta,(4)"""
function compute_land_surface_temperature(
    sensible_heat_flux::Float64,
    aerodynamic_resistance::Float64,
    air_temperature::Float64;
    density_air::Float64 = 1.225,
    specific_heat_capacity_air::Float64 = 1005.0,
)
    land_surface_temperature =
        (sensible_heat_flux * aerodynamic_resistance) /
        (density_air * specific_heat_capacity_air) + air_temperature
    return land_surface_temperature
end
