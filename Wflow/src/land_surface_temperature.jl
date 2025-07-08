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
Update LST model for a single timestep
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
        # net_radiation is now provided directly in atmospheric_forcing
        net_radiation = atmospheric_forcing.net_radiation[i]

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
                net_radiation,
                land_surface_temperature_model.variables.latent_heat_flux[i],
            )

        # Calculate aerodynamic resistance
        land_surface_temperature_model.variables.aerodynamic_resistance[i] =
            compute_aerodynamic_resistance(
                zm,
                zh,
                vegetation_parameters.canopy_height[i],
                k,
                atmospheric_forcing.wind_speed_2m[i],
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
    sensible_heat_flux = net_radiation - latent_heat_flux
    return sensible_heat_flux
end
""" 'aerodynamic resistance' :: ra = (ln(z/z0m) - psi_m) / (k^2 * u) """
#EQN: A10-13
function compute_aerodynamic_resistance(
    zm::Float64,
    zh::Float64,
    crop_height::Float64,
    von_karman_constant::Float64,
    wind_speed_2m::Float64,
)
    d = (2 / 3) * crop_height
    zom = 0.123 * crop_height
    zoh = 0.1 * zom

    # zm_eff = max(zm, d + 0.1)
    # zh_eff = max(zh, d + 0.1)

    numerator1 = (zm - d) / zom
    numerator2 = (zh - d) / zoh

    if numerator1 <= 0 || numerator2 <= 0
        @warn "Invalid log argument in aerodynamic resistance calc: crop_height=$crop_height, zm_eff=$zm, d=$d, zom=$zom, zoh=$zoh, numerator1=$numerator1, numerator2=$numerator2"
        return 1e6  # return a large resistance, essentially no mixing
    end

    a = log(numerator1) * log(numerator2)
    uz = max(wind_speed_2m, 0.5)
    b = von_karman_constant^2 * uz
    ra = a / b
    return ra
end
""" 'land surface temperature' :: Ts=(H ra) /(ρacp)+Ta,(4)"""
function compute_land_surface_temperature(
    sensible_heat_flux::Float64,
    aerodynamic_resistance::Float64,
    air_temperature::Float64,
)
    land_surface_temperature =
        (sensible_heat_flux * aerodynamic_resistance) / (density_air * cp) + air_temperature
    return land_surface_temperature
end