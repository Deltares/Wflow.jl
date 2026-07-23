# after Devi Purnamasari et al. 2025
# https://doi.org/10.5194/hess-29-1483-2025

const VON_KARMAN = 0.41

abstract type AbstractLandSurfaceTemperatureModel end
struct NoLandSurfaceTemperatureModel <: AbstractLandSurfaceTemperatureModel end

"Struct for storing land surface temperature model variables"
@with_kw struct LandSurfaceTemperatureVariables
    n::Int
    # Aerodynamic resistance (s/m)
    aerodynamic_resistance::Vector{Float64} = fill(MISSING_VALUE, n)
    # Latent heat flux (W/m2)
    latent_heat_flux::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sensible heat flux (W/m2)
    sensible_heat_flux::Vector{Float64} = fill(MISSING_VALUE, n)
    # Latent heat of vaporization (J/kg)
    latent_heat_of_vaporization::Vector{Float64} = fill(MISSING_VALUE, n)
    # Land surface temperature (K)
    land_surface_temperature::Vector{Float64} = fill(MISSING_VALUE, n)
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

"Update land surface temperature model for a single timestep."
function update_land_surface_temperature!(
    land_surface_temperature_model::LandSurfaceTemperatureModel,
    soil_model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    parameters::LandParameters,
    wind_measurement_height::Float64,
    dt::Float64,
)
    (; d0, z0m, z0h, skin_layer_height) = parameters
    n = length(land_surface_temperature_model.variables.land_surface_temperature)

    for i in 1:n
        # Use pre-calculated net radiation from forcing
        land_surface_temperature_model.variables.latent_heat_of_vaporization[i] =
            compute_latent_heat_of_vaporization(atmospheric_forcing.temperature[i])

        land_surface_temperature_model.variables.latent_heat_flux[i] =
            compute_latent_heat_flux(
                atmospheric_forcing.temperature[i],
                soil_model.variables.actevap[i],
            )

        # Calculate sensible heat flux
        land_surface_temperature_model.variables.sensible_heat_flux[i] =
            compute_sensible_heat_flux(
                atmospheric_forcing.net_radiation[i],
                land_surface_temperature_model.variables.latent_heat_flux[i],
            )

        # Calculate aerodynamic resistance using wind speed at canopy height
        land_surface_temperature_model.variables.aerodynamic_resistance[i] =
            compute_aerodynamic_resistance(
                atmospheric_forcing.wind_speed[i],
                wind_measurement_height,
                skin_layer_height[i],
                d0[i],
                z0m[i],
                z0h[i],
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
    parameters::LandParameters,
    wind_measurement_height::Float64,
    dt::Float64,
)
    return nothing
end

"Compute latent heat of vaporization"
function compute_latent_heat_of_vaporization(air_temperature::Float64)
    air_temperature_degc = from_SI(air_temperature, ABSOLUTE_DEGREES)
    return (2501.0 - 2.375 * air_temperature_degc) * 1000.0 # J/kg converted from kJ/kg
end

"Compute latent heat flux"
function compute_latent_heat_flux(
    air_temperature::Float64,
    actual_evapotranspiration::Float64,
)
    latent_heat_of_vaporization = compute_latent_heat_of_vaporization(air_temperature)
    latent_heat_flux =
        latent_heat_of_vaporization * WATER_DENSITY * actual_evapotranspiration
    return latent_heat_flux
end

"Compute sensible heat flux"
function compute_sensible_heat_flux(net_radiation::Float64, latent_heat_flux::Float64)
    # estimation of ground heat flux can be improved, currently estimated with a fixed
    # linear coefficient (0.1) applied to net_radiation
    ground_heat_flux = 0.1 * net_radiation
    sensible_heat_flux = net_radiation - latent_heat_flux - ground_heat_flux
    return sensible_heat_flux
end

"""
Compute aerodynamic resistance `ra` of the land surface "skin" (vegetation canopy, soil,
snow etc.) based on Thom's equation assuming neutral stability conditions.

The provided wind speed should be consistent with the aerodynamic roughness length for
momentum transfer `z0m` and zero-displacement height `d0` of the surface terrain (and not
based on "open-terrain" (short grass) as observed typically by meteorological stations). The
wind speed measurement height `z_measured` should be greater than the sum of `d0` and `z0m`
to convert the wind speed to the reference height `zm_ref` properly. The measurement heights
of wind and humidity are assumed to be equal, `z0h  is the aerodynamic roughness length for
heat transfer.
"""
function compute_aerodynamic_resistance(
    wind_speed::Float64,
    z_measured::Float64,
    skin_layer_height::Float64,
    d0::Float64,
    z0m::Float64,
    z0h::Float64,
)
    # set reference height (~2.0 m above surface skin layer height).
    zm_ref = round(skin_layer_height + 2.0)

    # wind speed at reference height and constrained to be greater than 0.5 m/s to consider
    # vapour exchange on the surface induced by air buoyancy and layer instability effects.
    min_wind_speed = 0.5
    wind_speed_ref = max(
        wind_speed * (log((zm_ref-d0) / z0m) / log((z_measured-d0) / z0m)),
        min_wind_speed,
    )

    # compute aerodynamic resistance based on Thom's equation.
    ra = (log((zm_ref-d0)/z0m)*log((zm_ref-d0)/z0h)) / (VON_KARMAN^2 * wind_speed_ref)

    return ra
end

"Compute land surface temperature"
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
