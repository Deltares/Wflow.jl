# after Devi Purnamasari et al. 2025 
# https://doi.org/10.5194/hess-29-1483-2025
#TODO: review comments that are too instructional, retain during dev.

abstract type AbstractLandSurfaceTemperatureModel end
struct NoLandSurfaceTemperatureModel <: AbstractLandSurfaceTemperatureModel end

#TODO: move to domain jl landparameters or veg parameters
const sigma = 5.67e-8  # Stefan-Boltzmann constant for blackbody radiation
const cp = 1005.0      # Specific heat capacity of air in J/kg/K
const density_air = 1.225     # Density of air in kg/m3 @ 15 degrees C
const density_water = 1000.0     # Density of water in kg/m3
const k = 0.41        # von Kármán constant
const zm = 2.0        # height of wind speed measurement = zh
const zh = 2.0        # height of air temperature measurement = zm

"Struct for storing LST model variables"
@with_kw struct LandSurfaceTemperatureVariables
    aerodynamic_resistance::Vector{Float64} = Float64[]      # Aerodynamic resistance (s/m)
    net_shortwave_radiation::Vector{Float64} = Float64[]     # Net shortwave radiation (W/m2)
    latent_heat_flux::Vector{Float64} = Float64[]      # Latent heat flux (W/m2)
    sensible_heat_flux::Vector{Float64} = Float64[]       # Sensible heat flux (W/m2)
    net_longwave_radiation::Vector{Float64} = Float64[]     # Net longwave radiation (W/m2)
    latent_heat_of_vaporization::Vector{Float64} = Float64[]  # Latent heat of vaporization (J/kg)
    net_radiation::Vector{Float64} = Float64[]    # Net radiation (W/m2)
    land_surface_temperature::Vector{Float64} = Float64[]     # Land surface temperature (K)
end

"Initialize Land Surface Temperature model variables"
function LandSurfaceTemperatureVariables(n::Int)
    return LandSurfaceTemperatureVariables(;
        aerodynamic_resistance = fill(MISSING_VALUE, n),
        net_shortwave_radiation = fill(MISSING_VALUE, n),
        latent_heat_flux = fill(MISSING_VALUE, n),
        sensible_heat_flux = fill(MISSING_VALUE, n),
        net_longwave_radiation = fill(MISSING_VALUE, n),
        latent_heat_of_vaporization = fill(MISSING_VALUE, n),
        net_radiation = fill(MISSING_VALUE, n),
        land_surface_temperature = fill(MISSING_VALUE, n),
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
    network::NetworkLand,
    vegetation_parameters::VegetationParameters,
    land_parameters::LandParameters,
    clock,
)
    n = length(land_surface_temperature_model.variables.land_surface_temperature)

    for i in 1:n
        land_surface_temperature_model.variables.net_shortwave_radiation[i] =
            compute_net_shortwave_radiation(
                land_parameters.albedo[i],
                atmospheric_forcing.shortwave_radiation_in[i],
            )

        land_surface_temperature_model.variables.net_longwave_radiation[i] =
            compute_net_longwave_radiation(
                atmospheric_forcing.temperature[i],
                atmospheric_forcing.shortwave_radiation_in[i],
                network.latitude[i],
                clock,
            )
        land_surface_temperature_model.variables.net_radiation[i] =
            land_surface_temperature_model.variables.net_shortwave_radiation[i] -
            land_surface_temperature_model.variables.net_longwave_radiation[i]

        land_surface_temperature_model.variables.latent_heat_of_vaporization[i] =
            compute_latent_heat_of_vaporization(atmospheric_forcing.temperature[i])

        land_surface_temperature_model.variables.latent_heat_flux[i] =
            compute_latent_heat_flux(
                atmospheric_forcing.temperature[i],
                soil_model.variables.actevap[i],
                clock,
            )

        # Calculate sensible heat flux
        land_surface_temperature_model.variables.sensible_heat_flux[i] =
            compute_sensible_heat_flux(
                land_surface_temperature_model.variables.net_radiation[i],
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
    network::NetworkLand,
    vegetation_parameters::VegetationParameters,
    land_parameters::LandParameters,
    clock,
)
    update_land_surface_temperature(
        land_surface_temperature_model,
        soil_model,
        atmospheric_forcing,
        network,
        vegetation_parameters,
        land_parameters,
        clock,
    )

    return nothing
end

function update!(
    model::NoLandSurfaceTemperatureModel,
    soil_model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    network::NetworkLand,
    vegetation_parameters::VegetationParameters,
    land_parameters::LandParameters,
    clock,
)
    return nothing
end

# wrapper methods
get_LandSurfaceTemperature(model::NoLandSurfaceTemperatureModel) = 0.0
get_LandSurfaceTemperature(model::AbstractLandSurfaceTemperatureModel) =
    model.variables.land_surface_temperature

""" 'net shortwave radiation' :: RSNet=(1−α)Rins(A4) """
#EQN: A4
function compute_net_shortwave_radiation(albedo::Float64, shortwave_radiation_in::Float64)
    net_shortwave_radiation = (1 - albedo) * shortwave_radiation_in
    return net_shortwave_radiation
end
#EQN: A8
function solar_declination(doy::Int)
    days_in_year = 365
    radians = sin((2π * doy / days_in_year) - 1.39)
    decl = 0.409 * radians
    return decl
end
#EQN: A8
function relative_distance(doy::Int)
    days_in_year = 365
    radians = cos(2π * doy / days_in_year)
    dist = radians * 0.033 + 1  # distance in AU
    return dist
end
#EQN: A8
function extraterrestrial_radiation(lat::Float64, clock)
    gsc = 118.08  # MJ/m²/day
    lat_rad = deg2rad(lat)
    doy = dayofyear(clock.time)
    decl = solar_declination(doy)
    dist = relative_distance(doy)
    sha = acos(-tan(lat_rad) * tan(decl))
    Ra =
        dist * gsc / π *
        (cos(lat_rad) * cos(decl) * sin(sha) + sha * sin(lat_rad) * sin(decl))
    return Ra
end
""" 'net longwave radiation' :: RLN=(σTa^4)(0.34−0.14−√ea)(1.35 (Rins/Rso) − 0.35) """
function compute_net_longwave_radiation(
    air_temperature::Float64,
    shortwave_radiation_in::Float64,
    lat::Float64,
    clock,
)
    sigma = 4.903e-9 # MJK-4m-2day-1 Steffan Boltzmann constant 
    a = sigma * (air_temperature + 273.15)^4
    ea = 0.611 * exp(17.27 * air_temperature / (237.3 + air_temperature)^2)
    b = 0.34 - 0.14 * sqrt(ea)
    Rso = extraterrestrial_radiation(lat, clock) * 0.75
    c = (1.35 * shortwave_radiation_in / Rso) - 0.35
    RLN = a * b * c
    return RLN
end
""" 'net radiation' :: Rn=Rins−Routs+Rinl+Rins, """
function compute_net_radiation(
    albedo::Float64,
    shortwave_radiation_in::Float64,
    lat::Float64,
    air_temperature::Float64,
    clock,
)
    net_shortwave_radiation =
        compute_net_shortwave_radiation(albedo, shortwave_radiation_in)
    net_longwave_radiation =
        compute_net_longwave_radiation(air_temperature, shortwave_radiation_in, lat, clock)
    net_radiation = net_shortwave_radiation - net_longwave_radiation
    return net_radiation
end
""" 'latent heat of vaporization' :: λ=2501−2.375Ta.(A1) """
function compute_latent_heat_of_vaporization(air_temperature::Float64)
    return 2501.0 - 2.375 * air_temperature
end
""" 'latent heat flux' :: LE=λ×ρwater×ET,(3)"""
function compute_latent_heat_flux(
    air_temperature::Float64,
    actual_evapotranspiration::Float64,
    clock,
)
    latent_heat_of_vaporization = compute_latent_heat_of_vaporization(air_temperature)
    # Convert actual_evapotranspiration from mm/Δt to m/s
    actual_evapotranspiration_ms = (actual_evapotranspiration / 1000.0) / tosecond(clock.dt)
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