# after Devi Purnamasari et al. 2025 
# https://doi.org/10.5194/hess-29-1483-2025
#TODO: review comments that are too instructional, retain during dev.

using Wflow
using Dates
import Wflow: AbstractLandModel, AtmosphericForcing, SbmSoilModel, Config
using NCDatasets
using Parameters

#abstract should define a common interface
abstract type AbstractLandSurfaceTemperatureModel end
#concrete type helps determine if LSTmodel is active or not
struct NoLandSurfaceTemperatureModel <: AbstractLandSurfaceTemperatureModel end

#TODO: move to domain jl landparameters or veg parameters
const sigma = 5.67e-8  # Stefan-Boltzmann constant for blackbody radiation
const cp = 1005.0      # Specific heat capacity of air in J/kg/K
const rho_a = 1.225     # Density of air in kg/m3 @ 15 degrees C
const rho_w = 1000.0     # Density of water in kg/m3
const k = 0.41        # von Kármán constant
#HYRAS input is 10m we correct to 2m in calculation for aerodynamic resistance
const zm = 2.0        # height of wind speed measurement = zh
const zh = 2.0        # height of air temperature measurement = zm
const dt = 86400.0     # timestep in seconds

"Struct for storing LST model variables"
@with_kw struct LandSurfaceTemperatureVariables
    # Input forcing (updated each timestep)
    RS_in::Vector{Float64}   # Downward shortwave radiation (W/m2)
    Ta::Vector{Float64}      # Air temperature (K)
    u2m::Vector{Float64}     # Wind speed at 2m height (m/s)
    albedo::Vector{Float64}  # Surface albedo
    #TODO: emissivity is a constant, not a variable
    # emissivity::Vector{Float64} # Surface emissivity

    #From LandHydrologySBM
    ET_a::Vector{Float64}    # Actual evapotranspiration from SBM (mm/Δt)

    # Calculated variables
    ra::Vector{Float64}      # Aerodynamic resistance (s/m)
    RSN::Vector{Float64}     # Net shortwave radiation (W/m2)
    LE::Vector{Float64}      # Latent heat flux (W/m2)
    H::Vector{Float64}       # Sensible heat flux (W/m2)
    RLN::Vector{Float64}     # Net longwave radiation (W/m2)
    lambda::Vector{Float64}  # Latent heat of vaporization (J/kg)
    Rnet::Vector{Float64}    # Net radiation (W/m2)
    LST::Vector{Float64}     # Land surface temperature (K)
end

"Initialize Land Surface Temperature model variables"
function LandSurfaceTemperatureVariables(n::Int)
    return LandSurfaceTemperatureVariables(;
        RS_in = fill(MISSING_VALUE, n),
        Ta = fill(MISSING_VALUE, n),
        u2m = fill(MISSING_VALUE, n),
        albedo = fill(MISSING_VALUE, n),
        #TODO: emissivity is a constant, not a variable
        # emissivity = fill(MISSING_VALUE, n),
        ET_a = fill(MISSING_VALUE, n),
        ra = fill(MISSING_VALUE, n),
        RSN = fill(MISSING_VALUE, n),
        LE = fill(MISSING_VALUE, n),
        H = fill(MISSING_VALUE, n),
        RLN = fill(MISSING_VALUE, n),
        lambda = fill(MISSING_VALUE, n),
        Rnet = fill(MISSING_VALUE, n),
        LST = fill(MISSING_VALUE, n),
    )
end

"Struct for storing LST model parameters"
@with_kw struct LandSurfaceTemperatureParameters
    # Static parameters (from netCDF or defaults)
    crop_height::Vector{Float64} # Crop height in m
    latitude::Vector{Float64}  # Latitude of each grid cell (degrees)
end

"Initialize LST model parameters"
function LandSurfaceTemperatureParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    lens = lens_input_parameter(config, "crop_height")
    crop_height =
        ncread(dataset, config, lens; sel = indices, defaults = 2.0, type = Float64)

    latitude = ncread(dataset, "lat"; sel = indices, type = Float64)

    return LandSurfaceTemperatureParameters(;
        crop_height = crop_height,
        latitude = latitude,
    )
end

@with_kw struct LandSurfaceTemperatureModel <: AbstractLandSurfaceTemperatureModel
    parameters::LandSurfaceTemperatureParameters
    variables::LandSurfaceTemperatureVariables
end

"""
Update actual evapotranspiration from SBM model
"""
function update_ET_a!(
    land_surface_temperature_model::LandSurfaceTemperatureModel,
    soil_model::SbmSoilModel,
)
    land_surface_temperature_model.variables.ET_a .= soil_model.variables.actevap
    return nothing
end

"""
Update forcing inputs for LST model
"""
function update_forcing!(
    land_surface_temperature_model::LandSurfaceTemperatureModel,
    atmospheric_forcing::AtmosphericForcing,
)
    land_surface_temperature_model.variables.Ta .= atmospheric_forcing.temperature

    land_surface_temperature_model.variables.RS_in .= atmospheric_forcing.RS_in

    land_surface_temperature_model.variables.albedo .= atmospheric_forcing.albedo

    land_surface_temperature_model.variables.u2m .= atmospheric_forcing.u2m

    return nothing
end

"""
Initialize LandSurfaceTemperatureModel with parameters from Gridded parameters file
"""
function LandSurfaceTemperatureModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    n = length(indices)

    params = LandSurfaceTemperatureParameters(dataset, config, indices)
    vars = LandSurfaceTemperatureVariables(n)

    return LandSurfaceTemperatureModel(; parameters = params, variables = vars)
end

"""
Update LST model for a single timestep
"""
function update_timestep_land_surface_temperature(
    land_surface_temperature_model::LandSurfaceTemperatureModel,
    soil_model::SbmSoilModel,
    current_time::DateTime,
    dt::Float64,
)
    n = length(land_surface_temperature_model.variables.LST)

    # Update ET_a from SBM model
    update_ET_a!(land_surface_temperature_model, soil_model)

    # Calculate day of year
    doy = dayofyear(current_time)

    # Calculate components for each grid cell
    for i in 1:n
        # Calculate net radiation
        land_surface_temperature_model.variables.RSN[i] = calculate_RSN(
            land_surface_temperature_model.variables.albedo[i],
            land_surface_temperature_model.variables.RS_in[i],
        )

        land_surface_temperature_model.variables.RLN[i] = calculate_RLN(
            land_surface_temperature_model.variables.Ta[i],
            land_surface_temperature_model.variables.RS_in[i],
            land_surface_temperature_model.parameters.latitude[i],
            doy,
        )
        land_surface_temperature_model.variables.Rnet[i] =
            land_surface_temperature_model.variables.RSN[i] -
            land_surface_temperature_model.variables.RLN[i]

        # Calculate latent heat flux using ET_a from SBM
        land_surface_temperature_model.variables.lambda[i] =
            calculate_lambda(land_surface_temperature_model.variables.Ta[i])
        land_surface_temperature_model.variables.LE[i] = calculate_LE(
            land_surface_temperature_model.variables.Ta[i],
            rho_w,
            land_surface_temperature_model.variables.ET_a[i],
            dt,
        )

        # Calculate sensible heat flux
        land_surface_temperature_model.variables.H[i] = calculate_H(
            land_surface_temperature_model.variables.Rnet[i],
            land_surface_temperature_model.variables.LE[i],
        )

        # Calculate LST
        land_surface_temperature_model.variables.LST[i] = calculate_LST(
            land_surface_temperature_model.variables.H[i],
            land_surface_temperature_model.variables.ra[i],
            rho_a,
            cp,
            land_surface_temperature_model.variables.Ta[i],
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
    current_time::DateTime,
    dt::Float64,
)
    # Update forcing inputs
    update_forcing!(land_surface_temperature_model, atmospheric_forcing)

    update_timestep_land_surface_temperature(
        land_surface_temperature_model,
        soil_model,
        current_time,
        dt,
    )

    return nothing
end

function update!(
    model::NoLandSurfaceTemperatureModel,
    soil_model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    current_time::DateTime,
    dt::Float64,
)
    return nothing
end

# wrapper methods
get_LandSurfaceTemperature(model::NoLandSurfaceTemperatureModel) = 0.0
get_LandSurfaceTemperature(model::AbstractLandSurfaceTemperatureModel) = model.variables.LST

""" 'net shortwave radiation' :: RSNet=(1−α)Rins(A4) """
#EQN: A4
function calculate_RSN(albedo::Float64, RS_in::Float64)
    RSN = (1 - albedo) * RS_in
    return RSN
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
function extraterrestrial_radiation(lat::Float64, doy::Int)
    gsc = 118.08  # MJ/m²/day
    lat_rad = deg2rad(lat)
    decl = solar_declination(doy)
    dist = relative_distance(doy)
    sha = acos(-tan(lat_rad) * tan(decl))
    Ra =
        dist * gsc / π *
        (cos(lat_rad) * cos(decl) * sin(sha) + sha * sin(lat_rad) * sin(decl))
    return Ra
end
""" 'net longwave radiation' :: RLN=(σTa^4)(0.34−0.14−√ea)(1.35 (Rins/Rso) − 0.35) """
function calculate_RLN(Ta::Float64, RS_in::Float64, lat::Float64, doy::Int)
    sigma = 4.903e-9 # MJK-4m-2day-1 Steffan Boltzmann constant 
    a = sigma * (Ta + 273.15)^4
    ea = 0.611 * exp(17.27 * Ta / (237.3 + Ta)^2)
    b = 0.34 - 0.14 * sqrt(ea)
    Rso = extraterrestrial_radiation(lat, doy) * 0.75
    c = (1.35 * RS_in / Rso) - 0.35
    RLN = a * b * c
    return RLN
end
""" 'net radiation' :: Rn=Rins−Routs+Rinl+Rins, """
function calculate_Rnet(
    albedo::Float64,
    RS_in::Float64,
    lat::Float64,
    doy::Int,
    Ta::Float64,
)
    RSN = calculate_RSN(albedo, RS_in)
    RLN = calculate_RLN(Ta, RS_in, lat, doy)
    Rnet = RSN - RLN
    return Rnet
end
""" 'latent heat of vaporization' :: λ=2501−2.375Ta.(A1) """
function calculate_lambda(Ta::Float64)
    return 2501.0 - 2.375 * Ta
end
""" 'latent heat flux' :: LE=λ×ρwater×ET,(3)"""
function calculate_LE(Ta::Float64, rho_water::Float64, ET_a::Float64, dt::Float64)
    lambda = calculate_lambda(Ta)
    # Convert ET_a from mm/Δt to m/s
    ET_a_ms = (ET_a / 1000.0) / dt
    LE = lambda * rho_water * ET_a_ms
    return LE
end
""" 'sensible heat flux' :: H  ≈ RNet - LE """
function calculate_H(Rnet::Float64, LE::Float64)
    H = Rnet - LE
    return H
end
#LIMITATION: this is only valid over short grass 
#u2 = uz * 4 . 87 / ln ( 67 . 8 * z − 5 . 42 ) See FAO Irrigation and Drainage Paper No. 56 https://www.researchgate.net/post/Is-that-possible-to-convert-wind-speed-measured-in-10-m-height-to-a-possible-2-m-height-wind-speed/59bf8c4796b7e406f877880c
""" 'wind speed' :: uz (Allen et al., 1998)"""
function calculate_uz(u10::Float64, v10::Float64)
    u10m = sqrt(u10^2 + v10^2)
    u2m = u10m * 4.87 / log((67.8 * 10 - 5.42))
    return u2m
end

""" 'aerodynamic resistance' :: ra = (ln(z/z0m) - psi_m) / (k^2 * u) """
#EQN: A10-13
function calculate_ra(zm::Float64, zh::Float64, hc::Float64, k::Float64, uz::Float64)
    d = (2 / 3) * hc      #zero-plane displacement height
    zom = 0.123 * hc     #roughness length for momentum
    zoh = 0.1 * zom      #roughness length for heat
    a = log((zm - d) / zom) * log((zh - d) / zoh)
    uz = max(uz, 0.5) # must ensure some mixing in low wind scenarios, set a floor of 0.5 m/s
    b = k^2 * uz
    ra = a / b
    return ra
end
""" 'land surface temperature' :: Ts=(H ra) /(ρacp)+Ta,(4)"""
function calculate_LST(H::Float64, Ra::Float64, rho_a::Float64, cp::Float64, Ta::Float64)
    LST = (H * Ra) / (rho_a * cp) + Ta
    return LST
end