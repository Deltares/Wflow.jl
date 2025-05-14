abstract type AbstractSnowModel end

"Struct for storing snow model variables"
@with_kw struct SnowVariables
    # Snow storage [mm]
    snow_storage::Vector{Float}
    # Liquid water content in the snow pack [mm]
    snow_water::Vector{Float}
    # Snow water equivalent (SWE) [mm]
    swe::Vector{Float}
    # Snow melt [mm Δt⁻¹]
    snow_melt::Vector{Float}
    # Runoff from snowpack [mm Δt⁻¹]
    runoff::Vector{Float}
end

"Initialize snow model variables"
function SnowVariables(n::Int)
    return SnowVariables(;
        snow_storage = fill(0.0, n),
        snow_water = fill(0.0, n),
        swe = fill(MISSING_VALUE, n),
        runoff = fill(MISSING_VALUE, n),
        snow_melt = fill(MISSING_VALUE, n),
    )
end

"Struct for storing snow model boundary conditions"
@with_kw struct SnowBC
    # Effective precipitation [mm Δt⁻¹]
    effective_precip::Vector{Float}
    # Snow precipitation [mm Δt⁻¹]
    snow_precip::Vector{Float}
    # Liquid precipitation [mm Δt⁻¹]
    liquid_precip::Vector{Float}
end

"Initialize snow model boundary conditions"
function SnowBC(n::Int)
    return SnowBC(;
        effective_precip = fill(MISSING_VALUE, n),
        snow_precip = fill(MISSING_VALUE, n),
        liquid_precip = fill(MISSING_VALUE, n),
    )
end

"Struct for storing snow HBV model parameters"
@with_kw struct SnowHbvParameters
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹]
    cfmax::Vector{Float}
    # Threshold temperature for snowfall [ᵒC]
    tt::Vector{Float}
    # Threshold temperature interval length [ᵒC]
    tti::Vector{Float}
    # Threshold temperature for snowmelt [ᵒC]
    ttm::Vector{Float}
    # Water holding capacity as fraction of current snow pack [-]
    whc::Vector{Float}
end

"Snow HBV model"
@with_kw struct SnowHbvModel <: AbstractSnowModel
    boundary_conditions::SnowBC
    parameters::SnowHbvParameters
    variables::SnowVariables
end

struct NoSnowModel <: AbstractSnowModel end

"Initialize snow HBV model parameters"
function SnowHbvParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    lens = lens_input_parameter(config, "snowpack__degree-day_coefficient")
    cfmax =
        ncread(dataset, config, lens; sel = indices, defaults = 3.75, type = Float) .*
        (dt / BASETIMESTEP)
    lens = lens_input_parameter(config, "atmosphere_air__snowfall_temperature_threshold")
    tt = ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float)

    lens = lens_input_parameter(config, "atmosphere_air__snowfall_temperature_interval")
    tti = ncread(dataset, config, lens; sel = indices, defaults = 1.0, type = Float)

    lens = lens_input_parameter(config, "snowpack__melting_temperature_threshold")
    ttm = ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float)

    lens = lens_input_parameter(config, "snowpack__liquid_water_holding_capacity")
    whc = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)
    snow_hbv_params =
        SnowHbvParameters(; cfmax = cfmax, tt = tt, tti = tti, ttm = ttm, whc = whc)
    return snow_hbv_params
end

"Initialize snow HBV model"
function SnowHbvModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    n = Int(length(indices))
    params = SnowHbvParameters(dataset, config, indices, dt)
    vars = SnowVariables(n)
    bc = SnowBC(n)
    model = SnowHbvModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

"Update boundary condition (effective precipitation provided by an interception model) of a snow model for a single timestep"
function update_boundary_conditions!(model::AbstractSnowModel, external_models::NamedTuple)
    (; effective_precip) = model.boundary_conditions
    (; interception) = external_models
    @. effective_precip =
        interception.variables.throughfall + interception.variables.stemflow
    return nothing
end

function update_boundary_conditions!(model::NoSnowModel, external_models::NamedTuple)
    return nothing
end

"Update snow HBV model for a single timestep"
function update!(model::SnowHbvModel, atmospheric_forcing::AtmosphericForcing)
    (; temperature) = atmospheric_forcing
    (; snow_storage, snow_water, swe, snow_melt, runoff) = model.variables
    (; effective_precip, snow_precip, liquid_precip) = model.boundary_conditions
    (; tt, tti, ttm, cfmax, whc) = model.parameters

    n = length(temperature)
    threaded_foreach(1:n; basesize = 1000) do i
        snow_precip[i], liquid_precip[i] =
            precipitation_hbv(effective_precip[i], temperature[i], tti[i], tt[i])
    end
    threaded_foreach(1:n; basesize = 1000) do i
        snow_storage[i], snow_water[i], swe[i], snow_melt[i], runoff[i] = snowpack_hbv(
            snow_storage[i],
            snow_water[i],
            snow_precip[i],
            liquid_precip[i],
            temperature[i],
            ttm[i],
            cfmax[i],
            whc[i],
        )
    end
    return nothing
end

function update!(model::NoSnowModel, atmospheric_forcing::AtmosphericForcing)
    return nothing
end

# wrapper methods
get_runoff(model::NoSnowModel) = 0.0
get_runoff(model::AbstractSnowModel) = model.variables.runoff
get_snow_storage(model::NoSnowModel) = 0.0
get_snow_storage(model::AbstractSnowModel) = model.variables.snow_storage
get_snow_water(model::NoSnowModel) = 0.0
get_snow_water(model::AbstractSnowModel) = model.variables.snow_water