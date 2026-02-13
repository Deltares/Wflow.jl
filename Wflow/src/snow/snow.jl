abstract type AbstractSnowModel end

"Struct for storing snow model variables"
@with_kw struct SnowVariables
    n::Int
    # Snow storage [mm]
    snow_storage::Vector{Float64} = zeros(n)
    # Liquid water content in the snow pack [mm]
    snow_water::Vector{Float64} = zeros(n)
    # Snow water equivalent (SWE) [mm]
    swe::Vector{Float64} = fill(MISSING_VALUE, n)
    # Snow melt [mm Δt⁻¹]
    snow_melt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Runoff from snowpack [mm Δt⁻¹]
    runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Lateral snow (SWE) transport from upstreams cells [mm Δt⁻¹]
    snow_in::Vector{Float64} = zeros(n)
    # Lateral snow (SWE) transport out of a cell [mm Δt⁻¹]
    snow_out::Vector{Float64} = zeros(n)
end

"Struct for storing snow model boundary conditions"
@with_kw struct SnowBC
    n::Int
    # Effective precipitation [mm Δt⁻¹]
    effective_precip::Vector{Float64} = fill(MISSING_VALUE, n)
    # Snow precipitation [mm Δt⁻¹]
    snow_precip::Vector{Float64} = fill(MISSING_VALUE, n)
    # Liquid precipitation [mm Δt⁻¹]
    liquid_precip::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing snow HBV model parameters"
@with_kw struct SnowHbvParameters
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹]
    cfmax::Vector{Float64}
    # Threshold temperature for snowfall [ᵒC]
    tt::Vector{Float64}
    # Threshold temperature interval length [ᵒC]
    tti::Vector{Float64}
    # Threshold temperature for snowmelt [ᵒC]
    ttm::Vector{Float64}
    # Water holding capacity as fraction of current snow pack [-]
    whc::Vector{Float64}
end

"Snow HBV model"
@with_kw struct SnowHbvModel <: AbstractSnowModel
    n::Int
    boundary_conditions::SnowBC = SnowBC(; n)
    parameters::SnowHbvParameters
    variables::SnowVariables = SnowVariables(; n)
end

struct NoSnowModel <: AbstractSnowModel
    n::Int
end

"Initialize snow HBV model parameters"
function SnowHbvParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    cfmax =
        ncread(
            dataset,
            config,
            "snowpack__degree_day_coefficient";
            sel = indices,
            defaults = 3.75,
            type = Float64,
        ) .* (dt / BASETIMESTEP)
    tt = ncread(
        dataset,
        config,
        "atmosphere_air__snowfall_temperature_threshold";
        sel = indices,
        defaults = 0.0,
        type = Float64,
    )
    tti = ncread(
        dataset,
        config,
        "atmosphere_air__snowfall_temperature_interval";
        sel = indices,
        defaults = 1.0,
        type = Float64,
    )
    ttm = ncread(
        dataset,
        config,
        "snowpack__melting_temperature_threshold";
        sel = indices,
        defaults = 0.0,
        type = Float64,
    )
    whc = ncread(
        dataset,
        config,
        "snowpack__liquid_water_holding_capacity";
        sel = indices,
        defaults = 0.1,
        type = Float64,
    )
    snow_hbv_params = SnowHbvParameters(; cfmax, tt, tti, ttm, whc)
    return snow_hbv_params
end

"Initialize snow HBV model"
function SnowHbvModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    n = length(indices)
    parameters = SnowHbvParameters(dataset, config, indices, dt)
    model = SnowHbvModel(; n, parameters)
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
get_runoff(model::NoSnowModel) = Zeros(model.n)
get_runoff(model::AbstractSnowModel) = model.variables.runoff
get_snow_storage(model::NoSnowModel) = Zeros(model.n)
get_snow_storage(model::AbstractSnowModel) = model.variables.snow_storage
get_snow_water(model::NoSnowModel) = Zeros(model.n)
get_snow_water(model::AbstractSnowModel) = model.variables.snow_water
get_snow_out(model::NoSnowModel) = Zeros(model.n)
get_snow_out(model::AbstractSnowModel) = model.variables.snow_out
get_snow_in(model::NoSnowModel) = Zeros(model.n)
get_snow_in(model::AbstractSnowModel) = model.variables.snow_in
