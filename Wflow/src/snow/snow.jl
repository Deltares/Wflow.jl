abstract type AbstractSnowModel end

"Struct for storing snow model variables"
@with_data_lookup struct SnowVariables
    n::Int
    # Snow storage [m]
    "snowpack_dry_snow__leq_depth"
    snow_storage::Vector{Float64} = zeros(n)
    # Liquid water content in the snow pack [m]
    "snowpack_liquid_water__depth"
    snow_water::Vector{Float64} = zeros(n)
    # Snow water equivalent (SWE) [m]
    "snowpack__leq_depth"
    swe::Vector{Float64} = fill(MISSING_VALUE, n)
    # Snow melt [m s⁻¹]
    "snowpack_meltwater__volume_flux"
    snow_melt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Runoff from snowpack [m s⁻¹]
    runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Lateral snow (SWE) transport from upstreams cells [m s⁻¹]
    snow_in::Vector{Float64} = zeros(n)
    # Lateral snow (SWE) transport out of a cell [m s⁻¹]
    snow_out::Vector{Float64} = zeros(n)
end

"Struct for storing snow model boundary conditions"
@kwdef struct SnowBC
    n::Int
    # Effective precipitation [m s⁻¹]
    effective_precip::Vector{Float64} = fill(MISSING_VALUE, n)
    # Snow precipitation [m s⁻¹]
    snow_precip::Vector{Float64} = fill(MISSING_VALUE, n)
    # Liquid precipitation [m s⁻¹]
    liquid_precip::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing snow HBV model parameters"
@with_data_lookup struct SnowHbvParameters
    # Degree-day factor [m K⁻¹ s⁻¹]
    "snowpack__degree_day_coefficient"
    cfmax::Vector{Float64}
    # Threshold temperature for snowfall [K]
    "atmosphere_air__snowfall_temperature_threshold"
    tt::Vector{Float64}
    # Threshold temperature interval length [K]
    "atmosphere_air__snowfall_temperature_interval"
    tti::Vector{Float64}
    # Threshold temperature for snowmelt [K]
    "snowpack__melting_temperature_threshold"
    ttm::Vector{Float64}
    # Water holding capacity as fraction of current snow pack [-]
    "snowpack__liquid_water_holding_capacity"
    whc::Vector{Float64}
end

"Snow HBV model"
@kwdef struct SnowHbvModel <: AbstractSnowModel
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
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
)
    cfmax = ncread(
        dataset,
        config,
        "snowpack__degree_day_coefficient",
        LandHydrologySBM;
        sel = indices,
    )
    tt = ncread(
        dataset,
        config,
        "atmosphere_air__snowfall_temperature_threshold",
        LandHydrologySBM;
        sel = indices,
    )
    tti = ncread(
        dataset,
        config,
        "atmosphere_air__snowfall_temperature_interval",
        LandHydrologySBM;
        sel = indices,
    )
    ttm = ncread(
        dataset,
        config,
        "snowpack__melting_temperature_threshold",
        LandHydrologySBM;
        sel = indices,
    )
    whc = ncread(
        dataset,
        config,
        "snowpack__liquid_water_holding_capacity",
        LandHydrologySBM;
        sel = indices,
    )
    snow_hbv_params = SnowHbvParameters(data_lookup; cfmax, tt, tti, ttm, whc)
    return snow_hbv_params
end

"Initialize snow HBV model"
function SnowHbvModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
)
    n = length(indices)
    parameters = SnowHbvParameters(dataset, config, indices; data_lookup)
    variables = SnowVariables(data_lookup; n)
    snow_model = SnowHbvModel(; n, parameters, variables)
    return snow_model
end

"Update boundary condition (effective precipitation provided by an interception model) of a snow model for a single timestep"
function update_bc_snow_model!(snow_model::AbstractSnowModel, external_models::NamedTuple)
    (; effective_precip) = snow_model.boundary_conditions
    (; interception) = external_models
    @. effective_precip =
        interception.variables.throughfall + interception.variables.stemflow
    return nothing
end

function update_bc_snow_model!(snow_model::NoSnowModel, ::NamedTuple)
    return nothing
end

"Update snow HBV model for a single timestep"
function update_snow_model!(
    snow_model::SnowHbvModel,
    atmospheric_forcing::AtmosphericForcing,
    dt::Float64,
)
    (; temperature) = atmospheric_forcing
    (; snow_storage, snow_water, swe, snow_melt, runoff) = snow_model.variables
    (; effective_precip, snow_precip, liquid_precip) = snow_model.boundary_conditions
    (; tt, tti, ttm, cfmax, whc) = snow_model.parameters

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
            dt,
        )
    end
    return nothing
end

function update_snow_model!(
    model::NoSnowModel,
    atmospheric_forcing::AtmosphericForcing,
    dt::Float64,
)
    return nothing
end

# wrapper methods
get_runoff(snow_model::NoSnowModel) = Zeros(snow_model.n)
get_runoff(snow_model::AbstractSnowModel) = snow_model.variables.runoff
get_snow_storage(snow_model::NoSnowModel) = Zeros(snow_model.n)
get_snow_storage(snow_model::AbstractSnowModel) = snow_model.variables.snow_storage
get_snow_water(snow_model::NoSnowModel) = Zeros(snow_model.n)
get_snow_water(snow_model::AbstractSnowModel) = snow_model.variables.snow_water
get_snow_out(snow_model::NoSnowModel) = Zeros(snow_model.n)
get_snow_out(snow_model::AbstractSnowModel) = snow_model.variables.snow_out
get_snow_in(snow_model::NoSnowModel) = Zeros(snow_model.n)
get_snow_in(snow_model::AbstractSnowModel) = snow_model.variables.snow_in
