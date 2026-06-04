abstract type AbstractSnowModel end

"Struct for storing snow model variables"
@with_kw struct SnowVariables
    n::Int
    # Snow storage [m]
    snow_storage::Vector{Float64} = zeros(n)
    # Liquid water content in the snow pack [m]
    snow_water::Vector{Float64} = zeros(n)
    # Snow water equivalent (SWE) [m]
    snow_water_equivalent::Vector{Float64} = fill(MISSING_VALUE, n)
    # Snow melt [m s⁻¹]
    snow_melt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Runoff from snowpack [m s⁻¹]
    runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Lateral snow (SWE) transport from upstreams cells [m s⁻¹]
    snow_in::Vector{Float64} = zeros(n)
    # Lateral snow (SWE) transport out of a cell [m s⁻¹]
    snow_out::Vector{Float64} = zeros(n)
end

"Struct for storing snow model boundary conditions"
@with_kw struct SnowBC
    n::Int
    # Effective precipitation [m s⁻¹]
    effective_precip::Vector{Float64} = fill(MISSING_VALUE, n)
    # Snow precipitation [m s⁻¹]
    snow_precip::Vector{Float64} = fill(MISSING_VALUE, n)
    # Liquid precipitation [m s⁻¹]
    liquid_precip::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing snow HBV model parameters"
@with_kw struct SnowHbvParameters
    # Degree-day factor [m K⁻¹ s⁻¹]
    degree_day_factor::Vector{Float64}
    # Threshold temperature for snowfall [K]
    temperature_threshold_snowfall::Vector{Float64}
    # Threshold temperature interval length [K]
    temperature_interval_snowfall::Vector{Float64}
    # Threshold temperature for snowmelt [K]
    temperature_threshold_melt::Vector{Float64}
    # Water holding capacity as fraction of current snow pack [-]
    water_holding_capacity::Vector{Float64}
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
)
    degree_day_factor = ncread(
        dataset,
        config,
        "snowpack__degree_day_coefficient",
        LandHydrologySBM;
        sel = indices,
    )
    temperature_threshold_snowfall = ncread(
        dataset,
        config,
        "atmosphere_air__snowfall_temperature_threshold",
        LandHydrologySBM;
        sel = indices,
    )
    temperature_interval_snowfall = ncread(
        dataset,
        config,
        "atmosphere_air__snowfall_temperature_interval",
        LandHydrologySBM;
        sel = indices,
    )
    temperature_threshold_melt = ncread(
        dataset,
        config,
        "snowpack__melting_temperature_threshold",
        LandHydrologySBM;
        sel = indices,
    )
    water_holding_capacity = ncread(
        dataset,
        config,
        "snowpack__liquid_water_holding_capacity",
        LandHydrologySBM;
        sel = indices,
    )
    snow_hbv_params = SnowHbvParameters(;
        degree_day_factor,
        temperature_threshold_snowfall,
        temperature_interval_snowfall,
        temperature_threshold_melt,
        water_holding_capacity,
    )
    return snow_hbv_params
end

"Initialize snow HBV model"
function SnowHbvModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = SnowHbvParameters(dataset, config, indices)
    snow_model = SnowHbvModel(; n, parameters)
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
    (; snow_storage, snow_water, snow_water_equivalent, snow_melt, runoff) =
        snow_model.variables
    (; effective_precip, snow_precip, liquid_precip) = snow_model.boundary_conditions
    (;
        temperature_threshold_snowfall,
        temperature_interval_snowfall,
        temperature_threshold_melt,
        degree_day_factor,
        water_holding_capacity,
    ) = snow_model.parameters

    n_cells = length(temperature)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        snow_precip[cell_idx], liquid_precip[cell_idx] = precipitation_hbv(
            effective_precip[cell_idx],
            temperature[cell_idx],
            temperature_interval_snowfall[cell_idx],
            temperature_threshold_snowfall[cell_idx],
        )
    end
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        snow_storage[cell_idx],
        snow_water[cell_idx],
        snow_water_equivalent[cell_idx],
        snow_melt[cell_idx],
        runoff[cell_idx] = snowpack_hbv(
            snow_storage[cell_idx],
            snow_water[cell_idx],
            snow_precip[cell_idx],
            liquid_precip[cell_idx],
            temperature[cell_idx],
            temperature_threshold_melt[cell_idx],
            degree_day_factor[cell_idx],
            water_holding_capacity[cell_idx],
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
