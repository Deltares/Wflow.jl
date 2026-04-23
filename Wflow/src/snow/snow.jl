abstract type AbstractSnowModel end

"Struct for storing snow model variables"
@with_kw struct SnowVariables
    n_land_cells::Int
    # Snow storage [mm]
    snow_storage::Vector{Float64} = zeros(n_land_cells)
    # Liquid water content in the snow pack [mm]
    snow_water::Vector{Float64} = zeros(n_land_cells)
    # Snow water equivalent (SWE) [mm]
    swe::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Snow melt [mm Δt⁻¹]
    snow_melt::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Runoff from snowpack [mm Δt⁻¹]
    runoff::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Lateral snow (SWE) transport from upstreams cells [mm Δt⁻¹]
    snow_in::Vector{Float64} = zeros(n_land_cells)
    # Lateral snow (SWE) transport out of a cell [mm Δt⁻¹]
    snow_out::Vector{Float64} = zeros(n_land_cells)
end

"Struct for storing snow model boundary conditions"
@with_kw struct SnowBC
    n_land_cells::Int
    # Effective precipitation [mm Δt⁻¹]
    effective_precip::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Snow precipitation [mm Δt⁻¹]
    snow_precip::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Liquid precipitation [mm Δt⁻¹]
    liquid_precip::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
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
    n_land_cells::Int
    boundary_conditions::SnowBC = SnowBC(; n_land_cells)
    parameters::SnowHbvParameters
    variables::SnowVariables = SnowVariables(; n_land_cells)
end

struct NoSnowModel <: AbstractSnowModel
    n_land_cells::Int
end

"Initialize snow HBV model parameters"
function SnowHbvParameters(
    dataset::NCDataset,
    config::Config,
    indices_land::Vector{CartesianIndex{2}},
    dt::Second,
)
    cfmax =
        ncread(
            dataset,
            config,
            "snowpack__degree_day_coefficient",
            LandHydrologySBM;
            sel=indices_land,
        ) .* (dt / BASETIMESTEP)
    tt = ncread(
        dataset,
        config,
        "atmosphere_air__snowfall_temperature_threshold",
        LandHydrologySBM;
        sel=indices_land,
    )
    tti = ncread(
        dataset,
        config,
        "atmosphere_air__snowfall_temperature_interval",
        LandHydrologySBM;
        sel=indices_land,
    )
    ttm = ncread(
        dataset,
        config,
        "snowpack__melting_temperature_threshold",
        LandHydrologySBM;
        sel=indices_land,
    )
    whc = ncread(
        dataset,
        config,
        "snowpack__liquid_water_holding_capacity",
        LandHydrologySBM;
        sel=indices_land,
    )
    snow_hbv_params = SnowHbvParameters(; cfmax, tt, tti, ttm, whc)
    return snow_hbv_params
end

"Initialize snow HBV model"
function SnowHbvModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
    dt::Second,
)
    n_land_cells = length(land_indices_2d)
    parameters = SnowHbvParameters(dataset, config, land_indices_2d, dt)
    snow_model = SnowHbvModel(; n_land_cells, parameters)
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
)
    (; boundary_conditions, parameters, variables, n_land_cells) = snow_model
    (; snow_storage, snow_water, swe, snow_melt, runoff) = variables
    (; effective_precip, snow_precip, liquid_precip) = boundary_conditions
    (; tt, tti, ttm, cfmax, whc) = parameters
    (; temperature) = atmospheric_forcing

    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        snow_precip[land_cell_idx], liquid_precip[land_cell_idx] = precipitation_hbv(
            effective_precip[land_cell_idx],
            temperature[land_cell_idx],
            tti[land_cell_idx],
            tt[land_cell_idx],
        )
    end
    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        snow_storage[land_cell_idx],
        snow_water[land_cell_idx],
        swe[land_cell_idx],
        snow_melt[land_cell_idx],
        runoff[land_cell_idx] = snowpack_hbv(
            snow_storage[land_cell_idx],
            snow_water[land_cell_idx],
            snow_precip[land_cell_idx],
            liquid_precip[land_cell_idx],
            temperature[land_cell_idx],
            ttm[land_cell_idx],
            cfmax[land_cell_idx],
            whc[land_cell_idx],
        )
    end
    return nothing
end

function update_snow_model!(
    snow_model::NoSnowModel,
    atmospheric_forcing::AtmosphericForcing,
)
    return nothing
end

# wrapper methods
get_runoff(snow_model::NoSnowModel) = Zeros(snow_model.n_land_cells)
get_runoff(snow_model::AbstractSnowModel) = snow_model.variables.runoff
get_snow_storage(snow_model::NoSnowModel) = Zeros(snow_model.n_land_cells)
get_snow_storage(snow_model::AbstractSnowModel) = snow_model.variables.snow_storage
get_snow_water(snow_model::NoSnowModel) = Zeros(snow_model.n_land_cells)
get_snow_water(snow_model::AbstractSnowModel) = snow_model.variables.snow_water
get_snow_out(snow_model::NoSnowModel) = Zeros(snow_model.n_land_cells)
get_snow_out(snow_model::AbstractSnowModel) = snow_model.variables.snow_out
get_snow_in(snow_model::NoSnowModel) = Zeros(snow_model.n_land_cells)
get_snow_in(snow_model::AbstractSnowModel) = snow_model.variables.snow_in
