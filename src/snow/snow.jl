abstract type AbstractSnowModel{T} end

"Struct for storing snow model variables"
@get_units @grid_loc @with_kw struct SnowVariables{T}
    # Snow storage [mm]
    snow_storage::Vector{T} | "mm"
    # Liquid water content in the snow pack [mm]
    snow_water::Vector{T} | "mm"
    # Snow water equivalent (SWE) [mm]
    swe::Vector{T} | "mm"
    # Snow melt [mm Δt⁻¹]
    snow_melt::Vector{T}
    # Runoff from snowpack [mm Δt⁻¹]
    runoff::Vector{T}
end

"Initialize snow model variables"
function SnowVariables(T::Type{<:AbstractFloat}, n::Int)
    return SnowVariables{T}(;
        snow_storage = fill(0.0, n),
        snow_water = fill(0.0, n),
        swe = fill(mv, n),
        runoff = fill(mv, n),
        snow_melt = fill(mv, n),
    )
end

"Struct for storing snow model boundary conditions"
@get_units @grid_loc @with_kw struct SnowBC{T}
    # Effective precipitation [mm Δt⁻¹]
    effective_precip::Vector{T}
    # Snow precipitation [mm Δt⁻¹]
    snow_precip::Vector{T}
    # Liquid precipitation [mm Δt⁻¹]
    liquid_precip::Vector{T}
end

"Initialize snow model boundary conditions"
function SnowBC(T::Type{<:AbstractFloat}, n::Int)
    return SnowBC{T}(;
        effective_precip = fill(mv, n),
        snow_precip = fill(mv, n),
        liquid_precip = fill(mv, n),
    )
end

"Struct for storing snow HBV model parameters"
@get_units @grid_loc @with_kw struct SnowHbvParameters{T}
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹]
    cfmax::Vector{T} | "mm ᵒC-1 dt-1"
    # Threshold temperature for snowfall [ᵒC]
    tt::Vector{T} | "ᵒC"
    # Threshold temperature interval length [ᵒC]
    tti::Vector{T} | "ᵒC"
    # Threshold temperature for snowmelt [ᵒC]
    ttm::Vector{T} | "ᵒC"
    # Water holding capacity as fraction of current snow pack [-]
    whc::Vector{T} | "-"
end

"Snow HBV model"
@with_kw struct SnowHbvModel{T} <: AbstractSnowModel{T}
    boundary_conditions::SnowBC{T}
    parameters::SnowHbvParameters{T}
    variables::SnowVariables{T}
end

struct NoSnowModel{T} <: AbstractSnowModel{T} end

"Initialize snow HBV model parameters"
function SnowHbvParameters(dataset, config, indices, dt)
    cfmax =
        ncread(
            dataset,
            config,
            "land.snow.parameters.cfmax";
            sel = indices,
            defaults = 3.75653,
            type = Float,
        ) .* (dt / basetimestep)
    tt = ncread(
        dataset,
        config,
        "land.snow.parameters.tt";
        sel = indices,
        defaults = 0.0,
        type = Float,
    )
    tti = ncread(
        dataset,
        config,
        "land.snow.parameters.tti";
        sel = indices,
        defaults = 1.0,
        type = Float,
    )
    ttm = ncread(
        dataset,
        config,
        "land.snow.parameters.ttm";
        sel = indices,
        defaults = 0.0,
        type = Float,
    )
    whc = ncread(
        dataset,
        config,
        "land.snow.parameters.whc";
        sel = indices,
        defaults = 0.1,
        type = Float,
    )
    snow_hbv_params =
        SnowHbvParameters(; cfmax = cfmax, tt = tt, tti = tti, ttm = ttm, whc = whc)
    return snow_hbv_params
end

"Initialize snow HBV model"
function SnowHbvModel(dataset, config, indices, dt)
    n = length(indices)
    params = SnowHbvParameters(dataset, config, indices, dt)
    vars = SnowVariables(Float, n)
    bc = SnowBC(Float, n)
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