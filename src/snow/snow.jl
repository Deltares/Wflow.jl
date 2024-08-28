
@get_units @with_kw struct SnowModelVars{T}
    # Snow storage [mm]
    snow_storage::Vector{T} | "mm"
    # Liquid water content in the snow pack [mm]
    snow_water::Vector{T} | "mm"
    # Snow water equivalent (SWE) [mm]
    swe::Vector{T} | "mm"
    # Runoff from snowpack [mm Δt⁻¹]
    runoff::Vector{T}
end

function snow_model_vars(n)
    vars = SnowModelVars(;
        snow_storage = fill(0.0, n),
        snow_water = fill(0.0, n),
        swe = fill(mv, n),
        runoff = fill(mv, n),
    )
    return vars
end

@get_units @with_kw struct SnowBC{T}
    # Effective precipitation [mm Δt⁻¹]
    effective_precip::Vector{T}
    # Snow precipitation [mm Δt⁻¹]
    snow_precip::Vector{T}
    # Liquid precipitation [mm Δt⁻¹]
    liquid_precip::Vector{T}
end

function snow_model_bc(n)
    bc = SnowBC(;
        effective_precip = fill(mv, n),
        snow_precip = fill(mv, n),
        liquid_precip = fill(mv, n),
    )
    return bc
end

@get_units @with_kw struct SnowHbvParameters{T}
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

abstract type AbstractSnowModel end

@get_units @with_kw struct SnowHbvModel{T} <: AbstractSnowModel
    boundary_conditions::SnowBC{T} | "-"
    parameters::SnowHbvParameters{T} | "-"
    variables::SnowModelVars{T} | "-"
end

struct NoSnowModel <: AbstractSnowModel end

function initialize_snow_hbv_params(nc, config, inds, dt)
    cfmax =
        ncread(
            nc,
            config,
            "vertical.snow.parameters.cfmax";
            sel = inds,
            defaults = 3.75653,
            type = Float,
        ) .* (dt / basetimestep)
    tt = ncread(
        nc,
        config,
        "vertical.snow.parameters.tt";
        sel = inds,
        defaults = 0.0,
        type = Float,
    )
    tti = ncread(
        nc,
        config,
        "vertical.snow.parameters.tti";
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    ttm = ncread(
        nc,
        config,
        "vertical.snow.parameters.ttm";
        sel = inds,
        defaults = 0.0,
        type = Float,
    )
    whc = ncread(
        nc,
        config,
        "vertical.snow.parameters.whc";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    snow_hbv_params =
        SnowHbvParameters(; cfmax = cfmax, tt = tt, tti = tti, ttm = ttm, whc = whc)
    return snow_hbv_params
end

function initialize_snow_hbv_model(nc, config, inds, dt)
    n = length(inds)
    params = initialize_snow_hbv_params(nc, config, inds, dt)
    vars = snow_model_vars(n)
    bc = snow_model_bc(n)
    model = SnowHbvModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

function update!(model::SnowHbvModel, atmospheric_forcing::AtmosphericForcing)
    (; temperature) = atmospheric_forcing
    (; snow_storage, snow_water, swe, runoff) = model.variables
    (; effective_precip, snow_precip, liquid_precip) = model.boundary_conditions
    (; tt, tti, ttm, cfmax, whc) = model.parameters

    n = length(temperature)
    threaded_foreach(1:n; basesize = 1000) do i
        snow_precip[i], liquid_precip[i] =
            precipitation_hbv(effective_precip[i], temperature[i], tti[i], tt[i])
    end
    threaded_foreach(1:n; basesize = 1000) do i
        snow_storage[i], snow_water[i], swe[i], runoff[i] = snowpack_hbv(
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
end

function update!(model::NoSnowModel, atmospheric_forcing::AtmosphericForcing)
    return nothing
end

get_runoff(model::NoSnowModel) = 0.0
get_runoff(model::AbstractSnowModel) = model.variables.runoff
get_snow_storage(model::NoSnowModel) = 0.0
get_snow_storage(model::AbstractSnowModel) = model.variables.snow_storage
get_snow_water(model::NoSnowModel) = 0.0
get_snow_water(model::AbstractSnowModel) = model.variables.snow_water