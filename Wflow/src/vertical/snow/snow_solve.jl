@enumx SnowSubtimestepType until_end no_liquid no_snow snow_saturated

function update_preamble!(
    model_vertical::NewSnowModel,
    land::LandHydrologySBM,
    i::Int,
)::Nothing
    (; integrators, p) = model_vertical
    (; cache, properties) = p
    (; snow_precip, liquid_precip, freeze_rate, melt_rate, snow_storage, snow_water) = cache
    (; tt, tti, ttm) = properties

    # Reset snow melt and runoff
    (; u) = integrators[i]
    u.snow_melt = 0
    u.runoff = 0

    # Obtain the sub time stepping initial condition
    u.snow_storage = snow_storage[i]
    u.snow_water = snow_water[i]

    # Precipitation rates
    temperature = land.atmospheric_forcing.temperature[i]
    effective_precip = cache.effective_precip[i]
    rainfrac = clamp(0.5 + (temperature - tt[i]) / tti[i], 0, 1)
    liquid_precip[i] = rainfrac * effective_precip
    snow_precip[i] = effective_precip - liquid_precip[i]

    # Freeze and melt rates
    cfmax = properties.cfmax[i]
    Δtemperature = temperature - ttm[i]
    melt_rate[i] = cfmax * max(0, Δtemperature)
    refreeze_coefficient = 0.05
    freeze_rate[i] = cfmax * refreeze_coefficient * max(0, -Δtemperature)
    return nothing
end

function set_instantaneous_rates!(
    du::SnowStateType,
    u::SnowStateType,
    p::NewSnowParameters,
    i::Int,
    progress::Float64,
)::Nothing
    (; cache, properties) = p
    (; snow_precip, liquid_precip, freeze_rate, melt_rate) = cache
    (; whc) = properties

    no_liquid = iszero(u.snow_water)
    no_snow = iszero(u.snow_storage)
    snow_saturated = (u.snow_water ≈ whc[i] * u.snow_storage)

    q_snow = snow_precip[i]
    q_rain = liquid_precip[i]
    q_freeze = no_liquid ? 0.0 : freeze_rate[i]
    q_melt = no_snow ? 0.0 : melt_rate[i]
    q_runoff = snow_saturated ? max(0.0, q_rain - q_freeze + q_melt) : 0.0

    du.snow_storage = q_snow + q_freeze - q_melt
    du.snow_water = q_rain - q_freeze + q_melt - q_runoff
    du.snow_melt = q_melt
    du.runoff = q_runoff
    return nothing
end

function set_sub_time_step!(model_vertical::NewSnowModel, i::Int)::SnowSubtimestepType.T
    (; p, integrators) = model_vertical
    whc = p.properties.whc[i]
    integrator = integrators[i]
    (; progress, u, du) = integrator

    # Until the end of the global time step
    dt_sub_end = 1 - progress

    # Until there is no liquid water left
    dt_sub_no_liquid = du.snow_water < 0 ? -u.snow_water / du.snow_water : Inf

    # Until there is no snow left
    dt_sub_no_snow = du.snow_storage < 0 ? -u.snow_storage / du.snow_storage : Inf

    # Until the snow saturation threshold has been reached
    dt_sub_saturation =
        -(u.snow_water - whc * u.snow_storage) / (du.snow_water - whc * du.snow_storage)
    dt_sub_saturation = (dt_sub_saturation > 0) ? dt_sub_saturation : Inf

    # Get the smallest dt_sub of the candidates
    # Candidates must be in the same order as the corresponding
    # sub time step type enumerator
    dt_sub_candidates =
        SVector(dt_sub_end, dt_sub_no_liquid, dt_sub_no_snow, dt_sub_saturation)
    dt_sub_idx = argmin(dt_sub_candidates)
    set_dt_sub!(integrator, dt_sub_candidates[dt_sub_idx])
    return SnowSubtimestepType.T(dt_sub_idx - 1)
end

function dt_sub_callback!(
    model_vertical::NewSnowModel,
    ::LandHydrologySBM,
    i::Int,
    dt_sub_type::SnowSubtimestepType.T,
)::Nothing
    (; integrators) = model_vertical
    integrator = integrators[i]
    (; u) = integrator

    if dt_sub_type == SnowSubtimestepType.until_end
        integrator.progress = 1.0
    elseif dt_sub_type == SnowSubtimestepType.no_liquid
        u.snow_water = 0
    elseif dt_sub_type == SnowSubtimestepType.no_snow
        u.snow_storage = 0
    end
    return nothing
end

function update_postamble!(
    model_vertical::NewSnowModel,
    ::LandHydrologySBM,
    i::Int,
)::Nothing
    (; p, integrators) = model_vertical
    (; cache) = p
    (; u) = integrators[i]

    # Copy integration result to cache
    cache.snow_storage[i] = u.snow_storage
    cache.snow_water[i] = u.snow_water
    cache.swe[i] = u.snow_storage + u.snow_water
    cache.snow_melt[i] = u.snow_melt
    cache.runoff[i] = u.runoff

    return nothing
end