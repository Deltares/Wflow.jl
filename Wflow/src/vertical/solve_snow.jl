function update_preamble!(model_vertical::NewSnowModel, model::Model, i::Int)::Nothing
    (; cache, properties) = model_vertical.p
    (; snow_precip, liquid_precip, freeze_rate, melt_rate) = cache
    (; tt, tti, ttm) = properties

    # Precipitation rates
    temperature = model.land.atmospheric_forcing.temperature[i]
    effective_precip = cache.effective_precip[i]
    rainfrac = clamp(0.5 + (temperature - tt[i]) / tti[i], 0, 1)
    liquid_precip[i] = rainfrac * effective_precip
    snow_precip[i] = effective_precip - liquid_precip

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

    no_liquid = (u.snow_water ≈ 0)
    no_snow = (u.snow_storage ≈ 0)
    snow_saturated = (u.snow_water ≈ whc[i] * u.snow_storage)

    q_snow = snow_precip[i]
    q_rain = liquid_precip[i]
    q_freeze = no_liquid ? 0.0 : freeze_rate[i]
    q_melt = no_snow ? 0.0 : melt_rate[i]
    q_runoff = snow_saturated ? max(0.0, q_rain - q_freeze + q_melt) : 0.0

    du.snow_storage = q_snow + q_freeze - q_melt
    du.snow_water = q_rain - q_freeze + q_melt - q_runoff
    du.cumulative_snow_melt = q_melt
    du.cumulative_runoff = q_runoff
    return nothing
end

function update_postamble!(model_vertical::NewSnowModel, model::Model, i::Int)::Nothing
    return nothing
end