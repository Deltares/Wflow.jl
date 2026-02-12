
"""
    snowpack_hbv(snow, snowwater, snow_precip, liquid_precip, temperature, ttm, cfmax, whc, dt; cfr = 0.05)

HBV type snowpack modeling using a temperature degree factor.
The refreezing efficiency factor `cfr` is set to 0.05.

# Arguments
- `snow` (snow storage)
- `snowwater` (liquid water content in the snow pack)
- `snow_precip` (snow precipitation)
- `liquid_precip` (liquid precipitation)
- `ttm` (melting threshold)
- `cfmax` (degree day factor, rate of snowmelt)
- `whc` (water holding capacity of snow)
- `cfr` refreeing efficiency constant in refreezing of liquied water in snow
- `dt` timestep

# Output
- `snow`
- `snowwater`
- `swe` (snow water equivalent)
- `runoff`
"""
function snowpack_hbv(
    snow_storage,
    snow_water,
    snow_precip,
    liquid_precip,
    temperature,
    ttm,
    cfmax,
    whc,
    dt;
    cfr = 0.05,
)
    if temperature > ttm
        # [m s⁻¹] = [m K⁻¹ s⁻¹] * ([K] - [K])
        potential_snow_melt = cfmax * (temperature - ttm)
        # [m s⁻¹] = min([m s⁻¹], [m] / [s])
        snow_melt = min(potential_snow_melt, snow_storage / dt)
        # [m] -= [m s⁻¹] * [s]
        snow_storage -= snow_melt * dt
        # [m] += [m s⁻¹] * [s]
        snow_water += snow_melt * dt
    else
        snow_melt = 0.0

        # [m s⁻¹] = [m K⁻¹ s⁻¹] * [-] * ([K] - [K])
        potential_refreezing = cfmax * cfr * (ttm - temperature)
        # [m] = min([m s⁻¹] * [s], [m])
        refreezing = min(potential_refreezing * dt, snow_water)
        # [m] += [m]
        snow_storage += refreezing
        # [m] -= [m]
        snow_water -= refreezing
    end

    # [m] += [m s⁻¹] * [s]
    snow_storage += snow_precip * dt # dry snow content
    # [m] += [m s⁻¹] * [s]
    snow_water += liquid_precip * dt # free water content in snow

    # [m] = [m] * [-]
    max_snow_water = snow_storage * whc  # max water in the snow

    if snow_water > max_snow_water
        # [m s⁻¹] = ([m] - [m]) / [s]
        runoff = (snow_water - max_snow_water) / dt
        snow_water = max_snow_water
    else
        runoff = 0.0
    end
    # [m] = [m] + [m]
    snow_water_equivalent = snow_water + snow_storage

    return snow_storage, snow_water, snow_water_equivalent, snow_melt, runoff
end

"""
    precipitation_hbv(precipitation, temperature, tti, tt; rfcf = 1.0, sfcf = 1.0)

HBV type precipitation routine to separate precipitation in snow precipitation and liquid precipitation.
All correction factors (RFCF and SFCF) are set to 1.

# Arguments
- `precipitation`
- `temperature` (threshold temperature for snowfall)
- `tti` (snowfall threshold interval length)
- `tt` (threshold temperature for snowfall)
- `rfcf` correction factor for liquid precipitation (rainfall)
- `sfcf` correction factor for snow precipitation (snowfall)

# Output
- `snow_precip`
- `liquid_precip`
"""
function precipitation_hbv(precipitation, temperature, tti, tt; rfcf = 1.0, sfcf = 1.0)
    # fraction of precipitation which falls as rain
    # [-]
    rainfrac = if iszero(tti)
        Float64(temperature > tt)
    else
        frac = (temperature - (tt - tti / 2.0)) / tti
        rainfrac = clamp(frac, 0.0, 1.0)
    end

    # fraction of precipitation which falls as snow
    snowfrac = 1.0 - rainfrac
    # different correction for liquid_precip and snow_precip
    # [m s⁻¹] = [-] * [-] * [m s⁻¹]
    snow_precip = snowfrac * sfcf * precipitation
    liquid_precip = rainfrac * rfcf * precipitation
    return snow_precip, liquid_precip
end
