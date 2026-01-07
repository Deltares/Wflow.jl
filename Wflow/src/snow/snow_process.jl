
"""
    snowpack_hbv(snow, snowwater, snow_precip, liquid_precip, temperature, ttm, cfmax, whc; cfr = 0.05)

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
    whc;
    cfr = 0.05,
)
    if temperature > ttm
        potential_snow_melt = cfmax * (temperature - ttm)
        snow_melt = min(potential_snow_melt, snow_storage)
        snow_storage -= snow_melt
        snow_water += snow_melt
    else
        snow_melt = 0.0
        potential_refreezing = cfmax * cfr * (ttm - temperature)
        refreezing = min(potential_refreezing, snow_water)
        snow_storage += refreezing
        snow_water -= refreezing
    end

    snow_storage += snow_precip
    snow_water += liquid_precip

    max_snow_water = snow_storage * whc
    runoff = max(snow_water - max_snow_water, 0)
    snow_water -= runoff
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
    rainfrac = if iszero(tti)
        Float64(temperature > tt)
    else
        frac = (temperature - (tt - tti / 2.0)) / tti
        min(frac, 1.0)
    end
    rainfrac = max(rainfrac, 0.0)

    # fraction of precipitation which falls as snow
    snowfrac = 1.0 - rainfrac
    # different correction for liquid_precip and snow_precip
    snow_precip = snowfrac * sfcf * precipitation  # snow_precip depth
    liquid_precip = rainfrac * rfcf * precipitation  # liquid_precip depth
    return snow_precip, liquid_precip
end
