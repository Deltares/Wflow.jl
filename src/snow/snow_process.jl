
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
    snow,
    snowwater,
    snow_precip,
    liquid_precip,
    temperature,
    ttm,
    cfmax,
    whc;
    cfr = 0.05,
)
    # potential snow melt, based on temperature
    potsnowmelt = temperature > ttm ? cfmax * (temperature - ttm) : 0.0
    # potential refreezing, based on temperature
    potrefreezing = temperature < ttm ? cfmax * cfr * (ttm - temperature) : 0.0
    # actual refreezing
    refreezing = temperature < ttm ? min(potrefreezing, snowwater) : 0.0

    # no landuse correction here
    snowmelt = min(potsnowmelt, snow)  # actual snow melt
    snow = snow + snow_precip + refreezing - snowmelt  # dry snow content
    snowwater = snowwater - refreezing  # free water content in snow
    maxsnowwater = snow * whc  # max water in the snow
    snowwater = snowwater + snowmelt + liquid_precip  # add all water and potentially supersaturate the snowpack
    runoff = max(snowwater - maxsnowwater, 0.0)  # rain + surpluss snowwater
    snowwater = snowwater - runoff
    swe = snowwater + snow # snow water equivalent

    return snow, snowwater, swe, snowmelt, runoff
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
        Float(temperature > tt)
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