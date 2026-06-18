
"""
    snowpack_hbv(snow_storage, snow_water, snow_precip, liquid_precip, temperature, temperature_threshold_melt, degree_day_factor, water_holding_capacity, dt; cfr = 0.05)

HBV type snowpack modeling using a temperature degree factor.
The refreezing efficiency factor `cfr` is set to 0.05.

# Arguments
- `snow_storage` (snow storage)
- `snow_water` (liquid water content in the snow pack)
- `snow_precip` (snow precipitation)
- `liquid_precip` (liquid precipitation)
- `temperature` (temperature)
- `temperature_threshold_melt` (melting threshold)
- `degree_day_factor` (degree day factor, rate of snowmelt)
- `water_holding_capacity` (water holding capacity of snow)
- `cfr` refreeing efficiency constant in refreezing of liquied water in snow
- `dt` timestep

# Output
- `snow`
- `snowwater`
- `snow_water_equivalent` (snow water equivalent)
- `runoff`
"""
function snowpack_hbv(
    snow_storage,
    snow_water,
    snow_precip,
    liquid_precip,
    temperature,
    temperature_threshold_melt,
    degree_day_factor,
    water_holding_capacity,
    dt;
    cfr = 0.05,
)
    if temperature > temperature_threshold_melt
        potential_snow_melt = degree_day_factor * (temperature - temperature_threshold_melt)
        snow_melt = min(potential_snow_melt, snow_storage / dt)
        snow_storage -= snow_melt * dt
        snow_water += snow_melt * dt
    else
        snow_melt = 0.0

        potential_refreezing =
            degree_day_factor * cfr * (temperature_threshold_melt - temperature)
        refreezing = min(potential_refreezing * dt, snow_water)
        snow_storage += refreezing
        snow_water -= refreezing
    end

    snow_storage += snow_precip * dt # dry snow content
    snow_water += liquid_precip * dt # free water content in snow

    max_snow_water = snow_storage * water_holding_capacity  # max water in the snow

    if snow_water > max_snow_water
        runoff = (snow_water - max_snow_water) / dt
        snow_water = max_snow_water
    else
        runoff = 0.0
    end
    snow_water_equivalent = snow_water + snow_storage

    return snow_storage, snow_water, snow_water_equivalent, snow_melt, runoff
end

"""
    precipitation_hbv(precipitation, temperature, temperature_interval_snowfall, temperature_threshold_snowfall; rfcf = 1.0, sfcf = 1.0)

HBV type precipitation routine to separate precipitation in snow precipitation and liquid precipitation.
All correction factors (RFCF and SFCF) are set to 1.

# Arguments
- `precipitation`
- `temperature` (threshold temperature for snowfall)
- `temperature_interval_snowfall` (snowfall threshold interval length)
- `temperature_threshold_snowfall` (threshold temperature for snowfall)
- `rfcf` correction factor for liquid precipitation (rainfall)
- `sfcf` correction factor for snow precipitation (snowfall)

# Output
- `snow_precip`
- `liquid_precip`
"""
function precipitation_hbv(
    precipitation,
    temperature,
    temperature_interval_snowfall,
    temperature_threshold_snowfall;
    rfcf = 1.0,
    sfcf = 1.0,
)
    # fraction of precipitation which falls as rain
    rainfrac = if iszero(temperature_interval_snowfall)
        Float64(temperature > temperature_threshold_snowfall)
    else
        frac =
            (
                temperature -
                (temperature_threshold_snowfall - temperature_interval_snowfall / 2.0)
            ) / temperature_interval_snowfall
        clamp(frac, 0.0, 1.0)
    end

    # fraction of precipitation which falls as snow
    snowfrac = 1.0 - rainfrac
    # different correction for liquid_precip and snow_precip
    snow_precip = snowfrac * sfcf * precipitation
    liquid_precip = rainfrac * rfcf * precipitation
    return snow_precip, liquid_precip
end
