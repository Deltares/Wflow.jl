"""
    glacier_hbv(glacier_frac, glacier_store, snow_storage, temperature, temperature_threshold_melt, degree_day_factor, snow_to_ice_fraction, maximum_snow_to_ice_rate)

HBV-light type of glacier modelling.
First, a fraction of the snowpack is converted into ice using the HBV-light
model (fraction between 0.001-0.005 per day).
Glacier melting is modelled using a temperature degree factor and only
occurs if the snow storage < 10 mm.

# Arguments
- `glacier_fraction` fraction covered by glaciers [-]
- `glacier_store` volume of the glacier [mm] w.e.
- `snow_storage` snow storage on top of glacier [mm]
- `temperature` air temperature [°C]
- `temperature_threshold_melt` temperature threshold for ice melting [°C]
- `degree_day_factor` ice degree-day factor in [mm/(°C/day)]
- `snow_to_ice_fraction` fraction of the snow turned into ice [-]
- `maximum_snow_to_ice_rate` maximum snow to glacier conversion rate

# Output
- `snow`
- `snow_to_glacier`
- `glacierstore`
- `glaciermelt`

"""
function glacier_hbv(
    glacier_fraction,
    glacier_store,
    snow_storage,
    temperature,
    temperature_threshold_melt,
    degree_day_factor,
    snow_to_ice_fraction,
    maximum_snow_to_ice_rate,
)

    # Fraction of the snow transformed into ice (HBV-light model)
    snow_to_glacier = if glacier_fraction > 0.0
        snow_to_ice_fraction * snow_storage
    else
        0.0
    end

    # Restrict snow_to_glacier conversion
    snow_to_glacier = min(snow_to_glacier, maximum_snow_to_ice_rate)

    snow_storage -= snow_to_glacier * glacier_fraction
    glacier_store += snow_to_glacier

    # Potential snow melt, based on temperature
    potential_melt =
        temperature > temperature_threshold_melt ?
        degree_day_factor * (temperature - temperature_threshold_melt) : 0.0

    # actual Glacier melt
    glacier_melt = snow_storage < 10.0 ? min(potential_melt, glacier_store) : 0.0
    glacier_store -= glacier_melt

    return snow_storage, snow_to_glacier, glacier_store, glacier_melt
end
