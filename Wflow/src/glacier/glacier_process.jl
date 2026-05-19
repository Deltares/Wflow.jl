"""
    glacier_hbv(glacier_frac, glacier_store, snow_storage, temperature, ttm, cfmax, g_sifrac, max_snow_to_glacier)

HBV-light type of glacier modelling.
First, a fraction of the snowpack is converted into ice using the HBV-light
model (fraction between 0.001-0.005 per day).
Glacier melting is modelled using a temperature degree factor and only
occurs if the snow storage < 10 mm.

# Arguments
- `glacier_frac` fraction covered by glaciers [-]
- `glacier_store` volume of the glacier [m] w.e.
- `snow_storage` snow storage on top of glacier [m]
- `temperature` air temperature [K]
- `ttm` temperature threshold for ice melting [K]
- `cfmax` ice degree-day factor in [m K⁻¹ s⁻¹]
- `g_sifrac` fraction of the snow turned into ice [s⁻¹]
- `max_snow_to_glacier` maximum snow to glacier conversion rate [m s⁻¹]

# Output
- `snow`
- `snow_to_glacier`
- `glacierstore`
- `glaciermelt`

"""
function glacier_hbv(
    glacier_frac,
    glacier_store,
    snow_storage,
    temperature,
    ttm,
    cfmax,
    g_sifrac,
    max_snow_to_glacier,
    dt,
)

    # Fraction of the snow transformed into ice (HBV-light model)
    snow_to_glacier = if glacier_frac > 0.0
        g_sifrac * snow_storage
    else
        0.0
    end

    # Restrict snow_to_glacier conversion
    snow_to_glacier = min(snow_to_glacier, max_snow_to_glacier)

    snow_storage -= snow_to_glacier * glacier_frac * dt
    glacier_store += snow_to_glacier * dt

    # Potential snow melt, based on temperature
    potential_melt = (temperature > ttm) ? cfmax * (temperature - ttm) : 0.0

    # actual Glacier melt
    glacier_melt = (snow_storage < 1e-2) ? min(potential_melt, glacier_store / dt) : 0.0
    glacier_store -= glacier_melt * dt

    return snow_storage, snow_to_glacier, glacier_store, glacier_melt
end
