"""
    glacier_hbv(glacier_frac, glacier_store, snow_storage, temperature, ttm, cfmax, g_sifrac, max_snow_to_glacier)

HBV-light type of glacier modelling.
First, a fraction of the snowpack is converted into ice using the HBV-light
model (fraction between 0.001-0.005 per day).
Glacier melting is modelled using a temperature degree factor and only
occurs if the snow storage < 10 mm.

# Arguments
- `glacier_frac` fraction covered by glaciers [-]
- `glacier_store` volume of the glacier [mm => m] w.e.
- `snow_storage` snow storage on top of glacier [mm => m]
- `temperature` air temperature [°C => K]
- `ttm` temperature threshold for ice melting [°C => K]
- `cfmax` ice degree-day factor in [mm °C⁻¹ day⁻¹ => m K⁻¹ s⁻¹]
- `g_sifrac` fraction of the snow turned into ice [-]
- `max_snow_to_glacier` maximum snow to glacier conversion rate

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
    # [m s⁻¹]
    snow_to_glacier = if glacier_frac > 0.0
        # [s⁻¹] * [m]
        g_sifrac * snow_storage
    else
        0.0
    end

    # Restrict snow_to_glacier conversion
    # [m s⁻¹] = min([m s⁻¹], [m s⁻¹])
    snow_to_glacier = min(snow_to_glacier, max_snow_to_glacier)

    # [m] = [m s⁻¹] * [-] * [s]
    snow_storage -= snow_to_glacier * glacier_frac * dt
    glacier_store += snow_to_glacier * dt

    # Potential snow melt, based on temperature
    # [m s⁻¹] = [m K⁻¹ s⁻¹] * [K]
    potential_melt = (temperature > ttm) ? cfmax * (temperature - ttm) : 0.0

    # actual Glacier melt
    # [m s⁻¹] = min([m s⁻¹], [m] / [s])
    glacier_melt = (snow_storage < 1e-2) ? min(potential_melt, glacier_store / dt) : 0.0
    # [m] -= [m s⁻¹] * [s]
    glacier_store -= glacier_melt * dt

    return snow_storage, snow_to_glacier, glacier_store, glacier_melt
end
