"""
    glacier_hbv(glacierfrac, glacierstore, snow_storage, temperature, tt, cfmax, g_sifrac, max_snow_to_glacier)

HBV-light type of glacier modelling.
First, a fraction of the snowpack is converted into ice using the HBV-light
model (fraction between 0.001-0.005 per day).
Glacier melting is modelled using a temperature degree factor and only
occurs if the snow storage < 10 mm.

# Arguments
- `glacierFrac` fraction covered by glaciers [-]
- `glacierstore` volume of the glacier [mm] w.e.
- `snow_storage` snow storage on top of glacier [mm]
- `temperature` air temperature [°C]
- `ttm` temperature threshold for ice melting [°C]
- `cfmax` ice degree-day factor in [mm/(°C/day)]
- `g_sifrac` fraction of the snow turned into ice [-]
- `max_snow_to_glacier` maximum snow to glacier conversion rate

# Output
- `snow`
- `snow_to_glacier`
- `glacierstore`
- `glaciermelt`

"""
function glacier_hbv(
    glacierfrac,
    glacierstore,
    snow,
    temperature,
    ttm,
    cfmax,
    g_sifrac,
    max_snow_to_glacier,
    dt,
)

    # Fraction of the snow transformed into ice (HBV-light model)
    # [m s⁻¹]
    snow_to_glacier = if glacierfrac > 0.0
        # [s⁻¹] * [m]
        g_sifrac * snow
    else
        0.0
    end

    # Restrict snow_to_glacier conversion
    snow_to_glacier = min(snow_to_glacier, max_snow_to_glacier)

    # [m] = [m s⁻¹] * [-] * [s]
    snow -= snow_to_glacier * glacierfrac * dt
    glacierstore += snow_to_glacier

    # Potential snow melt, based on temperature
    # [m s⁻¹] = [m K⁻¹ s⁻¹] * [K]
    potmelt = (temperature > ttm) ? cfmax * (temperature - ttm) : 0.0

    # actual Glacier melt
    # [m] = [m s⁻¹] * [s]
    glaciermelt = (snow < 10.0) ? min(potmelt * dt, glacierstore) : 0.0
    glacierstore -= glaciermelt

    return snow, snow_to_glacier, glacierstore, glaciermelt
end
