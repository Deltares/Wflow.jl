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
- `tt` temperature threshold for ice melting [°C]
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
    tt,
    cfmax,
    g_sifrac,
    max_snow_to_glacier,
)

    # Fraction of the snow transformed into ice (HBV-light model)
    snow_to_glacier = g_sifrac * snow
    snow_to_glacier = glacierfrac > 0.0 ? snow_to_glacier : 0.0

    # Restrict snow_to_glacier conversion
    snow_to_glacier = min(snow_to_glacier, max_snow_to_glacier)

    snow = snow - (snow_to_glacier * glacierfrac)
    glacierstore = glacierstore + snow_to_glacier

    # Potential snow melt, based on temperature
    potmelt = temperature > tt ? cfmax * (temperature - tt) : 0.0

    # actual Glacier melt
    glaciermelt = snow < 10.0 ? min(potmelt, glacierstore) : 0.0
    glacierstore = glacierstore - glaciermelt

    return snow, snow_to_glacier, glacierstore, glaciermelt
end
