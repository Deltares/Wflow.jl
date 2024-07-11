"""
    scurve(x, a, b, c)

Sigmoid "S"-shaped curve.

# Arguments
- `x::Real`: input
- `a::Real`: determines the centre level
- `b::Real`: determines the amplitude of the curve
- `c::Real`: determines the steepness or "stepwiseness" of the curve.
             The higher c the sharper the function. A negative c reverses the function.
"""
function scurve(x, a, b, c)
    s = one(x) / (b + exp(-c * (x - a)))
    return s
end

"""
    acttransp_unsat_sbm(rootingdepth, ustorelayerdepth, sumlayer, restpotevap, sum_actevapustore, c, usl, theta_s, theta_r, hb, ust::Bool = false)

Compute actual transpiration for unsaturated zone.
If `ust` is `true`, the whole unsaturated store is available for transpiration.

# Arguments
- `rootingdepth`
- `ustorelayerdepth`
- `sumlayer` (depth (z) of upper boundary unsaturated layer)
- `restpotevap` (remaining evaporation)
- `sum_actevapustore` (cumulative actual transpiration (more than one unsaturated layers))
- `c` (Brooks-Corey coefficient)
- `usl` (thickness of unsaturated zone)
- `theta_s`
- `theta_r`
- `hb` (air entry pressure)
- `ust`

# Output
- `ustorelayerdepth`
- `sum_actevapustore`
- `restpotevap`
"""
function acttransp_unsat_sbm(
    rootingdepth,
    ustorelayerdepth,
    sumlayer,
    restpotevap,
    sum_actevapustore,
    c,
    usl,
    theta_s,
    theta_r,
    hb,
    ust::Bool = false,
)

    # AvailCap is fraction of unsat zone containing roots
    if ust
        availcap = ustorelayerdepth * 0.99
    else
        if usl > 0.0
            availcap = min(1.0, max(0.0, (rootingdepth - sumlayer) / usl))
        else
            availcap = 0.0
        end
    end

    maxextr = availcap * ustorelayerdepth

    # Next step is to make use of the Feddes curve in order to decrease actevapustore when soil moisture values
    # occur above or below ideal plant growing conditions (see also Feddes et al., 1978). h1-h4 values are
    # actually negative, but all values are made positive for simplicity.
    h1 = hb  # cm (air entry pressure)
    h2 = 100.0  # cm (pF 2 for field capacity)
    h3 = 400.0  # cm (pF 3, critical pF value)
    h4 = 15849.0  # cm (pF 4.2, wilting point)

    # According to Brooks-Corey
    par_lambda = 2.0 / (c - 3.0)
    if usl > 0.0
        vwc = ustorelayerdepth / usl
    else
        vwc = 0.0
    end
    vwc = max(vwc, 0.0000001)
    head = hb / (pow(((vwc) / (theta_s - theta_r)), (1.0 / par_lambda)))  # Note that in the original formula, thetaR is extracted from vwc, but thetaR is not part of the numerical vwc calculation
    head = max(head, hb)

    # Transform h to a reduction coefficient value according to Feddes et al. (1978).
    # For now: no reduction for head < h2 until following improvement (todo):
    #       - reduction only applied to crops
    if head <= h1
        alpha = 1.0
    elseif head >= h4
        alpha = 0.0
    elseif (head < h2) && (head > h1)
        alpha = 1.0
    elseif (head > h3) && (head < h4)
        alpha = 1.0 - (head - h3) / (h4 - h3)
    else
        alpha = 1.0
    end

    actevapustore = (min(maxextr, restpotevap, ustorelayerdepth)) * alpha
    ustorelayerdepth = ustorelayerdepth - actevapustore
    restpotevap = restpotevap - actevapustore
    sum_actevapustore = actevapustore + sum_actevapustore

    return ustorelayerdepth, sum_actevapustore, restpotevap
end

"""
    infiltration(avail_forinfilt, pathfrac, cf_soil, tsoil, infiltcapsoil, infiltcappath, ustorecapacity, modelsnow::Bool, soilinfreduction::Bool)

Soil infiltration based on infiltration capacity soil `infiltcapsoil`, infiltration capacity compacted area
`infiltcappath` and capacity unsatured zone `ustorecapacity`. The soil infiltration capacity can be adjusted
in case the soil is frozen (`modelsnow` and `soilinfreduction` is `true`).

"""
function infiltration(
    avail_forinfilt,
    pathfrac,
    cf_soil,
    tsoil,
    infiltcapsoil,
    infiltcappath,
    ustorecapacity,
    modelsnow::Bool,
    soilinfreduction::Bool,
)
    # First determine if the soil infiltration capacity can deal with the amount of water
    # split between infiltration in undisturbed soil and compacted areas (paths)
    soilinf = avail_forinfilt * (1.0 - pathfrac)
    pathinf = avail_forinfilt * pathfrac
    if modelsnow && soilinfreduction
        bb = 1.0 / (1.0 - cf_soil)
        soilinfredu = scurve(tsoil, Float(0.0), bb, Float(8.0)) + cf_soil
    else
        soilinfredu = 1.0
    end
    max_infiltsoil = min(infiltcapsoil * soilinfredu, soilinf)
    max_infiltpath = min(infiltcappath * soilinfredu, pathinf)
    infiltsoilpath = min(max_infiltpath + max_infiltsoil, max(0.0, ustorecapacity))

    if max_infiltpath + max_infiltsoil > 0.0
        infiltsoil =
            max_infiltsoil *
            min(1.0, max(0.0, ustorecapacity) / (max_infiltpath + max_infiltsoil))
        infiltpath =
            max_infiltpath *
            min(1.0, max(0.0, ustorecapacity) / (max_infiltpath + max_infiltsoil))
    else
        infiltsoil = 0.0
        infiltpath = 0.0
    end

    infiltexcess = (soilinf - max_infiltsoil) + (pathinf - max_infiltpath)

    return infiltsoilpath, infiltsoil, infiltpath, soilinf, pathinf, infiltexcess
end

"""
    unsatzone_flow_layer(usd, kv_z, l_sat, c)

Assuming a unit head gradient, the transfer of water from an unsaturated store layer `usd` is controlled by the
vertical saturated hydraulic conductivity `kv_z` (bottom layer or water table), the effective saturation
degree of the layer (ratio `usd` and `l_sat`), and a Brooks-Corey power coefficient `c`.
"""
function unsatzone_flow_layer(usd, kv_z, l_sat, c)
    if usd <= 0.0
        return 0.0, 0.0
    end
    sum_ast = 0.0
    # first transfer soil water > maximum soil water capacity layer (iteration is not required because of steady theta (usd))
    st = kv_z * min(pow(usd / l_sat, c), 1.0)
    st_sat = max(0.0, usd - l_sat)
    usd -= min(st, st_sat)
    sum_ast = sum_ast + min(st, st_sat)
    ast = max(min(st - min(st, st_sat), usd), 0.0)
    # number of iterations (to reduce "overshooting") based on fixed maximum change in soil water per iteration step (0.2 mm / model timestep)
    its = Int(cld(ast, 0.2))
    for _ in 1:its
        st = (kv_z / its) * min(pow(usd / l_sat, c), 1.0)
        ast = min(st, usd)
        usd -= ast
        sum_ast += ast
    end

    return usd, sum_ast
end

"""
    unsatzone_flow_sbm(ustorelayerdepth, soilwatercapacity, satwaterdepth, kv_z, usl, theta_s, theta_r)

The transfer of water from the unsaturated store `ustorelayerdepth` to the saturated store `satwaterdepth`
is controlled by the vertical saturated hydraulic conductivity `kv_z` at the water table and the ratio between
`ustorelayerdepth` and the saturation deficit (`soilwatercapacity` minus `satwaterdepth`). This is the
original Topog_SBM vertical transfer formulation.

"""
function unsatzone_flow_sbm(
    ustorelayerdepth,
    soilwatercapacity,
    satwaterdepth,
    kv_z,
    usl,
    theta_s,
    theta_r,
)
    sd = soilwatercapacity - satwaterdepth
    if sd <= 0.00001
        ast = 0.0
    else
        st = kv_z * min(ustorelayerdepth, usl * (theta_s - theta_r)) / sd
        ast = min(st, ustorelayerdepth)
        ustorelayerdepth = ustorelayerdepth - ast
    end

    return ustorelayerdepth, ast
end

"""
    snowpack_hbv(snow, snowwater, precipitation, temperature, tti, tt, ttm, cfmax, whc)

HBV type snowpack modeling using a temperature degree factor.
All correction factors (RFCF and SFCF) are set to 1.
The refreezing efficiency factor is set to 0.05.

# Arguments
- `snow` (snow storage)
- `snowwater` (liquid water content in the snow pack)
- `precipitation` (throughfall + stemflow)
- `temperature`
- `tti` (snowfall threshold interval length)
- `tt` (threshold temperature for snowfall)
- `ttm` (melting threshold)
- `cfmax` (degree day factor, rate of snowmelt)
- `whc` (Water holding capacity of snow)
- `rfcf` correction factor for rainfall
- `sfcf` correction factor for snowfall
- `cfr` refreeing efficiency constant in refreezing of freewater in snow

# Output
- `snow`
- `snowwater`
- `snowmelt`
- `rainfall` (precipitation that occurs as rainfall)
- `snowfall` (precipitation that occurs as snowfall)
"""
function snowpack_hbv(
    snow,
    snowwater,
    precipitation,
    temperature,
    tti,
    tt,
    ttm,
    cfmax,
    whc;
    rfcf = 1.0,
    sfcf = 1.0,
    cfr = 0.05,
)

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
    # different correction for rainfall and snowfall
    snowfall = snowfrac * sfcf * precipitation  # snowfall depth
    rainfall = rainfrac * rfcf * precipitation  # rainfall depth
    # potential snow melt, based on temperature
    potsnowmelt = temperature > ttm ? cfmax * (temperature - ttm) : 0.0
    # potential refreezing, based on temperature
    potrefreezing = temperature < ttm ? cfmax * cfr * (ttm - temperature) : 0.0
    # actual refreezing
    refreezing = temperature < ttm ? min(potrefreezing, snowwater) : 0.0

    # no landuse correction here
    snowmelt = min(potsnowmelt, snow)  # actual snow melt
    snow = snow + snowfall + refreezing - snowmelt  # dry snow content
    snowwater = snowwater - refreezing  # free water content in snow
    maxsnowwater = snow * whc  # max water in the snow
    snowwater = snowwater + snowmelt + rainfall  # add all water and potentially supersaturate the snowpack
    rainfall = max(snowwater - maxsnowwater, 0.0)  # rain + surpluss snowwater
    snowwater = snowwater - rainfall

    return snow, snowwater, snowmelt, rainfall, snowfall
end

"""
    glacier_hbv(glacierfrac, glacierstore, snow, temperature, tt, cfmax, g_sifrac, dt)

HBV-light type of glacier modelling.
First, a fraction of the snowpack is converted into ice using the HBV-light
model (fraction between 0.001-0.005 per day).
Glacier melting is modelled using a temperature degree factor and only
occurs if the snow cover < 10 mm.

# Arguments
- `glacierFrac` fraction covered by glaciers [-]
- `glacierstore` volume of the glacier [mm] w.e.
- `snow` snow pack on top of glacier [mm]
- `temperature` air temperature [°C]
- `tt` temperature threshold for ice melting [°C]
- `cfmax` ice degree-day factor in [mm/(°C/day)]
- `g_sifrac` fraction of the snow turned into ice [-]
- `dt` model timestep [s]

# Output
- `snow`
- `snow2glacier`
- `glacierstore`
- `glaciermelt`

"""
function glacier_hbv(glacierfrac, glacierstore, snow, temperature, tt, cfmax, g_sifrac, dt)

    # Fraction of the snow transformed into ice (HBV-light model)
    snow2glacier = g_sifrac * snow
    snow2glacier = glacierfrac > 0.0 ? snow2glacier : 0.0

    # Max conversion to 8mm/day
    snow2glacier = min(snow2glacier, 8.0 * (dt / basetimestep))

    snow = snow - (snow2glacier * glacierfrac)
    glacierstore = glacierstore + snow2glacier

    # Potential snow melt, based on temperature
    potmelt = temperature > tt ? cfmax * (temperature - tt) : 0.0

    # actual Glacier melt
    glaciermelt = snow < 10.0 ? min(potmelt, glacierstore) : 0.0
    glacierstore = glacierstore - glaciermelt

    return snow, snow2glacier, glacierstore, glaciermelt
end
