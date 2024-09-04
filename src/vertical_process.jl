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
    rainfall_interception_gash(cmax, e_r, canopygapfraction, precipitation, canopystorage, maxevap)

Interception according to the Gash model (for daily timesteps). `cmax` is the maximum canopy storage and
`e_r` is the ratio of the average evaporation from the wet canopy and the average precipitation intensity on
a saturated canopy.
"""
function rainfall_interception_gash(
    cmax,
    e_r,
    canopygapfraction,
    precipitation,
    canopystorage,
    maxevap,
)
    # TODO: add other rainfall interception method (lui)
    # TODO: include subdaily Gash model
    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopygapfraction)
    pt = min(0.1 * canopygapfraction, 1.0 - canopygapfraction)
    pfrac = 1.0 - canopygapfraction - pt
    p_sat = (-cmax / e_r) * log(1.0 - min((e_r / pfrac), 1.0))
    p_sat = isinf(p_sat) ? 0.0 : p_sat

    # large storms P > P_sat
    largestorms = precipitation > p_sat

    iwet = largestorms ? (pfrac * p_sat) - cmax : precipitation * pfrac
    isat = largestorms ? (e_r) * (precipitation - p_sat) : 0.0
    idry = largestorms ? cmax : 0.0
    itrunc = 0.0

    stemflow = pt * precipitation

    throughfall = precipitation - iwet - idry - isat - itrunc - stemflow
    interception = iwet + idry + isat + itrunc

    # Now corect for area without any Interception (say open water Cmax -- zero)
    cmaxzero = cmax <= 0.0
    throughfall = cmaxzero ? precipitation : throughfall
    interception = cmaxzero ? 0.0 : interception
    stemflow = cmaxzero ? 0.0 : stemflow

    # Now corect for maximum potential evap
    canopy_drainage = interception > maxevap ? interception - maxevap : 0.0
    interception = min(interception, maxevap)

    # Add surpluss to the throughfall
    throughfall = throughfall + canopy_drainage

    return throughfall, interception, stemflow, canopystorage
end

"""
    rainfall_interception_modrut(precipitation, potential_evaporation, canopystorage, canopygapfraction, cmax)

Interception according to a modified Rutter model. The model is solved explicitly and there is no
drainage below `cmax`.
"""
function rainfall_interception_modrut(
    precipitation,
    potential_evaporation,
    canopystorage,
    canopygapfraction,
    cmax,
)

    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopygapfraction)
    pt = min(0.1 * canopygapfraction, 1.0 - canopygapfraction)

    # Amount of p that falls on the canopy
    precip_canopy = (1.0 - canopygapfraction - pt) * precipitation

    # Canopystorage cannot be larger than cmax, no gravity drainage below that. This check
    # is required because cmax can change over time
    canopy_drainage1 = canopystorage > cmax ? canopystorage - cmax : 0.0
    canopystorage = canopystorage - canopy_drainage1

    # Add the precipitation that falls on the canopy to the store
    canopystorage = canopystorage + precip_canopy

    # Now do the Evap, make sure the store does not get negative
    canopy_evap = min(canopystorage, potential_evaporation)
    canopystorage = canopystorage - canopy_evap

    # Amount of evap not used
    leftover = potential_evaporation - canopy_evap

    # Now drain the canopystorage again if needed...
    canopy_drainage2 = canopystorage > cmax ? canopystorage - cmax : 0.0
    canopystorage = canopystorage - canopy_drainage2

    # Calculate throughfall and stemflow
    throughfall = canopy_drainage1 + canopy_drainage2 + canopygapfraction * precipitation
    stemflow = precipitation * pt

    # Calculate interception, this is NET Interception
    netinterception = precipitation + canopy_drainage1 - throughfall - stemflow
    interception = canopy_evap

    return netinterception, throughfall, stemflow, leftover, interception, canopystorage
end

"""
    vwc_brooks_corey(h, hb, theta_s, theta_r, c)

Volumetric water content based on the Brooks-Corey soil hydraulic model.
"""
function vwc_brooks_corey(h, hb, theta_s, theta_r, c)
    if h < hb
        par_lambda = 2.0 / (c - 3.0)
        vwc = (theta_s - theta_r) * pow(hb / h, par_lambda) + theta_r
    else
        vwc = theta_s
    end
    return vwc
end

"""
    head_brooks_corey(vwc, theta_s, theta_r, c, hb)

Soil water pressure head based on the Brooks-Corey soil hydraulic model.
"""
function head_brooks_corey(vwc, theta_s, theta_r, c, hb)
    par_lambda = 2.0 / (c - 3.0)
    # Note that in the original formula, theta_r is extracted from vwc, but theta_r is not part of the numerical vwc calculation
    h = hb / (pow(((vwc) / (theta_s - theta_r)), (1.0 / par_lambda)))
    h = min(h, hb)
    return h
end

"""
    feddes_h3(h3_high, h3_low, tpot, Δt)

Return soil water pressure head `h3` of Feddes root water uptake reduction function.
"""
function feddes_h3(h3_high, h3_low, tpot, Δt)
    # value of h3 is a function of potential transpiration [mm/d]
    tpot_daily = tpot * (basetimestep / Δt)
    if (tpot_daily >= 0.0) && (tpot_daily <= 1.0)
        h3 = h3_low
    elseif (tpot_daily > 1.0) && (tpot_daily < 5.0)
        h3 = h3_high + ((h3_low - h3_high) * (5.0 - tpot_daily)) / (5.0 - 1.0)
    else
        h3 = h3_high
    end
    return h3
end

"""
    rwu_reduction_feddes(h, h1, h2, h3, h4, alpha_h1)

Root water uptake reduction factor based on Feddes.
"""
function rwu_reduction_feddes(h, h1, h2, h3, h4, alpha_h1)
    # root water uptake reduction coefficient alpha (see also Feddes et al., 1978)
    if alpha_h1 == 0.0
        if (h <= h4) || (h > h1)
            alpha = 0.0
        elseif (h > h2) && (h <= h1)
            alpha = (h - h1) / (h2 - h1)
        elseif (h >= h3) && (h <= h2)
            alpha = 1.0
        elseif (h >= h4) && (h < h3)
            alpha = (h - h4) / (h3 - h4)
        end
    else
        if h <= h4
            alpha = 0.0
        elseif h >= h3
            alpha = 1.0
        elseif (h >= h4) && (h < h3)
            alpha = (h - h4) / (h3 - h4)
        end
    end
    return alpha
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

    return infiltsoilpath,
    infiltsoil,
    infiltpath,
    soilinf,
    pathinf,
    infiltexcess,
    soilinfredu
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
