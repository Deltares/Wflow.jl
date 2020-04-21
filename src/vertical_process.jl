"""
sCurve function:

Input:
    - x input
    - c determines the steepness or "stepwiseness" of the curve.
      The higher c the sharper the function. A negative c reverses the function.
    - b determines the amplitude of the curve
    - a determines the centre level (default = 0)

Output:
    - result
"""
function scurve(x::Float64; a = 0.0, b = 1.0, c = 1.0)
    s = 1.0 / (b + exp(-c * (x - a)))
    return s
end

"""
Interception according to the Gash model (For daily timesteps).
"""
function rainfall_interception_gash(
    cmax::Float64,
    e_r::Float64,
    canopygapfraction::Float64,
    precipitation::Float64,
    canopystorage::Float64;
    maxevap = 9999.0,
)
    # TODO:  add other rainfall interception method (lui)
    # TODO: Include subdaily Gash model
    # TODO: add LAI variation in year
    # Hack for stemflow

    pt = 0.1 * canopygapfraction
    p_sat = max(
        0.0,
        (-cmax / e_r) * log(1.0 - (e_r / (1.0 - canopygapfraction - pt))),
    )
    # large storms P > P_sat
    largestorms = precipitation > p_sat

    iwet = largestorms ? ((1 - canopygapfraction - pt) * p_sat) - cmax :
        precipitation * (1 - canopygapfraction - pt)
    isat = largestorms ? (e_r) * (precipitation - p_sat) : 0.0
    idry = largestorms ? cmax : 0.0
    itrunc = 0

    stemflow = pt * precipitation

    throughfall = precipitation - iwet - idry - isat - itrunc - stemflow
    interception = iwet + idry + isat + itrunc

    # Now corect for area without any Interception (say open water Cmax -- zero)
    cmaxzero = cmax <= 0.0
    throughfall = cmaxzero ? precipitation : throughfall
    interception = cmaxzero ? 0.0 : interception
    stemflow = cmaxzero ? 0.0 : stemflow

    # Now corect for maximum potential evap
    overestimate = interception > maxevap ? interception - maxevap : 0.0
    interception = min(interception, maxevap)

    # Add surpluss to the thoughdfall
    throughfall = throughfall + overestimate

    return throughfall, interception, stemflow, canopystorage

end

"""
Actual transpiration function for unsaturated zone:

  if ust is True, all ustore is available for transpiration

Input:

    - RootingDepth, UStoreLayerDepth, sumLayer (depth (z) of upper boundary unsaturated layer),
      RestPotEvap (remaining evaporation), sumActEvapUStore (cumulative actual transpiration (more than one UStore layers))
      c (Brooks-Corey coefficient), L (thickness of unsaturated zone), thetaS, thetaR, hb (air entry pressure), ust

Output:

    - UStoreLayerDepth,  sumActEvapUStore, ActEvapUStore
"""
function acttransp_unsat_sbm(
    rootingdepth::Float64,
    ustorelayerdepth::Float64,
    sumlayer::Float64,
    restpotevap::Float64,
    sum_actevapustore::Float64,
    c::Float64,
    usl::Float64,
    θₛ::Float64,
    θᵣ::Float64,
    hb::Float64,
    ust = 0::Int,
)

    # AvailCap is fraction of unsat zone containing roots
    if ust >= 1
        availcap = ustorelayerdepth * 0.99
    else
        if usl > 0
            availcap = min(1.0, max(0.0, (rootingdepth - sumlayer) / usl))
        else
            availcap = 0.0
        end
    end

    maxextr = availcap * ustorelayerdepth

    # Next step is to make use of the Feddes curve in order to decrease ActEvapUstore when soil moisture values
    # occur above or below ideal plant growing conditions (see also Feddes et al., 1978). h1-h4 values are
    # actually negative, but all values are made positive for simplicity.
    h1 = hb  # cm (air entry pressure)
    h2 = 100.0  # cm (pF 2 for field capacity)
    h3 = 400.0  # cm (pF 3, critical pF value)
    h4 = 15849.0  # cm (pF 4.2, wilting point)

    # According to Brooks-Corey
    par_lambda = 2 / (c - 3)
    if usl > 0.0
        vwc = ustorelayerdepth / usl
    else
        vwc = 0.0
    end
    vwc = max(vwc, 0.0000001)
    head = hb / (((vwc) / (θₛ - θᵣ))^(1 / par_lambda))  # Note that in the original formula, thetaR is extracted from vwc, but thetaR is not part of the numerical vwc calculation
    head = max(head, hb)

    # Transform h to a reduction coefficient value according to Feddes et al. (1978).
    # For now: no reduction for head < h2 until following improvement (todo):
    #       - reduction only applied to crops
    if (head <= h1)
        alpha = 1
    elseif (head >= h4)
        alpha = 0
    elseif ((head < h2) & (head > h1))
        alpha = 1
    elseif ((head > h3) & (head < h4))
        alpha = 1 - (head - h3) / (h4 - h3)
    else
        alpha = 1
    end

    actevapustore = (min(maxextr, restpotevap, ustorelayerdepth)) * alpha
    ustorelayerdepth = ustorelayerdepth - actevapustore
    restpotevap = restpotevap - actevapustore
    sum_actevapustore = actevapustore + sum_actevapustore

    return ustorelayerdepth, sum_actevapustore, restpotevap
end


function infiltration(
    avail_forinfilt::Float64,
    pathfrac::Float64,
    cf_soil::Float64,
    tsoil::Float64,
    infiltcapsoil::Float64,
    infiltcappath::Float64,
    ustorecapacity::Float64,
    modelsnow::Bool,
    soilinfreduction::Bool,
)
    # First determine if the soil infiltration capacity can deal with the amount of water
    # split between infiltration in undisturbed soil and compacted areas (paths)
    soilinf = avail_forinfilt * (1 - pathfrac)
    pathinf = avail_forinfilt * pathfrac
    if modelsnow & soilinfreduction
        bb = 1.0 / (1.0 - cf_soil)
        soilinfredu = scurve(tsoil, a = 0.0, b = bb, c = 8.0) + cf_soil
    else
        soilinfredu = 1.0
    end
    max_infiltsoil = min(infiltcapsoil * soilinfredu, soilinf)
    max_infiltpath = min(infiltcappath * soilinfredu, pathinf)
    infiltsoilpath =
        min(max_infiltpath + max_infiltsoil, max(0.0, ustorecapacity))

    if max_infiltpath + max_infiltsoil > 0.0
        infiltsoil =
            max_infiltsoil * min(
                1.0,
                max(0.0, ustorecapacity) / (max_infiltpath + max_infiltsoil),
            )
        infiltpath =
            max_infiltpath * min(
                1.0,
                max(0.0, ustorecapacity) / (max_infiltpath + max_infiltsoil),
            )
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
    infiltexcess
end

function unsatzone_flow_layer(
    usd::Float64,
    kvfrac::Float64,
    kv::Float64,
    f::Float64,
    z::Float64,
    l_sat::Float64,
    c::Float64,
)

    sum_ast = 0
    #first transfer soil water > maximum soil water capacity layer (iteration is not required because of steady theta (usd))
    st = kvfrac * kv * exp(-f * z) * min((usd / l_sat)^c, 1.0)
    st_sat = max(0, usd - l_sat)
    usd = usd - min(st, st_sat)
    sum_ast = sum_ast + min(st, st_sat)
    ast = min(st - min(st, st_sat), usd)
    #number of iterations (to reduce "overshooting") based on fixed maximum change in soil water per iteration step (0.02 x maximum soil water capacity layer)
    its = max(Int(ceil(ast / (l_sat * 0.02))), 1)
    k = kvfrac * kv / its * exp(-f * z)
    for i = 1:its
        st = k * min((usd / l_sat)^c, 1.0)
        ast = min(st, usd)
        usd = usd - ast
        sum_ast = sum_ast + ast
    end

    return usd, sum_ast
end


function unsatzone_flow_sbm(
    ustorelayerdepth::Float64,
    soilwatercapacity::Float64,
    satwaterdepth::Float64,
    kvfrac::Float64,
    kv::Float64,
    usl::Float64,
    θₛ::Float64,
    θᵣ::Float64,
)

    sd = soilwatercapacity - satwaterdepth
    if sd <= 0.00001
        ast = 0.0
    else
        st = kvfrac * kv * min(ustorelayerdepth, usl * (θₛ - θᵣ)) / sd
        ast = min(st, ustorelayerdepth)
        ustorelayerdepth = ustorelayerdepth - ast
    end

    return ustorelayerdepth, ast

end


"""
snowpack(snow, snowwater, precipitation, temperature, tti, tt, ttm, cfmax, whc)

HBV type snowpack modeling using a temperature degree factor.
All correction factors (RFCF and SFCF) are set to 1.
The refreezing efficiency factor is set to 0.05.
"""
function snowpack_hbv(
    snow::Float64,
    snowwater::Float64,
    precipitation::Float64,
    temperature::Float64,
    tti::Float64,
    tt::Float64,
    ttm::Float64,
    cfmax::Float64,
    whc::Float64,
)

    rfcf = 1.0  # correction factor for rainfall
    cfr = 0.05  # refreeing efficiency constant in refreezing of freewater in snow
    sfcf = 1.0  # correction factor for snowfall

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
    # different correction for rainfall and snowfall
    precipitation =
        sfcf * snowfrac * precipitation + rfcf * rainfrac * precipitation

    snowfall = snowfrac * precipitation  #snowfall depth
    rainfall = rainfrac * precipitation  #rainfall depth
    # potential snow melt, based on temperature
    potsnowmelt = temperature > ttm ? cfmax * (temperature - ttm) : 0.0
    # potential refreezing, based on temperature
    potrefreezing = temperature < ttm ? cfmax * cfr * (ttm - temperature) : 0.0
    # actual refreezing
    refreezing = temperature < ttm ? min(potrefreezing, snowwater) : 0.0

    # no landuse correction here
    snowmelt = min(potsnowmelt, snow)  #actual snow melt
    snow = snow + snowfall + refreezing - snowmelt  #dry snow content
    snowwater = snowwater - refreezing  #free water content in snow
    maxsnowwater = snow * whc  # max water in the snow
    snowwater = snowwater + snowmelt + rainfall  # add all water and potentially supersaturate the snowpack
    rainfall = max(snowwater - maxsnowwater, 0.0)  # rain + surpluss snowwater
    snowwater = snowwater - rainfall

    return snow, snowwater, snowmelt, rainfall, snowfall
end
