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
