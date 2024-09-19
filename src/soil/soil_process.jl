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
    tpot_daily = tpot * (tosecond(basetimestep) / Δt)
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

function soil_fraction!(soil, runoff, glacier)
    (; canopygapfraction) = soil.parameters.vegetation_parameter_set
    (; soil_fraction) = soil.parameters
    (; waterfrac, riverfrac) = runoff.parameters
    glacier_fraction = get_glacier_fraction(glacier)

    @. soil_fraction =
        max(canopygapfraction - waterfrac - riverfrac - glacier_fraction, 0.0)
end