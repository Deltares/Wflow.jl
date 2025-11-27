"""
    infiltration(
        potential_infiltration,
        pathfrac,
        infiltcapsoil,
        infiltcappath,
        ustorecapacity,
        f_infiltration_reduction,
    )

Soil infiltration based on infiltration capacity soil `infiltcapsoil`, infiltration capacity
paved area `infiltcappath` and unsaturated store capacity `ustorecapacity`. The infiltration
capacity of the soil and paved area can be reduced with the infiltration reduction factor
`f_infiltration_reduction`.
"""
function infiltration(
    potential_infiltration,
    pathfrac,
    infiltcapsoil,
    infiltcappath,
    ustorecapacity,
    f_infiltration_reduction,
)
    # First determine if the soil infiltration capacity can deal with the amount of water
    # split between infiltration in undisturbed soil and paved areas (path).
    soilinf = potential_infiltration * (1.0 - pathfrac)
    pathinf = potential_infiltration * pathfrac

    max_infiltsoil = min(infiltcapsoil * f_infiltration_reduction, soilinf)
    max_infiltpath = min(infiltcappath * f_infiltration_reduction, pathinf)
    infiltsoilpath = min(max_infiltpath + max_infiltsoil, max(0.0, ustorecapacity))

    infiltexcess = (soilinf - max_infiltsoil) + (pathinf - max_infiltpath)

    return infiltsoilpath, infiltexcess
end

"""
    unsatzone_flow_layer(usd, kv_z, l_sat, c)

Assuming a unit head gradient, the transfer of water from an unsaturated store layer `usd`
is controlled by the vertical saturated hydraulic conductivity `kv_z` (bottom layer or water
table), the effective saturation degree of the layer (ratio `usd` and `l_sat`), and a
Brooks-Corey power coefficient `c`.
"""
function unsatzone_flow_layer(usd, kv_z, l_sat, c)
    if usd <= 0.0
        return 0.0, 0.0
    end
    # Excess soil water:
    # first transfer soil water > maximum soil water capacity layer (iteration is not
    # required because of steady theta (usd))
    st_sat = max(0.0, usd - l_sat)
    st = kv_z * min(pow(usd / l_sat, c), 1.0)
    sum_ast = min(st, st_sat)
    usd -= sum_ast

    # number of iterations (to reduce "overshooting") based on fixed maximum change in soil
    # water per iteration step (0.2 mm / model timestep)
    remainder = min(st - sum_ast, usd)
    its = Int(cld(remainder, 0.2))
    for _ in 1:its
        st = (kv_z / its) * min(pow(usd / l_sat, c), 1.0)
        ast = min(st, usd)
        usd -= ast
        sum_ast += ast
    end

    return usd, sum_ast
end

"""
    vwc_brooks_corey(h, hb, theta_s, theta_r, c)

Return volumetric water content based on the Brooks-Corey soil hydraulic model.
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

Return soil water pressure head based on the Brooks-Corey soil hydraulic model.
"""
function head_brooks_corey(vwc, theta_s, theta_r, c, hb)
    par_lambda = 2.0 / (c - 3.0)
    # Note that in the original formula, theta_r is extracted from vwc, but theta_r is not
    # part of the numerical vwc calculation
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
    tpot_daily = tpot * (tosecond(BASETIMESTEP) / Δt)
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
    soil_temperature(tsoil, w_soil, temperature))

Return the near surface soil temperature `tsoil` based on the near surface soil temperature
`tsoil_prev` at the previous timestep, and the difference between air `temperature` and near
surface soil temperature `tsoil_prev` at the previous timestep, weighted with the weighting
coefficient `w_soil` (Wigmosta et al., 2009).
"""
function soil_temperature(tsoil_prev, w_soil, temperature)
    tsoil = tsoil_prev + w_soil * (temperature - tsoil_prev)
    return tsoil
end

"""
    infiltration_reduction_factor(
        tsoil,
        cf_soil;
        modelsnow = false,
        soil_infiltration_reduction = false,
    )

When both `modelsnow` and `soil_infiltration_reduction` are `true` an infiltration reduction
factor `f_infiltration_reduction` is computed. The infiltration reduction factor is based on
the near surface soil temperature `tsoil`, parameter `cf_soil` and a s-curve to make a
smooth transition of `f_infiltration_reduction` as a function of `tsoil` and `cf_soil`.
Otherwise, `f_infiltration_reduction` is set to 1.0.
"""
function infiltration_reduction_factor(
    tsoil,
    cf_soil;
    modelsnow = false,
    soil_infiltration_reduction = false,
)
    if modelsnow && soil_infiltration_reduction
        bb = 1.0 / (1.0 - cf_soil)
        f_infiltration_reduction = scurve(tsoil, 0.0, bb, 8.0) + cf_soil
    else
        f_infiltration_reduction = 1.0
    end
    return f_infiltration_reduction
end

"Return soil evaporation from the unsaturated store"
function soil_evaporation_unsatured_store(
    potential_soilevaporation,
    ustorelayerdepth,
    ustorelayerthickness,
    n_unsatlayers,
    zi,
    theta_effective,
)
    if n_unsatlayers == 0
        soilevapunsat = 0.0
    elseif n_unsatlayers == 1
        # Check if groundwater level lies below the surface
        soilevapunsat =
            potential_soilevaporation * min(1.0, ustorelayerdepth / (zi * theta_effective))
    else
        # In case first layer contains no saturated storage
        soilevapunsat =
            potential_soilevaporation *
            min(1.0, ustorelayerdepth / (ustorelayerthickness * (theta_effective)))
    end
    return soilevapunsat
end

"Return soil evaporation from the saturated store"
function soil_evaporation_satured_store(
    potential_soilevaporation,
    n_unsatlayers,
    layerthickness,
    zi,
    theta_effective,
)
    if n_unsatlayers == 0 || n_unsatlayers == 1
        soilevapsat =
            potential_soilevaporation * min(1.0, (layerthickness - zi) / layerthickness)
        soilevapsat = min(soilevapsat, (layerthickness - zi) * theta_effective)
    else
        soilevapsat = 0.0
    end
    return soilevapsat
end

"Return actual infiltration rate for soil `actinfiltsoil` and paved area `actinfiltpath`"
function actual_infiltration_soil_path(
    potential_infiltration,
    actinfilt,
    pathfrac,
    infiltcapsoil,
    infiltcappath,
    f_infiltration_reduction,
)
    soilinf = potential_infiltration * (1.0 - pathfrac)
    pathinf = potential_infiltration * pathfrac
    if actinfilt > 0.0
        max_infiltsoil = min(infiltcapsoil * f_infiltration_reduction, soilinf)
        max_infiltpath = min(infiltcappath * f_infiltration_reduction, pathinf)

        actinfiltsoil = actinfilt * max_infiltsoil / (max_infiltpath + max_infiltsoil)
        actinfiltpath = actinfilt * max_infiltpath / (max_infiltpath + max_infiltsoil)

    else
        actinfiltsoil = 0.0
        actinfiltpath = 0.0
    end

    return actinfiltsoil, actinfiltpath
end
