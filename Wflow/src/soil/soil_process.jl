"""
    infiltration(
        potential_infiltration,
        compacted_soil_area_fraction,
        infiltration_capacity_soil,
        infiltration_capacity_compacted_soil,
        unsaturated_store_capacity,
        f_infiltration_reduction,
    )

Soil infiltration based on infiltration capacity soil `infiltration_capacity_soil`, infiltration capacity
paved area `infiltration_capacity_compacted_soil` and unsaturated store capacity `unsaturated_store_capacity`. The infiltration
capacity of the soil and paved area can be reduced with the infiltration reduction factor
`f_infiltration_reduction`.
"""
function infiltration(
    potential_infiltration,
    compacted_soil_area_fraction,
    infiltration_capacity_soil,
    infiltration_capacity_compacted_soil,
    unsaturated_store_capacity,
    f_infiltration_reduction,
    dt,
)
    # First determine if the soil infiltration capacity can deal with the amount of water
    # split between infiltration in undisturbed soil and paved areas (path).
    soilinf = potential_infiltration * (1.0 - compacted_soil_area_fraction)
    pathinf = potential_infiltration * compacted_soil_area_fraction

    max_infiltsoil = min(infiltration_capacity_soil * f_infiltration_reduction, soilinf)

    max_infiltpath =
        min(infiltration_capacity_compacted_soil * f_infiltration_reduction, pathinf)

    infiltration =
        min(max_infiltpath + max_infiltsoil, max(0.0, unsaturated_store_capacity / dt))

    infiltration_excess = (soilinf - max_infiltsoil) + (pathinf - max_infiltpath)

    return infiltration, infiltration_excess
end

"""
    unsatzone_flow_layer(usd, kv_z, l_sat, brooks_corey_exponent)

Assuming a unit head gradient, the transfer of water from an unsaturated store layer `usd`
is controlled by the vertical saturated hydraulic conductivity `kv_z` (bottom layer or water
table), the effective saturation degree of the layer (ratio `usd` and `l_sat`), and a
Brooks-Corey power coefficient `brooks_corey_exponent`.
"""
function unsatzone_flow_layer(
    unsaturated_layer_depth,
    kv_z,
    l_sat,
    brooks_corey_exponent,
    dt,
)
    if unsaturated_layer_depth <= 0.0
        return 0.0, 0.0
    end
    # Excess soil water:
    # first transfer soil water > maximum soil water capacity layer (iteration is not
    # required because of steady theta (unsaturated_layer_depth))

    st_sat = max(0.0, unsaturated_layer_depth - l_sat)

    st = kv_z * bounded_power(unsaturated_layer_depth / l_sat, brooks_corey_exponent)
    sum_ast = min(st, st_sat / dt)
    unsaturated_layer_depth -= sum_ast * dt

    # number of iterations (to reduce "overshooting") based on fixed maximum change in soil
    # water per iteration step (0.2 mm / model timestep)
    remainder = min((st - sum_ast) * dt, unsaturated_layer_depth)
    its = Int(cld(remainder, 2e-4))
    for _ in 1:its
        st =
            (kv_z / its) *
            bounded_power(unsaturated_layer_depth / l_sat, brooks_corey_exponent)
        st_max = unsaturated_layer_depth / dt

        if st < st_max
            unsaturated_layer_depth -= st * dt
            sum_ast += st
        else
            unsaturated_layer_depth = 0
            sum_ast += st_max
            break
        end
    end

    return unsaturated_layer_depth, sum_ast
end

"""
    vwc_brooks_corey(h, air_entry_pressure, theta_s, theta_r, brooks_corey_exponent)

Return volumetric water content based on the Brooks-Corey soil hydraulic model.
"""
function vwc_brooks_corey(h, air_entry_pressure, theta_s, theta_r, brooks_corey_exponent)
    return if h < air_entry_pressure
        par_lambda = 2.0 / (brooks_corey_exponent - 3.0)
        (theta_s - theta_r) * pow(air_entry_pressure / h, par_lambda) + theta_r
    else
        theta_s
    end
end

"""
    head_brooks_corey(volumetric_water_content, theta_s, theta_r, brooks_corey_exponent, air_entry_pressure)

Return soil water pressure head based on the Brooks-Corey soil hydraulic model.
"""
function head_brooks_corey(
    volumetric_water_content,
    theta_s,
    theta_r,
    brooks_corey_exponent,
    air_entry_pressure,
)
    par_lambda = 2 / (brooks_corey_exponent - 3.0)
    h = if par_lambda > 0
        # Note that in the original formula, theta_r is extracted from volumetric_water_content, but theta_r is not
        # part of the numerical volumetric_water_content calculation
        air_entry_pressure /
        pow(volumetric_water_content / (theta_s - theta_r), inv(par_lambda))
    else
        air_entry_pressure
    end
    return h
end

"""
    field_capacity(layer_thickness, n_layers, theta_s, theta_r, brooks_corey_exponent, air_entry_pressure)

Return water content at field capacity based on the Brooks-Corey soil hydraulic model.
"""
function field_capacity(
    layer_thickness,
    n_layers,
    theta_s,
    theta_r,
    brooks_corey_exponent,
    air_entry_pressure,
)
    theta_fc = 0.0
    total_depth = 0.0
    for i in 1:n_layers
        theta_fc +=
            vwc_brooks_corey(
                -1.0,
                air_entry_pressure,
                theta_s,
                theta_r,
                brooks_corey_exponent[i],
            ) * layer_thickness[i]
        total_depth += layer_thickness[i]
    end
    return theta_fc / total_depth
end

"""
    feddes_h3(h3_high, h3_low, tpot, Δt)

Return soil water pressure head `h3` of Feddes root water uptake reduction function.
"""
function feddes_h3(h3_high, h3_low, tpot_SI)
    # value of h3 is a function of potential transpiration [mm d⁻¹]
    tpot_daily = from_SI(tpot_SI, MM_PER_DAY)
    return if tpot_daily <= 1.0
        h3_low
    elseif tpot_daily < 5.0
        h3_low + (h3_high - h3_low) * (tpot_daily - 1.0) / (5.0 - 1.0)
    else
        h3_high
    end
end

"""
    rwu_reduction_feddes(h, h1, h2, h3, h4, alpha_h1)

Root water uptake reduction factor based on Feddes.
"""
function rwu_reduction_feddes(h, h1, h2, h3, h4, alpha_h1)
    # root water uptake reduction coefficient alpha (see also Feddes et al., 1978)
    return if h < h4
        0.0
    elseif h < h3
        (h - h4) / (h3 - h4)
    elseif iszero(alpha_h1)
        if h < h2
            1.0
        elseif h < h1
            (h1 - h) / (h1 - h2)
        else
            0.0
        end
    else
        1.0
    end
end

"""
    soil_temperature(tsoil_prev, w_soil, temperature))

Return the near surface soil temperature `soil_surface_temperature` based on the near surface soil temperature
`tsoil_prev` at the previous timestep, and the difference between air `temperature` and near
surface soil temperature `tsoil_prev` at the previous timestep, weighted with the weighting
coefficient `w_soil` (Wigmosta et al., 2009).
"""
function soil_temperature(tsoil_prev, w_soil, temperature)
    soil_surface_temperature = tsoil_prev + w_soil * (temperature - tsoil_prev)
    return soil_surface_temperature
end

"""
    infiltration_reduction_factor(
        soil_surface_temperature,
        cf_soil;
        modelsnow = false,
        soil_infiltration_reduction = false,
    )

When both `modelsnow` and `soil_infiltration_reduction` are `true` an infiltration reduction
factor `f_infiltration_reduction` is computed. The infiltration reduction factor is based on
the near surface soil temperature `soil_surface_temperature`, parameter `cf_soil` and a s-curve to make a
smooth transition of `f_infiltration_reduction` as a function of `soil_surface_temperature` and `cf_soil`.
Otherwise, `f_infiltration_reduction` is set to 1.0.
"""
function infiltration_reduction_factor(
    soil_surface_temperature,
    cf_soil;
    modelsnow = false,
    soil_infiltration_reduction = false,
)
    if modelsnow && soil_infiltration_reduction
        bb = 1.0 / (1.0 - cf_soil)
        f_infiltration_reduction =
            scurve(soil_surface_temperature, to_SI(0.0, ABSOLUTE_DEGREES), bb, 8.0) +
            cf_soil
    else
        f_infiltration_reduction = 1.0
    end
    return f_infiltration_reduction
end

"Return soil evaporation from the unsaturated store"
function soil_evaporation_unsaturated_store(
    potential_soilevaporation,
    unsaturated_layer_depth,
    unsaturated_layer_thickness,
    n_unsatlayers,
    water_table_depth,
    theta_effective,
)
    if n_unsatlayers == 0
        soilevapunsat = 0.0
    elseif n_unsatlayers == 1
        # Check if groundwater level lies below the surface
        soilevapunsat =
            potential_soilevaporation *
            min(1.0, unsaturated_layer_depth / (water_table_depth * theta_effective))
    else
        # In case first layer contains no saturated storage
        soilevapunsat =
            potential_soilevaporation * min(
                1.0,
                unsaturated_layer_depth / (unsaturated_layer_thickness * theta_effective),
            )
    end
    return soilevapunsat
end

"Return soil evaporation from the saturated store"
function soil_evaporation_saturated_store(
    potential_soilevaporation,
    n_unsatlayers,
    layerthickness,
    water_table_depth,
    theta_drainable,
    dt,
)
    if n_unsatlayers in (0, 1)
        soil_evaporation_saturated_zone =
            potential_soilevaporation *
            min(1.0, (layerthickness - water_table_depth) / layerthickness)
        soil_evaporation_saturated_zone = min(
            soil_evaporation_saturated_zone,
            (layerthickness - water_table_depth) * theta_drainable / dt,
        )
    else
        soil_evaporation_saturated_zone = 0.0
    end
    return soil_evaporation_saturated_zone
end

"Return actual infiltration rate for soil `actual_infiltration_soil` and paved area `actual_infiltration_compacted_soil`"
function actual_infiltration_soil_path(
    potential_infiltration,
    actual_infiltration,
    compacted_soil_area_fraction,
    infiltration_capacity_soil,
    infiltration_capacity_compacted_soil,
    f_infiltration_reduction,
)
    soilinf = potential_infiltration * (1.0 - compacted_soil_area_fraction)
    pathinf = potential_infiltration * compacted_soil_area_fraction
    if actual_infiltration > 0.0
        max_infiltsoil = min(infiltration_capacity_soil * f_infiltration_reduction, soilinf)
        max_infiltpath =
            min(infiltration_capacity_compacted_soil * f_infiltration_reduction, pathinf)

        actual_infiltration_soil =
            actual_infiltration * max_infiltsoil / (max_infiltpath + max_infiltsoil)
        actual_infiltration_compacted_soil =
            actual_infiltration * max_infiltpath / (max_infiltpath + max_infiltsoil)

    else
        actual_infiltration_soil = 0.0
        actual_infiltration_compacted_soil = 0.0
    end

    return actual_infiltration_soil, actual_infiltration_compacted_soil
end
