# ---- Scalar soil-physics functions ----

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
    feddes_h3(h3_high, h3_low, tpot_SI)

Return soil water pressure head `h3` of Feddes root water uptake reduction function.
`tpot_SI` is potential transpiration in m s⁻¹ (SI units).
"""
function feddes_h3(h3_high, h3_low, tpot_SI)
    # value of h3 is a function of potential transpiration [mm d⁻¹]
    tpot_daily = tpot_SI * 1000.0 * 86400.0
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
        # 273.15 K = 0.0 °C (freezing point in SI)
        f_infiltration_reduction =
            scurve(soil_surface_temperature, 273.15, bb, 8.0) + cf_soil
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

# ---- Per-cell mutating update functions ----

"Update total available water in the unsaturated zone of the SBM soil model for a single timestep"
function unsaturated_store_depth!(soil_model::SbmSoilModel)
    s = soil_model.variables.states
    d = soil_model.variables.diagnostic
    p = soil_model.parameters
    for i in eachindex(s.unsaturated_layer_depth)
        d.unsaturated_store_depth[i] =
            sum(@view s.unsaturated_layer_depth[i][1:p.number_of_layers[i]])
    end
    return nothing
end

"Update the infiltration reduction factor of the SBM soil model for a single timestep"
function infiltration_reduction_factor!(
    soil_model::SbmSoilModel;
    modelsnow = false,
    soil_infiltration_reduction = false,
)
    s = soil_model.variables.states
    im = soil_model.variables.intermediates
    p = soil_model.parameters

    n = length(s.soil_surface_temperature)
    threaded_foreach(1:n; basesize = 1000) do i
        im.f_infiltration_reduction[i] = infiltration_reduction_factor(
            s.soil_surface_temperature[i],
            p.cf_soil[i];
            modelsnow,
            soil_infiltration_reduction,
        )
    end
    return nothing
end

"""
    infiltration!(soil_model::SbmSoilModel, dt::Float64)

Update the infiltration rate `infiltration` and infiltration excess water rate
`infiltration_excess` of the SBM soil model for a single timestep.
"""
function infiltration!(soil_model::SbmSoilModel, dt::Float64)
    d = soil_model.variables.diagnostic
    im = soil_model.variables.intermediates
    f = soil_model.variables.fluxes
    p = soil_model.parameters
    (; water_flux_surface) = soil_model.boundary_conditions

    n = length(f.infiltration)
    threaded_foreach(1:n; basesize = 1000) do i
        f.infiltration[i], f.infiltration_excess[i] = infiltration(
            water_flux_surface[i],
            p.compacted_soil_area_fraction[i],
            p.infiltration_capacity_soil[i],
            p.infiltration_capacity_compacted_soil[i],
            d.unsaturated_store_capacity[i],
            im.f_infiltration_reduction[i],
            dt,
        )
    end
    return nothing
end

"""
    unsaturated_zone_flow!(soil_model::SbmSoilModel, dt::Float64)

Update unsaturated storage `unsaturated_layer_depth` and the `transfer` of water from the unsaturated
to the saturated store of the SBM soil model for a single timestep, based on the Brooks-Corey
approach.
"""
function unsaturated_zone_flow!(soil_model::SbmSoilModel, dt::Float64)
    s = soil_model.variables.states
    d = soil_model.variables.diagnostic
    f = soil_model.variables.fluxes
    p = soil_model.parameters

    n = length(f.transfer)
    threaded_foreach(1:n; basesize = 250) do i
        if d.n_unsatlayers[i] > 0
            # Brooks-Corey approach
            z = cumsum(d.unsaturated_layer_thickness[i])
            flow_rate = 0.0
            for m in 1:d.n_unsatlayers[i]
                l_sat = d.unsaturated_layer_thickness[i][m] * (p.theta_s[i] - p.theta_r[i])
                kv_z = hydraulic_conductivity_at_depth(
                    p.kv_profile,
                    p.vertical_hydraulic_conductivity_factor,
                    z[m],
                    i,
                    m,
                )
                unsaturated_layer_depth = if m == 1
                    s.unsaturated_layer_depth[i][m] + f.infiltration[i] * dt
                else
                    s.unsaturated_layer_depth[i][m] + flow_rate * dt
                end
                unsaturated_layer_depth, flow_rate = unsatzone_flow_layer(
                    unsaturated_layer_depth,
                    kv_z,
                    l_sat,
                    p.brooks_corey_exponent[i][m],
                    dt,
                )
                s.unsaturated_layer_depth[i] =
                    setindex(s.unsaturated_layer_depth[i], unsaturated_layer_depth, m)
            end
            f.transfer[i] = flow_rate
        else
            f.transfer[i] = 0.0
        end
    end
    return nothing
end

"""
    soil_evaporation!(soil_model::SbmSoilModel, dt::Float64)

Update soil evaporation from the saturated store `soil_evaporation_saturated_zone` and the total soil
evaporation from the unsaturated and saturated store `soil_evaporation` of the SBM soil model for a
single timestep. Also unsaturated storage `unsaturated_layer_depth` and the saturated store
`saturated_water_depth` are updated.
"""
function soil_evaporation!(soil_model::SbmSoilModel, dt::Float64)
    (; potential_soilevaporation) = soil_model.boundary_conditions
    s = soil_model.variables.states
    d = soil_model.variables.diagnostic
    f = soil_model.variables.fluxes
    p = soil_model.parameters

    n = length(potential_soilevaporation)
    threaded_foreach(1:n; basesize = 1000) do i
        potsoilevap = potential_soilevaporation[i]
        # First calculate the evaporation of unsaturated storage into the
        # atmosphere from the upper layer.
        soilevapunsat = soil_evaporation_unsaturated_store(
            potsoilevap,
            s.unsaturated_layer_depth[i][1],
            d.unsaturated_layer_thickness[i][1],
            d.n_unsatlayers[i],
            d.water_table_depth[i],
            p.theta_s[i] - p.theta_r[i],
        )
        # Ensure that the unsaturated evaporation rate does not exceed the
        # available unsaturated moisture
        soilevapunsat = min(soilevapunsat, s.unsaturated_layer_depth[i][1] / dt)
        # Update the additional atmospheric demand
        potsoilevap -= soilevapunsat
        s.unsaturated_layer_depth[i] = setindex(
            s.unsaturated_layer_depth[i],
            s.unsaturated_layer_depth[i][1] - soilevapunsat * dt,
            1,
        )
        theta_drainable = lower_bound_drainable_porosity(p.theta_s[i], p.theta_fc[i])
        soil_evaporation_saturated_zone = soil_evaporation_saturated_store(
            potsoilevap,
            d.n_unsatlayers[i],
            p.actual_layer_thickness[i][1],
            d.water_table_depth[i],
            theta_drainable,
            dt,
        )
        f.soil_evaporation_saturated_zone[i] = soil_evaporation_saturated_zone
        f.soil_evaporation[i] = soilevapunsat + soil_evaporation_saturated_zone
        d.drainable_water_depth[i] -= soil_evaporation_saturated_zone * dt
    end
    return nothing
end

"""
    transpiration!(soil_model::SbmSoilModel, dt::Float64)

Update total `transpiration`, transpiration from the unsaturated store `actual_evaporation_unsaturated_store` and
saturated store `actual_evaporation_saturated_zone` of the SBM soil model for a single timestep. Also unsaturated
storage `unsaturated_layer_depth` and the saturated store `saturated_water_depth` are updated.
"""
function transpiration!(soil_model::SbmSoilModel, dt::Float64)
    (; potential_transpiration) = soil_model.boundary_conditions
    s = soil_model.variables.states
    d = soil_model.variables.diagnostic
    im = soil_model.variables.intermediates
    f = soil_model.variables.fluxes
    p = soil_model.parameters

    rooting_depth = get_rootingdepth(soil_model)
    n = length(rooting_depth)

    threaded_foreach(1:n; basesize = 250) do i
        im.h3[i] = feddes_h3(p.h3_high[i], p.h3_low[i], potential_transpiration[i])

        # compute sum of root fraction in unsaturated soil layers and adapt root fraction
        # lowest unsaturated soil layer if water table depth intersects the unsaturated root
        # zone
        sum_rootfraction_unsat = 0.0
        rootfraction_unsat_lowest = 0.0
        for k in 1:d.n_unsatlayers[i]
            # the root fraction is valid for the root length in a soil layer, if zi decreases
            # the root length the root fraction needs to be adapted
            if k == d.n_unsatlayers[i] && d.water_table_depth[i] < rooting_depth[i]
                rootlength = min(
                    p.actual_layer_thickness[i][k],
                    rooting_depth[i] - p.cumulative_layer_depth[i][k],
                )
                rootfraction_unsat =
                    p.rootfraction[i][k] *
                    (d.unsaturated_layer_thickness[i][k] / rootlength)
            else
                rootfraction_unsat = p.rootfraction[i][k]
            end
            sum_rootfraction_unsat += rootfraction_unsat

            # rootfraction lowest unsaturated layer
            rootfraction_unsat_lowest = rootfraction_unsat
        end

        actevapustore = 0.0
        for k in 1:d.n_unsatlayers[i]
            # scale rootfraction soil layer unsaturated zone based on sum of rootfraction in
            # unsaturated zone
            if k < d.n_unsatlayers[i]
                rootfraction_unsat = p.rootfraction[i][k]
            else
                rootfraction_unsat = rootfraction_unsat_lowest
            end
            rootfraction_unsat_scaled =
                rooting_depth[i] > 0.0 ?
                max((1.0 / sum_rootfraction_unsat), 1.0) * rootfraction_unsat : 0.0
            volumetric_water_content = max(
                s.unsaturated_layer_depth[i][k] / d.unsaturated_layer_thickness[i][k],
                1e-7,
            )
            head = head_brooks_corey(
                volumetric_water_content,
                p.theta_s[i],
                p.theta_r[i],
                p.brooks_corey_exponent[i][k],
                p.air_entry_pressure[i],
            )
            alpha = rwu_reduction_feddes(
                head,
                p.h1[i],
                p.h2[i],
                im.h3[i],
                p.h4[i],
                p.alpha_h1[i],
            )
            availcap = min(
                1.0,
                max(
                    0.0,
                    (rooting_depth[i] - p.cumulative_layer_depth[i][k]) /
                    d.unsaturated_layer_thickness[i][k],
                ),
            )
            maxextr = s.unsaturated_layer_depth[i][k] * availcap / dt
            actevapustore_layer =
                min(alpha * rootfraction_unsat_scaled * potential_transpiration[i], maxextr)
            unsaturated_layer_depth =
                s.unsaturated_layer_depth[i][k] - actevapustore_layer * dt
            actevapustore += actevapustore_layer
            s.unsaturated_layer_depth[i] =
                setindex(s.unsaturated_layer_depth[i], unsaturated_layer_depth, k)
        end

        # transpiration from saturated store
        wetroots = scurve(
            d.water_table_depth[i],
            rooting_depth[i],
            1.0,
            p.wet_root_distribution_parameter[i],
        )
        alpha = rwu_reduction_feddes(
            Float64(0.0),
            p.h1[i],
            p.h2[i],
            im.h3[i],
            p.h4[i],
            p.alpha_h1[i],
        )
        restpottrans = potential_transpiration[i] - actevapustore
        actual_evaporation_saturated_zone =
            min(restpottrans * wetroots * alpha, d.drainable_water_depth[i] / dt)

        f.actual_evaporation_unsaturated_store[i] = actevapustore
        f.actual_evaporation_saturated_zone[i] = actual_evaporation_saturated_zone
        d.drainable_water_depth[i] -= actual_evaporation_saturated_zone * dt
        f.transpiration[i] = actevapustore + actual_evaporation_saturated_zone
    end
    return nothing
end

"""
    actual_infiltration!(soil_model::SbmSoilModel, dt::Float64)

Update the actual infiltration rate `actual_infiltration` of the SBM soil model for a single timestep.

A soil water balance check is performed. Unsaturated storage that exceeds the maximum
storage per unsaturated soil layer is transferred to the layer above (or surface), from the
bottom to the top unsaturated soil layer. The resulting excess water `ustoredepth_excess` is
subtracted from the infiltration rate `infiltration`.
"""
function actual_infiltration!(soil_model::SbmSoilModel, dt::Float64)
    s = soil_model.variables.states
    d = soil_model.variables.diagnostic
    f = soil_model.variables.fluxes
    p = soil_model.parameters

    n = length(f.actual_infiltration)
    threaded_foreach(1:n; basesize = 1000) do i
        # check soil moisture balance per layer
        ustoredepth_excess = 0.0
        for k in d.n_unsatlayers[i]:-1:1
            ustoredepth_excess = max(
                0.0,
                s.unsaturated_layer_depth[i][k] -
                d.unsaturated_layer_thickness[i][k] * (p.theta_s[i] - p.theta_r[i]),
            )
            s.unsaturated_layer_depth[i] = setindex(
                s.unsaturated_layer_depth[i],
                s.unsaturated_layer_depth[i][k] - ustoredepth_excess,
                k,
            )
            if k > 1
                s.unsaturated_layer_depth[i] = setindex(
                    s.unsaturated_layer_depth[i],
                    s.unsaturated_layer_depth[i][k - 1] + ustoredepth_excess,
                    k - 1,
                )
            end
        end
        f.actual_infiltration[i] = f.infiltration[i] - ustoredepth_excess / dt
    end
    return nothing
end

"""
    actual_infiltration_soil_path!(soil_model::SbmSoilModel)

Update the actual infiltration rate for soil `actual_infiltration_soil` and paved area
`actual_infiltration_compacted_soil` of the SBM soil model for a single timestep.
"""
function actual_infiltration_soil_path!(soil_model::SbmSoilModel)
    im = soil_model.variables.intermediates
    f = soil_model.variables.fluxes
    p = soil_model.parameters
    (; water_flux_surface) = soil_model.boundary_conditions

    n = length(water_flux_surface)
    threaded_foreach(1:n; basesize = 1000) do i
        f.actual_infiltration_soil[i], f.actual_infiltration_compacted_soil[i] =
            actual_infiltration_soil_path(
                water_flux_surface[i],
                f.actual_infiltration[i],
                p.compacted_soil_area_fraction[i],
                p.infiltration_capacity_soil[i],
                p.infiltration_capacity_compacted_soil[i],
                im.f_infiltration_reduction[i],
            )
    end
    return nothing
end

"""
    capillary_flux!(soil_model::SbmSoilModel, dt::Float64)

Update the capillary flux `actual_capillary_flux` of the SBM soil model for a single timestep.
"""
function capillary_flux!(soil_model::SbmSoilModel, dt::Float64)
    s = soil_model.variables.states
    d = soil_model.variables.diagnostic
    f = soil_model.variables.fluxes
    p = soil_model.parameters
    rooting_depth = get_rootingdepth(soil_model)

    n = length(rooting_depth)
    threaded_foreach(1:n; basesize = 1000) do i
        if d.n_unsatlayers[i] > 0
            ksat = hydraulic_conductivity_at_depth(
                p.kv_profile,
                p.vertical_hydraulic_conductivity_factor,
                d.water_table_depth[i],
                i,
                d.n_unsatlayers[i],
            )
            maxcapflux = max(
                0.0,
                min(
                    ksat,
                    f.actual_evaporation_unsaturated_store[i],
                    d.unsaturated_store_capacity[i] / dt,
                    d.drainable_water_depth[i] / dt,
                ),
            )

            capflux = if d.water_table_depth[i] > rooting_depth[i]
                maxcapflux * pow(
                    1.0 - min(d.water_table_depth[i], p.cap_hmax[i]) / (p.cap_hmax[i]),
                    p.cap_n[i],
                )
            else
                0.0
            end
            netcapflux = capflux
            actual_capillary_flux = 0.0
            for k in d.n_unsatlayers[i]:-1:1
                toadd = min(
                    netcapflux,
                    max(
                        (
                            d.unsaturated_layer_thickness[i][k] *
                            (p.theta_s[i] - p.theta_r[i]) -
                            s.unsaturated_layer_depth[i][k]
                        ) / dt,
                        0.0,
                    ),
                )
                s.unsaturated_layer_depth[i] = setindex(
                    s.unsaturated_layer_depth[i],
                    s.unsaturated_layer_depth[i][k] + toadd * dt,
                    k,
                )
                netcapflux -= toadd
                actual_capillary_flux += toadd
            end
            f.actual_capillary_flux[i] = actual_capillary_flux
        else
            f.actual_capillary_flux[i] = 0.0
        end
    end
    return nothing
end

"""
    leakage!(soil_model::SbmSoilModel, dt::Float64)

Update the actual leakage rate `actual_leakage` of the SBM soil model for a single timestep.
"""
function leakage!(soil_model::SbmSoilModel, dt::Float64)
    d = soil_model.variables.diagnostic
    f = soil_model.variables.fluxes
    p = soil_model.parameters

    n = length(f.actual_leakage)
    threaded_foreach(1:n; basesize = 1000) do i
        deepksat = hydraulic_conductivity_at_depth(
            p.kv_profile,
            p.vertical_hydraulic_conductivity_factor,
            p.soil_thickness[i],
            i,
            p.number_of_layers[i],
        )

        deeptransfer = min(d.drainable_water_depth[i] / dt, deepksat)
        f.actual_leakage[i] = max(0.0, min(p.maximum_leakage[i], deeptransfer))
    end
    return nothing
end

"""
    update_diagnostic_vars!(soil_model::SbmSoilModel)

Update diagnostic variables of `SbmSoilModel` that are critical for subsequent computations
and depend on state variables `saturated_water_depth` and `unsaturated_layer_depth`.
"""
function update_diagnostic_vars!(soil_model::SbmSoilModel)
    (; saturated_water_depth) = soil_model.variables.states
    (;
        water_table_depth,
        drainable_water_depth,
        unsaturated_layer_thickness,
        unsaturated_store_capacity,
        unsaturated_store_depth,
        total_soil_water_storage,
        n_unsatlayers,
    ) = soil_model.variables.diagnostic
    (;
        soil_thickness,
        theta_s,
        theta_r,
        theta_fc,
        soil_water_capacity,
        cumulative_layer_depth,
        actual_layer_thickness,
    ) = soil_model.parameters

    unsaturated_store_depth!(soil_model)
    @. water_table_depth =
        max(0.0, soil_thickness - saturated_water_depth / (theta_s - theta_r))
    @. drainable_water_depth =
        (soil_thickness - water_table_depth) *
        lower_bound_drainable_porosity(theta_s, theta_fc)
    @. unsaturated_store_capacity =
        soil_water_capacity - saturated_water_depth - unsaturated_store_depth
    @. unsaturated_layer_thickness = set_layerthickness(
        water_table_depth,
        cumulative_layer_depth,
        actual_layer_thickness,
    )
    @. n_unsatlayers = number_of_active_layers(unsaturated_layer_thickness)
    @. total_soil_water_storage = saturated_water_depth + unsaturated_store_depth
end

# wrapper method
get_rootingdepth(soil_model::SbmSoilModel) =
    soil_model.parameters.vegetation_parameter_set.rooting_depth

"""
    update_soil_water_flow!(
        soil_model::SbmSoilModel,
        atmospheric_forcing,
        dt::Float64;
        snow = nothing,
        snow__flag::Bool = false,
        soil_infiltration_reduction__flag::Bool = false,
    )

Update the SBM soil model (infiltration, unsaturated zone flow, soil evaporation and
transpiration, capillary flux and leakage) for a single timestep.
"""
function update_soil_water_flow!(
    soil_model::SbmSoilModel,
    dt::Float64;
    snow__flag::Bool = false,
    soil_infiltration_reduction__flag::Bool = false,
)
    (; water_flux_surface) = soil_model.boundary_conditions
    s = soil_model.variables.states
    d = soil_model.variables.diagnostic
    f = soil_model.variables.fluxes
    p = soil_model.parameters

    # mainly required for external state changes (e.g. through BMI)
    update_diagnostic_vars!(soil_model)
    # infiltration
    infiltration_reduction_factor!(
        soil_model;
        modelsnow = snow__flag,
        soil_infiltration_reduction = soil_infiltration_reduction__flag,
    )
    infiltration!(soil_model, dt)
    # unsaturated zone flow
    unsaturated_zone_flow!(soil_model, dt)
    # soil evaporation and transpiration
    soil_evaporation!(soil_model, dt)
    transpiration!(soil_model, dt)
    # actual infiltration and excess water
    actual_infiltration!(soil_model, dt)
    @. f.saturation_excess_water =
        (water_flux_surface - f.actual_infiltration) - f.infiltration_excess
    actual_infiltration_soil_path!(soil_model)
    @. f.excess_water_soil = max(
        water_flux_surface * (1.0 - p.compacted_soil_area_fraction) -
        f.actual_infiltration_soil,
        0.0,
    )
    @. f.excess_water_compacted_soil = max(
        water_flux_surface * p.compacted_soil_area_fraction -
        f.actual_infiltration_compacted_soil,
        0.0,
    )
    # recompute the unsaturated store and unsaturated_store_capacity (for capillary flux)
    unsaturated_store_depth!(soil_model)
    @. d.unsaturated_store_capacity =
        p.soil_water_capacity - s.saturated_water_depth - d.unsaturated_store_depth
    # capillary flux and leakage
    capillary_flux!(soil_model, dt)
    leakage!(soil_model, dt)
    # recharge rate to the saturated store
    @. f.recharge = (
        f.transfer - f.actual_capillary_flux - f.actual_leakage -
        f.actual_evaporation_saturated_zone - f.soil_evaporation_saturated_zone
    )
    # total actual evapotranspiration (soil + transpiration only; open water added by caller)
    f.actual_evapotranspiration .= f.soil_evaporation .+ f.transpiration
    return nothing
end
