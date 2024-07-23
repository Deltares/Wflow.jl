@get_units @with_kw struct SBM{IM, SM, GM, T}
    atmospheric_forcing::AtmosphericForcing | "-"
    veg_param_set::VegetationParameters | "-"
    interception_model::IM | "-"
    snow_model::SM | "-"
    glacier_model::GM | "-"
    soil_model::SoilSbmModel | "-"
    # Model time step [s]
    dt::T | "s"
end

function initialize_sbm(nc, config, riverfrac, inds)
    dt = Second(config.timestepsecs)
    n = length(inds)

    atmospheric_forcing = initialize_atmospheric_forcing(n)
    veg_param_set = initialize_vegetation_params(nc, config, inds)
    if dt >= Hour(23)
        interception_model =
            initialize_gash_interception_model(nc, config, inds, veg_param_set)
    else
        interception_model = initialize_rutter_interception_model(veg_param_set, n)
    end

    modelsnow = get(config.model, "snow", false)::Bool
    if modelsnow
        snow_model = initialize_snow_hbv_model(nc, config, inds, dt)
    else
        snow_model = NoSnowModel()
    end
    modelglacier = get(config.model, "glacier", false)::Bool
    if modelsnow && modelglacier
        glacier_bc = glacier_model_bc(snow_model.variables.snow)
        glacier_model = initialize_glacier_hbv_model(nc, config, inds, dt, glacier_bc)
    else
        glacier_model = NoGlacierModel()
    end
    soil_model = initialize_soil_sbm_model(nc, config, riverfrac, inds, dt)

    # TODO (part of refactor v1.0): simplify typeof arguments
    sbm =
        SBM{typeof(interception_model), typeof(snow_model), typeof(glacier_model), Float}(;
            atmospheric_forcing = atmospheric_forcing,
            veg_param_set = veg_param_set,
            interception_model = interception_model,
            snow_model = snow_model,
            glacier_model = glacier_model,
            soil_model = soil_model,
            dt = tosecond(dt),
        )
    return sbm
end

function update_until_snow(sbm::SBM, config)
    modelsnow = get(config.model, "snow", false)::Bool
    (; canopy_potevap, interception, throughfall, stemflow) =
        sbm.interception_model.variables
    (; pottrans) = sbm.soil_model.variables

    update(sbm.interception_model, sbm.atmospheric_forcing)
    @. pottrans = max(0.0, canopy_potevap - interception)

    if modelsnow
        (; effective_precip) = sbm.snow_model.boundary_conditions
        @. effective_precip = throughfall + stemflow
    end
    update(sbm.snow_model, sbm.atmospheric_forcing)

    update(sbm.glacier_model, sbm.atmospheric_forcing)
end

function update_until_recharge(sbm::SBM, config)

    # start dummy variables (should be generated from model reader and from Config.jl TOML)
    soilinfreduction = get(config.model, "soilinfreduction", false)::Bool
    modelglacier = get(config.model, "glacier", false)::Bool
    modelsnow = get(config.model, "snow", false)::Bool
    transfermethod = get(config.model, "transfermethod", false)::Bool
    ust = get(config.model, "whole_ust_available", false)::Bool # should be removed from optional setting and code?
    ksat_profile = get(config.input.vertical, "ksat_profile", "exponential")::String

    threaded_foreach(1:(sbm.n); basesize = 250) do i
        if modelsnow
            rainfallplusmelt = sbm.rainfallplusmelt[i]
            if modelglacier
                # Run Glacier module and add the snowpack on-top of it.
                # Estimate the fraction of snow turned into ice (HBV-light).
                # Estimate glacier melt.

                snow, _, glacierstore, glaciermelt = glacier_hbv(
                    sbm.glacierfrac[i],
                    sbm.glacierstore[i],
                    sbm.snow[i],
                    sbm.temperature[i],
                    sbm.g_tt[i],
                    sbm.g_cfmax[i],
                    sbm.g_sifrac[i],
                    Second(sbm.dt),
                )
                # Convert to mm per grid cell and add to snowmelt
                glaciermelt = glaciermelt * sbm.glacierfrac[i]
                rainfallplusmelt = rainfallplusmelt + glaciermelt
            end
        else
            rainfallplusmelt = sbm.stemflow[i] + sbm.throughfall[i]
        end

        avail_forinfilt = rainfallplusmelt
        ustoredepth = sum(@view sbm.ustorelayerdepth[i][1:sbm.nlayers[i]])

        runoff_river = min(1.0, sbm.riverfrac[i]) * avail_forinfilt
        runoff_land = min(1.0, sbm.waterfrac[i]) * avail_forinfilt
        avail_forinfilt = max(avail_forinfilt - runoff_river - runoff_land, 0.0)

        rootingdepth = min(sbm.soilthickness[i] * 0.99, sbm.rootingdepth[i])

        ae_openw_r = min(
            sbm.waterlevel_river[i] * sbm.riverfrac[i],
            sbm.riverfrac[i] * sbm.potential_evaporation[i],
        )
        ae_openw_l = min(
            sbm.waterlevel_land[i] * sbm.waterfrac[i],
            sbm.waterfrac[i] * sbm.potential_evaporation[i],
        )

        # evap available for soil evaporation
        soilevap_fraction = max(
            sbm.canopygapfraction[i] - sbm.riverfrac[i] - sbm.waterfrac[i] -
            sbm.glacierfrac[i],
            0.0,
        )
        potsoilevap = soilevap_fraction * sbm.potential_evaporation[i]

        # Calculate the initial capacity of the unsaturated store
        ustorecapacity = sbm.soilwatercapacity[i] - sbm.satwaterdepth[i] - ustoredepth

        # Calculate the infiltration flux into the soil column
        infiltsoilpath, infiltsoil, infiltpath, soilinf, pathinf, infiltexcess =
            infiltration(
                avail_forinfilt,
                sbm.pathfrac[i],
                sbm.cf_soil[i],
                sbm.tsoil[i],
                sbm.infiltcapsoil[i],
                sbm.infiltcappath[i],
                ustorecapacity,
                modelsnow,
                soilinfreduction,
            )

        usl, n_usl = set_layerthickness(sbm.zi[i], sbm.sumlayers[i], sbm.act_thickl[i])
        z = cumsum(usl)
        usld = sbm.ustorelayerdepth[i]

        ast = 0.0
        soilevapunsat = 0.0
        if n_usl > 0
            # Using the surface infiltration rate, calculate the flow rate between the
            # different soil layers that contain unsaturated storage assuming gravity
            # based flow only, estimate the gravity based flux rate to the saturated zone
            # (ast) and the updated unsaturated storage for each soil layer.
            if transfermethod && sbm.maxlayers == 1
                ustorelayerdepth = sbm.ustorelayerdepth[i][1] + infiltsoilpath
                kv_z = hydraulic_conductivity_at_depth(sbm, sbm.zi[i], i, 1, ksat_profile)
                ustorelayerdepth, ast = unsatzone_flow_sbm(
                    ustorelayerdepth,
                    sbm.soilwatercapacity[i],
                    sbm.satwaterdepth[i],
                    kv_z,
                    usl[1],
                    sbm.theta_s[i],
                    sbm.theta_r[i],
                )
                usld = setindex(usld, ustorelayerdepth, 1)
            else
                for m in 1:n_usl
                    l_sat = usl[m] * (sbm.theta_s[i] - sbm.theta_r[i])
                    kv_z = hydraulic_conductivity_at_depth(sbm, z[m], i, m, ksat_profile)
                    ustorelayerdepth = if m == 1
                        sbm.ustorelayerdepth[i][m] + infiltsoilpath
                    else
                        sbm.ustorelayerdepth[i][m] + ast
                    end
                    ustorelayerdepth, ast =
                        unsatzone_flow_layer(ustorelayerdepth, kv_z, l_sat, sbm.c[i][m])
                    usld = setindex(usld, ustorelayerdepth, m)
                end
            end

            # then evapotranspiration from layers
            # Calculate saturation deficit
            saturationdeficit = sbm.soilwatercapacity[i] - sbm.satwaterdepth[i]

            # First calculate the evaporation of unsaturated storage into the
            # atmosphere from the upper layer.
            if sbm.maxlayers == 1
                soilevapunsat =
                    potsoilevap * min(1.0, saturationdeficit / sbm.soilwatercapacity[i])
            else
                # In case only the most upper soil layer contains unsaturated storage
                if n_usl == 1
                    # Check if groundwater level lies below the surface
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (sbm.zi[i] * (sbm.theta_s[i] - sbm.theta_r[i])))
                else
                    # In case first layer contains no saturated storage
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (usl[1] * ((sbm.theta_s[i] - sbm.theta_r[i]))))
                end
            end
            # Ensure that the unsaturated evaporation rate does not exceed the
            # available unsaturated moisture
            soilevapunsat = min(soilevapunsat, usld[1])
            # Update the additional atmospheric demand
            potsoilevap = potsoilevap - soilevapunsat
            usld = setindex(usld, usld[1] - soilevapunsat, 1)
        end
        transfer = ast

        if sbm.maxlayers == 1
            soilevapsat = 0.0
        else
            if n_usl == 0 || n_usl == 1
                soilevapsat =
                    potsoilevap *
                    min(1.0, (sbm.act_thickl[i][1] - sbm.zi[i]) / sbm.act_thickl[i][1])
                soilevapsat = min(
                    soilevapsat,
                    (sbm.act_thickl[i][1] - sbm.zi[i]) * (sbm.theta_s[i] - sbm.theta_r[i]),
                )
            else
                soilevapsat = 0.0
            end
        end
        soilevap = soilevapunsat + soilevapsat
        satwaterdepth = sbm.satwaterdepth[i] - soilevapsat

        # transpiration from saturated store
        wetroots = scurve(sbm.zi[i], rootingdepth, Float(1.0), sbm.rootdistpar[i])
        actevapsat = min(sbm.pottrans[i] * wetroots, satwaterdepth)
        satwaterdepth = satwaterdepth - actevapsat
        restpottrans = sbm.pottrans[i] - actevapsat

        # actual transpiration from ustore
        actevapustore = 0.0
        for k in 1:n_usl
            ustorelayerdepth, actevapustore, restpottrans = acttransp_unsat_sbm(
                rootingdepth,
                usld[k],
                sbm.sumlayers[i][k],
                restpottrans,
                actevapustore,
                sbm.c[i][k],
                usl[k],
                sbm.theta_s[i],
                sbm.theta_r[i],
                sbm.hb[i],
                ust,
            )
            usld = setindex(usld, ustorelayerdepth, k)
        end

        # check soil moisture balance per layer
        du = 0.0
        for k in n_usl:-1:1
            du = max(0.0, usld[k] - usl[k] * (sbm.theta_s[i] - sbm.theta_r[i]))
            usld = setindex(usld, usld[k] - du, k)
            if k > 1
                usld = setindex(usld, usld[k - 1] + du, k - 1)
            end
        end

        actinfilt = infiltsoilpath - du
        excesswater = avail_forinfilt - infiltsoilpath - infiltexcess + du

        # Separation between compacted and non compacted areas (correction with the satflow du)
        # This is required for D-Emission/Delwaq
        if infiltsoil + infiltpath > 0.0
            actinfiltsoil = infiltsoil - du * infiltsoil / (infiltpath + infiltsoil)
            actinfiltpath = infiltpath - du * infiltpath / (infiltpath + infiltsoil)
        else
            actinfiltsoil = 0.0
            actinfiltpath = 0.0
        end
        excesswatersoil = max(soilinf - actinfiltsoil, 0.0)
        excesswaterpath = max(pathinf - actinfiltpath, 0.0)

        actcapflux = 0.0
        if n_usl > 0
            ksat = hydraulic_conductivity_at_depth(sbm, sbm.zi[i], i, n_usl, ksat_profile)
            ustorecapacity =
                sbm.soilwatercapacity[i] - satwaterdepth - sum(@view usld[1:sbm.nlayers[i]])
            maxcapflux = max(0.0, min(ksat, actevapustore, ustorecapacity, satwaterdepth))

            if sbm.zi[i] > rootingdepth
                capflux =
                    maxcapflux * pow(
                        1.0 - min(sbm.zi[i], sbm.cap_hmax[i]) / (sbm.cap_hmax[i]),
                        sbm.cap_n[i],
                    )
            else
                capflux = 0.0
            end

            netcapflux = capflux
            for k in n_usl:-1:1
                toadd = min(
                    netcapflux,
                    max(usl[k] * (sbm.theta_s[i] - sbm.theta_r[i]) - usld[k], 0.0),
                )
                usld = setindex(usld, usld[k] + toadd, k)
                netcapflux = netcapflux - toadd
                actcapflux = actcapflux + toadd
            end
        end
        deepksat = hydraulic_conductivity_at_depth(
            sbm,
            sbm.soilthickness[i],
            i,
            sbm.nlayers[i],
            ksat_profile,
        )
        deeptransfer = min(satwaterdepth, deepksat)
        actleakage = max(0.0, min(sbm.maxleakage[i], deeptransfer))

        # recharge (mm) for saturated zone
        recharge = (transfer - actcapflux - actleakage - actevapsat - soilevapsat)
        transpiration = actevapsat + actevapustore
        actevap = soilevap + transpiration + ae_openw_r + ae_openw_l + sbm.interception[i]

        # update the outputs and states
        sbm.n_unsatlayers[i] = n_usl
        sbm.net_runoff_river[i] = runoff_river - ae_openw_r
        sbm.avail_forinfilt[i] = avail_forinfilt
        sbm.actinfilt[i] = actinfilt
        sbm.infiltexcess[i] = infiltexcess
        sbm.recharge[i] = recharge
        sbm.transpiration[i] = transpiration
        sbm.soilevap[i] = soilevap
        sbm.soilevapsat[i] = soilevapsat
        sbm.ae_openw_r[i] = ae_openw_r
        sbm.ae_openw_l[i] = ae_openw_l
        sbm.runoff_land[i] = runoff_land
        sbm.runoff_river[i] = runoff_river
        sbm.actevapsat[i] = actevapsat
        sbm.actevap[i] = actevap
        sbm.ae_ustore[i] = actevapustore
        sbm.ustorelayerdepth[i] = usld
        sbm.transfer[i] = transfer
        sbm.actcapflux[i] = actcapflux
        sbm.actleakage[i] = actleakage
        sbm.actinfiltsoil[i] = actinfiltsoil
        sbm.actinfiltpath[i] = actinfiltpath
        sbm.excesswater[i] = excesswater
        sbm.excesswatersoil[i] = excesswatersoil
        sbm.excesswaterpath[i] = excesswaterpath
        sbm.rainfallplusmelt[i] = rainfallplusmelt
        sbm.infiltsoilpath[i] = infiltsoilpath
        sbm.satwaterdepth[i] = satwaterdepth
        if modelsnow
            if modelglacier
                sbm.snow[i] = snow
                sbm.glacierstore[i] = glacierstore
            end
        end
    end
end

function update_after_subsurfaceflow(sbm::SBM, zi, exfiltsatwater)
    threaded_foreach(1:(sbm.n); basesize = 1000) do i
        usl, n_usl = set_layerthickness(zi[i], sbm.sumlayers[i], sbm.act_thickl[i])
        # exfiltration from ustore
        usld = sbm.ustorelayerdepth[i]
        exfiltustore = 0.0
        for k in sbm.n_unsatlayers[i]:-1:1
            if k <= n_usl
                exfiltustore = max(0, usld[k] - usl[k] * (sbm.theta_s[i] - sbm.theta_r[i]))
            else
                exfiltustore = usld[k]
            end
            usld = setindex(usld, usld[k] - exfiltustore, k)
            if k > 1
                usld = setindex(usld, usld[k - 1] + exfiltustore, k - 1)
            end
        end

        ustoredepth = sum(@view usld[1:n_usl])

        runoff =
            exfiltustore +
            exfiltsatwater[i] +
            sbm.excesswater[i] +
            sbm.runoff_land[i] +
            sbm.infiltexcess[i]

        # volumetric water content per soil layer and root zone
        vwc = sbm.vwc[i]
        vwc_perc = sbm.vwc_perc[i]
        for k in 1:sbm.nlayers[i]
            if k <= n_usl
                vwc = setindex(
                    vwc,
                    (
                        usld[k] +
                        (sbm.act_thickl[i][k] - usl[k]) * (sbm.theta_s[i] - sbm.theta_r[i])
                    ) / sbm.act_thickl[i][k] + sbm.theta_r[i],
                    k,
                )
            else
                vwc = setindex(vwc, sbm.theta_s[i], k)
            end
            vwc_perc = setindex(vwc_perc, (vwc[k] / sbm.theta_s[i]) * 100.0, k)
        end

        rootstore_unsat = 0
        for k in 1:n_usl
            rootstore_unsat =
                rootstore_unsat +
                min(1.0, (max(0.0, sbm.rootingdepth[i] - sbm.sumlayers[i][k]) / usl[k])) *
                usld[k]
        end

        rootstore_sat =
            max(0.0, sbm.rootingdepth[i] - zi[i]) * (sbm.theta_s[i] - sbm.theta_r[i])
        rootstore = rootstore_sat + rootstore_unsat
        vwc_root = rootstore / sbm.rootingdepth[i] + sbm.theta_r[i]
        vwc_percroot = (vwc_root / sbm.theta_s[i]) * 100.0

        satwaterdepth = (sbm.soilthickness[i] - zi[i]) * (sbm.theta_s[i] - sbm.theta_r[i])

        # update the outputs and states
        sbm.n_unsatlayers[i] = n_usl
        sbm.ustorelayerdepth[i] = usld
        sbm.ustoredepth[i] = ustoredepth
        sbm.satwaterdepth[i] = satwaterdepth
        sbm.exfiltsatwater[i] = exfiltsatwater[i]
        sbm.exfiltustore[i] = exfiltustore
        sbm.runoff[i] = runoff
        sbm.net_runoff[i] = runoff - sbm.ae_openw_l[i]
        sbm.vwc[i] = vwc
        sbm.vwc_perc[i] = vwc_perc
        sbm.rootstore[i] = rootstore
        sbm.vwc_root[i] = vwc_root
        sbm.vwc_percroot[i] = vwc_percroot
        sbm.zi[i] = zi[i]
    end
end

"""
Update the total water storage per cell at the end of a timestep.

Takes the following parameters:
- sbm:
    The vertical concept (SBM struct)
- river_network:
    The indices of the river cells in relation to the active cells, i.e. model.network.index_river
- cell_xsize:
    Size in X direction of the cells acquired from model.network.land.xl
- cell_ysize:
    Size in Y direction of the cells acquired from model.network.land.yl
- river_routing:
    The river routing struct, i.e. model.lateral.river
- land_routing:
    The land routing struct, i.e. model.lateral.land
"""
function update_total_water_storage(
    sbm::SBM,
    river_network,
    cell_xsize,
    cell_ysize,
    river_routing,
    land_routing,
)
    # Get length active river cells
    nriv = length(river_network)

    # Set the total storage to zero
    fill!(sbm.total_storage, 0)

    # Burn the river routing values
    sbm.total_storage[river_network] = (
        (
            river_routing.h_av[1:nriv] .* river_routing.width[1:nriv] .*
            river_routing.dl[1:nriv]
        ) ./ (cell_xsize[river_network] .* cell_ysize[river_network]) * 1000 # Convert to mm
    )

    # Chunk the data for parallel computing
    threaded_foreach(1:(sbm.n); basesize = 1000) do i

        # Cumulate per vertical type
        # Maybe re-categorize in the future
        surface = (
            sbm.glacierstore[i] * sbm.glacierfrac[i] +
            sbm.snow[i] +
            sbm.snowwater[i] +
            sbm.canopystorage[i]
        )
        sub_surface = sbm.ustoredepth[i] + sbm.satwaterdepth[i]
        lateral = (
            land_routing.h_av[i] * (1 - sbm.riverfrac[i]) * 1000 # convert to mm
        )

        # Add everything to the total water storage
        sbm.total_storage[i] += (surface + sub_surface + lateral)
    end
end
