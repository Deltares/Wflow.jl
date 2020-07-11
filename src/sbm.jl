const mv = NaN

function update_until_snow(t, config)

    # # start dummy variables (should be generated from model reader and from Config.jl TOML)
    do_lai = !isequal(get(config.cyclic_parameters, "lai", nothing), nothing) ? true : false
    modelglacier = Bool(get(config.model, "modelglacier", 0))
    modelsnow = Bool(get(config.model, "modelsnow", 0))
    #potevap = 4.0
    #precipitation = 3.0
    #temperature = 10.0
    Δt = Second(config.input.timestepsecs)
    #basetimestep = Second(Day(1))
    # end dummpy variables

    for r in eachindex(t.nlayers)
        if do_lai
            cmax = t.sl[r] * t.lai[r] + t.swood[r]
            canopygapfraction = exp(-t.kext[r] * t.lai[r])
            ewet = (1.0 - exp(-t.kext[r] * t.lai[r])) * t.potevap[r]
            e_r = t.precipitation[r] > 0.0 ?
                min(0.25, ewet / max(0.0001, t.precipitation[r])) : 0.0
        end

        potevap = t.potevap[r] * t.et_reftopot[r]
        # should we include tempcor in SBM?
        # PotEvap = PotenEvap #??

        if Δt >= Hour(23)
            throughfall, interception, stemflow, canopystorage = rainfall_interception_gash(
                cmax,
                e_r,
                canopygapfraction,
                t.precipitation[r],
                t.canopystorage[r],
                maxevap = potevap,
            )
            pottrans_soil = max(0.0, potevap - interception) # now in mm
        else
            netinterception, throughfall, stemflow, leftover, interception, canopystorage =
                rainfall_interception_modrut(
                    t.precipitation[r],
                    potevap,
                    t.canopystorage[r],
                    canopygapfraction,
                    cmax,
                )
            pottrans_soil = max(0.0, leftover)  # now in mm
            interception = netinterception
        end

        if modelsnow
            tsoil = t.tsoil[r] + t.w_soil[r] * (t.temperature[r] - t.tsoil[r])
            snow, snowwater, snowmelt, rainfallplusmelt, snowfall = snowpack_hbv(
                t.snow[r],
                t.snowwater[r],
                throughfall + stemflow,
                t.temperature[r],
                t.tti[r],
                t.tt[r],
                t.ttm[r],
                t.cfmax[r],
                t.whc[r],
            )
        end

        # update the outputs and states
        t.cmax[r] = cmax
        t.canopygapfraction[r] = canopygapfraction
        t.canopystorage[r] = canopystorage
        t.interception[r] = interception
        t.stemflow[r] = stemflow
        t.throughfall[r] = throughfall
        t.pottrans_soil[r] = pottrans_soil
        if modelsnow
            t.snow[r] = snow
            t.snowwater[r] = snowwater
            t.tsoil[r] = tsoil
            t.rainfallplusmelt[r] = rainfallplusmelt
        end
    end
end

function update_until_recharge(t, config)

    # start dummy variables (should be generated from model reader and from Config.jl TOML)
    soilinfreduction = Bool(get(config.model, "soilinfreduction", 0))
    modelglacier = Bool(get(config.model, "modelglacier", 0))
    modelsnow = Bool(get(config.model, "modelsnow", 0))
    transfermethod = Bool(get(config.model, "transfermethod", 0))
    #potevap = 4.0
    #precipitation = 3.0
    #temperature = 10.0
    wl_land = 0.0 # from kinematic wave land
    wl_river = 0.10 # from kinematic river
    irsupply_mm = 0.0
    ust = Bool(get(config.model, "whole_ust_available", 0)) # should be removed from optional setting and code?
    Δt = Second(config.input.timestepsecs)
    #basetimestep = Second(Day(1))
    # end dummpy variables

    for r in eachindex(t.nlayers)
        if modelsnow
            rainfallplusmelt = t.rainfallplusmelt[r]
            if modelglacier
                # Run Glacier module and add the snowpack on-top of it.
                # Estimate the fraction of snow turned into ice (HBV-light).
                # Estimate glacier melt.

                snow, snow2glacier, glacierstore, glaciermelt = glacier_hbv(
                    t.glacierfrac[r],
                    t.glacierstore[r],
                    t.snow[r],
                    t.temperature[r],
                    t.g_tt[r],
                    t.g_cfmax[r],
                    t.g_sifrac[r],
                    Δt,
                    basetimestep,
                )
                # Convert to mm per grid cell and add to snowmelt
                glaciermelt = glaciermelt * t.glacierfrac[r]
                rainfallplusmelt = rainfallplusmelt + glaciermelt

            end
        else
            rainfallplusmelt = t.stemflow[r] + t.throughfall[r]
        end

        avail_forinfilt = rainfallplusmelt + irsupply_mm
        ustoredepth = sum(@view t.ustorelayerdepth[r][1:t.nlayers[r]])
        uStorecapacity = t.soilwatercapacity[r] - t.satwaterdepth[r] - ustoredepth

        runoff_river = min(1.0, t.riverfrac[r]) * avail_forinfilt
        runoff_land = min(1.0, t.waterfrac[r]) * avail_forinfilt
        avail_forinfilt = max(avail_forinfilt - runoff_river - runoff_land, 0.0)

        rootingdepth = min(t.soilthickness[r] * 0.99, t.rootingdepth[r])

        ae_openw_r =
            min(wl_river * 1000.0 * t.riverfrac[r], t.riverfrac[r] * t.pottrans_soil[r])
        ae_openw_l =
            min(wl_land * 1000.0 * t.waterfrac[r], t.waterfrac[r] * t.pottrans_soil[r])

        restevap = t.pottrans_soil[r] - ae_openw_r - ae_openw_l

        # evap available for soil evaporation and transpiration
        potsoilevap = restevap * t.canopygapfraction[r]
        pottrans = restevap * (1.0 - t.canopygapfraction[r])

        # Calculate the initial capacity of the unsaturated store
        ustorecapacity = t.soilwatercapacity[r] - t.satwaterdepth[r] - ustoredepth

        # Calculate the infiltration flux into the soil column
        infiltsoilpath, infiltsoil, infiltpath, soilinf, pathinf, infiltexcess =
            infiltration(
                avail_forinfilt,
                t.pathfrac[r],
                t.cf_soil[r],
                t.tsoil[r],
                t.infiltcapsoil[r],
                t.infiltcappath[r],
                ustorecapacity,
                modelsnow,
                soilinfreduction,
            )


        usl, n_usl = set_layerthickness(t.zi[r], t.sumlayers[r], t.act_thickl[r])
        z = cumsum(usl)
        usld = copy(t.ustorelayerdepth[r])

        ast = 0.0
        soilevapunsat = 0.0
        if n_usl > 0
            # Using the surface infiltration rate, calculate the flow rate between the
            # different soil layers that contain unsaturated storage assuming gravity
            # based flow only, estimate the gravity based flux rate to the saturated zone
            # (ast) and the updated unsaturated storage for each soil layer.
            if transfermethod == 1 && t.maxlayers[r] == 1
                ustorelayerdepth = t.ustorelayerdepth[r][1] + infiltsoilpath
                kv_z = t.kvfrac[r][1] * t.kv₀[r] * exp(-t.f[r] * t.zi[r])
                ustorelayerdepth, ast = unsatzone_flow_sbm(
                    ustorelayerdepth,
                    t.soilwatercapacity[r],
                    t.satwaterdepth[r],
                    kv_z,
                    usl[1],
                    t.θₛ[r],
                    t.θᵣ[r],
                )
                usld = setindex(usld, ustorelayerdepth, m)
            else
                for m = 1:n_usl
                    l_sat = usl[m] * (t.θₛ[r] - t.θᵣ[r])
                    kv_z = t.kvfrac[r][m] * t.kv₀[r] * exp(-t.f[r] * z[m])
                    ustorelayerdepth = m == 1 ? t.ustorelayerdepth[r][m] + infiltsoilpath :
                        t.ustorelayerdepth[r][m] + ast
                    ustorelayerdepth, ast =
                        unsatzone_flow_layer(ustorelayerdepth, kv_z, l_sat, t.c[r][m])
                    usld = setindex(usld, ustorelayerdepth, m)
                end
            end

            # then evapotranspiration from layers
            # Calculate saturation deficity
            saturationdeficit = t.soilwatercapacity[r] - t.satwaterdepth[r]

            # First calculate the evaporation of unsaturated storage into the
            # atmosphere from the upper layer.
            if t.maxlayers[r] == 1
                soilevapunsat =
                    potsoilevap * min(1.0, saturationdeficit / t.soilwatercapacity[r])
            else
                # In case only the most upper soil layer contains unsaturated storage
                if n_usl == 1
                    # Check if groundwater level lies below the surface
                    soilevapunsat =
                        potsoilevap * min(1.0, usld[1] / (t.zi[r] * (t.θₛ[r] - t.θᵣ[r])))
                else
                    # In case first layer contains no saturated storage
                    soilevapunsat =
                        potsoilevap * min(1.0, usld[1] / (usl[1] * ((t.θₛ[r] - t.θᵣ[r]))))
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

        if t.maxlayers[r] == 1
            soilevapsat = 0.0
        else
            # this check is an improvement compared to Python (only checked for n_usl == 1)
            if n_usl == 0 || n_usl == 1
                soilevapsat =
                    potsoilevap *
                    min(1.0, (t.act_thickl[r][1] - t.zi[r]) / t.act_thickl[r][1])
                soilevapsat =
                    min(soilevapsat, (t.act_thickl[r][1] - t.zi[r]) * (t.θₛ[r] - t.θᵣ[r]))
            else
                soilevapsat = 0.0
            end
        end
        soilevap = soilevapunsat + soilevapsat
        satwaterdepth = t.satwaterdepth[r] - soilevapsat

        # transpiration from saturated store
        wetroots = scurve(t.zi[r], a = rootingdepth, c = t.rootdistpar[r])
        actevapsat = min(pottrans * wetroots, satwaterdepth)
        satwaterdepth = satwaterdepth - actevapsat
        restpottrans = pottrans - actevapsat

        # actual transpiration from ustore
        actevapustore = 0.0
        for k = 1:n_usl
            ustorelayerdepth, actevapustore, restpottrans = acttransp_unsat_sbm(
                rootingdepth,
                usld[k],
                t.sumlayers[r][k],
                restpottrans,
                actevapustore,
                t.c[r][k],
                usl[k],
                t.θₛ[r],
                t.θᵣ[r],
                t.hb[r],
                ust,
            )
            usld = setindex(usld, ustorelayerdepth, k)
        end

        # check soil moisture balance per layer
        du = 0.0
        for k = n_usl:-1:1
            du = max(0.0, usld[k] - usl[k] * (t.θₛ[r] - t.θᵣ[r]))
            usld = setindex(usld, usld[k] - du, k)
            if k > 1
                usld = setindex(usld, usld[k-1] + du, k - 1)
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
            ksat = t.kvfrac[r][n_usl] * t.kv₀[r] * exp(-t.f[r] * t.zi[r])
            ustorecapacity =
                t.soilwatercapacity[r] - t.satwaterdepth[r] -
                sum(@view usld[1:t.nlayers[r]])
            maxcapflux =
                max(0.0, min(ksat, actevapustore, ustorecapacity, t.satwaterdepth[r]))

            if t.zi[r] > rootingdepth
                capfluxscale =
                    t.capscale[r] / (t.capscale[r] + t.zi[r] - rootingdepth) *
                    Float64(Δt.value) / Float64(basetimestep.value)
            else
                capfluxscale = 0.0
            end
            capflux = maxcapflux * capfluxscale

            netcapflux = capflux
            for k = n_usl:-1:1
                toadd = min(netcapflux, max(usl[k] * (t.θₛ[r] - t.θᵣ[r]) - usld[k], 0.0))
                usld = setindex(usld, usld[k] + toadd, k)
                netcapflux = netcapflux - toadd
                actcapflux = actcapflux + toadd
            end
        end
        deepksat = t.kv₀[r] * exp(-t.f[r] * t.soilthickness[r])
        deeptransfer = min(t.satwaterdepth[r], deepksat)
        actleakage = max(0.0, min(t.maxleakage[r], deeptransfer))

        # recharge (mm) for saturated zone
        recharge = (transfer - actcapflux - actleakage - actevapsat - soilevapsat)
        transpiration = actevapsat + actevapustore
        actevap = soilevap + transpiration + ae_openw_r + ae_openw_l

        # update the outputs and states
        t.actinfilt[r] = actinfilt
        t.infiltexcess[r] = infiltexcess
        t.recharge[r] = recharge
        t.transpiration[r] = transpiration
        t.soilevap[r] = soilevap
        t.ae_openw_r[r] = ae_openw_r
        t.ae_openw_l[r] = ae_openw_l
        t.runoff_land[r] = runoff_land
        t.runoff_river[r] = runoff_river
        t.actevapsat[r] = actevapsat
        t.actevap[r] = actevap
        t.ustorelayerdepth[r] = usld
        t.transfer[r] = transfer
        t.actinfiltsoil[r] = actinfiltsoil
        t.actinfiltpath[r] = actinfiltpath
        t.excesswater[r] = excesswater
        t.excesswatersoil[r] = excesswatersoil
        t.excesswaterpath[r] = excesswaterpath
        t.rainfallplusmelt[r] = rainfallplusmelt
    end
end

function update_after_lateralflow(t, zi, exfiltsatwater)

    for r in eachindex(t.nlayers)
        usl, n_usl = set_layerthickness(zi[r], t.sumlayers[r], t.act_thickl[r])
        # exfiltration from ustore
        usld = copy(t.ustorelayerdepth[r])
        exfiltustore = 0.0
        for k = n_usl:-1:1
            exfiltustore = max(0, usld[k] - usl[k] * (t.θₛ[r] - t.θᵣ[r]))
            usld = setindex(usld, usld[k] - exfiltustore, k)
            if k > 1
                usld = setindex(usld, usld[k-1] + exfiltustore, k - 1)
            end
        end

        ustoredepth = sum(@view usld[1:n_usl])

        runoff = max(
            exfiltustore +
            exfiltsatwater[r] +
            t.excesswater[r] +
            t.runoff_land[r] +
            t.infiltexcess[r] - t.ae_openw_l[r],
            0.0,
        )

        # volumetric water content per soil layer and root zone
        vwc = copy(t.vwc[r])
        vwc_perc = copy(t.vwc_perc[r])
        for k = 1:t.nlayers[r]
            if k <= n_usl
                vwc = setindex(
                    vwc,
                    (usld[k] + (t.act_thickl[r][k] - usl[k]) * (t.θₛ[r] - t.θᵣ[r])) /
                    usl[k] + t.θᵣ[r],
                    k,
                )
            else
                vwc = setindex(vwc, t.θₛ[r], k)
            end
            vwc_perc = setindex(vwc_perc, (vwc[k] / t.θₛ[r]) * 100.0, k)
        end

        rootstore_unsat = 0
        for k = 1:n_usl
            rootstore_unsat =
                rootstore_unsat +
                (max(0.0, t.rootingdepth[r] - t.sumlayers[r][k]) / usl[k]) * usld[k]
        end

        rootstore_sat = max(0.0, t.rootingdepth[r] - zi[r]) * (t.θₛ[r] - t.θᵣ[r])
        rootstore = rootstore_sat + rootstore_unsat
        vwc_root = rootstore / t.rootingdepth[r] + t.θᵣ[r]
        vwc_percroot = (vwc_root / t.θₛ[r]) * 100.0

        satwaterdepth = (t.soilthickness[r] - zi[r]) * (t.θₛ[r] - t.θᵣ[r])

        # update the outputs and states
        t.ustorelayerdepth[r] = usld
        t.ustoredepth[r] = ustoredepth
        t.satwaterdepth[r] = satwaterdepth
        t.exfiltsatwater[r] = exfiltsatwater[r]
        t.exfiltustore[r] = exfiltustore
        t.runoff[r] = runoff
        t.vwc[r] = vwc
        t.vwc_perc[r] = vwc_perc
        t.rootstore[r] = rootstore
        t.vwc_root[r] = vwc_root
        t.vwc_percroot[r] = vwc_percroot
        t.zi[r] = zi[r]
    end
end
