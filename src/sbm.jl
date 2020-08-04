Base.@kwdef struct SBM{T,N,M}
    # Maximum number of soil layers
    maxlayers::Int
    # number of cells
    n::Int
    # Number of soil layers
    nlayers::Vector{Int}
    # length of cells in y direction [m]
    yl::Vector{T}
    # length of cells in x direction [m]
    xl::Vector{T}
    # Fraction of river [-]
    riverfrac::Vector{T}
    # Saturated water content (porosity) [mm mm⁻¹]
    θₛ::Vector{T}
    # Residual water content [mm mm⁻¹]
    θᵣ::Vector{T}
    # Effictive porosity [mm mm⁻¹]
    θₑ::Vector{T}
    # Vertical hydraulic conductivity [mm Δt⁻¹] at soil surface
    kv₀::Vector{T}
    # Muliplication factor [-] applied to kv_z (vertical flow)
    kvfrac::Vector{SVector{N,T}}
    # Parameter [mm] controlling f
    m::Vector{T}
    # Air entry pressure [cm] of soil (Brooks-Corey)
    hb::Vector{T}
    # Soil thickness [mm]
    soilthickness::Vector{T}
    # Thickness of soil layers [mm]
    act_thickl::Vector{SVector{N,T}}
    # Cumulative sum of soil layers [mm], starting at soil surface (0)
    sumlayers::Vector{SVector{M,T}}
    # Infiltration capacity of the compacted areas [mm Δt⁻¹]
    infiltcappath::Vector{T}
    # Soil infiltration capacity [mm/Δt]
    infiltcapsoil::Vector{T}
    # Maximum leakage [mm/Δt] from saturated zone
    maxleakage::Vector{T}
    # Fraction of open water (excluding rivers) [-]
    waterfrac::Vector{T}
    # Fraction of compacted area  [-]
    pathfrac::Vector{T}
    # Vertical elevation [m]
    altitude::Vector{T}
    # Rooting depth [mm]
    rootingdepth::Vector{T}
    # Controls how roots are linked to water table [-]
    rootdistpar::Vector{T}
    # Parameter [mm] controlling capilary rise
    capscale::Vector{T}
    #Multiplication factor [-] to correct
    et_reftopot::Vector{T}
    # Brooks-Corey power coefﬁcient [-] for each soil layer
    c::Vector{SVector{N,T}}
    # Stemflow [mm]
    stemflow::Vector{T}
    # Throughfall [mm]
    throughfall::Vector{T}
    # A scaling parameter [mm⁻¹] (controls exponential decline of kv₀)
    f::Vector{T} = θₑ ./ m
    # Amount of water in the unsaturated store, per layer [mm]
    ustorelayerdepth::Vector{SVector{N,T}} = act_thickl .* 0.0
    # Saturated store [mm]
    satwaterdepth::Vector{T}
    # Pseudo-water table depth [mm] (top of the saturated zone)
    zi::Vector{T} = max.(0.0, soilthickness .- satwaterdepth ./ θₑ)
    # Soilwater capacity [mm]
    soilwatercapacity::Vector{T}
    # Canopy storage [mm]
    canopystorage::Vector{T} = fill(0.0, n)
    # Maximum canopy storage [mm]
    cmax::Vector{T}
    # Canopy gap fraction [-]
    canopygapfraction::Vector{T}
    # Gash interception model parameter, ratio of the average evaporation from the
    # wet canopy [mm Δt⁻¹] and the average precipitation intensity [mm Δt⁻¹] on a saturated canopy
    e_r::Vector{T}
    # Precipitation [mm]
    precipitation::Vector{T} = fill(mv, n)
    # Temperature [ᵒC]
    temperature::Vector{T} = fill(mv, n)
    # Potential evapotranspiration [mm]
    potevap::Vector{T} = fill(mv, n)
    # Potential transpiration, open water, river and soil evaporation (after subtracting interception from potevap)
    pottrans_soil::Vector{T} = fill(mv, n)
    # Transpiration [mm]
    transpiration::Vector{T} = fill(mv, n)
    # Actual evaporation from unsaturated store [mm]
    ae_ustore::Vector{T} = fill(mv, n)
    # Actual evaporation from saturated store [mm]
    ae_sat::Vector{T} = fill(mv, n)
    # Interception [mm]
    interception::Vector{T} = fill(mv, n)
    # Soil evaporation [mm]
    soilevap::Vector{T} = fill(mv, n)
    # Actual evaporation from saturated store (transpiration and soil evaporation) [mm]
    actevapsat::Vector{T} = fill(mv, n)
    # Total actual evapotranspiration [mm]
    actevap::Vector{T} = fill(mv, n)
    # Runoff from river based on riverfrac [mm]
    runoff_river::Vector{T} = fill(mv, n)
    # Runoff from land based on waterfrac [mm]
    runoff_land::Vector{T} = fill(mv, n)
    # Actual evaporation from open water (land) [mm]
    ae_openw_l::Vector{T} = fill(mv, n)
    # Actual evaporation from river [mm]
    ae_openw_r::Vector{T} = fill(mv, n)
    # Water available for infiltration [mm]
    avail_forinfilt::Vector{T} = fill(mv, n)
    # Actual infiltration into the unsaturated zone [mm]
    actinfilt::Vector{T} = fill(mv, n)
    # Actual infiltration non-compacted fraction [mm]
    actinfiltsoil::Vector{T} = fill(mv, n)
    # Actual infiltration compacted fraction [mm]
    actinfiltpath::Vector{T} = fill(mv, n)
    # Infiltration excess water [mm]
    infiltexcess::Vector{T} = fill(mv, n)
    # Water that cannot infiltrate due to saturated soil (saturation excess) [mm]
    excesswater::Vector{T} = fill(mv, n)
    # Water exfiltrating during saturation excess conditions [mm]
    exfiltsatwater::Vector{T} = fill(mv, n)
    # Water exfiltrating from unsaturated store because of change in water table [mm]
    exfiltustore::Vector{T} = fill(mv, n)
    # Excess water for non-compacted fraction [mm]
    excesswatersoil::Vector{T} = fill(mv, n)
    # Excess water for compacted fraction [mm]
    excesswaterpath::Vector{T} = fill(mv, n)
    # Total surface runoff from infiltration and saturation excess [mm]
    runoff::Vector{T} = fill(mv, n)
    # Volumetric water content [mm mm⁻¹] per soil layer (including θᵣ and saturated zone)
    vwc::Vector{SVector{N,T}}
    # Volumetric water content [%] per soil layer (including θᵣ and saturated zone)
    vwc_perc::Vector{SVector{N,T}}
    # Root water storage [mm] in unsaturated and saturated zone (excluding θᵣ)
    rootstore::Vector{T} = fill(mv, n)
    # Volumetric water content [mm mm⁻¹] in root zone (including θᵣ and saturated zone)
    vwc_root::Vector{T} = fill(mv, n)
    # Volumetric water content [%] in root zone (including θᵣ and saturated zone)
    vwc_percroot::Vector{T} = fill(mv, n)
    # Amount of available water in the unsaturated zone [mm]
    ustoredepth::Vector{T} = fill(mv, n)
    # Downward flux from unsaturated to saturated zone [mm]
    transfer::Vector{T} = fill(mv, n)
    # Capillary rise [mm]
    capflux::Vector{T} = fill(mv, n)
    # Net recharge to saturated store [mm]
    recharge::Vector{T} = fill(mv, n)
    ### Snow parameters ###
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹]
    cfmax::Vector{T} = fill(mv, n)
    # Threshold temperature for snowfall [ᵒC]
    tt::Vector{T} = fill(mv, n)
    # Threshold temperature interval length [ᵒC]
    tti::Vector{T} = fill(mv, n)
    # Threshold temperature for snowmelt [ᵒC]
    ttm::Vector{T} = fill(mv, n)
    # Water holding capacity as fraction of current snow pack [-]
    whc::Vector{T} = fill(mv, n)
    # Soil temperature smooth factor [-]
    w_soil::Vector{T} = fill(mv, n)
    # Controls soil infiltration reduction factor when soil is frozen [-]
    cf_soil::Vector{T} = fill(mv, n)
    # Snow storage [mm]
    snow::Vector{T} = fill(0.0, n)
    # Liquid water content in the snow pack [mm]
    snowwater::Vector{T} = fill(0.0, n)
    # Snow melt + precipitation as rainfall [mm]
    rainfallplusmelt::Vector{T} = fill(mv, n)
    # Top soil temperature [ᵒC]
    tsoil::Vector{T} = fill(10.0, n)
    ## Interception related to LAI climatology ###
    # Specific leaf storage [mm]
    sl::Vector{T} = fill(mv, n)
    # Storage woody part of vegetation [mm]
    swood::Vector{T} = fill(mv, n)
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    kext::Vector{T} = fill(mv, n)
    # Leaf area index [m² m⁻²]
    lai::Vector{T} = fill(mv, n)

    function SBM{T,N,M}(args...) where {T,N,M}
        equal_size_vectors(args)
        return new(args...)
    end
end

function update_until_snow(sbm::SBM, config)

    # # start dummy variables (should be generated from model reader and from Config.jl TOML)
    do_lai = haskey(config.cyclic.parameters, "lai")
    modelglacier = Bool(get(config.model, "modelglacier", 0))
    modelsnow = Bool(get(config.model, "modelsnow", 0))
    #potevap = 4.0
    #precipitation = 3.0
    #temperature = 10.0
    Δt = Second(config.timestepsecs)
    #basetimestep = Second(Day(1))
    # end dummpy variables
    for i = 1:sbm.n
        if do_lai
            cmax = sbm.sl[i] * sbm.lai[i] + sbm.swood[i]
            canopygapfraction = exp(-sbm.kext[i] * sbm.lai[i])
            ewet = (1.0 - exp(-sbm.kext[i] * sbm.lai[i])) * sbm.potevap[i]
            e_r = sbm.precipitation[i] > 0.0 ?
                min(0.25, ewet / max(0.0001, sbm.precipitation[i])) : 0.0
        end

        potevap = sbm.potevap[i] * sbm.et_reftopot[i]
        # should we include tempcor in SBM?
        # PotEvap = PotenEvap #??

        if Δt >= Hour(23)
            throughfall, interception, stemflow, canopystorage = rainfall_interception_gash(
                cmax,
                e_r,
                canopygapfraction,
                sbm.precipitation[i],
                sbm.canopystorage[i],
                maxevap = potevap,
            )
            pottrans_soil = max(0.0, potevap - interception) # now in mm
        else
            netinterception, throughfall, stemflow, leftover, interception, canopystorage =
                rainfall_interception_modrut(
                    sbm.precipitation[i],
                    potevap,
                    sbm.canopystorage[i],
                    canopygapfraction,
                    cmax,
                )
            pottrans_soil = max(0.0, leftover)  # now in mm
            interception = netinterception
        end

        if modelsnow
            tsoil = sbm.tsoil[i] + sbm.w_soil[i] * (sbm.temperature[i] - sbm.tsoil[i])
            snow, snowwater, snowmelt, rainfallplusmelt, snowfall = snowpack_hbv(
                sbm.snow[i],
                sbm.snowwater[i],
                throughfall + stemflow,
                sbm.temperature[i],
                sbm.tti[i],
                sbm.tt[i],
                sbm.ttm[i],
                sbm.cfmax[i],
                sbm.whc[i],
            )
        end

        # update the outputs and states
        sbm.cmax[i] = cmax
        sbm.canopygapfraction[i] = canopygapfraction
        sbm.canopystorage[i] = canopystorage
        sbm.interception[i] = interception
        sbm.stemflow[i] = stemflow
        sbm.throughfall[i] = throughfall
        sbm.pottrans_soil[i] = pottrans_soil
        if modelsnow
            sbm.snow[i] = snow
            sbm.snowwater[i] = snowwater
            sbm.tsoil[i] = tsoil
            sbm.rainfallplusmelt[i] = rainfallplusmelt
        end
    end
end

function update_until_recharge(sbm::SBM, config)

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
    Δt = Second(config.timestepsecs)
    #basetimestep = Second(Day(1))
    # end dummpy variables

    for i = 1:sbm.n
        if modelsnow
            rainfallplusmelt = sbm.rainfallplusmelt[i]
            if modelglacier
                # Run Glacier module and add the snowpack on-top of it.
                # Estimate the fraction of snow turned into ice (HBV-light).
                # Estimate glacier melt.

                snow, snow2glacier, glacierstore, glaciermelt = glacier_hbv(
                    sbm.glacierfrac[i],
                    sbm.glacierstore[i],
                    sbm.snow[i],
                    sbm.temperature[i],
                    sbm.g_tt[i],
                    sbm.g_cfmax[i],
                    sbm.g_sifrac[i],
                    Δt,
                )
                # Convert to mm per grid cell and add to snowmelt
                glaciermelt = glaciermelt * sbm.glacierfrac[i]
                rainfallplusmelt = rainfallplusmelt + glaciermelt

            end
        else
            rainfallplusmelt = sbm.stemflow[i] + sbm.throughfall[i]
        end

        avail_forinfilt = rainfallplusmelt + irsupply_mm
        ustoredepth = sum(@view sbm.ustorelayerdepth[i][1:sbm.nlayers[i]])
        uStorecapacity = sbm.soilwatercapacity[i] - sbm.satwaterdepth[i] - ustoredepth

        runoff_river = min(1.0, sbm.riverfrac[i]) * avail_forinfilt
        runoff_land = min(1.0, sbm.waterfrac[i]) * avail_forinfilt
        avail_forinfilt = max(avail_forinfilt - runoff_river - runoff_land, 0.0)

        rootingdepth = min(sbm.soilthickness[i] * 0.99, sbm.rootingdepth[i])

        ae_openw_r = min(
            wl_river * 1000.0 * sbm.riverfrac[i],
            sbm.riverfrac[i] * sbm.pottrans_soil[i],
        )
        ae_openw_l = min(
            wl_land * 1000.0 * sbm.waterfrac[i],
            sbm.waterfrac[i] * sbm.pottrans_soil[i],
        )

        restevap = sbm.pottrans_soil[i] - ae_openw_r - ae_openw_l

        # evap available for soil evaporation and transpiration
        potsoilevap = restevap * sbm.canopygapfraction[i]
        pottrans = restevap * (1.0 - sbm.canopygapfraction[i])

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
        usld = copy(sbm.ustorelayerdepth[i])

        ast = 0.0
        soilevapunsat = 0.0
        if n_usl > 0
            # Using the surface infiltration rate, calculate the flow rate between the
            # different soil layers that contain unsaturated storage assuming gravity
            # based flow only, estimate the gravity based flux rate to the saturated zone
            # (ast) and the updated unsaturated storage for each soil layer.
            if transfermethod == 1 && sbm.maxlayers == 1
                ustorelayerdepth = sbm.ustorelayerdepth[i][1] + infiltsoilpath
                kv_z = sbm.kvfrac[i][1] * sbm.kv₀[i] * exp(-sbm.f[i] * sbm.zi[i])
                ustorelayerdepth, ast = unsatzone_flow_sbm(
                    ustorelayerdepth,
                    sbm.soilwatercapacity[i],
                    sbm.satwaterdepth[i],
                    kv_z,
                    usl[1],
                    t.θₛ[i],
                    t.θᵣ[i],
                )
                usld = setindex(usld, ustorelayerdepth, m)
            else
                for m = 1:n_usl
                    l_sat = usl[m] * (sbm.θₛ[i] - sbm.θᵣ[i])
                    kv_z = sbm.kvfrac[i][m] * sbm.kv₀[i] * exp(-sbm.f[i] * z[m])
                    ustorelayerdepth =
                        m == 1 ? sbm.ustorelayerdepth[i][m] + infiltsoilpath :
                        sbm.ustorelayerdepth[i][m] + ast
                    ustorelayerdepth, ast =
                        unsatzone_flow_layer(ustorelayerdepth, kv_z, l_sat, sbm.c[i][m])
                    usld = setindex(usld, ustorelayerdepth, m)
                end
            end

            # then evapotranspiration from layers
            # Calculate saturation deficity
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
                        min(1.0, usld[1] / (sbm.zi[i] * (sbm.θₛ[i] - sbm.θᵣ[i])))
                else
                    # In case first layer contains no saturated storage
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (usl[1] * ((sbm.θₛ[i] - sbm.θᵣ[i]))))
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
            # this check is an improvement compared to Python (only checked for n_usl == 1)
            if n_usl == 0 || n_usl == 1
                soilevapsat =
                    potsoilevap *
                    min(1.0, (sbm.act_thickl[i][1] - sbm.zi[i]) / sbm.act_thickl[i][1])
                soilevapsat = min(
                    soilevapsat,
                    (sbm.act_thickl[i][1] - sbm.zi[i]) * (sbm.θₛ[i] - sbm.θᵣ[i]),
                )
            else
                soilevapsat = 0.0
            end
        end
        soilevap = soilevapunsat + soilevapsat
        satwaterdepth = sbm.satwaterdepth[i] - soilevapsat

        # transpiration from saturated store
        wetroots = scurve(sbm.zi[i], a = rootingdepth, c = sbm.rootdistpar[i])
        actevapsat = min(pottrans * wetroots, satwaterdepth)
        satwaterdepth = satwaterdepth - actevapsat
        restpottrans = pottrans - actevapsat

        # actual transpiration from ustore
        actevapustore = 0.0
        for k = 1:n_usl
            ustorelayerdepth, actevapustore, restpottrans = acttransp_unsat_sbm(
                rootingdepth,
                usld[k],
                sbm.sumlayers[i][k],
                restpottrans,
                actevapustore,
                sbm.c[i][k],
                usl[k],
                sbm.θₛ[i],
                sbm.θᵣ[i],
                sbm.hb[i],
                ust,
            )
            usld = setindex(usld, ustorelayerdepth, k)
        end

        # check soil moisture balance per layer
        du = 0.0
        for k = n_usl:-1:1
            du = max(0.0, usld[k] - usl[k] * (sbm.θₛ[i] - sbm.θᵣ[i]))
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
            ksat = sbm.kvfrac[i][n_usl] * sbm.kv₀[i] * exp(-sbm.f[i] * sbm.zi[i])
            ustorecapacity =
                sbm.soilwatercapacity[i] - sbm.satwaterdepth[i] -
                sum(@view usld[1:sbm.nlayers[i]])
            maxcapflux =
                max(0.0, min(ksat, actevapustore, ustorecapacity, sbm.satwaterdepth[i]))

            if sbm.zi[i] > rootingdepth
                capfluxscale =
                    sbm.capscale[i] / (sbm.capscale[i] + sbm.zi[i] - rootingdepth) *
                    Float64(Δt.value) / Float64(basetimestep.value)
            else
                capfluxscale = 0.0
            end
            capflux = maxcapflux * capfluxscale

            netcapflux = capflux
            for k = n_usl:-1:1
                toadd =
                    min(netcapflux, max(usl[k] * (sbm.θₛ[i] - sbm.θᵣ[i]) - usld[k], 0.0))
                usld = setindex(usld, usld[k] + toadd, k)
                netcapflux = netcapflux - toadd
                actcapflux = actcapflux + toadd
            end
        end
        deepksat = sbm.kv₀[i] * exp(-sbm.f[i] * sbm.soilthickness[i])
        deeptransfer = min(sbm.satwaterdepth[i], deepksat)
        actleakage = max(0.0, min(sbm.maxleakage[i], deeptransfer))

        # recharge (mm) for saturated zone
        recharge = (transfer - actcapflux - actleakage - actevapsat - soilevapsat)
        transpiration = actevapsat + actevapustore
        actevap = soilevap + transpiration + ae_openw_r + ae_openw_l

        # update the outputs and states
        sbm.actinfilt[i] = actinfilt
        sbm.infiltexcess[i] = infiltexcess
        sbm.recharge[i] = recharge
        sbm.transpiration[i] = transpiration
        sbm.soilevap[i] = soilevap
        sbm.ae_openw_r[i] = ae_openw_r
        sbm.ae_openw_l[i] = ae_openw_l
        sbm.runoff_land[i] = runoff_land
        sbm.runoff_river[i] = runoff_river
        sbm.actevapsat[i] = actevapsat
        sbm.actevap[i] = actevap
        sbm.ustorelayerdepth[i] = usld
        sbm.transfer[i] = transfer
        sbm.actinfiltsoil[i] = actinfiltsoil
        sbm.actinfiltpath[i] = actinfiltpath
        sbm.excesswater[i] = excesswater
        sbm.excesswatersoil[i] = excesswatersoil
        sbm.excesswaterpath[i] = excesswaterpath
        sbm.rainfallplusmelt[i] = rainfallplusmelt
    end
end

function update_after_lateralflow(sbm::SBM, zi, exfiltsatwater)

    for i = 1:sbm.n
        usl, n_usl = set_layerthickness(zi[i], sbm.sumlayers[i], sbm.act_thickl[i])
        # exfiltration from ustore
        usld = copy(sbm.ustorelayerdepth[i])
        exfiltustore = 0.0
        for k = n_usl:-1:1
            exfiltustore = max(0, usld[k] - usl[k] * (sbm.θₛ[i] - sbm.θᵣ[i]))
            usld = setindex(usld, usld[k] - exfiltustore, k)
            if k > 1
                usld = setindex(usld, usld[k-1] + exfiltustore, k - 1)
            end
        end

        ustoredepth = sum(@view usld[1:n_usl])

        runoff = max(
            exfiltustore +
            exfiltsatwater[i] +
            sbm.excesswater[i] +
            sbm.runoff_land[i] +
            sbm.infiltexcess[i] - sbm.ae_openw_l[i],
            0.0,
        )

        # volumetric water content per soil layer and root zone
        vwc = copy(sbm.vwc[i])
        vwc_perc = copy(sbm.vwc_perc[i])
        for k = 1:sbm.nlayers[i]
            if k <= n_usl
                vwc = setindex(
                    vwc,
                    (usld[k] + (sbm.act_thickl[i][k] - usl[k]) * (sbm.θₛ[i] - sbm.θᵣ[i])) / usl[k] + sbm.θᵣ[i],
                    k,
                )
            else
                vwc = setindex(vwc, sbm.θₛ[i], k)
            end
            vwc_perc = setindex(vwc_perc, (vwc[k] / sbm.θₛ[i]) * 100.0, k)
        end

        rootstore_unsat = 0
        for k = 1:n_usl
            rootstore_unsat =
                rootstore_unsat +
                (max(0.0, sbm.rootingdepth[i] - sbm.sumlayers[i][k]) / usl[k]) * usld[k]
        end

        rootstore_sat = max(0.0, sbm.rootingdepth[i] - zi[i]) * (sbm.θₛ[i] - sbm.θᵣ[i])
        rootstore = rootstore_sat + rootstore_unsat
        vwc_root = rootstore / sbm.rootingdepth[i] + sbm.θᵣ[i]
        vwc_percroot = (vwc_root / sbm.θₛ[i]) * 100.0

        satwaterdepth = (sbm.soilthickness[i] - zi[i]) * (sbm.θₛ[i] - sbm.θᵣ[i])

        # update the outputs and states
        sbm.ustorelayerdepth[i] = usld
        sbm.ustoredepth[i] = ustoredepth
        sbm.satwaterdepth[i] = satwaterdepth
        sbm.exfiltsatwater[i] = exfiltsatwater[i]
        sbm.exfiltustore[i] = exfiltustore
        sbm.runoff[i] = runoff
        sbm.vwc[i] = vwc
        sbm.vwc_perc[i] = vwc_perc
        sbm.rootstore[i] = rootstore
        sbm.vwc_root[i] = vwc_root
        sbm.vwc_percroot[i] = vwc_percroot
        sbm.zi[i] = zi[i]
    end
end
