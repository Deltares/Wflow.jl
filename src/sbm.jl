const mv = NaN

Base.@kwdef struct SBM{T,N,M}
    maxlayers::Int              # Maximum number of soil layers
    nlayers::Int                # Number of soil layers
    cfmax::T                    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹]
    tt::T                       # Threshold temperature for snowfall [ᵒC]
    ttm::T                      # Threshold temperature for snowmelt [ᵒC]
    tti::T                      # Threshold temperature interval length [ᵒC]
    whc::T                      # Water holding capacity as fraction of current snow pack [-]
    cf_soil::T                  # Controls soil infiltration reduction factor when soil is frozen [-]
    w_soil::T                   # Soil temperature smooth factor [-]
    soilthickness::T            # Soil thickness [mm]
    infiltcapsoil::T            # Soil infiltration capacity [mm/Δt]
    infiltcappath::T            # Infiltration capacity of the compacted areas [mm Δt⁻¹]
    pathfrac::T                 # Fraction of compacted area  [-]
    waterfrac::T                # Fraction of open water (excluding rivers) [-]
    riverfrac::T                # Fraction of river [-]
    θₛ::T                       # Saturated water content (porosity) [mm mm⁻¹]
    θᵣ::T                       # Residual water content [mm mm⁻¹]
    hb::T                       # Air entry pressure [cm] of soil (Brooks-Corey)
    kv₀::T                      # Vertical hydraulic conductivity [mm Δt⁻¹] at soil surface
    kvfrac::SVector{N,T}        # Muliplication factor [-] applied to kv_z (vertical flow)
    maxleakage::T               # Maximum leakage [mm/Δt] from saturated zone
    c::SVector{N,T}             # Brooks-Corey power coefﬁcient [-] for each soil layer
    m::T                        # Parameter [mm] controlling f
    f::T = (θₛ - θᵣ) / m        # A scaling parameter [mm⁻¹] (controls exponential decline of kv₀)
    capscale::T                 # Parameter [mm] controlling capilary rise
    rootdistpar::T              # Controls how roots are linked to water table [-]
    rootingdepth::T             # Rooting depth [mm]
    lai::T                      # Leaf area index [m² m⁻²]
    sl::T                       # Specific leaf storage [mm]
    kext::T                     # Extinction coefficient [-] (to calculate canopy gap fraction)
    swood::T                    # Storage woody part of vegetation [mm]
    cmax::T                     # Maximum canopy storage [mm]
    canopygapfraction::T        # Canopy gap fraction [-]
    e_r::T                      # Gash interception model parameter, ratio of the average evaporation from the wet canopy [mm Δt⁻¹] and the average precipitation intensity [mm Δt⁻¹] on a saturated canopy
    et_reftopot::T              # Multiplication factor [-] to correct
    altitude::T                 # Vertical elevation [m]
    precipitation::T = mv       # Precipitation [mm]
    temperature::T = mv         # Temperature [ᵒC]
    potevap::T = mv             # Potential evapotranspiration [mm]
    pottrans_soil::T = mv       # Potential transpiration, open water, river and soil evaporation (after subtracting interception from potevap)
    transpiration::T = mv       # Transpiration [mm]
    ae_ustore::T = mv           # Actual evaporation from unsaturated store [mm]
    ae_sat::T = mv              # Actual evaporation from saturated store [mm]
    interception::T = mv        # Interception [mm]
    soilevap::T = mv            # Soil evaporation [mm]
    actevapsat::T = mv          # Actual evaporation from saturated store (transpiration and soil evaporation) [mm]
    actevap::T = mv             # Total actual evaporation (transpiration + soil evapation + open water evaporation) [mm]
    ae_openw_l::T = mv          # Actual evaporation from open water (land) [mm]
    ae_openw_r::T = mv          # Actual evaporation from river [mm]
    avail_forinfilt::T = mv     # Water available for infiltration [mm]
    actinfilt::T = mv           # Actual infiltration into the unsaturated zone [mm]
    actinfiltsoil::T = mv       # Actual infiltration non-compacted fraction [mm]
    actinfiltpath::T = mv       # Actual infiltration compacted fraction [mm]
    infiltexcess::T = mv        # Infiltration excess water [mm]
    excesswater::T = mv         # Water that cannot infiltrate due to saturated soil (saturation excess) [mm]
    exfiltsatwater::T = mv      # Water exfiltrating during saturation excess conditions [mm]
    exfiltustore::T = mv        # Water exfiltrating from unsaturated store because of change in water table [mm]
    excesswatersoil::T = mv     # Excess water for non-compacted fraction [mm]
    excesswaterpath::T = mv     # Excess water for compacted fraction [mm]
    runoff::T = mv              # Total surface runoff from infiltration and saturation excess [mm]
    act_thickl::SVector{N,T}    # Thickness of soil layers [mm]
    ustorelayerdepth::SVector{N,T} = fill(0.0, SVector{N,T}) .* act_thickl # Amount of water in the unsaturated store, per layer [mm]
    vwc::SVector{N,T} = fill(mv, SVector{N,T})              # Volumetric water content [mm mm⁻¹] per soil layer (including θᵣ and saturated zone)
    vwc_perc::SVector{N,T} = fill(mv, SVector{N,T})         # Volumetric water content [%] per soil layer (including θᵣ and saturated zone)
    rootstore::T = mv           # Root water storage [mm] in unsaturated and saturated zone (excluding θᵣ)
    vwc_root::T = mv            # Volumetric water content [mm mm⁻¹] in root zone (including θᵣ and saturated zone)
    vwc_percroot::T = mv        # Volumetric water content [%] in root zone (including θᵣ and saturated zone)
    sumlayers::SVector{M,T}     # Cumulative sum of soil layers [mm], starting at soil surface (0)
    ustoredepth::T = mv         # Amount of available water in the unsaturated zone [mm]
    transfer::T = mv            # Downward flux from unsaturated to saturated zone [mm]
    capflux::T = mv             # Capilary rise [mm]
    recharge::T = mv            # Net recharge to saturated store [mm]
    soilwatercapacity::T = soilthickness * (θₛ - θᵣ)            # Soilwater capacity [mm]
    satwaterdepth::T = 0.85 * soilwatercapacity                 # Saturated store [mm]
    zi::T = max(0.0, soilthickness - satwaterdepth / (θₛ - θᵣ)) # Pseudo-water table depth [mm] (top of the saturated zone)
    snow::T = 0.0               # Snow storage [mm]
    snowwater::T = 0.0          # Liquid water content in the snow pack [mm]
    tsoil::T = 10.0             # Top soil temperature [ᵒC]
    canopystorage::T = 0.0      # Canopy storage [mm]
end

"""
    statenames(::Type{SBM})

Returns Array{Symbol,1} for extracting model state fields from SBM struct.
"""
function statenames(::Type{SBM})

    # depends on ini file settings (optional: glaciers, snow, irrigation)
    states = [
        :satwaterdepth,
        :snow,
        :tsoil,
        :ustorelayerdepth,
        :snowwater,
        :canopystorage,
    ]
    # TODO: (warm) states read from netcdf file or cold state (reinit=1, setting in ini file)

end

function update_before_lateralflow(sbm)

    # start dummy variables (should be generated from model reader and from Config.jl TOML)
    do_lai = true
    glacierfrac = false
    modelsnow = true
    soilinfreduction = false
    transfermethod = 0
    potevap = 4.0
    precipitation = 3.0
    temperature = 10.0
    wl_land = 0.0 # from kinematic wave land
    wl_river = 0.10 # from kinematic river
    irsupply_mm = 0.0
    ust = false
    Δt = Second(Day(1))
    basetimestep = Second(Day(1))
    # end dummpy variables

    # TODO iterate over typedtable?
    for v in eachindex(sbm.nlayers)
        if do_lai
            cmax = sbm.sl[v] * sbm.lai[v] + sbm.swood[v]
            canopygapfraction = exp(-sbm.kext[v] * sbm.lai[v])
            ewet = (1.0 - exp(-sbm.kext[v] * sbm.lai[v])) * potevap
            e_r =
                precipitation > 0.0 ? min(0.25, ewet / max(0.0001, precipitation)) :
                0.0
        end

        potevap = potevap * sbm.et_reftopot[v]
        # should we include tempcor in SBM?
        # PotEvap = PotenEvap #??

        if Δt >= Hour(23)
            throughfall, interception, stemflow, canopystorage =
                rainfall_interception_gash(
                    cmax,
                    e_r,
                    canopygapfraction,
                    precipitation,
                    sbm.canopystorage[v],
                    maxevap = potevap,
                )
            pottrans_soil = max(0.0, potevap - interception) # now in mm
        else
            netinterception,
            throughfall,
            stemflow,
            leftover,
            interception,
            canopystorage = rainfall_interception_modrut(
                precipitation,
                potevap,
                sbm.canopystorage[v],
                canopygapfraction,
                cmax,
            )
            pottrans_soil = max(0.0, leftover)  # now in mm
            interception = netinterception
        end

        eff_precipitation = throughfall + stemflow

        if modelsnow
            tsoil = sbm.tsoil[v] + sbm.w_soil[v] * (temperature - sbm.tsoil[v])
            snow, snowwater, snowmelt, rainfallplusmelt, snowfall = snowpack_hbv(
                sbm.snow[v],
                sbm.snowwater[v],
                eff_precipitation,
                temperature,
                sbm.tti[v],
                sbm.tt[v],
                sbm.ttm[v],
                sbm.cfmax[v],
                sbm.whc[v],
            )
            if glacierfrac
                # Run Glacier module and add the snowpack on-top of it.
                # Estimate the fraction of snow turned into ice (HBV-light).
                # Estimate glacier melt.

                snow, snow2glacier, glacierstore, glaciermelt = glacier_hbv(
                    sbm.glacierfrac[v],
                    sbm.glacierstore[v],
                    sbm.snow[v],
                    temperature,
                    sbm.g_tt[v],
                    sbm.g_cfmax[v],
                    sbm.g_sifrac[v],
                    Δt,
                    basetimestep,
                )
                # Convert to mm per grid cell and add to snowmelt
                glaciermelt = glaciermelt * sbm.glacierfrac[v]
                rainfallplusmelt = rainfallplusmelt + glaciermelt

            end
        else
            rainfallplusmelt = eff_precipitation
        end

        avail_forinfilt = rainfallplusmelt + irsupply_mm
        ustoredepth = sum(@view sbm.ustorelayerdepth[v][1:sbm.nlayers[v]])
        uStorecapacity = sbm.soilwatercapacity[v] - sbm.satwaterdepth[v] - ustoredepth

        runoff_river = min(1.0, sbm.riverfrac[v]) * avail_forinfilt
        runoff_land = min(1.0, sbm.waterfrac[v]) * avail_forinfilt
        avail_forinfilt = max(avail_forinfilt - runoff_river - runoff_land, 0.0)

        rootingdepth = min(sbm.soilthickness[v] * 0.99, sbm.rootingdepth[v])

        ae_openw_r =
            min(wl_river * 1000.0 * sbm.riverfrac[v], sbm.riverfrac[v] * pottrans_soil)
        ae_openw_l =
            min(wl_land * 1000.0 * sbm.waterfrac[v], sbm.waterfrac[v] * pottrans_soil)

        restevap = pottrans_soil - ae_openw_r - ae_openw_l

        # evap available for soil evaporation and transpiration
        potsoilevap = restevap * canopygapfraction
        pottrans = restevap * (1.0 - canopygapfraction)

        # Calculate the initial capacity of the unsaturated store
        ustorecapacity = sbm.soilwatercapacity[v] - sbm.satwaterdepth[v] - ustoredepth

        # Calculate the infiltration flux into the soil column
        infiltsoilpath, infiltsoil, infiltpath, soilinf, pathinf, infiltexcess =
            infiltration(
                avail_forinfilt,
                sbm.pathfrac[v],
                sbm.cf_soil[v],
                tsoil,
                sbm.infiltcapsoil[v],
                sbm.infiltcappath[v],
                ustorecapacity,
                modelsnow,
                soilinfreduction,
            )


        usl, n_usl = set_layerthickness(sbm.zi[v], sbm.sumlayers[v], sbm.act_thickl[v])
        z = cumsum(usl)
        usld = copy(sbm.ustorelayerdepth[v])

        ast = 0.0
        soilevapunsat = 0.0
        if n_usl > 0
            # Using the surface infiltration rate, calculate the flow rate between the
            # different soil layers that contain unsaturated storage assuming gravity
            # based flow only, estimate the gravity based flux rate to the saturated zone
            # (ast) and the updated unsaturated storage for each soil layer.
            if transfermethod == 1 && sbm.maxlayers[v] == 1
                ustorelayerdepth = sbm.ustorelayerdepth[v][1] + infiltsoilpath
                kv_z = sbm.kvfrac[v][1] * sbm.kv₀[v] * exp(-sbm.f[v] * sbm.zi[v])
                ustorelayerdepth, ast = unsatzone_flow_sbm(
                    ustorelayerdepth,
                    sbm.soilwatercapacity[v],
                    sbm.satwaterdepth[v],
                    kv_z,
                    usl[1],
                    sbm.θₛ[v],
                    sbm.θᵣ[v],
                )
                usld = setindex(usld, ustorelayerdepth, m)
            else
                for m = 1:n_usl
                    l_sat = usl[m] * (sbm.θₛ[v] - sbm.θᵣ[v])
                    kv_z = sbm.kvfrac[v][m] * sbm.kv₀[v] * exp(-sbm.f[v] * z[m])
                    ustorelayerdepth =
                        m == 1 ? sbm.ustorelayerdepth[v][m] + infiltsoilpath :
                        sbm.ustorelayerdepth[v][m] + ast
                    ustorelayerdepth, ast = unsatzone_flow_layer(
                        ustorelayerdepth,
                        kv_z,
                        l_sat,
                        sbm.c[v][m],
                    )
                    usld = setindex(usld, ustorelayerdepth, m)
                end
            end

            # then evapotranspiration from layers
            # Calculate saturation deficity
            saturationdeficit = sbm.soilwatercapacity[v] - sbm.satwaterdepth[v]

            # First calculate the evaporation of unsaturated storage into the
            # atmosphere from the upper layer.
            if sbm.maxlayers[v] == 1
                soilevapunsat =
                    potsoilevap *
                    min(1.0, saturationdeficit / sbm.soilwatercapacity[v])
            else
                # In case only the most upper soil layer contains unsaturated storage
                if n_usl == 1
                    # Check if groundwater level lies below the surface
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (sbm.zi[v] * (sbm.θₛ[v] - sbm.θᵣ[v])))
                else
                    # In case first layer contains no saturated storage
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (usl[1] * ((sbm.θₛ[v] - sbm.θᵣ[v]))))
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

        if sbm.maxlayers[v] == 1
            soilevapsat = 0.0
        else
            # this check is an improvement compared to Python (only checked for n_usl == 1)
            if n_usl == 0 || n_usl == 1
                soilevapsat =
                    potsoilevap *
                    min(1.0, (sbm.act_thickl[v][1] - sbm.zi[v]) / sbm.act_thickl[v][1])
                soilevapsat = min(
                    soilevapsat,
                    (sbm.act_thickl[v][1] - sbm.zi[v]) * (sbm.θₛ[v] - sbm.θᵣ[v]),
                )
            else
                soilevapsat = 0.0
            end
        end
        soilevap = soilevapunsat + soilevapsat
        satwaterdepth = sbm.satwaterdepth[v] - soilevapsat

        # transpiration from saturated store
        wetroots = scurve(sbm.zi[v], a = rootingdepth, c = sbm.rootdistpar[v])
        actevapsat = min(pottrans * wetroots, satwaterdepth)
        satwaterdepth = satwaterdepth - actevapsat
        restpottrans = pottrans - actevapsat

        # actual transpiration from ustore
        actevapustore = 0.0
        for k = 1:n_usl
            ustorelayerdepth, actevapustore, restpottrans = acttransp_unsat_sbm(
                rootingdepth,
                usld[k],
                sbm.sumlayers[v][k],
                restpottrans,
                actevapustore,
                sbm.c[v][k],
                usl[k],
                sbm.θₛ[v],
                sbm.θᵣ[v],
                sbm.hb[v],
                ust,
            )
            usld = setindex(usld, ustorelayerdepth, k)
        end

        # check soil moisture balance per layer
        du = 0.0
        for k = n_usl:-1:1
            du = max(0.0, usld[k] - usl[k] * (sbm.θₛ[v] - sbm.θᵣ[v]))
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
            ksat = sbm.kvfrac[v][n_usl] * sbm.kv₀[v] * exp(-sbm.f[v] * sbm.zi[v])
            ustorecapacity =
                sbm.soilwatercapacity[v] -
                sbm.satwaterdepth[v] - sum(@view usld[1:sbm.nlayers[v]])
            maxcapflux = max(
                0.0,
                min(ksat, actevapustore, ustorecapacity, sbm.satwaterdepth[v]),
            )

            if sbm.zi[v] > rootingdepth
                capfluxscale =
                    sbm.capscale[v] / (sbm.capscale[v] + sbm.zi[v] - rootingdepth) *
                    Float64(Δt.value) / Float64(basetimestep.value)
            else
                capfluxscale = 0.0
            end
            capflux = maxcapflux * capfluxscale

            netcapflux = capflux
            for k = n_usl:-1:1
                toadd =
                    min(netcapflux, max(usl[k] * (sbm.θₛ[v] - sbm.θᵣ[v]) - usld[k], 0.0))
                usld = setindex(usld, usld[k] + toadd, k)
                netcapflux = netcapflux - toadd
                actcapflux = actcapflux + toadd
            end
        end
        deepksat = sbm.kv₀[v] * exp(-sbm.f[v] * sbm.soilthickness[v])
        deeptransfer = min(sbm.satwaterdepth[v], deepksat)
        actleakage = max(0.0, min(sbm.maxleakage[v], deeptransfer))

        # recharge (mm) for saturated zone
        recharge = (transfer - actcapflux - actleakage - actevapsat - soilevapsat)
        transpiration = actevapsat + actevapustore
        actevap = soilevap + transpiration + ae_openw_r + ae_openw_l

        # update the outputs and states
        sbm.cmax[v] = cmax
        sbm.canopygapfraction[v] = canopygapfraction
        sbm.canopystorage[v] = canopystorage
        sbm.snow[v] = snow
        sbm.snowwater[v] = snowwater
        sbm.tsoil[v] = tsoil
        sbm.actinfilt[v] = actinfilt
        sbm.recharge[v] = recharge
        sbm.transpiration[v] = transpiration
        sbm.soilevap[v] = soilevap
        sbm.interception[v] = interception
        sbm.ae_openw_r[v] = ae_openw_r
        sbm.ae_openw_l[v] = ae_openw_l
        sbm.actevapsat[v] = actevapsat
        sbm.actevap[v] = actevap
        sbm.ustorelayerdepth[v] = usld
        sbm.transfer[v] = transfer
        sbm.actinfiltsoil[v] = actinfiltsoil
        sbm.actinfiltpath[v] = actinfiltpath
        sbm.excesswater[v] = excesswater
        sbm.excesswatersoil[v] = excesswatersoil
        sbm.excesswaterpath[v] = excesswaterpath
    end
end

function update_after_lateralflow(sbm, zi, exfiltsatwater)

    for v in eachindex(sbm.nlayers)
        usl, n_usl = set_layerthickness(zi[v], sbm.sumlayers[v], sbm.act_thickl[v])
        # exfiltration from ustore
        usld = copy(sbm.ustorelayerdepth[v])
        exfiltustore = 0.0
        for k = n_usl:-1:1
            exfiltustore = max(0, usld[k] - usl[k] * (sbm.θₛ[v] - sbm.θᵣ[v]))
            usld = setindex(usld, usld[k] - exfiltustore, k)
            if k > 1
                usld = setindex(usld, usld[k - 1] + exfiltustore, k - 1)
            end
        end

        ustoredepth = sum(@view usld[1:n_usl])

        runoff = max(
            exfiltustore +
            exfiltsatwater[v] +
            sbm.excesswater[v] +
            sbm.infiltexcess[v] +
            sbm.ae_openw_l[v],
            0.0,
        )

        # volumetric water content per soil layer and root zone
        vwc = copy(sbm.vwc[v])
        vwc_perc = copy(sbm.vwc_perc[v])
        for k = 1:sbm.nlayers[v]
            if k <= n_usl
                vwc = setindex(
                    vwc,
                    (usld[k] + (sbm.act_thickl[v][k] - usl[k]) * (sbm.θₛ[v] - sbm.θᵣ[v])) / usl[k] + sbm.θᵣ[v],
                    k,
                )
            else
                vwc = setindex(vwc, sbm.θₛ[v], k)
            end
            vwc_perc = setindex(vwc_perc, (vwc[k] / sbm.θₛ[v]) * 100.0, k)
        end

        rootstore_unsat = 0
        for k = 1:n_usl
            rootstore_unsat =
                rootstore_unsat +
                (max(0.0, sbm.rootingdepth[v] - sbm.sumlayers[v][k]) / usl[k]) * usld[k]
        end

        rootstore_sat = max(0.0, sbm.rootingdepth[v] - zi[v]) * (sbm.θₛ[v] - sbm.θᵣ[v])
        rootstore = rootstore_sat + rootstore_unsat
        vwc_root = rootstore / sbm.rootingdepth[v] + sbm.θᵣ[v]
        vwc_percroot = (vwc_root / sbm.θₛ[v]) * 100.0

        satwaterdepth = (sbm.soilthickness[v] - zi[v]) * (sbm.θₛ[v] - sbm.θᵣ[v])

        # update the outputs and states
        sbm.ustorelayerdepth[v] = usld
        sbm.ustoredepth[v] = ustoredepth
        sbm.satwaterdepth[v] = satwaterdepth
        sbm.exfiltsatwater[v] = exfiltsatwater[v]
        sbm.exfiltustore[v] = exfiltustore
        sbm.runoff[v] = runoff
        sbm.vwc[v] = vwc
        sbm.vwc_perc[v] = vwc_perc
        sbm.rootstore[v] = rootstore
        sbm.vwc_root[v] = vwc_root
        sbm.vwc_percroot[v] = vwc_percroot
        sbm.zi[v] = zi[v]
    end
end
