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

function update_before_lateralflow(t)

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

    for r in eachindex(t.nlayers)
        if do_lai
            cmax = t.sl[r] * t.lai[r] + t.swood[r]
            canopygapfraction = exp(-t.kext[r] * t.lai[r])
            ewet = (1.0 - exp(-t.kext[r] * t.lai[r])) * potevap
            e_r =
                precipitation > 0.0 ? min(0.25, ewet / max(0.0001, precipitation)) :
                0.0
        end

        potevap = potevap * t.et_reftopot[r]
        # should we include tempcor in SBM?
        # PotEvap = PotenEvap #??

        if Δt >= Hour(23)
            throughfall, interception, stemflow, canopystorage =
                rainfall_interception_gash(
                    cmax,
                    e_r,
                    canopygapfraction,
                    precipitation,
                    t.canopystorage[r],
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
                t.canopystorage[r],
                canopygapfraction,
                cmax,
            )
            pottrans_soil = max(0.0, leftover)  # now in mm
            interception = netinterception
        end

        eff_precipitation = throughfall + stemflow

        if modelsnow
            tsoil = t.tsoil[r] + t.w_soil[r] * (temperature - t.tsoil[r])
            snow, snowwater, snowmelt, rainfallplusmelt, snowfall = snowpack_hbv(
                t.snow[r],
                t.snowwater[r],
                eff_precipitation,
                temperature,
                t.tti[r],
                t.tt[r],
                t.ttm[r],
                t.cfmax[r],
                t.whc[r],
            )
            if glacierfrac
                # Run Glacier module and add the snowpack on-top of it.
                # Estimate the fraction of snow turned into ice (HBV-light).
                # Estimate glacier melt.

                snow, snow2glacier, glacierstore, glaciermelt = glacier_hbv(
                    t.glacierfrac[r],
                    t.glacierstore[r],
                    t.snow[r],
                    temperature,
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
            rainfallplusmelt = eff_precipitation
        end

        avail_forinfilt = rainfallplusmelt + irsupply_mm
        ustoredepth = sum(@view t.ustorelayerdepth[r][1:t.nlayers[r]])
        uStorecapacity = t.soilwatercapacity[r] - t.satwaterdepth[r] - ustoredepth

        runoff_river = min(1.0, t.riverfrac[r]) * avail_forinfilt
        runoff_land = min(1.0, t.waterfrac[r]) * avail_forinfilt
        avail_forinfilt = max(avail_forinfilt - runoff_river - runoff_land, 0.0)

        rootingdepth = min(t.soilthickness[r] * 0.99, t.rootingdepth[r])

        ae_openw_r =
            min(wl_river * 1000.0 * t.riverfrac[r], t.riverfrac[r] * pottrans_soil)
        ae_openw_l =
            min(wl_land * 1000.0 * t.waterfrac[r], t.waterfrac[r] * pottrans_soil)

        restevap = pottrans_soil - ae_openw_r - ae_openw_l

        # evap available for soil evaporation and transpiration
        potsoilevap = restevap * canopygapfraction
        pottrans = restevap * (1.0 - canopygapfraction)

        # Calculate the initial capacity of the unsaturated store
        ustorecapacity = t.soilwatercapacity[r] - t.satwaterdepth[r] - ustoredepth

        # Calculate the infiltration flux into the soil column
        infiltsoilpath, infiltsoil, infiltpath, soilinf, pathinf, infiltexcess =
            infiltration(
                avail_forinfilt,
                t.pathfrac[r],
                t.cf_soil[r],
                tsoil,
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
                    ustorelayerdepth =
                        m == 1 ? t.ustorelayerdepth[r][m] + infiltsoilpath :
                        t.ustorelayerdepth[r][m] + ast
                    ustorelayerdepth, ast = unsatzone_flow_layer(
                        ustorelayerdepth,
                        kv_z,
                        l_sat,
                        t.c[r][m],
                    )
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
                    potsoilevap *
                    min(1.0, saturationdeficit / t.soilwatercapacity[r])
            else
                # In case only the most upper soil layer contains unsaturated storage
                if n_usl == 1
                    # Check if groundwater level lies below the surface
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (t.zi[r] * (t.θₛ[r] - t.θᵣ[r])))
                else
                    # In case first layer contains no saturated storage
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (usl[1] * ((t.θₛ[r] - t.θᵣ[r]))))
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
                soilevapsat = min(
                    soilevapsat,
                    (t.act_thickl[r][1] - t.zi[r]) * (t.θₛ[r] - t.θᵣ[r]),
                )
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
            ksat = t.kvfrac[r][n_usl] * t.kv₀[r] * exp(-t.f[r] * t.zi[r])
            ustorecapacity =
                t.soilwatercapacity[r] -
                t.satwaterdepth[r] - sum(@view usld[1:t.nlayers[r]])
            maxcapflux = max(
                0.0,
                min(ksat, actevapustore, ustorecapacity, t.satwaterdepth[r]),
            )

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
                toadd =
                    min(netcapflux, max(usl[k] * (t.θₛ[r] - t.θᵣ[r]) - usld[k], 0.0))
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
        t.cmax[r] = cmax
        t.canopygapfraction[r] = canopygapfraction
        t.canopystorage[r] = canopystorage
        t.snow[r] = snow
        t.snowwater[r] = snowwater
        t.tsoil[r] = tsoil
        t.actinfilt[r] = actinfilt
        t.recharge[r] = recharge
        t.transpiration[r] = transpiration
        t.soilevap[r] = soilevap
        t.interception[r] = interception
        t.ae_openw_r[r] = ae_openw_r
        t.ae_openw_l[r] = ae_openw_l
        t.actevapsat[r] = actevapsat
        t.actevap[r] = actevap
        t.ustorelayerdepth[r] = usld
        t.transfer[r] = transfer
        t.actinfiltsoil[r] = actinfiltsoil
        t.actinfiltpath[r] = actinfiltpath
        t.excesswater[r] = excesswater
        t.excesswatersoil[r] = excesswatersoil
        t.excesswaterpath[r] = excesswaterpath
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
                usld = setindex(usld, usld[k - 1] + exfiltustore, k - 1)
            end
        end

        ustoredepth = sum(@view usld[1:n_usl])

        runoff = max(
            exfiltustore +
            exfiltsatwater[r] +
            t.excesswater[r] +
            t.infiltexcess[r] +
            t.ae_openw_l[r],
            0.0,
        )

        # volumetric water content per soil layer and root zone
        vwc = copy(t.vwc[r])
        vwc_perc = copy(t.vwc_perc[r])
        for k = 1:t.nlayers[r]
            if k <= n_usl
                vwc = setindex(
                    vwc,
                    (usld[k] + (t.act_thickl[r][k] - usl[k]) * (t.θₛ[r] - t.θᵣ[r])) / usl[k] + t.θᵣ[r],
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
