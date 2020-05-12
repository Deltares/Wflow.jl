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
    #TODO: (warm) states read from netcdf file or cold state (reinit=1, setting in ini file)

end

function update_before_lateralflow(sbm::SBM)

    #start dummy variables (should be generated from model reader and from Config.jl TOML)
    do_lai = true
    glacierfrac = false
    modelsnow = true
    soilinfreduction = false
    transfermethod = 0
    potevap = 4.0
    precipitation = 3.0
    temperature = 10.0
    wl_land = 0.0 #from kinematic wave land
    wl_river = 0.10 #from kinematic river
    irsupply_mm = 0.0
    ust = false
    Δt = Second(Day(1))
    basetimestep = Second(Day(1))
    #end dummpy variables

    if do_lai
        cmax = sbm.sl * sbm.lai + sbm.swood
        canopygapfraction = exp(-sbm.kext * sbm.lai)
        ewet = (1.0 - exp(-sbm.kext * sbm.lai)) * potevap
        e_r =
            precipitation > 0.0 ? min(0.25, ewet / max(0.0001, precipitation)) :
            0.0
    end

    potevap = potevap * sbm.et_reftopot
    # should we include tempcor in SBM?
    # PotEvap = PotenEvap #??

    if Δt >= Hour(23)
        throughfall, interception, stemflow, canopystorage =
            rainfall_interception_gash(
                cmax,
                e_r,
                canopygapfraction,
                precipitation,
                sbm.canopystorage,
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
            sbm.canopystorage,
            canopygapfraction,
            cmax,
        )
        pottrans_soil = max(0.0, leftover)  # now in mm
        interception = netinterception
    end

    eff_precipitation = throughfall + stemflow

    if modelsnow
        tsoil = sbm.tsoil + sbm.w_soil * (temperature - sbm.tsoil)
        snow, snowwater, snowmelt, rainfallplusmelt, snowfall = snowpack_hbv(
            sbm.snow,
            sbm.snowwater,
            eff_precipitation,
            temperature,
            sbm.tti,
            sbm.tt,
            sbm.ttm,
            sbm.cfmax,
            sbm.whc,
        )
        if glacierfrac
            # Run Glacier module and add the snowpack on-top of it.
            # Estimate the fraction of snow turned into ice (HBV-light).
            # Estimate glacier melt.

            snow, snow2glacier, glacierstore, glaciermelt = glacier_hbv(
                sbm.glacierfrac,
                sbm.glacierstore,
                sbm.snow,
                temperature,
                sbm.g_tt,
                sbm.g_cfmax,
                sbm.g_sifrac,
                Δt,
                basetimestep,
            )
            # Convert to mm per grid cell and add to snowmelt
            glaciermelt = glaciermelt * sbm.glacierfrac
            rainfallplusmelt = rainfallplusmelt + glaciermelt

        end
    else
        rainfallplusmelt = eff_precipitation
    end

    avail_forinfilt = rainfallplusmelt + irsupply_mm
    ustoredepth = sum(sbm.ustorelayerdepth[1:sbm.nlayers])
    uStorecapacity = sbm.soilwatercapacity - sbm.satwaterdepth - ustoredepth

    runoff_river = min(1.0, sbm.riverfrac) * avail_forinfilt
    runoff_land = min(1.0, sbm.waterfrac) * avail_forinfilt
    avail_forinfilt = max(avail_forinfilt - runoff_river - runoff_land, 0.0)

    rootingdepth = min(sbm.soilthickness * 0.99, sbm.rootingdepth)

    ae_openw_r =
        min(wl_river * 1000.0 * sbm.riverfrac, sbm.riverfrac * pottrans_soil)
    ae_openw_l =
        min(wl_land * 1000.0 * sbm.waterfrac, sbm.waterfrac * pottrans_soil)

    restevap = pottrans_soil - ae_openw_r - ae_openw_l

    # evap available for soil evaporation and transpiration
    potsoilevap = restevap * canopygapfraction
    pottrans = restevap * (1.0 - canopygapfraction)

    # Calculate the initial capacity of the unsaturated store
    ustorecapacity = sbm.soilwatercapacity - sbm.satwaterdepth - ustoredepth

    # Calculate the infiltration flux into the soil column
    infiltsoilpath, infiltsoil, infiltpath, soilinf, pathinf, infiltexcess =
        infiltration(
            avail_forinfilt,
            sbm.pathfrac,
            sbm.cf_soil,
            tsoil,
            sbm.infiltcapsoil,
            sbm.infiltcappath,
            ustorecapacity,
            modelsnow,
            soilinfreduction,
        )

    usl, n_usl = set_layerthickness(sbm.zi, sbm.sumlayers, sbm.act_thickl)
    z = cumsum(usl)

    usld = copy(sbm.ustorelayerdepth)
    # Using the surface infiltration rate, calculate the flow rate between the
    # different soil layers that contain unsaturated storage assuming gravity
    # based flow only, estimate the gravity based flux rate to the saturated zone
    # (ast) and the updated unsaturated storage for each soil layer.
    if transfermethod == 1 && sbm.maxlayers == 1
        ustorelayerdepth = sbm.ustorelayerdepth[1] + infiltsoilpath
        kv_z = sbm.kvfrac[1] * sbm.kv₀ * exp(-sbm.f * sbm.zi)
        ustorelayerdepth, ast = unsatzone_flow_sbm(
            ustorelayerdepth,
            sbm.soilwatercapacity,
            sbm.satwaterdepth,
            kv_z,
            usl[1],
            sbm.θₛ,
            sbm.θᵣ,
        )
        usld = setindex(usld, ustorelayerdepth, m)
    else
        for m = 1:n_usl
            l_sat = usl[m] * (sbm.θₛ - sbm.θᵣ)
            kv_z = sbm.kvfrac[m] * sbm.kv₀ * exp(-sbm.f * z[m])
            ustorelayerdepth =
                m == 1 ? sbm.ustorelayerdepth[m] + infiltsoilpath :
                sbm.ustorelayerdepth[m] + ast
            ustorelayerdepth, ast =
                unsatzone_flow_layer(ustorelayerdepth, kv_z, l_sat, sbm.c[m])
            usld = setindex(usld, ustorelayerdepth, m)
        end
    end

    transfer = ast

    # then evapotranspiration from layers
    # Calculate saturation deficity
    saturationdeficit = sbm.soilwatercapacity - sbm.satwaterdepth

    # First calculate the evaporation of unsaturated storage into the
    # atmosphere from the upper layer.
    if sbm.maxlayers == 1
        soilevapunsat =
            potsoilevap * min(1.0, saturationdeficit / sbm.soilwatercapacity)
    else
        #In case only the most upper soil layer contains unsaturated storage
        if n_usl == 1
            # Check if groundwater level lies below the surface
            soilevapunsat = sbm.zi > 0.0 ?
                potsoilevap * min(1.0, usld[1] / (sbm.zi * (sbm.θₛ - sbm.θᵣ))) :
                0.0
        else
            # In case first layer contains no saturated storage
            soilevapunsat =
                potsoilevap * min(1.0, usld[1] / (usl[1] * ((sbm.θₛ - sbm.θᵣ))))
        end
    end

    # Ensure that the unsaturated evaporation rate does not exceed the
    # available unsaturated moisture
    soilevapunsat = min(soilevapunsat, usld[1])
    # Update the additional atmospheric demand
    potsoilevap = potsoilevap - soilevapunsat
    usld = setindex(usld, usld[1] - soilevapunsat, 1)

    if sbm.maxlayers == 1
        soilevapsat = 0.0
    else
        # this check is an improvement compared to Python (only checked for n_usl == 1)
        if n_usl == 0 || n_usl == 1
            soilevapsat = potsoilevap * min(1.0, (sbm.act_thickl[1] - sbm.zi) / sbm.act_thickl[1])
            soilevapsat =
                min(soilevapsat, (sbm.act_thickl[1] - sbm.zi) * (sbm.θₛ - sbm.θᵣ))
        else
            soilevapsat = 0.0
        end
    end
    soilevap = soilevapunsat + soilevapsat
    satwaterdepth = sbm.satwaterdepth - soilevapsat

    # transpiration from saturated store
    wetroots = scurve(sbm.zi, a = rootingdepth, c = sbm.rootdistpar)
    actevapsat = min(pottrans * wetroots, satwaterdepth)
    satwaterdepth = satwaterdepth - actevapsat
    restpottrans = pottrans - actevapsat

    # actual transpiration from ustore
    actevapustore = 0.0
    for k = 1:n_usl
        ustorelayerdepth, actevapustore, restpottrans = acttransp_unsat_sbm(
            rootingdepth,
            usld[k],
            sbm.sumlayers[k],
            restpottrans,
            actevapustore,
            sbm.c[k],
            usl[k],
            sbm.θₛ,
            sbm.θᵣ,
            sbm.hb,
            ust,
        )
        usld = setindex(usld, ustorelayerdepth, k)
    end

    #check soil moisture balance per layer
    du = 0.0
    for k = n_usl:-1:1
        du = max(0.0, usld[k] - usl[k] * (sbm.θₛ - sbm.θᵣ))
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

    ksat = sbm.kvfrac[n_usl] * sbm.kv₀ * exp(-sbm.f * sbm.zi)
    ustorecapacity = sbm.soilwatercapacity - sbm.satwaterdepth - sum(usld[1:sbm.nlayers])
    maxcapflux =
        max(0.0, min(ksat, actevapustore, ustorecapacity, sbm.satwaterdepth))

    if sbm.zi > rootingdepth
        capfluxscale =
            sbm.capscale / (sbm.capscale + sbm.zi - rootingdepth) *
            Float64(Δt.value) / Float64(basetimestep.value)
    else
        capfluxscale = 0.0
    end
    capflux = maxcapflux * capfluxscale

    netcapflux = capflux
    actcapflux = 0.0
    for k = n_usl:-1:1
        toadd = min(netcapflux, max(usl[k] * (sbm.θₛ - sbm.θᵣ) - usld[k], 0.0))
        usld = setindex(usld, usld[k] + toadd, k)
        netcapflux = netcapflux - toadd
        actcapflux = actcapflux + toadd
    end
    deepksat = sbm.kv₀ * exp(-sbm.f * sbm.soilthickness)
    deeptransfer = min(sbm.satwaterdepth, deepksat)
    actleakage = max(0.0, min(sbm.maxleakage, deeptransfer))

    # recharge (mm) for saturated zone
    recharge = (transfer - actcapflux - actleakage - actevapsat - soilevapsat)
    transpiration = actevapsat + actevapustore
    actevap = soilevap + transpiration + ae_openw_r + ae_openw_l

    return setproperties(sbm, (
        cmax = cmax,
        canopygapfraction = canopygapfraction,
        canopystorage = canopystorage,
        snow = snow,
        snowwater = snowwater,
        tsoil = tsoil,
        actinfilt = actinfilt,
        recharge = recharge,
        transpiration = transpiration,
        soilevap = soilevap,
        interception = interception,
        ae_openw_r = ae_openw_r,
        ae_openw_l = ae_openw_l,
        actevapsat = actevapsat,
        actevap = actevap,
        ustorelayerdepth = usld,
        transfer = transfer,
        actinfiltsoil = actinfiltsoil,
        actinfiltpath = actinfiltpath,
        excesswater = excesswater,
        excesswatersoil = excesswatersoil,
        excesswaterpath = excesswaterpath))
end

function update_after_lateralflow(sbm::SBM, zi, exfiltsatwater)

    usl, n_usl = set_layerthickness(zi, sbm.sumlayers, sbm.act_thickl)
    # exfiltration from ustore
    usld = copy(sbm.ustorelayerdepth)
    exfiltustore = 0.0
    for k = n_usl:-1:1
        exfiltustore = max(0, usld[k] - usl[k] * (sbm.θₛ - sbm.θᵣ))
        usld = setindex(usld, usld[k] - exfiltustore, k)
        if k > 1
            usld = setindex(usld, usld[k-1] + exfiltustore, k - 1)
        end
    end

    ustoredepth = sum(usld[1:n_usl])

    runoff = max(
        exfiltustore +
        exfiltsatwater +
        sbm.excesswater +
        sbm.infiltexcess +
        sbm.ae_openw_l,
        0.0,
    )

    # volumetric water content per soil layer and root zone
    vwc = copy(sbm.vwc)
    vwc_perc = copy(sbm.vwc_perc)
    for k = 1:sbm.nlayers
        if k <= n_usl
            vwc = setindex(
                vwc,
                (usld[k] + (sbm.act_thickl[k] - usl[k]) * (sbm.θₛ - sbm.θᵣ)) / usl[k] + sbm.θᵣ,
                k,
            )
        else
            vwc = setindex(vwc, sbm.θₛ, k)
        end
        vwc_perc = setindex(vwc_perc, (vwc[k] / sbm.θₛ) * 100.0, k)
    end

    rootstore_unsat = 0
    for k = 1:n_usl
        rootstore_unsat =
            rootstore_unsat +
            (max(0.0, sbm.rootingdepth - sbm.sumlayers[k]) / usl[k]) * usld[k]
    end

    rootstore_sat = max(0.0, sbm.rootingdepth - zi) * (sbm.θₛ - sbm.θᵣ)
    rootstore = rootstore_sat + rootstore_unsat
    vwc_root = rootstore / sbm.rootingdepth + sbm.θᵣ
    vwc_percroot = (vwc_root / sbm.θₛ) * 100.0

    satwaterdepth =  (sbm.soilthickness - zi) * (sbm.θₛ - sbm.θᵣ)

    return setproperties(sbm, (
        ustorelayerdepth = usld,
        ustoredepth = ustoredepth,
        satwaterdepth = satwaterdepth,
        exfiltsatwater = exfiltsatwater,
        exfiltustore = exfiltustore,
        runoff = runoff,
        vwc = vwc,
        vwc_perc = vwc_perc,
        rootstore = rootstore,
        vwc_root = vwc_root,
        vwc_percroot = vwc_percroot,
        zi = zi))
end
