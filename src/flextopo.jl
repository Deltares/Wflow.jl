@get_units @with_kw struct FLEXTOPO{T,N}
    # Model time step [s]
    Δt::T | "s"
    # Number of classes
    nclass::Int | "-"
    # Number of cells
    n::Int | "-"
    #dictionary with all possible functions for each store
    dic_function::Dict
    #current class
    kclass::Vector{Int64}
    classes::Vector{String}
    select_snow::Vector{String}
    select_interception::Vector{String}
    select_hortonponding::Vector{String}
    select_hortonrunoff::Vector{String}
    select_rootzone::Vector{String}
    select_fast::Vector{String}
    select_slow::Vector{String}

    #fraction of each class
    hrufrac::Vector{SVector{N,T}} | "-"
    ## PARAMETERS
    ##SNOW
    # Correction factor for precipitation [-]
    pcorr::Vector{T} | "-"
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹]
    cfmax::Vector{T} | "mm ᵒC-1 Δt-1"
    # Threshold temperature for snowfall [ᵒC]
    tt::Vector{T} | "ᵒC"
    # Threshold temperature interval length [ᵒC]
    tti::Vector{T} | "ᵒC"
    # Threshold temperature for snowmelt [ᵒC]
    ttm::Vector{T} | "ᵒC"
    # Water holding capacity as fraction of current snow pack [-]
    whc::Vector{T} | "-"
    # Refreezing efficiency constant in refreezing of freewater in snow [-]
    cfr::Vector{T} | "-"
    # Correction factor for precipitation [-]
    rfcf::Vector{T} | "-"
    # Correction factor for snowfall [-]
    sfcf::Vector{T} | "-"
    ## GLACIER
    # Threshold temperature for snowfall above glacier [ᵒC]
    g_tt::Vector{T} | "ᵒC"
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹] for glacier
    g_cfmax::Vector{T} | "mm ᵒC-1 Δt-1"
    # Fraction of the snowpack on top of the glacier converted into ice [Δt⁻¹]
    g_sifrac::Vector{T} | "Δt-1"
    # Water within the glacier [mm]
    glacierstore::Vector{T} | "mm"
    # Fraction covered by a glacier [-]
    glacierfrac::Vector{T} | "-"
    ##INTERCEPTION
    # Maximum interception storage (in forested and non-forested areas) [mm]
    imax::Vector{SVector{N,T}} | "mm"
    # Evap correction [-]
    ecorr::Vector{T} | "-"
    ##HORTON
    # Maximum storage capacity in the hortonian ponding storage [mm]
    shmax::Vector{SVector{N,T}} | "mm"
    #recession coefficient of the hortonian runoff storage [Δt-1]
    khf::Vector{SVector{N,T}} | "Δt-1"
    #maximum modelled accumulated frost resulting in shmin [ᵒC Δt]
    facc0::Vector{SVector{N,T}} | "ᵒC"
    #minimum modelled accumulated frost resulting in shmax [ᵒC Δt]
    facc1::Vector{SVector{N,T}} | "ᵒC"
    #exponent for the decline of infiltration capacity [-]
    fdec::Vector{SVector{N,T}} | "-"
    #maximum infiltration capacity from horton ponding [mm Δt-1]
    fmax::Vector{SVector{N,T}} | "mm Δt-1"
    #minimum storage capacity in horton ponding (relative to shmax) [-]
    shmin::Vector{SVector{N,T}} | "-"
    #melt coefficient for melt of frozen topsoil [-]
    kmf::Vector{SVector{N,T}} | "-"
    ##ROOTZONE
    # maximum root-zone storage capacity [mm]
    srmax::Vector{SVector{N,T}} | "mm"
    # Fraction of root zone storage below which actual evaporation is potential evaporation [-]
    lp::Vector{SVector{N,T}} | "-"
    # Exponent in soil runoff generation equation [-]
    beta::Vector{SVector{N,T}} | "-"
    # maximum percolation rate [mm Δt⁻¹]
    perc::Vector{SVector{N,T}} | "mm Δt-1"
    # maximum capillary rise rate [mm Δt⁻¹]
    cap::Vector{SVector{N,T}} | "mm Δt-1"
    #FAST
    # Exponent for non linear recession [-]
    alfa::Vector{SVector{N,T}} | "-"
    #recession coefficient of fast storage [Δt-1]
    kf::Vector{SVector{N,T}} | "Δt-1"
    # fraction of qrootzone to slowstorage (1-ds to faststorage)
    ds::Vector{SVector{N,T}} | "-"
    # SLOW
    #recession coefficient of slow storage [Δt-1]
    ks::Vector{T} | "Δt-1"

    ## STATES
    ##SNOW
    # Snow water equivalent [mm]
    snow::Vector{T} | "mm"
    # Liquid water content in the snow pack [mm]
    snowwater::Vector{T} | "mm"
    # Interception storage [mm]
    interceptionstorage::Vector{SVector{N,T}} | "mm"
    # Storage in the hortonian ponding reservoir [mm]
    hortonpondingstorage::Vector{SVector{N,T}} | "mm"
    # Storage in the hortonian runoff generation [mm]
    hortonrunoffstorage::Vector{SVector{N,T}} | "mm"
    # Storage in the root-zone [mm]
    rootzonestorage::Vector{SVector{N,T}} | "mm"
    # Storage in the root-zone relative to maximum root-zone storage capacity [mm]
    srootzone_over_srmax::Vector{SVector{N,T}} | "mm"
    # Storage in the fast store [mm]
    faststorage::Vector{SVector{N,T}} | "mm"
    # Storage in the slow reservoir (for qcapillary calc) [mm]
    slowstorage::Vector{T} | "mm"
    #states previous time step to calc water balance [mm]
    states_::Vector{SVector{N,T}} | "mm"
    #states previous time step to calc water balance combined based on perc class. [mm]
    states_m::Vector{T} | "mm"
    #states averaged over classes
    interceptionstorage_m::Vector{T} | "mm"
    hortonpondingstorage_m::Vector{T} | "mm"
    hortonrunoffstorage_m::Vector{T} | "mm"
    srootzone_m::Vector{T} | "mm"
    faststorage_m::Vector{T} | "mm"
    srootzone_over_srmax_m::Vector{T} | "mm"

    ## FLUXES
    #SNOW
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{T} | "mm Δt-1"
    # Temperature [ᵒC]
    temperature::Vector{T} | "ᵒC"
    # Potential evapotranspiration [mm Δt⁻¹]
    potential_evaporation::Vector{T} | "mm Δt-1"
    # Potential evapotranspiration corrected [mm Δt⁻¹]
    epotcorr::Vector{T} | "mm Δt-1"
    # Precipitation corrected [mm Δt⁻¹]
    precipcorr::Vector{T} | "mm Δt-1"
    # Snow melt + precipitation as rainfall [mm]
    rainfallplusmelt::Vector{T} | "mm Δt-1"
    # Snowfall [mm]
    snowfall::Vector{T} | "mm Δt-1"
    # Snowmelt [mm]
    snowmelt::Vector{T} | "mm Δt-1"
    #INTERCEPTION
    # Potential soil evaporation [mm Δt⁻¹]
    potsoilevap::Vector{SVector{N,T}} | "mm Δt-1"
    # Evaporation from interception storage [mm Δt⁻¹]
    intevap::Vector{SVector{N,T}} | "mm Δt-1"
    #effective precipitation [mm Δt⁻¹]
    precipeffective::Vector{SVector{N,T}} | "mm Δt-1"
    # Evaporation from interception sum classes [mm Δt⁻¹]
    intevap_m::Vector{T} | "mm Δt-1"
    #HORTONPONDING
    # Evaporation from the hortonion ponding storage [-]
    hortonevap::Vector{SVector{N,T}} | "mm Δt-1"
    # Flux from the hortonian ponding storage to the hortonian runoff storage [mm Δt⁻¹]
    qhortonpond::Vector{SVector{N,T}} | "mm Δt-1"
    # Flux from the hortonian ponding storage to the root zone storage [mm Δt⁻¹]
    qhortonrootzone::Vector{SVector{N,T}} | "mm Δt-1"
    # modeled accumulated frost [ᵒC Δt]
    facc::Vector{SVector{N,T}} | "ᵒC Δt"
    # Evaporation from the hortonian sum classes [mm Δt⁻¹]
    hortonevap_m::Vector{T} | "mm Δt-1"
    #HORTONRUNOFF
    # Flux from the hortonian runoff storage [mm Δt⁻¹]
    qhortonrun::Vector{SVector{N,T}} | "mm Δt-1"
    #ROOTZONE
    # Evaporation from the root-zone storage [mm Δt⁻¹]
    rootevap::Vector{SVector{N,T}} | "mm Δt-1"
    # Flux from the root-zone storage [mm Δt⁻¹]
    qrootzone::Vector{SVector{N,T}} | "mm Δt-1"
    # Pref. recharge to fast storage [mm Δt⁻¹]
    qrootzonefast::Vector{SVector{N,T}} | "mm Δt-1"
    # Pref. recharge to slow storage sum classes [mm Δt⁻¹]
    qrootzoneslow_m::Vector{T} | "mm Δt-1"
    # Capillary flux from the slow to the root-zone storage [mm Δt⁻¹]
    qcapillary::Vector{SVector{N,T}} | "mm Δt-1"
    # Capillary flux from the slow to the root-zone storage sum classes [mm Δt⁻¹]
    qcapillary_m::Vector{T} | "mm Δt-1"
    # Percolation flux from the root-zone to the slow storage [mm Δt⁻¹]
    qpercolation::Vector{SVector{N,T}} | "mm Δt-1"
    # Percolation flux from the root-zone to the slow storage sum classes [mm Δt⁻¹]
    qpercolation_m::Vector{T} | "mm Δt-1"
    # Evaporation from the root-zone storage, interception and hortonian [mm Δt⁻¹]
    actevap::Vector{SVector{N,T}} | "mm Δt-1"
    # Evaporation from the root-zone storage, interception and hortonian sum classes [mm Δt⁻¹]
    actevap_m::Vector{T} | "mm Δt-1"
    # Evaporation from the root-zone storage sum classes [mm Δt⁻¹]
    rootevap_m::Vector{T} | "mm Δt-1"
    #FAST
    # runoff from fast reservoir [mm Δt⁻¹]
    qfast::Vector{SVector{N,T}} | "mm Δt-1"
    #SLOW
    # runoff from slow reservoir [mm Δt⁻¹]
    qslow::Vector{T} | "mm Δt-1"
    #Total [mm Δt⁻¹]
    runoff::Vector{T} | "mm Δt-1"
    # fast runoff sum classes [mm Δt⁻¹]
    qfast_tot = Vector{T} | "mm Δt-1"


    ## WATERBALANCES
    #water balance snow store
    wb_snow::Vector{T} | "mm Δt-1"
    #water balance interception storage [mm Δt⁻¹]
    wb_interception::Vector{SVector{N,T}} | "mm Δt-1"
    #water balance hortonian ponding storage [mm Δt⁻¹]
    wb_hortonponding::Vector{SVector{N,T}} | "mm Δt-1"
    #water balance hortonian runoff storage [mm Δt⁻¹]
    wb_hortonrunoff::Vector{SVector{N,T}} | "mm Δt-1"
    #water balance root-zone storage [mm Δt⁻¹]
    wb_rootzone::Vector{SVector{N,T}} | "mm Δt-1"
    #water balance fast storage [mm Δt⁻¹]
    wb_fast::Vector{SVector{N,T}} | "mm Δt-1"
    #water balance slow storage [mm Δt⁻¹]
    wb_slow::Vector{T} | "mm Δt-1"
    #total water balance [mm Δt⁻¹]
    wb_tot::Vector{T} | "mm Δt-1"

end

statevars(::FLEXTOPO) = (
    :snow,
    :snowwater,
    :interceptionstorage,
    :hortonpondingstorage,
    :hortonrunoffstorage,
    :rootzonestorage,
    :faststorage,
    :slowstorage,
)

function common_snow_hbv(flextopo::FLEXTOPO)
    for i = 1:flextopo.n
        # precip correction
        precipcorr = flextopo.precipitation[i] * flextopo.pcorr[i]

        #hbv snow
        snow, snowwater, snowmelt, rainfallplusmelt, snowfall = snowpack_hbv(
            flextopo.snow[i],
            flextopo.snowwater[i],
            precipcorr,
            flextopo.temperature[i],
            flextopo.tti[i],
            flextopo.tt[i],
            flextopo.ttm[i],
            flextopo.cfmax[i],
            flextopo.whc[i],
        )

        #wb
        wb_snow =
            precipcorr - rainfallplusmelt - snow + flextopo.snow[i] - snowwater +
            flextopo.snowwater[i]

        #update stores
        # states_ is sum of states of previous time steps to compute WB at the end (per class); states_m is for combined stores
        flextopo.states_m[i] = flextopo.snow[i] + flextopo.snowwater[i]
        flextopo.snowfall[i] = snowfall
        flextopo.snowmelt[i] = snowmelt
        flextopo.snow[i] = snow
        flextopo.snowwater[i] = snowwater
        flextopo.rainfallplusmelt[i] = rainfallplusmelt
        flextopo.precipcorr[i] = precipcorr
        flextopo.wb_snow[i] = wb_snow
    end
end

function common_snow_no_storage(flextopo::FLEXTOPO)
    for i = 1:flextopo.n
        snow = 0.0
        snowwater = 0.0
        snowfall = 0.0
        snowmelt = 0.0
        # precip correction
        precipcorr = flextopo.precipitation[i] * flextopo.pcorr[i]
        rainfallplusmelt = precipcorr

        wb_snow = precipcorr - rainfallplusmelt - snow + flextopo.snow[i]

        #update stores
        # states_ is sum of states of previous time steps to compute WB at the end (per class); states_m is for combined stores
        flextopo.states_m[i] = flextopo.snow[i] + flextopo.snowwater[i]
        flextopo.snowfall[i] = snowfall
        flextopo.snowmelt[i] = snowmelt
        flextopo.snow[i] = snow
        flextopo.snowwater[i] = snowwater
        flextopo.rainfallplusmelt[i] = rainfallplusmelt
        flextopo.wb_snow[i] = wb_snow
        flextopo.precipcorr[i] = precipcorr
    end
end

function common_glaciers(flextopo::FLEXTOPO, config)
    modelglacier = get(config.model, "glacier", false)::Bool
    for i = 1:flextopo.n
        if modelglacier
            # Run Glacier module and add the snowpack on-top of it.
            # Estimate the fraction of snow turned into ice (HBV-light).
            # Estimate glacier melt.

            flextopo.snow[i], snow2glacier, flextopo.glacierstore[i], glaciermelt =
                glacier_hbv(
                    flextopo.glacierfrac[i],
                    flextopo.glacierstore[i],
                    flextopo.snow[i],
                    flextopo.temperature[i],
                    flextopo.g_tt[i],
                    flextopo.g_cfmax[i],
                    flextopo.g_sifrac[i],
                    Second(flextopo.Δt),
                )
            # Convert to mm per grid cell and add to snowmelt
            glaciermelt = glaciermelt * flextopo.glacierfrac[i]
            rainfallplusmelt = flextopo.rainfallplusmelt[i] + glaciermelt
        else
            rainfallplusmelt = flextopo.rainfallplusmelt[i]
        end
    end
end


function interception_no_storage(flextopo::FLEXTOPO)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        intevap = 0.0
        precipeffective = max(flextopo.rainfallplusmelt[i], 0.0)
        interceptionstorage = 0.0

        # correction for potential evaporation
        epotcorr = flextopo.ecorr[i] * flextopo.potential_evaporation[i]
        restevap = max(0.0, epotcorr - intevap)

        #wb interceptionstore
        wb_interception =
            flextopo.rainfallplusmelt[i] - intevap - precipeffective - interceptionstorage +
            flextopo.interceptionstorage[i][k]

        #update stores
        flextopo.states_[i] =
            setindex(flextopo.states_[i], flextopo.interceptionstorage[i][k], k)
        flextopo.interceptionstorage[i] =
            setindex(flextopo.interceptionstorage[i], interceptionstorage, k)
        flextopo.intevap[i] = setindex(flextopo.intevap[i], intevap, k)
        flextopo.precipeffective[i] =
            setindex(flextopo.precipeffective[i], precipeffective, k)
        flextopo.potsoilevap[i] = setindex(flextopo.potsoilevap[i], restevap, k)
        flextopo.epotcorr[i] = epotcorr
        flextopo.wb_interception[i] =
            setindex(flextopo.wb_interception[i], wb_interception, k)

        #average storage over classes
        flextopo.interceptionstorage_m[i] =
            sum(flextopo.interceptionstorage[i] .* flextopo.hrufrac[i])
        flextopo.intevap_m[i] = sum(flextopo.intevap[i] .* flextopo.hrufrac[i])

    end
end

function interception_overflow(flextopo::FLEXTOPO)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        # rainfall added to interception store interceptionstorage
        interception = min(
            flextopo.rainfallplusmelt[i],
            flextopo.imax[i][k] - flextopo.interceptionstorage[i][k],
        )
        # current interception storage
        interceptionstorage = flextopo.interceptionstorage[i][k] + interception
        precipeffective = flextopo.rainfallplusmelt[i] - interception

        # correction for potential evaporation
        epotcorr = flextopo.ecorr[i] * flextopo.potential_evaporation[i]

        # evaporation from interception storage
        intevap = min(interceptionstorage, epotcorr)
        interceptionstorage = interceptionstorage - intevap
        restevap = max(0.0, epotcorr - intevap)

        #wb interceptionstore
        wb_interception =
            flextopo.rainfallplusmelt[i] - intevap - precipeffective - interceptionstorage +
            flextopo.interceptionstorage[i][k]

        #update stores
        flextopo.states_[i] =
            setindex(flextopo.states_[i], flextopo.interceptionstorage[i][k], k)
        flextopo.interceptionstorage[i] =
            setindex(flextopo.interceptionstorage[i], interceptionstorage, k)
        flextopo.intevap[i] = setindex(flextopo.intevap[i], intevap, k)
        flextopo.precipeffective[i] =
            setindex(flextopo.precipeffective[i], precipeffective, k)
        flextopo.potsoilevap[i] = setindex(flextopo.potsoilevap[i], restevap, k)
        flextopo.epotcorr[i] = epotcorr
        flextopo.wb_interception[i] =
            setindex(flextopo.wb_interception[i], wb_interception, k)

        #average storage over classes
        flextopo.interceptionstorage_m[i] =
            sum(flextopo.interceptionstorage[i] .* flextopo.hrufrac[i])
        flextopo.intevap_m[i] = sum(flextopo.intevap[i] .* flextopo.hrufrac[i])

    end
end


function hortonponding_no_storage(flextopo::FLEXTOPO)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        hortonevap = 0.0
        qhortonpond = 0.0
        qhortonrootzone = max(flextopo.precipeffective[i][k], 0.0)
        hortonpondingstorage = 0.0

        wb_hortonponding =
            flextopo.precipeffective[i][k] - hortonevap - qhortonpond - qhortonrootzone -
            hortonpondingstorage + flextopo.hortonpondingstorage[i][k]

        #update with setindex
        flextopo.states_[i] = setindex(
            flextopo.states_[i],
            flextopo.states_[i][k] + flextopo.hortonpondingstorage[i][k],
            k,
        )
        flextopo.qhortonpond[i] = setindex(flextopo.qhortonpond[i], qhortonpond, k)
        flextopo.qhortonrootzone[i] =
            setindex(flextopo.qhortonrootzone[i], qhortonrootzone, k)
        flextopo.hortonpondingstorage[i] =
            setindex(flextopo.hortonpondingstorage[i], hortonpondingstorage, k)
        flextopo.hortonevap[i] = setindex(flextopo.hortonevap[i], hortonevap, k)
        flextopo.wb_hortonponding[i] =
            setindex(flextopo.wb_hortonponding[i], wb_hortonponding, k)

        #average storage over classes
        flextopo.hortonpondingstorage_m[i] =
            sum(flextopo.hortonpondingstorage[i] .* flextopo.hrufrac[i])
        flextopo.hortonevap_m[i] = sum(flextopo.hortonevap[i] .* flextopo.hrufrac[i])
    end
end

function hortonponding(flextopo::FLEXTOPO)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n

        #calculate reduction of shmax due to frost
        facc = min(flextopo.facc[i][k] + flextopo.temperature[i] * flextopo.kmf[i][k], 0.0)
        ft = min(
            max(
                facc / (flextopo.facc1[i][k] - flextopo.facc0[i][k]) -
                flextopo.facc0[i][k] / (flextopo.facc1[i][k] - flextopo.facc0[i][k]),
                flextopo.shmin[i][k],
            ),
            1.0,
        )
        shmax_frost = ft * flextopo.shmax[i][k]

        #add effective precipitation from interception to the horton ponding storage.
        hortonpondingstorage =
            flextopo.hortonpondingstorage[i][k] + flextopo.precipeffective[i][k]
        # if soil is filled until max capacity, additional water runs of directly
        directrunoff_h = max(hortonpondingstorage - shmax_frost, 0.0)
        #update hortonpondingstorage
        hortonpondingstorage = hortonpondingstorage - directrunoff_h
        #net water which infiltrates in horton (careful, shmax_frost can decrease from one timestep to another)
        netin_hortonpond = max(flextopo.precipeffective[i][k] - directrunoff_h, 0.0)

        #evaporation from the horton ponding storage
        hortonevap =
            hortonpondingstorage > shmax_frost * flextopo.lp[i][k] ?
            min(hortonpondingstorage, flextopo.potsoilevap[i][k]) :
            min(
                flextopo.potsoilevap[i][k] *
                (hortonpondingstorage / (shmax_frost * flextopo.lp[i][k])),
            )
        #update storage
        hortonpondingstorage = hortonpondingstorage - hortonevap
        #update restevap
        restevap = max(0.0, flextopo.potsoilevap[i][k] - hortonevap)

        #excess water from beta function of netin_hortonpond (due to rouding hortonpondingstorage could be slightly larger than shmax_frost, make sure it is always less than 1)
        qhorton_in =
            netin_hortonpond * (
                1.0 - pow(
                    (1.0 - min(hortonpondingstorage / shmax_frost, 1.0)),
                    flextopo.beta[i][k],
                )
            )
        #update storage
        hortonpondingstorage = hortonpondingstorage - qhorton_in

        #total water out of horton ponding consists of directrunoff and excess water from beta function
        qhortonpond = directrunoff_h + qhorton_in

        #infiltration to root-zone
        qhortonrootzone =
            min(hortonpondingstorage / shmax_frost, 1.0) > 0.0 ?
            min(
                flextopo.fmax[i][k] * exp(
                    -flextopo.fdec[i][k] *
                    (1.0 - min(hortonpondingstorage / shmax_frost, 1.0)),
                ),
                hortonpondingstorage,
            ) : 0.0
        hortonpondingstorage = hortonpondingstorage - qhortonrootzone

        #water balance
        wb_hortonponding =
            flextopo.precipeffective[i][k] - hortonevap - qhortonpond - qhortonrootzone -
            hortonpondingstorage + flextopo.hortonpondingstorage[i][k]

        #update states
        flextopo.states_[i] = setindex(
            flextopo.states_[i],
            flextopo.states_[i][k] + flextopo.hortonpondingstorage[i][k],
            k,
        )
        flextopo.potsoilevap[i] = setindex(flextopo.potsoilevap[i], restevap, k)
        flextopo.qhortonpond[i] = setindex(flextopo.qhortonpond[i], qhortonpond, k)
        flextopo.qhortonrootzone[i] =
            setindex(flextopo.qhortonrootzone[i], qhortonrootzone, k)
        flextopo.hortonpondingstorage[i] =
            setindex(flextopo.hortonpondingstorage[i], hortonpondingstorage, k)
        flextopo.hortonevap[i] = setindex(flextopo.hortonevap[i], hortonevap, k)
        flextopo.wb_hortonponding[i] =
            setindex(flextopo.wb_hortonponding[i], wb_hortonponding, k)

        #average storage over classes
        flextopo.hortonpondingstorage_m[i] =
            sum(flextopo.hortonpondingstorage[i] .* flextopo.hrufrac[i])
        flextopo.hortonevap_m[i] = sum(flextopo.hortonevap[i] .* flextopo.hrufrac[i])
    end
end

function hortonrunoff_no_storage(flextopo::FLEXTOPO)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        qhortonrun = 0.0
        hortonrunoffstorage = 0.0

        wb_hortonrunoff =
            flextopo.qhortonpond[i][k] - qhortonrun - hortonrunoffstorage +
            flextopo.hortonrunoffstorage[i][k]

        #update with setindex
        flextopo.states_[i] = setindex(
            flextopo.states_[i],
            flextopo.states_[i][k] + flextopo.hortonrunoffstorage[i][k],
            k,
        )
        flextopo.qhortonrun[i] = setindex(flextopo.qhortonrun[i], qhortonrun, k)
        flextopo.hortonrunoffstorage[i] =
            setindex(flextopo.hortonrunoffstorage[i], hortonrunoffstorage, k)
        flextopo.wb_hortonrunoff[i] =
            setindex(flextopo.wb_hortonrunoff[i], wb_hortonrunoff, k)

        #average storage over classes
        flextopo.hortonrunoffstorage_m[i] =
            sum(flextopo.hortonrunoffstorage[i] .* flextopo.hrufrac[i])
    end
end

function hortonrunoff(flextopo::FLEXTOPO)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n

        qhortonrun = min(
            flextopo.hortonrunoffstorage[i][k],
            flextopo.hortonrunoffstorage[i][k] * flextopo.khf[i][k],
        )
        hortonrunoffstorage =
            flextopo.hortonrunoffstorage[i][k] + flextopo.qhortonpond[i][k] - qhortonrun

        wb_hortonrunoff =
            flextopo.qhortonpond[i][k] - qhortonrun - hortonrunoffstorage +
            flextopo.hortonrunoffstorage[i][k]

        #update with setindex
        flextopo.states_[i] = setindex(
            flextopo.states_[i],
            flextopo.states_[i][k] + flextopo.hortonrunoffstorage[i][k],
            k,
        )
        flextopo.qhortonrun[i] = setindex(flextopo.qhortonrun[i], qhortonrun, k)
        flextopo.hortonrunoffstorage[i] =
            setindex(flextopo.hortonrunoffstorage[i], hortonrunoffstorage, k)
        flextopo.wb_hortonrunoff[i] =
            setindex(flextopo.wb_hortonrunoff[i], wb_hortonrunoff, k)

        #average storage over classes
        flextopo.hortonrunoffstorage_m[i] =
            sum(flextopo.hortonrunoffstorage[i] .* flextopo.hrufrac[i])

    end
end

function rootzone_no_storage(flextopo::FLEXTOPO)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        rootevap = 0.0
        qpercolation = 0.0
        qcapillary = 0.0
        qrootzone =
            max(flextopo.qhortonrootzone[i][k] + flextopo.rootzonestorage[i][k], 0.0) #if store not empty initial conditions
        rootzonestorage = 0.0

        #compute water balance rootzonestorage storage
        wb_rootzone =
            flextopo.qhortonrootzone[i][k] - rootevap - qrootzone - qpercolation +
            qcapillary - rootzonestorage + flextopo.rootzonestorage[i][k]

        #update states and fluxes with setindex
        flextopo.states_[i] = setindex(
            flextopo.states_[i],
            flextopo.states_[i][k] + flextopo.rootzonestorage[i][k],
            k,
        )
        flextopo.qrootzone[i] = setindex(flextopo.qrootzone[i], qrootzone, k)
        flextopo.rootzonestorage[i] =
            setindex(flextopo.rootzonestorage[i], rootzonestorage, k)
        flextopo.srootzone_over_srmax[i] = setindex(
            flextopo.srootzone_over_srmax[i],
            flextopo.rootzonestorage[i][k] / flextopo.srmax[i][k],
            k,
        )
        flextopo.qcapillary[i] = setindex(flextopo.qcapillary[i], qcapillary, k)
        flextopo.qpercolation[i] = setindex(flextopo.qpercolation[i], qpercolation, k)
        flextopo.rootevap[i] = setindex(flextopo.rootevap[i], rootevap, k)
        flextopo.actevap[i] = setindex(
            flextopo.actevap[i],
            flextopo.rootevap[i][k] + flextopo.intevap[i][k] + flextopo.hortonevap[i][k],
            k,
        )
        flextopo.wb_rootzone[i] = setindex(flextopo.wb_rootzone[i], wb_rootzone, k)

        #average storage over classes
        flextopo.srootzone_m[i] = sum(flextopo.rootzonestorage[i] .* flextopo.hrufrac[i])
        flextopo.srootzone_over_srmax_m[i] =
            sum(flextopo.srootzone_over_srmax[i] .* flextopo.hrufrac[i])
        flextopo.rootevap_m[i] = sum(flextopo.rootevap[i] .* flextopo.hrufrac[i])
        flextopo.actevap_m[i] = sum(flextopo.actevap[i] .* flextopo.hrufrac[i])
    end
end

function rootzone_storage(flextopo::FLEXTOPO)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n

        #added water to the root-zone. NB: if no horton storages: qhortonrootzone is in fact Pe!! (effective precip).
        rootzonestorage = flextopo.rootzonestorage[i][k] + flextopo.qhortonrootzone[i][k]
        # if soil is filled until max capacity, additional water runs of directly
        directrunoff = max(rootzonestorage - flextopo.srmax[i][k], 0.0)
        #update rootzonestorage
        rootzonestorage = rootzonestorage - directrunoff
        #net water which infiltrates in root-zone
        netin_rootzone = flextopo.qhortonrootzone[i][k] - directrunoff

        #evaporation from the rootzone
        rootevap =
            rootzonestorage > flextopo.srmax[i][k] * flextopo.lp[i][k] ?
            min(rootzonestorage, flextopo.potsoilevap[i][k]) :
            min(
                flextopo.potsoilevap[i][k] *
                (rootzonestorage / (flextopo.srmax[i][k] * flextopo.lp[i][k])),
            )
        #update storage
        rootzonestorage = rootzonestorage - rootevap

        #excess water from beta function of netin_rootzone
        qrootzone_in =
            netin_rootzone *
            (1.0 - pow(1.0 - rootzonestorage / flextopo.srmax[i][k], flextopo.beta[i][k]))
        #update storage
        rootzonestorage = rootzonestorage - qrootzone_in

        #total water out of root-zone consists of directrunoff and excess water from beta function
        qrootzone = directrunoff + qrootzone_in

        #percolation
        qpercolation = flextopo.perc[i][k] * rootzonestorage / flextopo.srmax[i][k]
        rootzonestorage = rootzonestorage - qpercolation

        #capillary rise
        qcapillary = min(
            flextopo.cap[i][k] * (1.0 - rootzonestorage / flextopo.srmax[i][k]),
            flextopo.slowstorage[i],
        )
        rootzonestorage = rootzonestorage + qcapillary

        #compute water balance rootzonestorage storage
        wb_rootzone =
            flextopo.qhortonrootzone[i][k] - rootevap - qrootzone - qpercolation +
            qcapillary - rootzonestorage + flextopo.rootzonestorage[i][k]

        #update states and fluxes with setindex
        flextopo.states_[i] = setindex(
            flextopo.states_[i],
            flextopo.states_[i][k] + flextopo.rootzonestorage[i][k],
            k,
        )
        flextopo.qrootzone[i] = setindex(flextopo.qrootzone[i], qrootzone, k)
        flextopo.rootzonestorage[i] =
            setindex(flextopo.rootzonestorage[i], rootzonestorage, k)
        flextopo.srootzone_over_srmax[i] = setindex(
            flextopo.srootzone_over_srmax[i],
            flextopo.rootzonestorage[i][k] / flextopo.srmax[i][k],
            k,
        )
        flextopo.qcapillary[i] = setindex(flextopo.qcapillary[i], qcapillary, k)
        flextopo.qpercolation[i] = setindex(flextopo.qpercolation[i], qpercolation, k)
        flextopo.rootevap[i] = setindex(flextopo.rootevap[i], rootevap, k)
        flextopo.actevap[i] = setindex(
            flextopo.actevap[i],
            flextopo.rootevap[i][k] + flextopo.intevap[i][k] + flextopo.hortonevap[i][k],
            k,
        )
        flextopo.wb_rootzone[i] = setindex(flextopo.wb_rootzone[i], wb_rootzone, k)

        #average storage over classes
        flextopo.srootzone_m[i] = sum(flextopo.rootzonestorage[i] .* flextopo.hrufrac[i])
        flextopo.srootzone_over_srmax_m[i] =
            sum(flextopo.srootzone_over_srmax[i] .* flextopo.hrufrac[i])
        flextopo.rootevap_m[i] = sum(flextopo.rootevap[i] .* flextopo.hrufrac[i])
        flextopo.actevap_m[i] = sum(flextopo.actevap[i] .* flextopo.hrufrac[i])
    end
end


function fast_no_storage(flextopo::FLEXTOPO)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        #make sure ds does not exceed 1
        ds = flextopo.ds[i][k] > 1.0 ? 1.0 : flextopo.ds[i][k]
        #calc inflow to faststorage
        qrootzonefast = flextopo.qrootzone[i][k] * (1.0 - ds)
        qfast = qrootzonefast + flextopo.faststorage[i][k] #if store not empty initial conditions, make sure to empty
        faststorage = 0.0

        #compute wb faststorage
        wb_fast = qrootzonefast - qfast - faststorage + flextopo.faststorage[i][k]

        #update with setindex
        flextopo.states_[i] = setindex(
            flextopo.states_[i],
            flextopo.states_[i][k] + flextopo.faststorage[i][k],
            k,
        )
        flextopo.qrootzonefast[i] = setindex(flextopo.qrootzonefast[i], qrootzonefast, k)
        flextopo.qfast[i] = setindex(flextopo.qfast[i], qfast, k)
        flextopo.faststorage[i] = setindex(flextopo.faststorage[i], faststorage, k)
        flextopo.wb_fast[i] = setindex(flextopo.wb_fast[i], wb_fast, k)

        #average storage over classes
        flextopo.faststorage_m[i] = sum(flextopo.faststorage[i] .* flextopo.hrufrac[i])
    end
end

function fast_storage(flextopo::FLEXTOPO)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n

        #make sure ds does not exceed 1
        ds = flextopo.ds[i][k] > 1.0 ? 1.0 : flextopo.ds[i][k]

        #split part of the outflow from the root-zone to the fast runoff (and part as preferential recharge to the slow reservoir)
        qrootzonefast = flextopo.qrootzone[i][k] * (1.0 - ds)

        #fast runoff
        qfast = min(
            flextopo.faststorage[i][k],
            pow(flextopo.faststorage[i][k], flextopo.alfa[i][k]) * flextopo.kf[i][k],
        )
        #update store
        faststorage = flextopo.faststorage[i][k] + qrootzonefast - qfast

        #compute wb faststorage
        wb_fast = qrootzonefast - qfast - faststorage + flextopo.faststorage[i][k]

        #update with setindex
        flextopo.states_[i] = setindex(
            flextopo.states_[i],
            flextopo.states_[i][k] + flextopo.faststorage[i][k],
            k,
        )
        flextopo.qrootzonefast[i] = setindex(flextopo.qrootzonefast[i], qrootzonefast, k)
        flextopo.qfast[i] = setindex(flextopo.qfast[i], qfast, k)
        flextopo.faststorage[i] = setindex(flextopo.faststorage[i], faststorage, k)
        flextopo.wb_fast[i] = setindex(flextopo.wb_fast[i], wb_fast, k)

        #average storage over classes
        flextopo.faststorage_m[i] = sum(flextopo.faststorage[i] .* flextopo.hrufrac[i])
    end
end

function slow_no_storage(flextopo::FLEXTOPO)
    for i = 1:flextopo.n
        #make sure ds does not exceed 1
        pref_recharge = flextopo.qrootzone[i] .* min.(flextopo.ds[i], 1.0)
        pref_recharge_sum_classes = sum(pref_recharge .* flextopo.hrufrac[i])

        Qcap_sum_classes = sum(flextopo.qcapillary[i] .* flextopo.hrufrac[i])
        Qperc_sum_classes = sum(flextopo.qpercolation[i] .* flextopo.hrufrac[i])

        Qsin = pref_recharge_sum_classes + Qperc_sum_classes - Qcap_sum_classes
        qslow = Qsin + flextopo.slowstorage[i] # if at start store is not empty
        slowstorage = 0.0

        qfast_tot =
            sum(flextopo.qfast[i] .* flextopo.hrufrac[i]) +
            sum(flextopo.qhortonrun[i] .* flextopo.hrufrac[i])
        runoff = max(0.0, qslow + qfast_tot)

        wb_slow = Qsin - qslow - slowstorage + flextopo.slowstorage[i]

        #update
        flextopo.qfast_tot[i] = qfast_tot
        flextopo.states_m[i] =
            flextopo.states_m[i] +
            sum(flextopo.states_[i] .* flextopo.hrufrac[i]) +
            flextopo.slowstorage[i]
        flextopo.runoff[i] = runoff
        flextopo.qrootzoneslow_m[i] = pref_recharge_sum_classes
        flextopo.qpercolation_m[i] = Qperc_sum_classes
        flextopo.qcapillary_m[i] = Qcap_sum_classes
        flextopo.qslow[i] = qslow
        flextopo.slowstorage[i] = slowstorage
        flextopo.wb_slow[i] = wb_slow
    end
end

function common_slow_storage(flextopo::FLEXTOPO)
    for i = 1:flextopo.n
        #make sure ds does not exceed 1
        pref_recharge = flextopo.qrootzone[i] .* min.(flextopo.ds[i], 1.0)
        pref_recharge_sum_classes = sum(pref_recharge .* flextopo.hrufrac[i])

        Qcap_sum_classes = sum(flextopo.qcapillary[i] .* flextopo.hrufrac[i])
        Qperc_sum_classes = sum(flextopo.qpercolation[i] .* flextopo.hrufrac[i])

        Qsin = pref_recharge_sum_classes + Qperc_sum_classes - Qcap_sum_classes
        slowstorage = flextopo.slowstorage[i] + Qsin

        qslow = min(flextopo.slowstorage[i], flextopo.slowstorage[i] * flextopo.ks[i])
        slowstorage = slowstorage - qslow

        qfast_tot =
            sum(flextopo.qfast[i] .* flextopo.hrufrac[i]) +
            sum(flextopo.qhortonrun[i] .* flextopo.hrufrac[i])
        runoff = max(0.0, qslow + qfast_tot)

        wb_slow = Qsin - qslow - slowstorage + flextopo.slowstorage[i]

        #update
        flextopo.qfast_tot[i] = qfast_tot
        flextopo.states_m[i] =
            flextopo.states_m[i] +
            sum(flextopo.states_[i] .* flextopo.hrufrac[i]) +
            flextopo.slowstorage[i]
        flextopo.runoff[i] = runoff
        flextopo.qrootzoneslow_m[i] = pref_recharge_sum_classes
        flextopo.qpercolation_m[i] = Qperc_sum_classes
        flextopo.qcapillary_m[i] = Qcap_sum_classes
        flextopo.qslow[i] = qslow
        flextopo.slowstorage[i] = slowstorage
        flextopo.wb_slow[i] = wb_slow
    end
end

function watbal(flextopo::FLEXTOPO)
    for i = 1:flextopo.n
        states =
            flextopo.snow[i] +
            flextopo.snowwater[i] +
            sum(flextopo.interceptionstorage[i] .* flextopo.hrufrac[i]) +
            sum(flextopo.hortonpondingstorage[i] .* flextopo.hrufrac[i]) +
            sum(flextopo.hortonrunoffstorage[i] .* flextopo.hrufrac[i]) +
            sum(flextopo.rootzonestorage[i] .* flextopo.hrufrac[i]) +
            sum(flextopo.faststorage[i] .* flextopo.hrufrac[i]) +
            sum(flextopo.slowstorage[i] .* flextopo.hrufrac[i])
        wb_tot =
            flextopo.precipitation[i] - (
                sum(flextopo.intevap[i] .* flextopo.hrufrac[i]) +
                sum(flextopo.hortonevap[i] .* flextopo.hrufrac[i]) +
                sum(flextopo.rootevap[i] .* flextopo.hrufrac[i])
            ) - flextopo.runoff[i] - states + flextopo.states_m[i]
        #update wb
        flextopo.wb_tot[i] = wb_tot
    end
end
