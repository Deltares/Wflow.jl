Base.@kwdef struct HBV{T}
    n::Int                      # Number of cells
    yl::Vector{T}               # Length of cells in y direction [m]
    xl::Vector{T}               # Length of cells in x direction [m]
    altitude::Vector{T}         # Vertical elevation [m]
    fc::Vector{T}               # Field capacity [mm]
    betaseepage::Vector{T}      # Exponent in soil runoff generation equation [-]
    lp::Vector{T}               # Fraction of field capacity below which actual evaporation=potential evaporation
    threshold::Vector{T}        # Threshold soilwater storage above which AE=PE
    k4::Vector{T}               # Recession constant baseflow
    kquickflow::Vector{T}       
    suz::Vector{T}              # Level over which k0 is used
    k0::Vector{T}               
    khq::Vector{T}              # Recession rate at flow hq
    hq::Vector{T}               # High flow rate hq for which recession rate of upper reservoir is known
    alphanl::Vector{T}          # Measure of non-linearity of upper reservoir
    perc::Vector{T}             # Percolation from upper to lower zone [mm/day]
    cfr::Vector{T}              # Refreezing efficiency constant in refreezing of freewater in snow
    pcorr::Vector{T}            # Correction factor for precipitation
    rfcf::Vector{T}             # Correction factor for rainfall
    sfcf::Vector{T}             # Vorrection factor for snowfall
    cflux::Vector{T}            # Maximum capillary rise from runoff response routine to soil moisture routine
    icf::Vector{T}              # Maximum interception storage (in forested and non-forested areas)
    cevpf::Vector{T}            # Correction factor for potential evaporation
    epf::Vector{T}              # Exponent of correction factor for evaporation on days with precipitation
    ecorr::Vector{T}            # Evap correction
    tti::Vector{T}              # Critical temperature for snowmelt and refreezing [ᵒC]
    tt::Vector{T}               # Defines interval in which precipitation falls as rainfall and snowfall [ᵒC]
    ttm::Vector{T}              # Threshold temperature for snowmelt [ᵒC]
    cfmax::Vector{T}            # Meltconstant in temperature-index [-]
    whc::Vector{T}              # Fraction of snow volume that can store water [-]
    precipitation::Vector{T} = fill(mv, n)          # Precipitation [mm]
    temperature::Vector{T} = fill(mv, n)           # Temperature [ᵒC]
    potential_evaporation::Vector{T} = fill(mv, n)                # Potential evapotranspiration [mm]
    potsoilevap::Vector{T} = fill(mv, n)            # Potential soil evaporation [mm]
    soilevap::Vector{T} = fill(mv, n)               # Soil evaporation [mm] 
    intevap::Vector{T} = fill(mv, n)                # Evaporation from interception storage [mm]
    actevap::Vector{T} = fill(mv, n)                # Total actual evapotranspiration (intevap + soilevap) [mm] 
    interceptionstorage::Vector{T} = fill(mv, n)    # Actual interception storage [mm]
    snowwater::Vector{T} = fill(mv, n)              # Available free water in snow [mm]
    snow::Vector{T} = fill(mv, n)                   # Snow pack [mm]
    rainfallplusmelt::Vector{T} = fill(mv, n)       # Snow melt + precipitation as rainfall [mm]
    soilmoisture::Vector{T} = fill(mv, n)           # Actual soil moisture [mm]
    directrunoff::Vector{T} = fill(mv, n)           # Direct runoff to upper zone [mm]
    hbv_seepage::Vector{T} = fill(mv, n)            # Recharge to upper zone [mm]
    in_upperzone::Vector{T} = fill(mv, n)           # Water inflow into upper zone [mm]
    upperzonestorage::Vector{T} = fill(mv, n)       # Water content of the upper zone [mm]
    quickflow::Vector{T} = fill(mv, n)              # Specific runoff (quickflow part) [mm]
    real_quickflow::Vector{T} = fill(mv, n)         # specific runoff (quickflow), if K upper zone is precalculated [mm]
    percolation::Vector{T} = fill(mv, n)            # Actual percolation to the lower zone [mm]
    capflux::Vector{T} = fill(mv, n)                # Capilary rise [mm]
    lowerzonestorage::Vector{T} = fill(mv, n)       # Water content of the lower zone [mm]
    baseflow::Vector{T} = fill(mv, n)               # Specific runoff (baseflow part) per cell [mm]
    runoff::Vector{T} = fill(mv, n)                 # Total specific runoff per cell [mm]
end

statevars(::HBV) = (:soilmoisture, :snow, :snowwater, :upperzonestorage, :lowerzonestorage, :interceptionstorage,)

function update_until_snow(hbv::HBV, config)

    for i = 1:hbv.n
        precipitation = hbv.precipitation[i] * hbv.pcorr[i]
        # fraction of precipitation which falls as rain
        rainfrac = if iszero(hbv.tti[i])
            Float64(hbv.temperature[i] > hbv.tt[i])
        else
            frac = (hbv.temperature[i] - (hbv.tt[i] - hbv.tti[i] / 2.0)) / hbv.tti[i]
            min(frac, 1.0)
        end
        rainfrac = max(rainfrac, 0.0)

        # fraction of precipitation which falls as snow
        snowfrac = 1.0 - rainfrac
        # different correction for rainfall and snowfall
        precipitation = hbv.sfcf[i] * snowfrac * precipitation + hbv.rfcf[i] * rainfrac * precipitation

        # interception
        interception = min(precipitation, hbv.icf[i] - hbv.interceptionstorage[i])
        # current interception storage
        interceptionstorage = hbv.interceptionstorage[i] + interception
        precipitation = precipitation - interception

        # correction for potential evaporation on wet days
        potevap = exp(-hbv.epf[i] * precipitation) * hbv.ecorr[i] * hbv.potential_evaporation[i]
        # correct per landuse
        potevap = hbv.cevpf[i] * potevap  

        # evaporation from interception storage
        intevap = min(interceptionstorage, potevap)
        interceptionstorage = interceptionstorage - intevap
        restevap = max(0.0, potevap - intevap)

        snow, snowwater, snowmelt, rainfallplusmelt, snowfall = snowpack_hbv(
            hbv.snow[i],
            hbv.snowwater[i],
            precipitation,
            hbv.temperature[i],
            hbv.tti[i],
            hbv.tt[i],
            hbv.ttm[i],
            hbv.cfmax[i],
            hbv.whc[i],
        )

        # update the outputs and states
        hbv.rainfallplusmelt[i] = rainfallplusmelt
        hbv.snowwater[i] = snowwater
        hbv.snow[i] = snow
        hbv.interceptionstorage[i] = interceptionstorage
        hbv.potsoilevap[i] = restevap
        hbv.intevap[i] = intevap    
    end
end

function update_after_snow(hbv::HBV, config)

    modelglacier = get(config.model, "glacier", false)
    set_kquickflow = get(config.model, "set_kquickflow", false)
    external_qbase = get(config.model, "external_qbase", false)

    for i = 1:hbv.n
        if modelglacier
            # Run Glacier module and add the snowpack on-top of it.
            # Estimate the fraction of snow turned into ice (HBV-light).
            # Estimate glacier melt.

            drysnow, snow2glacier, glacierstore, glaciermelt = glacier_hbv(
                hbv.glacierfrac[i],
                hbv.glacierstore[i],
                hbv.drysnow[i],
                hbv.temperature[i],
                hbv.g_tt[i],
                hbv.g_cfmax[i],
                hbv.g_sifrac[i],
                Δt,
            )
            # Convert to mm per grid cell and add to snowmelt
            glaciermelt = glaciermelt * hbv.glacierfrac[i]
            rainfallplusmelt = hbv.rainfallplusmelt[i] + glaciermelt
        else
            rainfallplusmelt = hbv.rainfallplusmelt[i]        
        end

        soilmoisture = hbv.soilmoisture[i] + rainfallplusmelt
        # if soil is filled to capacity: abundant water runs of directly
        directrunoff = max(soilmoisture - hbv.fc[i], 0.0)
        soilmoisture = soilmoisture - directrunoff
        # net water which infiltrates into soil
        netinsoil = rainfallplusmelt - directrunoff
        
        # soil evapotranspiration
        soilevap = soilmoisture > hbv.threshold[i] ? min(soilmoisture, hbv.potsoilevap[i]) : min(hbv.potsoilevap[i] * (soilmoisture/hbv.threshold[i]))
        # evaporation from soil moisture storage
        soilmoisture = soilmoisture - soilevap
        # sum of evaporation components (IntEvap+SoilEvap)
        actevap = hbv.intevap[i] + soilevap
        # runoff water from soil
        hbv_seepage = pow(min(soilmoisture/hbv.fc[i],1.0), hbv.betaseepage[i]) * netinsoil
        soilmoisture = soilmoisture - hbv_seepage
        # correction for extremely wet periods: soil is filled to capacity
        back_tosoil = min(hbv.fc[i] - soilmoisture, directrunoff)
        directrunoff = directrunoff - back_tosoil
        soilmoisture = soilmoisture + back_tosoil
        # total water available for runoff
        in_upperzone = directrunoff + hbv_seepage
        ### calculations for upper zone ###
        upperzonestorage = hbv.upperzonestorage[i] + in_upperzone
        percolation = min(hbv.perc[i], upperzonestorage - in_upperzone / 2.0)
        upperzonestorage = upperzonestorage - percolation
        # capillary flux flowing back to soil  
        capflux = hbv.cflux[i] * (( hbv.fc[i] - soilmoisture) / hbv.fc[i])
        capflux = min(hbv.fc[i] - soilmoisture, capflux)
        upperzonestorage = upperzonestorage - capflux
        soilmoisture = soilmoisture + capflux

        real_quickflow = 0.0
        if set_kquickflow == false

            if percolation < hbv.perc[i]
                quickflow = 0.0
            else
                quickflow = min(hbv.kquickflow[i] * pow((upperzonestorage - min(in_upperzone / 2.0, upperzonestorage)), 1.0 + hbv.alphanl[i]), upperzonestorage)
            
            upperzonestorage = percolation < hbv.perc[i] ? upperzonestorage : max(upperzonestorage - quickflow, 0.0)
            end
        else
            quickflow = hbv.kquickflow[i] * upperzonestorage
            real_quickflow = max(0.0, hbv.k0[i] * (upperzonestorage - hbv.suz[i]))
            upperzonestorage = upperzonestorage - quickflow - real_quickflow

        end

        ### calculations for lower zone ###
        lowerzonestorage = hbv.lowerzonestorage[i] + percolation
        # baseflow in mm/timestep
        baseflow = min(lowerzonestorage, hbv.k4[i] * lowerzonestorage)
        lowerzonestorage = lowerzonestorage - baseflow

        if external_qbase
            directrunoffstorage = quickflow + seepage + real_quickflow
        else
            directrunoffstorage = quickflow + baseflow + real_quickflow
        end

        runoff = max(0.0, directrunoffstorage)
        
        # update the outputs and states
        hbv.runoff[i] = runoff
        hbv.rainfallplusmelt[i] = rainfallplusmelt
        hbv.soilmoisture[i] = soilmoisture
        hbv.soilevap[i] = soilevap
        hbv.actevap[i] = actevap
        hbv.quickflow[i] = quickflow
        hbv.real_quickflow[i] = real_quickflow
        hbv.upperzonestorage[i] = upperzonestorage
        hbv.lowerzonestorage[i] = lowerzonestorage
        hbv.baseflow[i] = baseflow
        hbv.hbv_seepage[i] = hbv_seepage
        hbv.percolation[i] = percolation
        hbv.capflux[i] = capflux
        hbv.in_upperzone[i] = in_upperzone
    end

end