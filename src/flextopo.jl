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
    selectSw::Vector{String}
    selectSi::Vector{String}
    selectSh::Vector{String}
    selectShf::Vector{String}
    selectSr::Vector{String}
    selectSf::Vector{String}
    selectSs::Vector{String} 

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
    # Fraction of the snowpack on top of the glacier converted into ice [-]
    g_sifrac::Vector{T} | "-"           
    # Water within the glacier [mm]
    glacierstore::Vector{T} | "mm"      
    # Fraction covered by a glacier [-]
    glacierfrac::Vector{T} | "-"        
    ##INTERCEPTION
    # Maximum interception storage (in forested and non-forested areas) [mm]
    # imax::Vector{T} | "mm"           
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
    #lag time of the hortonian runoff storage [Δt]
    # Thf::Vector{T} | "Δt"            
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
    # fraction of Qr to Ss (1-ds to Sf)
    ds::Vector{SVector{N,T}} | "-"
    # SLOW
    #recession coefficient of slow storage [Δt-1]
    ks::Vector{T} | "Δt-1"            

    ## STATES
    ##SNOW
    # Snow water equivalent [mm]
    Sw::Vector{T} | "mm"
    # Liquid water content in the snow pack [mm]
    Sww::Vector{T} | "mm"
    # Interception storage [mm]
    Si::Vector{SVector{N,T}} | "mm"  
    # Storage in the hortonian ponding reservoir [mm]
    Sh::Vector{SVector{N,T}} | "mm"
    # Storage in the hortonian runoff generation [mm]
    Shf::Vector{SVector{N,T}} | "mm"  
    # Storage in the root-zone [mm]
    Sr::Vector{SVector{N,T}} | "mm" 
    # Storage in the root-zone relative to maximum root-zone storage capacity [mm]
    Sr_over_srmax::Vector{SVector{N,T}} | "mm" 
    # Storage in the fast store [mm]
    Sf::Vector{SVector{N,T}} | "mm" 
    # Storage in the slow reservoir (for Qcap calc) [mm]
    Ss::Vector{T} | "mm"  
    #states previous time step to calc water balance [mm]
    states_::Vector{SVector{N,T}} | "mm"  
    #states previous time step to calc water balance combined based on perc class. [mm]
    states_m::Vector{T} | "mm"  
    #states averaged over classes
    # Sw_m::Vector{T} | "mm"  
    # Sww_m::Vector{T} | "mm"  
    Si_m::Vector{T} | "mm"  
    Sh_m::Vector{T} | "mm"  
    Shf_m::Vector{T} | "mm"  
    Sr_m::Vector{T} | "mm"  
    Sf_m::Vector{T} | "mm"  
    Sr_over_srmax_m::Vector{T} | "mm"  

    ## FLUXES
    #SNOW
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{T} | "mm Δt-1"
    # Temperature [ᵒC]
    temperature::Vector{T} | "ᵒC"       
    # Potential evapotranspiration [mm Δt⁻¹]
    potential_evaporation::Vector{T} | "mm Δt-1"    
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
    Ei::Vector{SVector{N,T}} | "mm Δt-1"                      
    #effective precipitation [mm Δt⁻¹]
    Pe::Vector{SVector{N,T}} | "mm Δt-1"       
    # Evaporation from interception sum classes [mm Δt⁻¹]
    Ei_m::Vector{T} | "mm Δt-1"                  
    #HORTONPONDING
    # Evaporation from the hortonion ponding storage [-]
    Eh::Vector{SVector{N,T}} | "mm Δt-1"           
    # Flux from the hortonian ponding storage to the hortonian runoff storage [mm Δt⁻¹]
    Qh::Vector{SVector{N,T}} | "mm Δt-1"                  
    # Flux from the hortonian ponding storage to the root zone storage [mm Δt⁻¹]
    Qhr::Vector{SVector{N,T}} | "mm Δt-1"                 
    # modeled accumulated frost [ᵒC Δt]
    facc::Vector{SVector{N,T}} | "ᵒC Δt"
    # Evaporation from the hortonian sum classes [mm Δt⁻¹]
    Eh_m::Vector{T} | "mm Δt-1"                  
    #HORTONRUNOFF
    # Flux from the hortonian runoff storage [mm Δt⁻¹]
    Qhf::Vector{SVector{N,T}} | "mm Δt-1"       
    #ROOTZONE
    # Evaporation from the root-zone storage [mm Δt⁻¹]
    Er::Vector{SVector{N,T}} | "mm Δt-1"                  
    # Flux from the root-zone storage [mm Δt⁻¹]
    Qr::Vector{SVector{N,T}} | "mm Δt-1"                  
    # Pref. recharge to fast storage [mm Δt⁻¹]
    Qrf::Vector{SVector{N,T}} | "mm Δt-1"                  
    # Pref. recharge to slow storage sum classes [mm Δt⁻¹]
    Qrs_m::Vector{T} | "mm Δt-1"                  
    # Capillary flux from the slow to the root-zone storage [mm Δt⁻¹]
    Qcap::Vector{SVector{N,T}} | "mm Δt-1"        
    # Capillary flux from the slow to the root-zone storage sum classes [mm Δt⁻¹]
    Qcap_m::Vector{T} | "mm Δt-1"                  
    # Percolation flux from the root-zone to the slow storage [mm Δt⁻¹]
    Qperc::Vector{SVector{N,T}} | "mm Δt-1"                 
    # Percolation flux from the root-zone to the slow storage sum classes [mm Δt⁻¹]
    Qperc_m::Vector{T} | "mm Δt-1"                 
    # Evaporation from the root-zone storage, interception and hortonian [mm Δt⁻¹]
    Ea::Vector{SVector{N,T}} | "mm Δt-1"                  
    # Evaporation from the root-zone storage, interception and hortonian sum classes [mm Δt⁻¹]
    Ea_m::Vector{T} | "mm Δt-1"               
    # Evaporation from the root-zone storage sum classes [mm Δt⁻¹]
    Er_m::Vector{T} | "mm Δt-1"                     
    #FAST 
    # runoff from fast reservoir [mm Δt⁻¹]
    Qf::Vector{SVector{N,T}} | "mm Δt-1"             
    #SLOW
    # runoff from slow reservoir [mm Δt⁻¹]
    Qs::Vector{T} | "mm Δt-1"                  
    #Total [mm Δt⁻¹]
    runoff::Vector{T} | "mm Δt-1"                  
    # fast runoff sum classes [mm Δt⁻¹]
    Qftotal = Vector{T} | "mm Δt-1"                  


    ## WATERBALANCES
    #water balance snow store
    wbSw::Vector{T} | "mm Δt-1" 
    #water balance interception storage [mm Δt⁻¹]
    wbSi::Vector{SVector{N,T}} | "mm Δt-1" 
    #water balance hortonian ponding storage [mm Δt⁻¹]
    wbSh::Vector{SVector{N,T}} | "mm Δt-1" 
    #water balance hortonian runoff storage [mm Δt⁻¹]
    wbShf::Vector{SVector{N,T}} | "mm Δt-1" 
    #water balance root-zone storage [mm Δt⁻¹]
    wbSr::Vector{SVector{N,T}} | "mm Δt-1" 
    #water balance fast storage [mm Δt⁻¹]
    wbSf::Vector{SVector{N,T}} | "mm Δt-1" 
    #water balance slow storage [mm Δt⁻¹]
    wbSs::Vector{T} | "mm Δt-1" 
    #total water balance [mm Δt⁻¹]
    wbtot::Vector{T} | "mm Δt-1" 
    

    # function FLEXTOPO{T,N}(args...) where {T,N}
    #     equal_size_vectors(args)
    #     return new(args...)
    # end
end

statevars(::FLEXTOPO) = (
    :Sw,
    :Sww,
    :Si,
    :Sh,
    :Shf,
    :Sr,
    :Sf,
    :Ss
)

function common_snow_hbv(flextopo::FLEXTOPO, config)
    for i = 1:flextopo.n
        # precip correction 
        precipitation = flextopo.precipitation[i] * flextopo.pcorr[i]

        #hbv snow 
        Sw, Sww, snowmelt, rainfallplusmelt, snowfall = snowpack_hbv(
            flextopo.Sw[i],
            flextopo.Sww[i],
            precipitation,
            flextopo.temperature[i],
            flextopo.tti[i],
            flextopo.tt[i],
            flextopo.ttm[i],
            flextopo.cfmax[i],
            flextopo.whc[i],
            )
        
        #update stores
        # states_ is sum of states of previous time steps to compute WB at the end (per class); states_m is for combined stores
        flextopo.states_m[i] = flextopo.Sw[i] + flextopo.Sww[i]
        flextopo.snowfall[i] =  snowfall
        flextopo.snowmelt[i] = snowmelt
        flextopo.Sw[i] = Sw
        flextopo.Sww[i] = Sww
        flextopo.rainfallplusmelt[i] = rainfallplusmelt
        #TODO add wbSw wbSww
    end
end

function common_snow_no_storage(flextopo::FLEXTOPO, config)
    for i = 1:flextopo.n
        Sw = 0.0
        Sww = 0.0
        snowfall = 0.0
        snowmelt = 0.0
        # precip correction 
        precipitation = flextopo.precipitation[i] * flextopo.pcorr[i] 
        rainfallplusmelt = precipitation

        wbSw = precipitation - rainfallplusmelt - Sw + flextopo.Sw[i]

        #update stores 
        # states_ is sum of states of previous time steps to compute WB at the end (per class); states_m is for combined stores
        flextopo.states_m[i] = flextopo.Sw[i] + flextopo.Sww[i]
        flextopo.snowfall[i] =  snowfall
        flextopo.snowmelt[i] = snowmelt
        flextopo.Sw[i] = Sw
        flextopo.Sww[i] = Sww
        flextopo.rainfallplusmelt[i] = rainfallplusmelt
        flextopo.wbSw[i] = wbSw
        #TODO add wbSww
    end
end

function common_glaciers(flextopo::FLEXTOPO, config)
    modelglacier = get(config.model, "glacier", false)::Bool
    for i = 1:flextopo.n
        if modelglacier
            # Run Glacier module and add the snowpack on-top of it.
            # Estimate the fraction of snow turned into ice (HBV-light).
            # Estimate glacier melt.

            flextopo.Sw[i], snow2glacier, flextopo.glacierstore[i], glaciermelt = glacier_hbv(
                flextopo.glacierfrac[i],
                flextopo.glacierstore[i],
                flextopo.Sw[i],
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
        

function interception_no_storage(flextopo::FLEXTOPO, config)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        Ei = 0.0
        Pe = max(flextopo.rainfallplusmelt[i], 0)
        Si = 0.0

        # correction for potential evaporation 
        Ep_corr = flextopo.ecorr[i] * flextopo.potential_evaporation[i]
        restevap = max(0.0, Ep_corr - Ei)

        #wb interceptionstore 
        wbSi = flextopo.rainfallplusmelt[i] - Ei - Pe - Si + flextopo.Si[i][k]

        #update stores
        flextopo.states_[i] = setindex(flextopo.states_[i], flextopo.Si[i][k], k)
        flextopo.Si[i] = setindex(flextopo.Si[i], Si, k)
        flextopo.Ei[i] = setindex(flextopo.Ei[i], Ei, k)
        flextopo.Pe[i] = setindex(flextopo.Pe[i], Pe, k)
        flextopo.potsoilevap[i] = setindex(flextopo.potsoilevap[i], restevap, k)
        flextopo.wbSi[i] = setindex(flextopo.wbSi[i], wbSi, k)
        
        #average storage over classes
        flextopo.Si_m[i] = sum(flextopo.Si[i] .* flextopo.hrufrac[i])
        flextopo.Ei_m[i] = sum(flextopo.Ei[i] .* flextopo.hrufrac[i])

    end
end

function interception_overflow(flextopo::FLEXTOPO, config)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        # rainfall added to interception store Si
        interception = min(flextopo.rainfallplusmelt[i], flextopo.imax[i][k] - flextopo.Si[i][k])
        # current interception storage
        Si = flextopo.Si[i][k] + interception
        Pe = flextopo.rainfallplusmelt[i] - interception

        # correction for potential evaporation 
        Ep_corr = flextopo.ecorr[i] * flextopo.potential_evaporation[i]

        # evaporation from interception storage
        Ei = min(Si, Ep_corr)
        Si = Si - Ei
        restevap = max(0.0, Ep_corr - Ei)

        #wb interceptionstore 
        wbSi = flextopo.rainfallplusmelt[i] - Ei - Pe - Si + flextopo.Si[i][k]

        #update stores
        flextopo.states_[i] = setindex(flextopo.states_[i], flextopo.Si[i][k], k)
        flextopo.Si[i] = setindex(flextopo.Si[i], Si, k)
        flextopo.Ei[i] = setindex(flextopo.Ei[i], Ei, k)
        flextopo.Pe[i] = setindex(flextopo.Pe[i], Pe, k)
        flextopo.potsoilevap[i] = setindex(flextopo.potsoilevap[i], restevap, k)
        flextopo.wbSi[i] = setindex(flextopo.wbSi[i], wbSi, k)
        
        #average storage over classes
        flextopo.Si_m[i] = sum(flextopo.Si[i] .* flextopo.hrufrac[i])
        flextopo.Ei_m[i] = sum(flextopo.Ei[i] .* flextopo.hrufrac[i])

    end
end


function hortonponding_no_storage(flextopo::FLEXTOPO, config)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        Eh = 0.0
        Qh = 0.0
        Qhr = max(flextopo.Pe[i][k], 0)
        Sh = 0.0

        wbSh = flextopo.Pe[i][k] - Eh - Qh - Qhr - Sh + flextopo.Sh[i][k]

        #update with setindex
        flextopo.states_[i] = setindex(flextopo.states_[i], flextopo.states_[i][k]  + flextopo.Sh[i][k], k)
        flextopo.Qh[i] = setindex(flextopo.Qh[i], Qh, k)
        flextopo.Qhr[i] = setindex(flextopo.Qhr[i], Qhr, k)
        flextopo.Sh[i] = setindex(flextopo.Sh[i], Sh, k)
        flextopo.Eh[i] = setindex(flextopo.Eh[i], Eh, k)
        flextopo.wbSh[i] = setindex(flextopo.wbSh[i], wbSh, k)

        #average storage over classes
        flextopo.Sh_m[i] = sum(flextopo.Sh[i] .* flextopo.hrufrac[i])
        flextopo.Eh_m[i] = sum(flextopo.Eh[i] .* flextopo.hrufrac[i])
    end
end

function hortonponding(flextopo::FLEXTOPO, config)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n

        #calculate reduction of shmax due to frost
        facc = min(flextopo.facc[i][k] + flextopo.temperature[i] * flextopo.kmf[i][k], 0.0) #TODO check units kmf, facc
        ft = min(max(facc / (flextopo.facc1[i][k] - flextopo.facc0[i][k]) - flextopo.facc0[i][k] / (flextopo.facc1[i][k] - flextopo.facc0[i][k]), flextopo.shmin[i][k]), 1.0)
        shmax_frost = ft * flextopo.shmax[i][k]

        #add effective precipitation from interception to the horton ponding storage.  
        Sh = flextopo.Sh[i][k] + flextopo.Pe[i][k]
        # if soil is filled until max capacity, additional water runs of directly
        directrunoff_h = max(Sh - shmax_frost, 0.0)
        #update Sh
        Sh = Sh - directrunoff_h
        #net water which infiltrates in horton (careful, shmax_frost can decrease from one timestep to another)
        netinSh = max(flextopo.Pe[i][k] - directrunoff_h, 0)

        #evaporation from the horton ponding storage
        Eh = Sh > shmax_frost * flextopo.lp[i][k] ? min(Sh, flextopo.potsoilevap[i][k]) :
            min(flextopo.potsoilevap[i][k] * (Sh / (shmax_frost * flextopo.lp[i][k])))
        #update storage
        Sh = Sh - Eh
        #update restevap 
        restevap = max(0.0, flextopo.potsoilevap[i][k] - Eh)

        #excess water from beta function of netinSh (due to rouding Sh could be slightly larger than shmax_frost, make sure it is always less than 1)
        Qh1 = netinSh * (1 - pow((1 - min(Sh / shmax_frost, 1)), flextopo.beta[i][k]))
        #update storage
        Sh = Sh - Qh1

        #total water out of horton ponding consists of directrunoff and excess water from beta function
        Qh = directrunoff_h + Qh1

        #infiltration to root-zone 
        # Qhr = max(0.0, flextopo.fmax[i][k] * exp(-flextopo.fdec[i][k] * (1 - Sh/shmax_frost)))
        Qhr = min(Sh / shmax_frost, 1) > 0 ? min(flextopo.fmax[i][k] * exp(-flextopo.fdec[i][k] * (1 - min(Sh / shmax_frost, 1))), Sh) : 0
        Sh = Sh - Qhr

        #water balance
        wbSh = flextopo.Pe[i][k] - Eh - Qh - Qhr - Sh + flextopo.Sh[i][k]

        #update states
        flextopo.states_[i] = setindex(flextopo.states_[i], flextopo.states_[i][k]  + flextopo.Sh[i][k], k)
        flextopo.potsoilevap[i] = setindex(flextopo.potsoilevap[i], restevap, k)
        flextopo.Qh[i] = setindex(flextopo.Qh[i], Qh, k)
        flextopo.Qhr[i] = setindex(flextopo.Qhr[i], Qhr, k)
        flextopo.Sh[i] = setindex(flextopo.Sh[i], Sh, k)
        flextopo.Eh[i] = setindex(flextopo.Eh[i], Eh, k)
        flextopo.wbSh[i] = setindex(flextopo.wbSh[i], wbSh, k)

        #average storage over classes
        flextopo.Sh_m[i] = sum(flextopo.Sh[i] .* flextopo.hrufrac[i])
        flextopo.Eh_m[i] = sum(flextopo.Eh[i] .* flextopo.hrufrac[i])
    end
end

function hortonrunoff_no_storage(flextopo::FLEXTOPO, config)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        Qhf = 0
        Shf = 0

        wbShf = flextopo.Qh[i][k] - Qhf - Shf + flextopo.Shf[i][k]

        #update with setindex
        flextopo.states_[i] = setindex(flextopo.states_[i], flextopo.states_[i][k]  + flextopo.Shf[i][k], k)
        flextopo.Qhf[i] = setindex(flextopo.Qhf[i], Qhf, k)
        flextopo.Shf[i] = setindex(flextopo.Shf[i], Shf, k)
        flextopo.wbShf[i] = setindex(flextopo.wbShf[i], wbShf, k)

        #average storage over classes
        flextopo.Shf_m[i] = sum(flextopo.Shf[i] .* flextopo.hrufrac[i])
    end
end

function hortonrunoff(flextopo::FLEXTOPO, config)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n

        Qhf = min(flextopo.Shf[i][k], flextopo.Shf[i][k] * flextopo.khf[i][k]) 
        Shf = flextopo.Shf[i][k] + flextopo.Qh[i][k] - Qhf

        wbShf = flextopo.Qh[i][k] - Qhf - Shf + flextopo.Shf[i][k]

        #update with setindex
        flextopo.states_[i] = setindex(flextopo.states_[i], flextopo.states_[i][k]  + flextopo.Shf[i][k], k)
        flextopo.Qhf[i] = setindex(flextopo.Qhf[i], Qhf, k)
        flextopo.Shf[i] = setindex(flextopo.Shf[i], Shf, k)
        flextopo.wbShf[i] = setindex(flextopo.wbShf[i], wbShf, k)

        #average storage over classes
        flextopo.Shf_m[i] = sum(flextopo.Shf[i] .* flextopo.hrufrac[i])
        
    end
end

function rootzone_no_storage(flextopo::FLEXTOPO, config)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        Er = 0
        Qperc = 0
        Qcap = 0
        Qr = max(flextopo.Qhr[i][k] + flextopo.Sr[i][k], 0) #if store not empty initial conditions
        Sr = 0

        #compute water balance Sr storage
        wbSr = flextopo.Qhr[i][k] - Er - Qr - Qperc + Qcap - Sr + flextopo.Sr[i][k]

        #update states and fluxes with setindex
        flextopo.states_[i] = setindex(flextopo.states_[i], flextopo.states_[i][k]  + flextopo.Sr[i][k], k)
        flextopo.Qr[i] = setindex(flextopo.Qr[i], Qr, k)
        flextopo.Sr[i] = setindex(flextopo.Sr[i], Sr, k)
        flextopo.Sr_over_srmax[i] = setindex(flextopo.Sr_over_srmax[i], flextopo.Sr[i][k] / flextopo.srmax[i][k], k)
        flextopo.Qcap[i] =  setindex(flextopo.Qcap[i], Qcap, k)
        flextopo.Qperc[i] = setindex(flextopo.Qperc[i], Qperc, k)
        flextopo.Er[i] = setindex(flextopo.Er[i], Er, k)
        flextopo.Ea[i] = setindex(flextopo.Ea[i], flextopo.Er[i][k] + flextopo.Ei[i][k] + flextopo.Eh[i][k], k)
        flextopo.wbSr[i] = setindex(flextopo.wbSr[i], wbSr, k)

        #average storage over classes
        flextopo.Sr_m[i] = sum(flextopo.Sr[i] .* flextopo.hrufrac[i])
        flextopo.Sr_over_srmax_m[i] = sum(flextopo.Sr_over_srmax[i] .* flextopo.hrufrac[i])
        flextopo.Er_m[i] = sum(flextopo.Er[i] .* flextopo.hrufrac[i])
        flextopo.Ea_m[i] = sum(flextopo.Ea[i] .* flextopo.hrufrac[i])
    end
end

function rootzone_storage(flextopo::FLEXTOPO, config)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        # Sr = flextopo.Sr[i] + Qhr > flextopo.srmax[i] ? flextopo.srmax[i] : flextopo.Sr[i] + Qhr
        
        #added water to the root-zone. NB: if no horton storages: Qhr is in fact Pe!! (effective precip). 
        Sr = flextopo.Sr[i][k] + flextopo.Qhr[i][k]
        # if soil is filled until max capacity, additional water runs of directly
        directrunoff = max(Sr - flextopo.srmax[i][k], 0)
        #update Sr
        Sr = Sr - directrunoff
        #net water which infiltrates in root-zone
        netinSr = flextopo.Qhr[i][k] - directrunoff

        #evaporation from the rootzone
        Er = Sr > flextopo.srmax[i][k] * flextopo.lp[i][k] ? min(Sr, flextopo.potsoilevap[i][k]) :
            min(flextopo.potsoilevap[i][k] * (Sr / (flextopo.srmax[i][k] * flextopo.lp[i][k])))
        #update storage
        Sr = Sr - Er

        #excess water from beta function of netinSr
        Qr1 = netinSr * (1 - pow(1 - Sr / flextopo.srmax[i][k], flextopo.beta[i][k]))
        #update storage
        Sr = Sr - Qr1

        #total water out of root-zone consists of directrunoff and excess water from beta function
        Qr = directrunoff + Qr1

        #percolation 
        Qperc = flextopo.perc[i][k] * Sr / flextopo.srmax[i][k]
        Sr = Sr - Qperc

        #TODO? eventueel niet een volgorde voor de verschillende processen opleggen - eerst Er, dan Qr, dan Perc, maar alles naar rato laten gebeuren zoals in code python. 

        #capillary rise 
        Qcap = min(flextopo.cap[i][k] * (1 - Sr / flextopo.srmax[i][k]), flextopo.Ss[i])
        Sr = Sr + Qcap

        #compute water balance Sr storage
        wbSr = flextopo.Qhr[i][k] - Er - Qr - Qperc + Qcap - Sr + flextopo.Sr[i][k]

        #update states and fluxes with setindex
        flextopo.states_[i] = setindex(flextopo.states_[i], flextopo.states_[i][k]  + flextopo.Sr[i][k], k)
        flextopo.Qr[i] = setindex(flextopo.Qr[i], Qr, k)
        flextopo.Sr[i] = setindex(flextopo.Sr[i], Sr, k)
        flextopo.Sr_over_srmax[i] = setindex(flextopo.Sr_over_srmax[i], flextopo.Sr[i][k] / flextopo.srmax[i][k], k)
        flextopo.Qcap[i] =  setindex(flextopo.Qcap[i], Qcap, k)
        flextopo.Qperc[i] = setindex(flextopo.Qperc[i], Qperc, k)
        flextopo.Er[i] = setindex(flextopo.Er[i], Er, k)
        flextopo.Ea[i] = setindex(flextopo.Ea[i], flextopo.Er[i][k] + flextopo.Ei[i][k] + flextopo.Eh[i][k], k)
        flextopo.wbSr[i] = setindex(flextopo.wbSr[i], wbSr, k)

        #average storage over classes
        flextopo.Sr_m[i] = sum(flextopo.Sr[i] .* flextopo.hrufrac[i])
        flextopo.Sr_over_srmax_m[i] = sum(flextopo.Sr_over_srmax[i] .* flextopo.hrufrac[i])
        flextopo.Er_m[i] = sum(flextopo.Er[i] .* flextopo.hrufrac[i])
        flextopo.Ea_m[i] = sum(flextopo.Ea[i] .* flextopo.hrufrac[i])
    end
end


function fast_no_storage(flextopo::FLEXTOPO, config)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        #make sure ds does not exceed 1 
        ds = flextopo.ds[i][k] > 1 ? 1 : flextopo.ds[i][k]
        #calc inflow to Sf
        Qrf = flextopo.Qr[i][k] * (1 - ds)
        Qf = Qrf + flextopo.Sf[i][k] #if store not empty initial conditions, make sure to empty
        Sf = 0.0

        #compute wb Sf
        wbSf = Qrf - Qf - Sf + flextopo.Sf[i][k]

        #update with setindex
        flextopo.states_[i] = setindex(flextopo.states_[i], flextopo.states_[i][k]  + flextopo.Sf[i][k], k)
        flextopo.Qrf[i] = setindex(flextopo.Qrf[i], Qrf, k)
        flextopo.Qf[i] = setindex(flextopo.Qf[i], Qf, k)
        flextopo.Sf[i] = setindex(flextopo.Sf[i], Sf, k)
        flextopo.wbSf[i] = setindex(flextopo.wbSf[i], wbSf, k)

        #average storage over classes
        flextopo.Sf_m[i] = sum(flextopo.Sf[i] .* flextopo.hrufrac[i])
    end
end

function fast_storage(flextopo::FLEXTOPO, config)
    k = flextopo.kclass[1]
    for i = 1:flextopo.n
        #TODO check to add convolution time lag? 
        
        #make sure ds does not exceed 1 
        ds = flextopo.ds[i][k] > 1 ? 1 : flextopo.ds[i][k]
        
        #split part of the outflow from the root-zone to the fast runoff (and part as preferential recharge to the slow reservoir) 
        Qrf = flextopo.Qr[i][k] * (1 - ds)

        #fast runoff 
        Qf = min(flextopo.Sf[i][k], pow(flextopo.Sf[i][k], flextopo.alfa[i][k]) * flextopo.kf[i][k])
        #update store
        Sf = flextopo.Sf[i][k] + Qrf - Qf

        #compute wb Sf
        wbSf = Qrf - Qf - Sf + flextopo.Sf[i][k]

        #update with setindex
        flextopo.states_[i] = setindex(flextopo.states_[i], flextopo.states_[i][k]  + flextopo.Sf[i][k], k)
        flextopo.Qrf[i] = setindex(flextopo.Qrf[i], Qrf, k)
        flextopo.Qf[i] = setindex(flextopo.Qf[i], Qf, k)
        flextopo.Sf[i] = setindex(flextopo.Sf[i], Sf, k)
        flextopo.wbSf[i] = setindex(flextopo.wbSf[i], wbSf, k)

        #average storage over classes
        flextopo.Sf_m[i] = sum(flextopo.Sf[i] .* flextopo.hrufrac[i])
    end
end

function slow_no_storage(flextopo::FLEXTOPO, config)
    for i = 1:flextopo.n
        #make sure ds does not exceed 1 
        pref_recharge = flextopo.Qr[i] .* min.(flextopo.ds[i], 1)
        pref_recharge_sum_classes = sum(pref_recharge .* flextopo.hrufrac[i])

        Qcap_sum_classes = sum(flextopo.Qcap[i] .* flextopo.hrufrac[i])
        Qperc_sum_classes = sum(flextopo.Qperc[i] .* flextopo.hrufrac[i])

        Qsin = pref_recharge_sum_classes + Qperc_sum_classes - Qcap_sum_classes
        Qs = Qsin + flextopo.Ss[i] # if at start store is not empty 
        Ss = 0.0

        Qftotal = sum(flextopo.Qf[i] .* flextopo.hrufrac[i]) + sum(flextopo.Qhf[i] .* flextopo.hrufrac[i])
        runoff =  max(0, Qs + Qftotal)
        
        wbSs = Qsin - Qs - Ss + flextopo.Ss[i]

        #update
        flextopo.Qftotal[i] = Qftotal
        flextopo.states_m[i] = flextopo.states_m[i] + sum(flextopo.states_[i] .* flextopo.hrufrac[i])  + flextopo.Ss[i]
        flextopo.runoff[i] = runoff
        flextopo.Qrs_m[i] = pref_recharge_sum_classes
        flextopo.Qperc_m[i] = Qperc_sum_classes
        flextopo.Qcap_m[i] = Qcap_sum_classes
        flextopo.Qs[i] = Qs
        flextopo.Ss[i] = Ss
        flextopo.wbSs[i] = wbSs
    end
end

function common_slow_storage(flextopo::FLEXTOPO, config)
    for i = 1:flextopo.n
        #make sure ds does not exceed 1 
        pref_recharge = flextopo.Qr[i] .* min.(flextopo.ds[i], 1)
        pref_recharge_sum_classes = sum(pref_recharge .* flextopo.hrufrac[i])

        Qcap_sum_classes = sum(flextopo.Qcap[i] .* flextopo.hrufrac[i])
        Qperc_sum_classes = sum(flextopo.Qperc[i] .* flextopo.hrufrac[i])

        Qsin = pref_recharge_sum_classes + Qperc_sum_classes - Qcap_sum_classes
        Ss = flextopo.Ss[i] + Qsin
        
        Qs = min(flextopo.Ss[i], flextopo.Ss[i] * flextopo.ks[i]) 
        Ss = Ss - Qs

        Qftotal = sum(flextopo.Qf[i] .* flextopo.hrufrac[i]) + sum(flextopo.Qhf[i] .* flextopo.hrufrac[i])
        runoff =  max(0, Qs + Qftotal)

        wbSs = Qsin - Qs - Ss + flextopo.Ss[i]

        #update
        flextopo.Qftotal[i] = Qftotal
        flextopo.states_m[i] = flextopo.states_m[i] + sum(flextopo.states_[i] .* flextopo.hrufrac[i])  + flextopo.Ss[i]
        flextopo.runoff[i] = runoff
        flextopo.Qrs_m[i] = pref_recharge_sum_classes
        flextopo.Qperc_m[i] = Qperc_sum_classes
        flextopo.Qcap_m[i] = Qcap_sum_classes
        flextopo.Qs[i] = Qs
        flextopo.Ss[i] = Ss
        flextopo.wbSs[i] = wbSs
    end
end

function watbal(flextopo::FLEXTOPO, config)
    for i = 1:flextopo.n
        states = flextopo.Sw[i] + flextopo.Sww[i] + sum(flextopo.Si[i] .* flextopo.hrufrac[i]) + sum(flextopo.Sh[i] .* flextopo.hrufrac[i]) + sum(flextopo.Shf[i] .* flextopo.hrufrac[i]) + sum(flextopo.Sr[i] .* flextopo.hrufrac[i]) + sum(flextopo.Sf[i] .* flextopo.hrufrac[i]) + sum(flextopo.Ss[i] .* flextopo.hrufrac[i])
        wbtot = flextopo.precipitation[i] - (sum(flextopo.Ei[i] .* flextopo.hrufrac[i]) + sum(flextopo.Eh[i] .* flextopo.hrufrac[i]) + sum(flextopo.Er[i] .* flextopo.hrufrac[i])) - flextopo.runoff[i] - states + flextopo.states_m[i]
        #update wb
        flextopo.wbtot[i] = wbtot
    end
end
