"""
    initialize_flextopo_model(config::Config)

Initial part of the FlexTopo model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_flextopo_model(config::Config)
    # unpack the paths to the NetCDF files
    static_path = input_path(config, config.input.path_static)

    reader = prepare_reader(config)
    clock = Clock(config, reader)
    Δt = clock.Δt

    reinit = get(config.model, "reinit", true)::Bool
    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_pits = get(config.model, "pits", false)::Bool

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool

    classes = get(config.model, "classes", "")
    # println(classes)
    nclass = length(classes)
    @show classes nclass
    kclass = [1]

    # dictionary of available functions for each store
    dic_function =  Dict{String, Function}(
        "snow_hbv" => snow_hbv,
        "snow_no_storage" => snow_no_storage,
        "interception_overflow" => interception_overflow,
        "interception_no_storage" => interception_no_storage,
        "hortonponding" => hortonponding,
        "hortonponding_no_storage" => hortonponding_no_storage,
        "hortonrunoff" => hortonrunoff,
        "hortonrunoff_no_storage" => hortonrunoff_no_storage,
        "rootzone_storage" => rootzone_storage,
        "rootzone_no_storage" => rootzone_no_storage,
        "fast_no_storage" => fast_no_storage,
        "fast_storage" => fast_storage,
        "slow_no_storage" => slow_no_storage,
        "common_slow_storage" => common_slow_storage,
    )

    selectSw = get(config.model, "selectSw", ["snow_hbv"])
    selectSi = get(config.model, "selectSi", ["interception_overflow"])
    selectSh = get(config.model, "selectSh", ["hortonponding_no_storage"])
    selectShf = get(config.model, "selectShf", ["hortonrunoff_no_storage"])
    selectSr = get(config.model, "selectSr", ["rootzone_storage"])
    selectSf = get(config.model, "selectSf", ["fast_storage"])
    selectSs = get(config.model, "selectSs", ["common_slow_storage"])

    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, param(config, "input.subcatchment"); allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)

    
    # # glacier parameters
    # g_tt = ncread(
    #     nc,
    #     param(config, "input.vertical.g_tt", nothing);
    #     sel = inds,
    #     defaults = 0.0,
    #     type = Float,
    #     fill = 0.0,
    # )
    # g_cfmax =
    #     ncread(
    #         nc,
    #         param(config, "input.vertical.g_cfmax", nothing);
    #         sel = inds,
    #         defaults = 3.0,
    #         type = Float,
    #         fill = 0.0,
    #     ) .* (Δt / basetimestep)
    # g_sifrac = ncread(
    #     nc,
    #     param(config, "input.vertical.g_sifrac", nothing);
    #     sel = inds,
    #     defaults = 0.001,
    #     type = Float,
    #     fill = 0.0,
    # )
    # glacierfrac = ncread(
    #     nc,
    #     param(config, "input.vertical.glacierfrac", nothing);
    #     sel = inds,
    #     defaults = 0.0,
    #     type = Float,
    #     fill = 0.0,
    # )
    # glacierstore = ncread(
    #     nc,
    #     param(config, "input.vertical.glacierstore", nothing);
    #     sel = inds,
    #     defaults = 5500.0,
    #     type = Float,
    #     fill = 0.0,
    # )

    #parameters which are not class specific
    pcorr = ncread(
        nc,
        param(config, "input.vertical.pcorr", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    ecorr = ncread(
        nc,
        param(config, "input.vertical.ecorr", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    rfcf = ncread(
        nc,
        param(config, "input.vertical.rfcf", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    sfcf = ncread(
        nc,
        param(config, "input.vertical.sfcf", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    ks = ncread(
        nc,
        param(config, "input.vertical.ks", nothing);
        sel = inds,
        defaults = 0.006,
        type = Float,
    ) .* (Δt / basetimestep)

    # #initialize parameters that differ per class
    hrufrac = zeros(Float, nclass, n)
    cfmax = zeros(Float, nclass, n)
    tt = zeros(Float, nclass, n)
    tti = zeros(Float, nclass, n)
    ttm = zeros(Float, nclass, n)
    whc = zeros(Float, nclass, n)
    cfr = zeros(Float, nclass, n)

    imax = zeros(Float, nclass, n)

    shmax = zeros(Float, nclass, n)
    khf = zeros(Float, nclass, n)
    facc0 = zeros(Float, nclass, n)
    facc1 = zeros(Float, nclass, n)
    fdec = zeros(Float, nclass, n)
    fmax = zeros(Float, nclass, n)
    shmin = zeros(Float, nclass, n)
    kmf = zeros(Float, nclass, n)

    srmax = zeros(Float, nclass, n)
    beta = zeros(Float, nclass, n)
    lp = zeros(Float, nclass, n)
    perc = zeros(Float, nclass, n)
    cap = zeros(Float, nclass, n)

    kf = zeros(Float, nclass, n)
    alfa = zeros(Float, nclass, n)

    ds = zeros(Float, nclass, n)


    #fill matrix parameter
    for (k, class) in enumerate(classes)
        hrufrac_k = ncread(
            nc,
            param(config, "input.vertical.hrufrac"*class, nothing);
            sel = inds,
            defaults = 1.0/length(classes),
            type = Float,
            # dimname = :classes,
        )
        hrufrac[k,:] = hrufrac_k

        cfmax_k =
            ncread(
                nc,
                param(config, "input.vertical.cfmax"*class, nothing);
                sel = inds,
                defaults = 3.75653,
                type = Float,
            ) .* (Δt / basetimestep)
        cfmax[k,:] = cfmax_k

        tt_k = ncread(
            nc,
            param(config, "input.vertical.tt"*class, nothing);
            sel = inds,
            defaults = -1.41934,
            type = Float,
        )
        tt[k,:] = tt_k

        tti_k = ncread(
            nc,
            param(config, "input.vertical.tti"*class, nothing);
            sel = inds,
            defaults = 1.0,
            type = Float,
        )
        tti[k,:] = tti_k

        ttm_k = ncread(
            nc,
            param(config, "input.vertical.ttm"*class, nothing);
            sel = inds,
            defaults = -1.41934,
            type = Float,
        )
        ttm[k,:] = ttm_k
        
        whc_k = ncread(
            nc,
            param(config, "input.vertical.whc"*class, nothing);
            sel = inds,
            defaults = 0.1,
            type = Float,
        )
        whc[k,:] = whc_k

        cfr_k = ncread(
            nc,
            param(config, "input.vertical.cfr"*class, nothing);
            sel = inds,
            defaults = 0.05,
            type = Float,
        )
        cfr[k,:] = cfr_k

        imax_k = ncread(
            nc,
            param(config, "input.vertical.imax"*class, nothing);
            sel = inds,
            defaults = 3.0,
            type = Float,
            # dimname = :classes,
        )
        imax[k,:] = imax_k

        shmax_k = ncread(
            nc,
            param(config, "input.vertical.shmax"*class, nothing);
            sel = inds,
            defaults = 30.0,
            type = Float,
        )
        shmax[k,:] = shmax_k

        khf_k = ncread(
            nc,
            param(config, "input.vertical.khf"*class, nothing);
            sel = inds,
            defaults = 0.5,
            type = Float,
        ) .* (Δt / basetimestep)
        khf[k,:] = khf_k

        facc0_k = ncread(
            nc,
            param(config, "input.vertical.facc0"*class, nothing);
            sel = inds,
            defaults = -3.0,
            type = Float,
        )
        facc0[k,:] = facc0_k

        facc1_k = ncread(
            nc,
            param(config, "input.vertical.facc1"*class, nothing);
            sel = inds,
            defaults = 0.0,
            type = Float,
        )
        facc1[k,:] = facc1_k

        fdec_k = ncread(
            nc,
            param(config, "input.vertical.fdec"*class, nothing);
            sel = inds,
            defaults = 0.2,
            type = Float,
        )
        fdec[k,:] = fdec_k

        fmax_k = ncread(
            nc,
            param(config, "input.vertical.fmax"*class, nothing);
            sel = inds,
            defaults = 2.0,
            type = Float,
        ) .* (Δt / basetimestep)
        fmax[k,:] = fmax_k

        shmin_k = ncread(
            nc,
            param(config, "input.vertical.shmin"*class, nothing);
            sel = inds,
            defaults = 0.2,
            type = Float,
        )
        shmin[k,:] = shmin_k

        kmf_k = ncread(
            nc,
            param(config, "input.vertical.kmf"*class, nothing);
            sel = inds,
            defaults = 1.0,
            type = Float,
        )
        kmf[k,:] = kmf_k

        srmax_k = ncread(
            nc,
            param(config, "input.vertical.srmax"*class, nothing);
            sel = inds,
            defaults = 260.0,
            type = Float,
        )
        srmax[k,:] = srmax_k

        beta_k = ncread(
            nc,
            param(config, "input.vertical.beta"*class, nothing);
            sel = inds,
            defaults = 0.3,
            type = Float,
        )
        beta[k,:] = beta_k

        lp_k = ncread(
            nc,
            param(config, "input.vertical.lp"*class, nothing);
            sel = inds,
            defaults = 0.3,
            type = Float,
        )
        lp[k,:] = lp_k

        perc_k = ncread(
            nc,
            param(config, "input.vertical.perc"*class, nothing);
            sel = inds,
            defaults = 0.30,
            type = Float,
        ) .* (Δt / basetimestep)
        perc[k,:] = perc_k
        
        cap_k = ncread(
            nc,
            param(config, "input.vertical.cap"*class, nothing);
            sel = inds,
            defaults = 0.20,
            type = Float,
        ) .* (Δt / basetimestep)
        cap[k,:] = cap_k

        kf_k = ncread(
            nc,
            param(config, "input.vertical.kf"*class, nothing);
            sel = inds,
            defaults = 0.1,
            type = Float,
        ) .* (Δt / basetimestep)
        kf[k,:] = kf_k

        alfa_k = ncread(
            nc, 
            param(config, "input.vertical.alfa"*class, nothing);
            sel = inds,
            defaults = 1.3,
            type = Float,
        )
        alfa[k,:] = alfa_k

        ds_k = ncread(
            nc,
            param(config, "input.vertical.ds"*class, nothing);
            sel = inds,
            defaults = 0.2,
            type = Float,
        )
        ds[k,:] = ds_k

    end

    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))


    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)


    # dummy_ = zeros(Float, nclass, n)
    # dummy = svectorscopy(dummy_, Val{nclass}())

    Sw = zeros(Float, nclass, n)
    Sww = zeros(Float, nclass, n)
    Si = zeros(Float, nclass, n)
    Sh = zeros(Float, nclass, n)
    Shf = zeros(Float, nclass, n)
    Sr = zeros(Float, nclass, n)
    Sf = zeros(Float, nclass, n)
    Ss = zeros(Float, nclass, n)
    Sr_over_srmax = zeros(Float, nclass, n)   

    states_  = fill(mv, nclass, n)

    rainfallplusmelt = fill(mv, nclass, n)
    snowfall = fill(mv, nclass, n)
    snowmelt = fill(mv, nclass, n)
    potsoilevap = fill(mv, nclass, n)
    Ei = fill(mv, nclass, n)
    Pe = fill(mv, nclass, n)
    Eh = fill(mv, nclass, n)
    Er = fill(mv, nclass, n)
    Qh = fill(mv, nclass, n)
    Qhr = fill(mv, nclass, n)
    facc = zeros(Float, nclass, n)
    Qhf = fill(mv, nclass, n)
    Qr = fill(mv, nclass, n)
    Qrf = fill(mv, nclass, n)
    # Qrs = fill(mv, nclass, n)
    Qf = fill(mv, nclass, n)
    Ea = fill(mv, nclass, n)
    # directrunoff = fill(mv, nclass, n)
    Qperc = fill(mv, nclass, n)
    Qcap = fill(mv, nclass, n)

    wbSw = fill(mv, nclass, n)
    wbSi = fill(mv, nclass, n)
    wbSh = fill(mv, nclass, n)
    wbShf = fill(mv, nclass, n)
    wbSr = fill(mv, nclass, n)
    wbSf = fill(mv, nclass, n)    
    


    # store = FLEXTOPO{Float, nclass}(
    # store = InterceptionStorage{Float}(
    flextopo = FLEXTOPO{Float, nclass}(
        Δt = Float(tosecond(Δt)),
        nclass = nclass,
        n = n,
        dic_function = dic_function,
        kclass = kclass,
        classes = classes,
        selectSw = selectSw,
        selectSi = selectSi,
        selectSh = selectSh,
        selectShf = selectShf,
        selectSr = selectSr,
        selectSf = selectSf,
        selectSs = selectSs,

        hrufrac = svectorscopy(hrufrac, Val{nclass}()),
        pcorr = pcorr,
        ecorr = ecorr,
        # snow
        tti = svectorscopy(tti, Val{nclass}()),
        tt = svectorscopy(tt, Val{nclass}()),
        ttm = svectorscopy(ttm, Val{nclass}()),
        cfmax = svectorscopy(cfmax, Val{nclass}()),
        whc = svectorscopy(whc, Val{nclass}()),
        cfr = svectorscopy(cfr, Val{nclass}()),
        rfcf = rfcf,
        sfcf = sfcf,
        #interception 
        # imax = imax,
        imax = svectorscopy(imax, Val{nclass}()),
        #horton
        shmax = svectorscopy(shmax, Val{nclass}()),
        khf = svectorscopy(khf, Val{nclass}()),
        facc0 = svectorscopy(facc0, Val{nclass}()),
        facc1 = svectorscopy(facc1, Val{nclass}()),
        fdec = svectorscopy(fdec, Val{nclass}()),
        fmax = svectorscopy(fmax, Val{nclass}()),
        shmin = svectorscopy(shmin, Val{nclass}()),
        kmf = svectorscopy(kmf, Val{nclass}()),
        #root zone
        srmax = svectorscopy(srmax, Val{nclass}()),
        lp = svectorscopy(lp, Val{nclass}()),
        beta = svectorscopy(beta, Val{nclass}()),
        perc = svectorscopy(perc, Val{nclass}()),
        cap = svectorscopy(cap, Val{nclass}()),
        #fast
        ds = svectorscopy(ds, Val{nclass}()),
        alfa = svectorscopy(alfa, Val{nclass}()),
        kf = svectorscopy(kf, Val{nclass}()),
        #slow
        ks = ks,

        # kquickflow = set_kquickflow ? kquickflow :
        #              pow.(khq, 1.0 .+ alphanl) .* pow.(hq, -alphanl),
        
        # # glacier parameters
        # g_tt = g_tt,
        # g_sifrac = g_sifrac,
        # g_cfmax = g_cfmax,
        # glacierstore = glacierstore,
        # glacierfrac = glacierfrac,
        
        # # default (cold) states:
        # Sw = zeros(Float, n), 
        # Sw = zero(dummy),
        Sw = svectorscopy(Sw, Val{nclass}()), 
        Sww = svectorscopy(Sww, Val{nclass}()), 
        # Si = zero(dummy),
        Si = svectorscopy(Si, Val{nclass}()), 
        Sh = svectorscopy(Sh, Val{nclass}()), 
        Shf = svectorscopy(Shf, Val{nclass}()), 
        # Sr = copy(srmax),
        # Sf = 0.2 .* srmax,
        Sr = svectorscopy(srmax, Val{nclass}()), 
        Sf = 0.2 .* svectorscopy(srmax, Val{nclass}()), 
        # Sr = svectorscopy(Sr, Val{nclass}()), 
        # Sf = svectorscopy(Sf, Val{nclass}()), 
        Sr_over_srmax = svectorscopy(Sr_over_srmax, Val{nclass}()), 
        Ss = 1.0 ./ (3.0 .* ks),
        #states previous time step
        states_m = fill(mv, n),
        states_ = svectorscopy(states_, Val{nclass}()), 
        #states averaged over all classes
        Sw_m = zeros(Float, n),
        Sww_m = zeros(Float, n),
        Si_m = zeros(Float, n),
        Sh_m = zeros(Float, n),
        Shf_m = zeros(Float, n),
        Sr_m = zeros(Float, n),
        Sf_m = zeros(Float, n),
        Sr_over_srmax_m = zeros(Float, n),

        # variables:
        precipitation = fill(mv, n),
        temperature = fill(mv, n),
        potential_evaporation = fill(mv, n),
        rainfallplusmelt = svectorscopy(rainfallplusmelt, Val{nclass}()), 
        snowfall = svectorscopy(snowfall, Val{nclass}()), 
        snowmelt = svectorscopy(snowmelt, Val{nclass}()), 
        potsoilevap = svectorscopy(potsoilevap, Val{nclass}()), 
        # Pe = fill(mv, n),
        Pe = svectorscopy(Pe, Val{nclass}()), 
        # Pe = zero(dummy),
        # Ei = fill(mv, n),
        Ei = svectorscopy(Ei, Val{nclass}()), 
        # Ei = zero(dummy),
        Eh = svectorscopy(Eh, Val{nclass}()), 
        Er = svectorscopy(Er, Val{nclass}()), 
        Qh = svectorscopy(Qh, Val{nclass}()), 
        Qhr = svectorscopy(Qhr, Val{nclass}()), 
        facc = svectorscopy(facc, Val{nclass}()), 
        Qhf = svectorscopy(Qhf, Val{nclass}()), 
        Qr = svectorscopy(Qr, Val{nclass}()), 
        Qrf = svectorscopy(Qrf, Val{nclass}()), 
        # Qrs = svectorscopy(Qrs, Val{nclass}()), 
        Qrs_m = fill(mv, n),
        Qf = svectorscopy(Qf, Val{nclass}()), 
        Ea = svectorscopy(Ea, Val{nclass}()), 
        # directrunoff = svectorscopy(potsoilevap, Val{nclass}()), 
        # hbv_seepage = svectorscopy(potsoilevap, Val{nclass}()), 
        Qperc = svectorscopy(Qperc, Val{nclass}()), 
        Qcap = svectorscopy(Qcap, Val{nclass}()), 
        Qperc_m = fill(mv, n),
        Qcap_m = fill(mv, n),
        # combined for the classes
        Ea_m = fill(mv, n),
        Ei_m = fill(mv, n),
        Eh_m = fill(mv, n),
        Er_m = fill(mv, n),
        Qs = fill(mv, n),
        Qftotal = fill(mv, n),
        runoff = fill(mv, n),


        wbSw = svectorscopy(wbSw, Val{nclass}()), 
        wbSi = svectorscopy(wbSi, Val{nclass}()), 
        wbSh = svectorscopy(wbSh, Val{nclass}()), 
        wbShf = svectorscopy(wbShf, Val{nclass}()), 
        wbSr = svectorscopy(wbSr, Val{nclass}()), 
        wbSf = svectorscopy(wbSf, Val{nclass}()), 
        wbSs = fill(mv, n), 
        wbtot = fill(mv, n),  
    )

    modelsize_2d = size(subcatch_2d)
    river_2d = ncread(nc, param(config, "input.river_location"); type = Bool, fill = false)
    river = river_2d[inds]
    riverwidth_2d =
        ncread(nc, param(config, "input.lateral.river.width"); type = Float, fill = 0)
    riverwidth = riverwidth_2d[inds]
    riverlength_2d =
        ncread(nc, param(config, "input.lateral.river.length"); type = Float, fill = 0)
    riverlength = riverlength_2d[inds]

    inds_riv, rev_inds_riv = active_indices(river_2d, 0)
    nriv = length(inds_riv)

    # reservoirs
    pits = zeros(Bool, modelsize_2d)
    if do_reservoirs
        reservoirs, resindex, reservoir, pits =
            initialize_simple_reservoir(config, nc, inds_riv, nriv, pits, tosecond(Δt))
    else
        reservoir = ()
        reservoirs = nothing
        resindex = fill(0, nriv)    
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits = initialize_natural_lake(
            config,
            nc,
            inds_riv,
            nriv,
            pits,
            tosecond(Δt),
        )
    else
        lake = ()
        lakes = nothing
        lakeindex = fill(0, nriv)
    end

    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d = ncread(nc, param(config, "input.pits"); type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, inds)
    end

    βₗ = ncread(nc, param(config, "input.lateral.land.slope"); sel = inds, type = Float)
    clamp!(βₗ, 0.00001, Inf)

    dl = map(detdrainlength, ldd, xl, yl)
    dw = (xl .* yl) ./ dl
    olf = initialize_surfaceflow_land(
        nc,
        config,
        inds;
        sl = βₗ,
        dl = dl,
        width = map(det_surfacewidth, dw, riverwidth, river),
        wb_pit = pits[inds],
        iterate = kinwave_it,
        tstep = kw_land_tstep,
        Δt = Δt,
    )

    graph = flowgraph(ldd, inds, pcr_dir)

    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]

    ldd_riv = ldd_2d[inds_riv]
    if do_pits
        ldd_riv = set_pit_ldd(pits_2d, ldd_riv, inds_riv)
    end
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, βₗ, n)

    rf = initialize_surfaceflow_river(
        nc,
        config,
        inds_riv;
        dl = riverlength,
        width = riverwidth,
        wb_pit = pits[inds_riv],
        reservoir_index = resindex,
        reservoir = reservoirs,
        lake_index = lakeindex,
        lake = lakes,
        river = river,
        iterate = kinwave_it,
        tstep = kw_river_tstep,
        Δt = Δt,
    )

    # setup subdomains for the land and river kinematic wave domain, if nthreads = 1
    # subdomain is equal to the complete domain
    toposort = topological_sort_by_dfs(graph)
    toposort_riv = topological_sort_by_dfs(graph_riv)
    index_pit_land = findall(x -> x == 5, ldd)
    index_pit_river = findall(x -> x == 5, ldd_riv)
    subbas_order, indices_subbas, topo_subbas =
        kinwave_set_subdomains(config, graph, toposort, index_pit_land)
    subriv_order, indices_subriv, topo_subriv =
        kinwave_set_subdomains(config, graph_riv, toposort_riv, index_pit_river)

    modelmap = (vertical = flextopo, lateral = (land = olf, river = rf))
    # modelmap = (vertical = store, lateral = (land = olf, river = rf))
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
    )
    writer = prepare_writer(
        config, 
        reader, 
        modelmap, 
        indices_reverse, 
        x_nc, 
        y_nc, 
        nc,
        #TODO: add:
        # nclass = flextopo.nclass,
        )
    close(nc)

    # for each domain save:
    # - the directed acyclic graph (graph),
    # - the traversion order (order),
    # - upstream_nodes,
    # - subdomains for the kinematic wave domains for parallel execution (execution order of
    #   subbasins (subdomain_order), traversion order per subbasin (topo_subdomain) and
    #   Vector indices per subbasin matching the traversion order of the complete domain
    #   (indices_subdomain)) 
    # - the indices that map it back to the two dimensional grid (indices)

    # for the land domain the x and y length [m] of the grid cells are stored
    # for reservoirs and lakes indices information is available from the initialization
    # functions
    land = (
        graph = graph,
        upstream_nodes = filter_upsteam_nodes(graph, olf.wb_pit),
        subdomain_order = subbas_order,
        topo_subdomain = topo_subbas,
        indices_subdomain = indices_subbas,
        order = toposort,
        indices = inds,
        reverse_indices = rev_inds,
        xl = xl,
        yl = yl,
    )
    river = (
        graph = graph_riv,
        upstream_nodes = filter_upsteam_nodes(graph_riv, rf.wb_pit),
        subdomain_order = subriv_order,
        topo_subdomain = topo_subriv,
        indices_subdomain = indices_subriv,
        order = toposort_riv,
        indices = inds_riv,
        reverse_indices = rev_inds_riv,
    )

    model = Model(
        config,
        (; land, river, reservoir, lake, index_river, frac_toriver),
        (subsurface = nothing, land = olf, river = rf),
        # store,
        flextopo,
        clock,
        reader,
        writer,
        FlextopoModel(),
    )

    # read and set states in model object if reinit=true
    if reinit == false
        instate_path = input_path(config, config.state.path_input)
        state_ncnames = ncnames(config.state)
        set_states(instate_path, model, state_ncnames; type = Float)
        # update kinematic wave volume for river and land domain
        @unpack lateral = model
        lateral.land.volume .= lateral.land.h .* lateral.land.width .* lateral.land.dl
        lateral.river.volume .= lateral.river.h .* lateral.river.width .* lateral.river.dl

        if do_lakes
            # storage must be re-initialized after loading the state with the current
            # waterlevel otherwise the storage will be based on the initial water level
            lakes.storage .=
                initialize_storage(lakes.storfunc, lakes.area, lakes.waterlevel, lakes.sh)
        end
    end

    return model
end

function update(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:FlextopoModel}
    @unpack lateral, vertical, network, clock, config = model

    inds_riv = network.index_river

    for (k, class) in enumerate(vertical.classes)
        vertical.kclass[1] = k
        
        #SNOW
        # snow_no_storage_k(vertical, config)
        vertical.dic_function[vertical.selectSw[k]](vertical, config)

        # TODO add mass wasting and glaciers?

        #INTERCEPTION
        vertical.dic_function[vertical.selectSi[k]](vertical, config)
        # interception_overflow(vertical, config)

        #HORTON
        vertical.dic_function[vertical.selectSh[k]](vertical, config)
        vertical.dic_function[vertical.selectShf[k]](vertical, config)

        #ROOT-ZONE
        vertical.dic_function[vertical.selectSr[k]](vertical, config)

        #FAST
        vertical.dic_function[vertical.selectSf[k]](vertical, config)

    end

    #COMMON SLOW
    vertical.dic_function[vertical.selectSs[1]](vertical, config)

    # WAT BAL TODO
    watbal(vertical, config)

    
    # # lateral snow transport
    # if get(config.model, "masswasting", false)::Bool
    #     lateral_snow_transport!(
    #         vertical.snow,
    #         vertical.snowwater,
    #         lateral.land.sl,
    #         network.land,
    #     )
    # end


    surface_routing(model)

    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
