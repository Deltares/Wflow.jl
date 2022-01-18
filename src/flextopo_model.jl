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
    nclass = length(classes)
    @show classes nclass
    kclass = [1]    #needed to initialize

    # dictionary of available functions for each store
    dic_function = Dict{String,Function}(
        "common_snow_hbv" => common_snow_hbv,
        "common_snow_no_storage" => common_snow_no_storage,
        "common_glaciers" => common_glaciers,
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

    selectSw = get(config.model, "selectSw", ["common_snow_hbv"])
    selectSi = get(config.model, "selectSi", ["interception_overflow"])
    selectSh = get(config.model, "selectSh", ["hortonponding_no_storage"])
    selectShf = get(config.model, "selectShf", ["hortonrunoff_no_storage"])
    selectSr = get(config.model, "selectSr", ["rootzone_storage"])
    selectSf = get(config.model, "selectSf", ["fast_storage"])
    selectSs = get(config.model, "selectSs", ["common_slow_storage"])

    nc = NCDataset(static_path)

    subcatch_2d =
        ncread(nc, config.input, "subcatchment"; optional = false, allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)


    # glacier parameters
    g_tt = ncread(
        nc,
        config.input,
        "vertical.g_tt";
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    g_cfmax =
        ncread(
            nc,
            config.input,
            "vertical.g_cfmax";
            sel = inds,
            defaults = 3.0,
            type = Float,
            fill = 0.0,
        ) .* (Δt / basetimestep)

    g_sifrac = ncread(
        nc,
        config.input,
        "vertical.g_sifrac";
        sel = inds,
        defaults = 0.001,
        type = Float,
        fill = 0.0,
    )
    glacierfrac = ncread(
        nc,
        config.input,
        "vertical.glacierfrac";
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    glacierstore = ncread(
        nc,
        config.input,
        "vertical.glacierstore";
        sel = inds,
        defaults = 5500.0,
        type = Float,
        fill = 0.0,
    )

    #snow param (single class)
    cfmax =
        ncread(
            nc,
            config.input,
            "vertical.cfmax";
            sel = inds,
            defaults = 3.75653,
            type = Float,
        ) .* (Δt / basetimestep)

    tt = ncread(
        nc,
        config.input,
        "vertical.tt";
        sel = inds,
        defaults = -1.41934,
        type = Float,
    )

    tti = ncread(nc, config.input, "vertical.tti"; sel = inds, defaults = 1.0, type = Float)

    ttm = ncread(
        nc,
        config.input,
        "vertical.ttm";
        sel = inds,
        defaults = -1.41934,
        type = Float,
    )

    whc = ncread(nc, config.input, "vertical.whc"; sel = inds, defaults = 0.1, type = Float)

    cfr =
        ncread(nc, config.input, "vertical.cfr"; sel = inds, defaults = 0.05, type = Float)

    #parameters which are not class specific
    pcorr =
        ncread(nc, config.input, "vertical.pcorr"; sel = inds, defaults = 1.0, type = Float)
    ecorr =
        ncread(nc, config.input, "vertical.ecorr"; sel = inds, defaults = 1.0, type = Float)
    rfcf =
        ncread(nc, config.input, "vertical.rfcf"; sel = inds, defaults = 1.0, type = Float)
    sfcf =
        ncread(nc, config.input, "vertical.sfcf"; sel = inds, defaults = 1.0, type = Float)
    ks =
        ncread(
            nc,
            config.input,
            "vertical.ks";
            sel = inds,
            defaults = 0.006,
            type = Float,
        ) .* (Δt / basetimestep)

    #initialize parameters that differ per class
    hrufrac = ncread(
        nc,
        config.input,
        "vertical.hrufrac";
        sel = inds,
        defaults = 1.0 / length(classes),
        type = Float,
        dimname = :classes,
    )
    imax = ncread(
        nc,
        config.input,
        "vertical.imax";
        sel = inds,
        defaults = 3.0,
        type = Float,
        dimname = :classes,
    )
    shmax = ncread(
        nc,
        config.input,
        "vertical.shmax";
        sel = inds,
        defaults = 30.0,
        type = Float,
        dimname = :classes,
    )
    khf =
        ncread(
            nc,
            config.input,
            "vertical.khf";
            sel = inds,
            defaults = 0.5,
            type = Float,
            dimname = :classes,
        ) .* (Δt / basetimestep)
    facc0 = ncread(
        nc,
        config.input,
        "vertical.facc0";
        sel = inds,
        defaults = -3.0,
        type = Float,
        dimname = :classes,
    )
    facc1 = ncread(
        nc,
        config.input,
        "vertical.facc1";
        sel = inds,
        defaults = 0.0,
        type = Float,
        dimname = :classes,
    )
    fdec = ncread(
        nc,
        config.input,
        "vertical.fdec";
        sel = inds,
        defaults = 0.2,
        type = Float,
        dimname = :classes,
    )
    fmax =
        ncread(
            nc,
            config.input,
            "vertical.fmax";
            sel = inds,
            defaults = 2.0,
            type = Float,
            dimname = :classes,
        ) .* (Δt / basetimestep)
    shmin = ncread(
        nc,
        config.input,
        "vertical.shmin";
        sel = inds,
        defaults = 0.2,
        type = Float,
        dimname = :classes,
    )
    kmf = ncread(
        nc,
        config.input,
        "vertical.kmf";
        sel = inds,
        defaults = 1.0,
        type = Float,
        dimname = :classes,
    )
    srmax = ncread(
        nc,
        config.input,
        "vertical.srmax";
        sel = inds,
        defaults = 260.0,
        type = Float,
        dimname = :classes,
    )
    beta = ncread(
        nc,
        config.input,
        "vertical.beta";
        sel = inds,
        defaults = 0.3,
        type = Float,
        dimname = :classes,
    )
    lp = ncread(
        nc,
        config.input,
        "vertical.lp";
        sel = inds,
        defaults = 0.3,
        type = Float,
        dimname = :classes,
    )
    perc =
        ncread(
            nc,
            config.input,
            "vertical.perc";
            sel = inds,
            defaults = 0.30,
            type = Float,
            dimname = :classes,
        ) .* (Δt / basetimestep)
    cap =
        ncread(
            nc,
            config.input,
            "vertical.cap";
            sel = inds,
            defaults = 0.20,
            type = Float,
            dimname = :classes,
        ) .* (Δt / basetimestep)
    kf =
        ncread(
            nc,
            config.input,
            "vertical.kf";
            sel = inds,
            defaults = 0.1,
            type = Float,
            dimname = :classes,
        ) .* (Δt / basetimestep)
    alfa = ncread(
        nc,
        config.input,
        "vertical.alfa";
        sel = inds,
        defaults = 1.3,
        type = Float,
        dimname = :classes,
    )
    ds = ncread(
        nc,
        config.input,
        "vertical.ds";
        sel = inds,
        defaults = 0.2,
        type = Float,
        dimname = :classes,
    )

    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))


    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)


    # dummy_ = zeros(Float, nclass, n)
    # dummy = svectorscopy(dummy_, Val{nclass}())

    # Sw = zeros(Float, nclass, n)
    # Sww = zeros(Float, nclass, n)
    Si = zeros(Float, nclass, n)
    Sh = zeros(Float, nclass, n)
    Shf = zeros(Float, nclass, n)
    Sr = zeros(Float, nclass, n)
    Sf = zeros(Float, nclass, n)
    Ss = zeros(Float, nclass, n)
    Sr_over_srmax = zeros(Float, nclass, n)

    states_ = fill(mv, nclass, n)

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

    wbSi = fill(mv, nclass, n)
    wbSh = fill(mv, nclass, n)
    wbShf = fill(mv, nclass, n)
    wbSr = fill(mv, nclass, n)
    wbSf = fill(mv, nclass, n)



    # store = FLEXTOPO{Float, nclass}(
    # store = InterceptionStorage{Float}(
    flextopo = FLEXTOPO{Float,nclass}(
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
        #glaciers
        g_tt = g_tt,
        g_cfmax = g_cfmax,
        g_sifrac = g_sifrac,
        glacierfrac = glacierfrac,
        glacierstore = glacierstore,
        # snow
        tti = tti,
        tt = tt,
        ttm = ttm,
        cfmax = cfmax,
        whc = whc,
        cfr = cfr,
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

        # # default (cold) states:
        Sw = zeros(Float, n),
        Sww = zeros(Float, n),
        # Sw = zero(dummy),
        # Sw = svectorscopy(Sw, Val{nclass}()), 
        # Sww = svectorscopy(Sww, Val{nclass}()), 
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
        # Sw_m = zeros(Float, n),
        # Sww_m = zeros(Float, n),
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
        rainfallplusmelt = fill(mv, n),
        snowfall = fill(mv, n),
        snowmelt = fill(mv, n),
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
        wbSw = fill(mv, n),
        wbSi = svectorscopy(wbSi, Val{nclass}()),
        wbSh = svectorscopy(wbSh, Val{nclass}()),
        wbShf = svectorscopy(wbShf, Val{nclass}()),
        wbSr = svectorscopy(wbSr, Val{nclass}()),
        wbSf = svectorscopy(wbSf, Val{nclass}()),
        wbSs = fill(mv, n),
        wbtot = fill(mv, n),
    )

    modelsize_2d = size(subcatch_2d)
    river_2d = ncread(
        nc,
        config.input,
        "river_location";
        optional = false,
        type = Bool,
        fill = false,
    )
    river = river_2d[inds]
    riverwidth_2d = ncread(
        nc,
        config.input,
        "lateral.river.width";
        optional = false,
        type = Float,
        fill = 0,
    )
    riverwidth = riverwidth_2d[inds]
    riverlength_2d = ncread(
        nc,
        config.input,
        "lateral.river.length";
        optional = false,
        type = Float,
        fill = 0,
    )
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
        lakes, lakeindex, lake, pits =
            initialize_natural_lake(config, nc, inds_riv, nriv, pits, tosecond(Δt))
    else
        lake = ()
        lakes = nothing
        lakeindex = fill(0, nriv)
    end

    ldd_2d = ncread(nc, config.input, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d =
            ncread(nc, config.input, "pits"; optional = false, type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, inds)
    end

    βₗ = ncread(
        nc,
        config.input,
        "lateral.land.slope";
        optional = false,
        sel = inds,
        type = Float,
    )
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
        extra_dim = (name = "classes", value = flextopo.classes),
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
        set_states(instate_path, model, state_ncnames; type = Float, dimname = :classes)
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

    #COMMON SNOW
    vertical.dic_function[vertical.selectSw[1]](vertical, config)

    #lateral snow transport
    if get(config.model, "masswasting", false)::Bool
        lateral_snow_transport!(vertical.Sw, vertical.Sww, lateral.land.sl, network.land)
    end

    if get(config.model, "glacier", false)::Bool
        common_glaciers(vertical, config)
    end

    for (k, class) in enumerate(vertical.classes)
        vertical.kclass[1] = k

        #INTERCEPTION
        vertical.dic_function[vertical.selectSi[k]](vertical, config)

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

    # WAT BAL 
    watbal(vertical, config)

    surface_routing(model)

    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
