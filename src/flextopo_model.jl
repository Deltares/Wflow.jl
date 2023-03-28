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
    @info "Classes are set to names `$(join(classes,", "))` and have size `$nclass`."
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

    select_snow = get(config.model, "select_snow", ["common_snow_hbv"])
    select_interception =
        get(config.model, "select_interception", ["interception_overflow"])
    select_hortonponding =
        get(config.model, "select_hortonponding", ["hortonponding_no_storage"])
    select_hortonrunoff =
        get(config.model, "select_hortonrunoff", ["hortonrunoff_no_storage"])
    select_rootzone = get(config.model, "select_rootzone", ["rootzone_storage"])
    select_fast = get(config.model, "select_fast", ["fast_storage"])
    select_slow = get(config.model, "select_slow", ["common_slow_storage"])

    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, config, "subcatchment"; optional = false, allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)


    # glacier parameters
    g_tt = ncread(
        nc,
        config,
        "vertical.g_tt";
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    g_cfmax =
        ncread(
            nc,
            config,
            "vertical.g_cfmax";
            sel = inds,
            defaults = 3.0,
            type = Float,
            fill = 0.0,
        ) .* (Δt / basetimestep)

    g_sifrac =
        ncread(
            nc,
            config,
            "vertical.g_sifrac";
            sel = inds,
            defaults = 0.001,
            type = Float,
            fill = 0.0,
        ) .* (Δt / basetimestep)
    glacierfrac = ncread(
        nc,
        config,
        "vertical.glacierfrac";
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    glacierstore = ncread(
        nc,
        config,
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
            config,
            "vertical.cfmax";
            sel = inds,
            defaults = 3.75653,
            type = Float,
        ) .* (Δt / basetimestep)

    tt = ncread(nc, config, "vertical.tt"; sel = inds, defaults = -1.41934, type = Float)

    tti = ncread(nc, config, "vertical.tti"; sel = inds, defaults = 1.0, type = Float)

    ttm = ncread(nc, config, "vertical.ttm"; sel = inds, defaults = -1.41934, type = Float)

    whc = ncread(nc, config, "vertical.whc"; sel = inds, defaults = 0.1, type = Float)

    cfr = ncread(nc, config, "vertical.cfr"; sel = inds, defaults = 0.05, type = Float)

    #parameters which are not class specific
    pcorr = ncread(nc, config, "vertical.pcorr"; sel = inds, defaults = 1.0, type = Float)
    ecorr = ncread(nc, config, "vertical.ecorr"; sel = inds, defaults = 1.0, type = Float)
    rfcf = ncread(nc, config, "vertical.rfcf"; sel = inds, defaults = 1.0, type = Float)
    sfcf = ncread(nc, config, "vertical.sfcf"; sel = inds, defaults = 1.0, type = Float)
    ks =
        ncread(nc, config, "vertical.ks"; sel = inds, defaults = 0.006, type = Float) .*
        (Δt / basetimestep)

    #initialize parameters that differ per class
    hrufrac = ncread(
        nc,
        config,
        "vertical.hrufrac";
        sel = inds,
        defaults = 1.0 / length(classes),
        type = Float,
        dimname = :classes,
    )
    imax = ncread(
        nc,
        config,
        "vertical.imax";
        sel = inds,
        defaults = 3.0,
        type = Float,
        dimname = :classes,
    )
    shmax = ncread(
        nc,
        config,
        "vertical.shmax";
        sel = inds,
        defaults = 30.0,
        type = Float,
        dimname = :classes,
    )
    khf =
        ncread(
            nc,
            config,
            "vertical.khf";
            sel = inds,
            defaults = 0.5,
            type = Float,
            dimname = :classes,
        ) .* (Δt / basetimestep)
    facc0 = ncread(
        nc,
        config,
        "vertical.facc0";
        sel = inds,
        defaults = -3.0,
        type = Float,
        dimname = :classes,
    )
    facc1 = ncread(
        nc,
        config,
        "vertical.facc1";
        sel = inds,
        defaults = 0.0,
        type = Float,
        dimname = :classes,
    )
    fdec = ncread(
        nc,
        config,
        "vertical.fdec";
        sel = inds,
        defaults = 0.2,
        type = Float,
        dimname = :classes,
    )
    fmax =
        ncread(
            nc,
            config,
            "vertical.fmax";
            sel = inds,
            defaults = 2.0,
            type = Float,
            dimname = :classes,
        ) .* (Δt / basetimestep)
    shmin = ncread(
        nc,
        config,
        "vertical.shmin";
        sel = inds,
        defaults = 0.2,
        type = Float,
        dimname = :classes,
    )
    kmf = ncread(
        nc,
        config,
        "vertical.kmf";
        sel = inds,
        defaults = 1.0,
        type = Float,
        dimname = :classes,
    )
    srmax = ncread(
        nc,
        config,
        "vertical.srmax";
        sel = inds,
        defaults = 260.0,
        type = Float,
        dimname = :classes,
    )
    beta = ncread(
        nc,
        config,
        "vertical.beta";
        sel = inds,
        defaults = 0.3,
        type = Float,
        dimname = :classes,
    )
    lp = ncread(
        nc,
        config,
        "vertical.lp";
        sel = inds,
        defaults = 0.3,
        type = Float,
        dimname = :classes,
    )
    perc =
        ncread(
            nc,
            config,
            "vertical.perc";
            sel = inds,
            defaults = 0.30,
            type = Float,
            dimname = :classes,
        ) .* (Δt / basetimestep)
    cap =
        ncread(
            nc,
            config,
            "vertical.cap";
            sel = inds,
            defaults = 0.20,
            type = Float,
            dimname = :classes,
        ) .* (Δt / basetimestep)
    kf =
        ncread(
            nc,
            config,
            "vertical.kf";
            sel = inds,
            defaults = 0.1,
            type = Float,
            dimname = :classes,
        ) .* (Δt / basetimestep)
    alfa = ncread(
        nc,
        config,
        "vertical.alfa";
        sel = inds,
        defaults = 1.0,
        type = Float,
        dimname = :classes,
    )
    ds = ncread(
        nc,
        config,
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


    interceptionstorage = zeros(Float, nclass, n)
    hortonpondingstorage = zeros(Float, nclass, n)
    hortonrunoffstorage = zeros(Float, nclass, n)
    srootzone_over_srmax = zeros(Float, nclass, n) .+ 1.0

    states_ = fill(mv, nclass, n)

    potsoilevap = fill(mv, nclass, n)
    intevap = fill(mv, nclass, n)
    precipeffective = fill(mv, nclass, n)
    hortonevap = fill(mv, nclass, n)
    rootevap = fill(mv, nclass, n)
    qhortonpond = fill(mv, nclass, n)
    qhortonrootzone = fill(mv, nclass, n)
    facc = zeros(Float, nclass, n)
    qhortonrun = fill(mv, nclass, n)
    qrootzone = fill(mv, nclass, n)
    qrootzonefast = fill(mv, nclass, n)
    qfast = fill(mv, nclass, n)
    actevap = fill(mv, nclass, n)
    qpercolation = fill(mv, nclass, n)
    qcapillary = fill(mv, nclass, n)

    wb_interception = fill(mv, nclass, n)
    wb_hortonponding = fill(mv, nclass, n)
    wb_hortonrunoff = fill(mv, nclass, n)
    wb_rootzone = fill(mv, nclass, n)
    wb_fast = fill(mv, nclass, n)


    flextopo = FLEXTOPO{Float,nclass}(
        Δt = Float(tosecond(Δt)),
        nclass = nclass,
        n = n,
        dic_function = dic_function,
        kclass = kclass,
        classes = classes,
        select_snow = select_snow,
        select_interception = select_interception,
        select_hortonponding = select_hortonponding,
        select_hortonrunoff = select_hortonrunoff,
        select_rootzone = select_rootzone,
        select_fast = select_fast,
        select_slow = select_slow,
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

        # # default (cold) states:
        snow = zeros(Float, n),
        snowwater = zeros(Float, n),
        interceptionstorage = svectorscopy(interceptionstorage, Val{nclass}()),
        hortonpondingstorage = svectorscopy(hortonpondingstorage, Val{nclass}()),
        hortonrunoffstorage = svectorscopy(hortonrunoffstorage, Val{nclass}()),
        rootzonestorage = svectorscopy(srmax, Val{nclass}()),
        faststorage = 0.0 .* svectorscopy(srmax, Val{nclass}()),
        srootzone_over_srmax = svectorscopy(srootzone_over_srmax, Val{nclass}()),
        slowstorage = zeros(Float, n) .+ 30.0,
        #states previous time step
        states_m = fill(mv, n),
        states_ = svectorscopy(states_, Val{nclass}()),
        #states averaged over all classes
        interceptionstorage_m = zeros(Float, n),
        hortonpondingstorage_m = zeros(Float, n),
        hortonrunoffstorage_m = zeros(Float, n),
        srootzone_m = zeros(Float, n),
        faststorage_m = zeros(Float, n),
        srootzone_over_srmax_m = zeros(Float, n),

        # variables:
        precipitation = fill(mv, n),
        temperature = fill(mv, n),
        potential_evaporation = fill(mv, n),
        epotcorr = fill(mv, n),
        precipcorr = fill(mv, n),
        rainfallplusmelt = fill(mv, n),
        snowfall = fill(mv, n),
        snowmelt = fill(mv, n),
        potsoilevap = svectorscopy(potsoilevap, Val{nclass}()),
        precipeffective = svectorscopy(precipeffective, Val{nclass}()),
        intevap = svectorscopy(intevap, Val{nclass}()),
        hortonevap = svectorscopy(hortonevap, Val{nclass}()),
        rootevap = svectorscopy(rootevap, Val{nclass}()),
        qhortonpond = svectorscopy(qhortonpond, Val{nclass}()),
        qhortonrootzone = svectorscopy(qhortonrootzone, Val{nclass}()),
        facc = svectorscopy(facc, Val{nclass}()),
        qhortonrun = svectorscopy(qhortonrun, Val{nclass}()),
        qrootzone = svectorscopy(qrootzone, Val{nclass}()),
        qrootzonefast = svectorscopy(qrootzonefast, Val{nclass}()),
        qrootzoneslow_m = fill(mv, n),
        qfast = svectorscopy(qfast, Val{nclass}()),
        actevap = svectorscopy(actevap, Val{nclass}()),
        qpercolation = svectorscopy(qpercolation, Val{nclass}()),
        qcapillary = svectorscopy(qcapillary, Val{nclass}()),
        qpercolation_m = fill(mv, n),
        qcapillary_m = fill(mv, n),
        # combined for the classes
        actevap_m = fill(mv, n),
        intevap_m = fill(mv, n),
        hortonevap_m = fill(mv, n),
        rootevap_m = fill(mv, n),
        qslow = fill(mv, n),
        qfast_tot = fill(mv, n),
        runoff = fill(mv, n),
        wb_snow = fill(mv, n),
        wb_interception = svectorscopy(wb_interception, Val{nclass}()),
        wb_hortonponding = svectorscopy(wb_hortonponding, Val{nclass}()),
        wb_hortonrunoff = svectorscopy(wb_hortonrunoff, Val{nclass}()),
        wb_rootzone = svectorscopy(wb_rootzone, Val{nclass}()),
        wb_fast = svectorscopy(wb_fast, Val{nclass}()),
        wb_slow = fill(mv, n),
        wb_tot = fill(mv, n),
    )

    modelsize_2d = size(subcatch_2d)
    river_2d =
        ncread(nc, config, "river_location"; optional = false, type = Bool, fill = false)
    river = river_2d[inds]
    riverwidth_2d =
        ncread(nc, config, "lateral.river.width"; optional = false, type = Float, fill = 0)
    riverwidth = riverwidth_2d[inds]
    riverlength_2d =
        ncread(nc, config, "lateral.river.length"; optional = false, type = Float, fill = 0)
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

    ldd_2d = ncread(nc, config, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d = ncread(nc, config, "pits"; optional = false, type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, inds)
    end

    βₗ =
        ncread(nc, config, "lateral.land.slope"; optional = false, sel = inds, type = Float)
    clamp!(βₗ, 0.00001, Inf)

    graph = flowgraph(ldd, inds, pcr_dir)
    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, βₗ, n)

    dl = map(detdrainlength, ldd, xl, yl)
    dw = (xl .* yl) ./ dl
    olf = initialize_surfaceflow_land(
        nc,
        config,
        inds;
        sl = βₗ,
        dl = dl,
        width = map(det_surfacewidth, dw, riverwidth, river),
        frac_toriver,
        iterate = kinwave_it,
        tstep = kw_land_tstep,
        Δt = Δt,
    )

    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]

    ldd_riv = ldd_2d[inds_riv]
    if do_pits
        ldd_riv = set_pit_ldd(pits_2d, ldd_riv, inds_riv)
    end
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    rf = initialize_surfaceflow_river(
        nc,
        config,
        inds_riv;
        dl = riverlength,
        width = riverwidth,
        reservoir_index = resindex,
        reservoir = reservoirs,
        lake_index = lakeindex,
        lake = lakes,
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
    streamorder = stream_order(graph, toposort)
    min_streamorder_land = get(config.model, "min_streamorder_land", 5)
    subbas_order, indices_subbas, topo_subbas = kinwave_set_subdomains(
        graph,
        toposort,
        index_pit_land,
        streamorder,
        min_streamorder_land,
    )
    min_streamorder_river = get(config.model, "min_streamorder_river", 6)
    subriv_order, indices_subriv, topo_subriv = kinwave_set_subdomains(
        graph_riv,
        toposort_riv,
        index_pit_river,
        streamorder[index_river],
        min_streamorder_river,
    )

    modelmap = (vertical = flextopo, lateral = (land = olf, river = rf))
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
    )

    # if fews_run = true, then classes are specified as Float64 index so FEWS can import
    # NetCDF output (3D data/SVector).
    fews_run = get(config, "fews_run", false)::Bool
    classes = fews_run ? Float64.(collect(1:length(flextopo.classes))) : flextopo.classes

    writer = prepare_writer(
        config,
        modelmap,
        indices_reverse,
        x_nc,
        y_nc,
        nc,
        extra_dim = (name = "classes", value = classes),
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
        upstream_nodes = filter_upsteam_nodes(graph, pits[inds]),
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
        upstream_nodes = filter_upsteam_nodes(graph_riv, pits[inds_riv]),
        subdomain_order = subriv_order,
        topo_subdomain = topo_subriv,
        indices_subdomain = indices_subriv,
        order = toposort_riv,
        indices = inds_riv,
        reverse_indices = rev_inds_riv,
    )

    model = Model(
        config,
        (; land, river, reservoir, lake, index_river),
        (subsurface = nothing, land = olf, river = rf),
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
    vertical.dic_function[vertical.select_snow[1]](vertical)

    #lateral snow transport
    if get(config.model, "masswasting", false)::Bool
        lateral_snow_transport!(
            vertical.snow,
            vertical.snowwater,
            lateral.land.sl,
            network.land,
        )
    end

    if get(config.model, "glacier", false)::Bool
        common_glaciers(vertical, config)
    end

    for (k, class) in enumerate(vertical.classes)
        vertical.kclass[1] = k

        #INTERCEPTION
        vertical.dic_function[vertical.select_interception[k]](vertical)

        #HORTON
        vertical.dic_function[vertical.select_hortonponding[k]](vertical)
        vertical.dic_function[vertical.select_hortonrunoff[k]](vertical)

        #ROOT-ZONE
        vertical.dic_function[vertical.select_rootzone[k]](vertical)

        #FAST
        vertical.dic_function[vertical.select_fast[k]](vertical)

    end

    #COMMON SLOW
    vertical.dic_function[vertical.select_slow[1]](vertical)

    # WAT BAL
    watbal(vertical)

    surface_routing(model)

    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
