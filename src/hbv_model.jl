"""
    initialize_hbv_model(config::Config)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_hbv_model(config::Config)
    # unpack the paths to the NetCDF files
    tomldir = dirname(config)
    static_path = joinpath(tomldir, config.input.path_static)
    dynamic_path = joinpath(tomldir, config.input.path_forcing)

    reader = prepare_reader(dynamic_path, static_path, config)
    clock = Clock(config, reader)
    Δt = clock.Δt

    reinit = get(config.model, "reinit", true)::Bool
    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_pits = get(config.model, "pits", false)::Bool
    set_kquickflow = get(config.model, "set_kquickflow", false)::Bool

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool

    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, param(config, "input.subcatchment"); allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)

    cfmax =
        ncread(
            nc,
            param(config, "input.vertical.cfmax", nothing);
            sel = inds,
            defaults = 3.75653,
            type = Float,
        ) .* (Δt / basetimestep)
    tt = ncread(
        nc,
        param(config, "input.vertical.tt", nothing);
        sel = inds,
        defaults = -1.41934,
        type = Float,
    )
    tti = ncread(
        nc,
        param(config, "input.vertical.tti", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    ttm = ncread(
        nc,
        param(config, "input.vertical.ttm", nothing);
        sel = inds,
        defaults = -1.41934,
        type = Float,
    )
    whc = ncread(
        nc,
        param(config, "input.vertical.whc", nothing);
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    # glacier parameters
    g_tt = ncread(
        nc,
        param(config, "input.vertical.g_tt", nothing);
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    g_cfmax =
        ncread(
            nc,
            param(config, "input.vertical.g_cfmax", nothing);
            sel = inds,
            defaults = 3.0,
            type = Float,
            fill = 0.0,
        ) .* (Δt / basetimestep)
    g_sifrac = ncread(
        nc,
        param(config, "input.vertical.g_sifrac", nothing);
        sel = inds,
        defaults = 0.001,
        type = Float,
        fill = 0.0,
    )
    glacierfrac = ncread(
        nc,
        param(config, "input.vertical.glacierfrac", nothing);
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    glacierstore = ncread(
        nc,
        param(config, "input.vertical.glacierstore", nothing);
        sel = inds,
        defaults = 5500.0,
        type = Float,
        fill = 0.0,
    )
    fc = ncread(
        nc,
        param(config, "input.vertical.fc", nothing);
        sel = inds,
        defaults = 260.0,
        type = Float,
    )
    betaseepage = ncread(
        nc,
        param(config, "input.vertical.betaseepage", nothing);
        sel = inds,
        defaults = 1.8,
        type = Float,
    )
    lp = ncread(
        nc,
        param(config, "input.vertical.lp", nothing);
        sel = inds,
        defaults = 0.53,
        type = Float,
    )
    k4 =
        ncread(
            nc,
            param(config, "input.vertical.k4", nothing);
            sel = inds,
            defaults = 0.02307,
            type = Float,
        ) .* (Δt / basetimestep)
    kquickflow =
        ncread(
            nc,
            param(config, "input.vertical.kquickflow", nothing);
            sel = inds,
            defaults = 0.09880,
            type = Float,
        ) .* (Δt / basetimestep)
    suz = ncread(
        nc,
        param(config, "input.vertical.suz", nothing);
        sel = inds,
        defaults = 100.0,
        type = Float,
    )
    k0 =
        ncread(
            nc,
            param(config, "input.vertical.k0", nothing);
            sel = inds,
            defaults = 0.30,
            type = Float,
        ) .* (Δt / basetimestep)
    khq =
        ncread(
            nc,
            param(config, "input.vertical.khq", nothing);
            sel = inds,
            defaults = 0.09880,
            type = Float,
        ) .* (Δt / basetimestep)
    hq =
        ncread(
            nc,
            param(config, "input.vertical.hq", nothing);
            sel = inds,
            defaults = 3.27,
            type = Float,
        ) .* (Δt / basetimestep)
    alphanl = ncread(
        nc,
        param(config, "input.vertical.alphanl", nothing);
        sel = inds,
        defaults = 1.1,
        type = Float,
    )
    perc =
        ncread(
            nc,
            param(config, "input.vertical.perc", nothing);
            sel = inds,
            defaults = 0.4,
            type = Float,
        ) .* (Δt / basetimestep)
    cfr = ncread(
        nc,
        param(config, "input.vertical.cfr", nothing);
        sel = inds,
        defaults = 0.05,
        type = Float,
    )
    pcorr = ncread(
        nc,
        param(config, "input.vertical.pcorr", nothing);
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
    cflux =
        ncread(
            nc,
            param(config, "input.vertical.cflux", nothing);
            sel = inds,
            defaults = 2.0,
            type = Float,
        ) .* (Δt / basetimestep)
    icf = ncread(
        nc,
        param(config, "input.vertical.icf", nothing);
        sel = inds,
        defaults = 2.0,
        type = Float,
    )
    cevpf = ncread(
        nc,
        param(config, "input.vertical.cevpf", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    epf = ncread(
        nc,
        param(config, "input.vertical.epf", nothing);
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

    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))


    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)

    threshold = fc .* lp

    hbv = HBV{Float}(
        Δt = Float(tosecond(Δt)),
        n = n,
        fc = fc,
        betaseepage = betaseepage,
        lp = lp,
        threshold = threshold,
        k4 = k4,
        kquickflow = set_kquickflow ? kquickflow :
                     pow.(khq, 1.0 .+ alphanl) .* pow.(hq, -alphanl),
        suz = suz,
        k0 = k0,
        khq = khq,
        hq = hq,
        alphanl = alphanl,
        perc = perc,
        cfr = cfr,
        pcorr = pcorr,
        rfcf = rfcf,
        sfcf = sfcf,
        cflux = cflux,
        icf = icf,
        cevpf = cevpf,
        epf = epf,
        ecorr = ecorr,
        tti = tti,
        tt = tt,
        ttm = ttm,
        cfmax = cfmax,
        whc = whc,
        # glacier parameters
        g_tt = g_tt,
        g_sifrac = g_sifrac,
        g_cfmax = g_cfmax,
        glacierstore = glacierstore,
        glacierfrac = glacierfrac,
        # default (cold) states:
        interceptionstorage = zeros(Float, n),
        snow = zeros(Float, n),
        snowwater = zeros(Float, n),
        soilmoisture = copy(fc),
        upperzonestorage = 0.2 .* fc,
        lowerzonestorage = 1.0 ./ (3.0 .* k4),
        # variables:
        precipitation = fill(mv, n),
        temperature = fill(mv, n),
        potential_evaporation = fill(mv, n),
        potsoilevap = fill(mv, n),
        soilevap = fill(mv, n),
        intevap = fill(mv, n),
        actevap = fill(mv, n),
        rainfallplusmelt = fill(mv, n),
        directrunoff = fill(mv, n),
        hbv_seepage = fill(mv, n),
        in_upperzone = fill(mv, n),
        quickflow = fill(mv, n),
        real_quickflow = fill(mv, n),
        percolation = fill(mv, n),
        capflux = fill(mv, n),
        baseflow = fill(mv, n),
        runoff = fill(mv, n),
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
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits = initialize_natural_lake(
            config,
            static_path,
            nc,
            inds_riv,
            nriv,
            pits,
            tosecond(Δt),
        )
    else
        lake = ()
    end

    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d = ncread(nc, param(config, "input.pits"); type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, inds)
    end

    βₗ = ncread(nc, param(config, "input.lateral.land.slope"); sel = inds, type = Float)
    clamp!(βₗ, 0.00001, Inf)
    dl = fill(mv, n)
    dw = fill(mv, n)

    for i = 1:n
        dl[i] = detdrainlength(ldd[i], xl[i], yl[i])
        dw[i] = detdrainwidth(ldd[i], xl[i], yl[i])
    end

    n_land = ncread(
        nc,
        param(config, "input.lateral.land.n", nothing);
        sel = inds,
        defaults = 0.072,
        type = Float,
    )

    alpha_pow = Float((2.0 / 3.0) * 0.6)
    β = Float(0.6)
    olf = SurfaceFlow(
        β = β,
        sl = βₗ,
        n = n_land,
        dl = dl,
        q = zeros(Float, n),
        qin = zeros(Float, n),
        q_av = zeros(Float, n),
        qlat = zeros(Float, n),
        inwater = zeros(Float, n),
        volume = zeros(Float, n),
        h = zeros(Float, n),
        h_av = zeros(Float, n),
        h_bankfull = zeros(Float, n),
        Δt = Float(tosecond(Δt)),
        its = kw_land_tstep > 0 ? Int(cld(tosecond(Δt), kw_land_tstep)) : kw_land_tstep,
        width = dw,
        wb_pit = pits[inds],
        alpha_pow = alpha_pow,
        alpha_term = fill(mv, n),
        α = fill(mv, n),
        cel = zeros(Float, n),
        to_river = zeros(Float, n),
        rivercells = fill(false, n),
        reservoir_index = fill(0, n),
        lake_index = fill(0, n),
        reservoir = nothing,
        lake = nothing,
        kinwave_it = kinwave_it,
    )

    graph = flowgraph(ldd, inds, pcr_dir)

    riverslope =
        ncread(nc, param(config, "input.lateral.river.slope"); sel = inds_riv, type = Float)
    clamp!(riverslope, 0.00001, Inf)
    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]
    n_river = ncread(
        nc,
        param(config, "input.lateral.river.n", nothing);
        sel = inds_riv,
        defaults = 0.036,
        type = Float,
    )
    h_bankfull = ncread(
        nc,
        param(config, "input.lateral.river.h_bankfull", nothing);
        sel = inds_riv,
        defaults = 1.0,
        type = Float,
    )
    ldd_riv = ldd_2d[inds_riv]
    if do_pits
        ldd_riv = set_pit_ldd(pits_2d, ldd_riv, inds_riv)
    end
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, βₗ, n)

    rf = SurfaceFlow(
        β = β,
        sl = riverslope,
        n = n_river,
        dl = riverlength,
        q = zeros(Float, nriv),
        qin = zeros(Float, nriv),
        q_av = zeros(Float, nriv),
        qlat = zeros(Float, nriv),
        inwater = zeros(Float, nriv),
        volume = zeros(Float, nriv),
        h = zeros(Float, nriv),
        h_av = zeros(Float, nriv),
        h_bankfull = h_bankfull,
        Δt = Float(tosecond(Δt)),
        its = kw_river_tstep > 0 ? ceil(Int(tosecond(Δt) / kw_river_tstep)) :
              kw_river_tstep,
        width = riverwidth,
        wb_pit = pits[inds_riv],
        alpha_pow = alpha_pow,
        alpha_term = fill(mv, nriv),
        α = fill(mv, nriv),
        cel = zeros(Float, nriv),
        to_river = zeros(Float, nriv),
        reservoir_index = do_reservoirs ? resindex : fill(0, nriv),
        lake_index = do_lakes ? lakeindex : fill(0, nriv),
        reservoir = do_reservoirs ? reservoirs : nothing,
        lake = do_lakes ? lakes : nothing,
        rivercells = river,
        kinwave_it = kinwave_it,
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

    modelmap = (vertical = hbv, lateral = (land = olf, river = rf))
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
    )
    writer = prepare_writer(config, reader, modelmap, indices_reverse, x_nc, y_nc, nc)
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
        (land = olf, river = rf),
        hbv,
        clock,
        reader,
        writer,
    )

    # read and set states in model object if reinit=true
    if reinit == false
        instate_path = joinpath(tomldir, config.state.path_input)
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

    # make sure the forcing is already loaded
    # it's fine to run twice, and may help catching errors earlier
    update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end
    return model
end

function update(model::Model{N,L,V,R,W}) where {N,L,V<:HBV,R,W}
    @unpack lateral, vertical, network, clock, config = model

    inds_riv = network.index_river
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool

    update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end

    # vertical hbv concept is updated until snow state, after that (optional)
    # snow transport is possible
    update_until_snow(vertical, config)

    # lateral snow transport
    if get(config.model, "masswasting", false)::Bool
        lateral_snow_transport!(
            vertical.snow,
            vertical.snowwater,
            lateral.land.sl,
            network.land,
        )
    end

    # update vertical hbv concept
    update_after_snow(vertical, config)

    # determine lateral inflow for overland flow based on vertical runoff [mm] from vertical
    # hbv concept
    lateral.land.inwater .=
        (vertical.runoff .* network.land.xl .* network.land.yl .* 0.001) ./ lateral.land.Δt
    lateral.land.qlat .= lateral.land.inwater ./ lateral.land.dl

    # run kinematic wave for overland flow
    update(lateral.land, network.land, frac_toriver = network.frac_toriver)

    # determine lateral inflow (from overland flow) for river flow
    lateral.river.inwater .= copy(lateral.land.to_river[inds_riv])
    lateral.river.qlat .= lateral.river.inwater ./ lateral.river.dl

    # run kinematic wave for river flow
    # check if reservoirs or lakes are defined, the inflow from overland flow is required
    if !isnothing(lateral.river.reservoir) || !isnothing(lateral.river.lake)
        update(
            lateral.river,
            network.river,
            inflow_wb = lateral.land.q_av[inds_riv],
            doy = dayofyear(clock.time),
        )
    else
        update(lateral.river, network.river, doy = dayofyear(clock.time))
    end

    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
