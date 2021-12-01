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

    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, param(config, "input.subcatchment"); allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)

    # cfmax =
    #     ncread(
    #         nc,
    #         param(config, "input.vertical.cfmax", nothing);
    #         sel = inds,
    #         defaults = 3.75653,
    #         type = Float,
    #     ) .* (Δt / basetimestep)
    # tt = ncread(
    #     nc,
    #     param(config, "input.vertical.tt", nothing);
    #     sel = inds,
    #     defaults = -1.41934,
    #     type = Float,
    # )
    # tti = ncread(
    #     nc,
    #     param(config, "input.vertical.tti", nothing);
    #     sel = inds,
    #     defaults = 1.0,
    #     type = Float,
    # )
    # ttm = ncread(
    #     nc,
    #     param(config, "input.vertical.ttm", nothing);
    #     sel = inds,
    #     defaults = -1.41934,
    #     type = Float,
    # )
    # whc = ncread(
    #     nc,
    #     param(config, "input.vertical.whc", nothing);
    #     sel = inds,
    #     defaults = 0.1,
    #     type = Float,
    # )
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
    # fc = ncread(
    #     nc,
    #     param(config, "input.vertical.fc", nothing);
    #     sel = inds,
    #     defaults = 260.0,
    #     type = Float,
    # )
    # betaseepage = ncread(
    #     nc,
    #     param(config, "input.vertical.betaseepage", nothing);
    #     sel = inds,
    #     defaults = 1.8,
    #     type = Float,
    # )
    # lp = ncread(
    #     nc,
    #     param(config, "input.vertical.lp", nothing);
    #     sel = inds,
    #     defaults = 0.53,
    #     type = Float,
    # )
    # k4 =
    #     ncread(
    #         nc,
    #         param(config, "input.vertical.k4", nothing);
    #         sel = inds,
    #         defaults = 0.02307,
    #         type = Float,
    #     ) .* (Δt / basetimestep)
    # kquickflow =
    #     ncread(
    #         nc,
    #         param(config, "input.vertical.kquickflow", nothing);
    #         sel = inds,
    #         defaults = 0.09880,
    #         type = Float,
    #     ) .* (Δt / basetimestep)
    # suz = ncread(
    #     nc,
    #     param(config, "input.vertical.suz", nothing);
    #     sel = inds,
    #     defaults = 100.0,
    #     type = Float,
    # )
    # k0 =
    #     ncread(
    #         nc,
    #         param(config, "input.vertical.k0", nothing);
    #         sel = inds,
    #         defaults = 0.30,
    #         type = Float,
    #     ) .* (Δt / basetimestep)
    # khq =
    #     ncread(
    #         nc,
    #         param(config, "input.vertical.khq", nothing);
    #         sel = inds,
    #         defaults = 0.09880,
    #         type = Float,
    #     ) .* (Δt / basetimestep)
    # hq =
    #     ncread(
    #         nc,
    #         param(config, "input.vertical.hq", nothing);
    #         sel = inds,
    #         defaults = 3.27,
    #         type = Float,
    #     ) .* (Δt / basetimestep)
    # alphanl = ncread(
    #     nc,
    #     param(config, "input.vertical.alphanl", nothing);
    #     sel = inds,
    #     defaults = 1.1,
    #     type = Float,
    # )
    # perc =
    #     ncread(
    #         nc,
    #         param(config, "input.vertical.perc", nothing);
    #         sel = inds,
    #         defaults = 0.4,
    #         type = Float,
    #     ) .* (Δt / basetimestep)
    # cfr = ncread(
    #     nc,
    #     param(config, "input.vertical.cfr", nothing);
    #     sel = inds,
    #     defaults = 0.05,
    #     type = Float,
    # )
    pcorr = ncread(
        nc,
        param(config, "input.vertical.pcorr", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    # rfcf = ncread(
    #     nc,
    #     param(config, "input.vertical.rfcf", nothing);
    #     sel = inds,
    #     defaults = 1.0,
    #     type = Float,
    # )
    # sfcf = ncread(
    #     nc,
    #     param(config, "input.vertical.sfcf", nothing);
    #     sel = inds,
    #     defaults = 1.0,
    #     type = Float,
    # )
    # cflux =
    #     ncread(
    #         nc,
    #         param(config, "input.vertical.cflux", nothing);
    #         sel = inds,
    #         defaults = 2.0,
    #         type = Float,
    #     ) .* (Δt / basetimestep)
    icf = ncread(
        nc,
        param(config, "input.vertical.icf", nothing);
        sel = inds,
        defaults = 2.0,
        type = Float,
        # dimname = :classes,
    )
    # if size(icf, 1) != nclass
    #     parname = param(config, "input.vertical.icf")
    #     size1 = size(icf, 1)
    #     error("$parname needs a class dimension of size $nclass, but is $size1")
    # end

    cevpf = ncread(
        nc,
        param(config, "input.vertical.cevpf", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    # epf = ncread(
    #     nc,
    #     param(config, "input.vertical.epf", nothing);
    #     sel = inds,
    #     defaults = 1.0,
    #     type = Float,
    # )
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

    # threshold = fc .* lp
    #,nclass
    flextopo = FLEXTOPO{Float, nclass}(
        Δt = Float(tosecond(Δt)),
        nclass = nclass,
        n = n,
        # fc = fc,
        # betaseepage = betaseepage,
        # lp = lp,
        # threshold = threshold,
        # k4 = k4,
        # kquickflow = set_kquickflow ? kquickflow :
        #              pow.(khq, 1.0 .+ alphanl) .* pow.(hq, -alphanl),
        # suz = suz,
        # k0 = k0,
        # khq = khq,
        # hq = hq,
        # alphanl = alphanl,
        # perc = perc,
        # cfr = cfr,
        pcorr = pcorr,
        # rfcf = rfcf,
        # sfcf = sfcf,
        # cflux = cflux,
        icf = icf,
        # icf = svectorscopy(icf, Val{nclass}()),
        cevpf = cevpf,
        # epf = epf,
        ecorr = ecorr,
        # tti = tti,
        # tt = tt,
        # ttm = ttm,
        # cfmax = cfmax,
        # whc = whc,

        # # glacier parameters
        # g_tt = g_tt,
        # g_sifrac = g_sifrac,
        # g_cfmax = g_cfmax,
        # glacierstore = glacierstore,
        # glacierfrac = glacierfrac,
        
        # # default (cold) states:
        interceptionstorage = zeros(Float, n),
        # snow = zeros(Float, n),
        # snowwater = zeros(Float, n),
        # soilmoisture = copy(fc),
        # upperzonestorage = 0.2 .* fc,
        # lowerzonestorage = 1.0 ./ (3.0 .* k4),

        # variables:
        precipitation = fill(mv, n),
        temperature = fill(mv, n),
        potential_evaporation = fill(mv, n),
        potsoilevap = fill(mv, n),
        # soilevap = fill(mv, n),
        intevap = fill(mv, n),
        # actevap = fill(mv, n),
        # rainfallplusmelt = fill(mv, n),
        # directrunoff = fill(mv, n),
        # hbv_seepage = fill(mv, n),
        # in_upperzone = fill(mv, n),
        # quickflow = fill(mv, n),
        # real_quickflow = fill(mv, n),
        # percolation = fill(mv, n),
        # capflux = fill(mv, n),
        # baseflow = fill(mv, n),
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

    # vertical hbv concept is updated until snow state, after that (optional)
    # snow transport is possible
    update_until_snow(vertical, config)

    # # lateral snow transport
    # if get(config.model, "masswasting", false)::Bool
    #     lateral_snow_transport!(
    #         vertical.snow,
    #         vertical.snowwater,
    #         lateral.land.sl,
    #         network.land,
    #     )
    # end

    # # update vertical hbv concept
    # update_after_snow(vertical, config)

    surface_routing(model)

    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
