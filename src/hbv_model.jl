"""
    initialize_hbv_model(config::Config)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_hbv_model(config::Config)
    # unpack the paths to the NetCDF files
    tomldir = dirname(config)
    static_path = joinpath(tomldir, config.input.path_static)
    cyclic_path = joinpath(tomldir, config.input.path_static)
    dynamic_path = joinpath(tomldir, config.input.path_forcing)
    instate_path = joinpath(tomldir, config.state.path_input)
    output_path = joinpath(tomldir, config.output.path)

    Δt = Second(config.timestepsecs)

    sizeinmetres = get(config.model, "sizeinmetres", false)
    reinit = get(config.model, "reinit", true)
    do_reservoirs = get(config.model, "reservoirs", false)
    do_lakes = get(config.model, "lakes", false)
    do_pits = get(config.model, "pits", false)
    set_kquickflow = get(config.model, "set_kquickflow", false)

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)

    nc = NCDataset(static_path)
    dims = dimnames(nc[param(config, "input.subcatchment")])

    # There is no need to permute the dimensions of the data, since the active indices are
    # correctly calculated in both ways.
    # The dimension order only needs to be known for interpreting the LDD directions
    # and creating the coordinate maps.
    dims_xy = dims[2] in ("y", "lat")
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
            type = Float64,
        ) .* (Δt / basetimestep)
    tt = ncread(
        nc,
        param(config, "input.vertical.tt", nothing);
        sel = inds,
        defaults = -1.41934,
        type = Float64,
    )
    tti = ncread(
        nc,
        param(config, "input.vertical.tti", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
    )
    ttm = ncread(
        nc,
        param(config, "input.vertical.ttm", nothing);
        sel = inds,
        defaults = -1.41934,
        type = Float64,
    )
    whc = ncread(
        nc,
        param(config, "input.vertical.whc", nothing);
        sel = inds,
        defaults = 0.1,
        type = Float64,
    )
    # glacier parameters
    g_tt = ncread(
        nc,
        param(config, "input.vertical.g_tt", nothing);
        sel = inds,
        defaults = 0.0,
        type = Float64,
        fill = 0.0,
    )
    g_cfmax =
        ncread(
            nc,
            param(config, "input.vertical.g_cfmax", nothing);
            sel = inds,
            defaults = 3.0,
            type = Float64,
            fill = 0.0,
        ) .* (Δt / basetimestep)
    g_sifrac = ncread(
        nc,
        param(config, "input.vertical.g_sifrac", nothing);
        sel = inds,
        defaults = 0.001,
        type = Float64,
        fill = 0.0,
    )
    glacierfrac = ncread(
        nc,
        param(config, "input.vertical.glacierfrac", nothing);
        sel = inds,
        defaults = 0.0,
        type = Float64,
        fill = 0.0,
    )
    glacierstore = ncread(
        nc,
        param(config, "input.vertical.glacierstore", nothing);
        sel = inds,
        defaults = 5500.0,
        type = Float64,
        fill = 0.0,
    )
    fc = ncread(
        nc,
        param(config, "input.vertical.fc", nothing);
        sel = inds,
        defaults = 260.0,
        type = Float64,
    )
    betaseepage = ncread(
        nc,
        param(config, "input.vertical.betaseepage", nothing);
        sel = inds,
        defaults = 1.8,
        type = Float64,
    )
    lp = ncread(
        nc,
        param(config, "input.vertical.lp", nothing);
        sel = inds,
        defaults = 0.53,
        type = Float64,
    )
    k4 =
        ncread(
            nc,
            param(config, "input.vertical.k4", nothing);
            sel = inds,
            defaults = 0.02307,
            type = Float64,
        ) .* (Δt / basetimestep)
    kquickflow =
        ncread(
            nc,
            param(config, "input.vertical.kquickflow", nothing);
            sel = inds,
            defaults = 0.09880,
            type = Float64,
        ) .* (Δt / basetimestep)
    suz = ncread(
        nc,
        param(config, "input.vertical.suz", nothing);
        sel = inds,
        defaults = 100.0,
        type = Float64,
    )
    k0 =
        ncread(
            nc,
            param(config, "input.vertical.k0", nothing);
            sel = inds,
            defaults = 0.30,
            type = Float64,
        ) .* (Δt / basetimestep)
    khq =
        ncread(
            nc,
            param(config, "input.vertical.khq", nothing);
            sel = inds,
            defaults = 0.09880,
            type = Float64,
        ) .* (Δt / basetimestep)
    hq =
        ncread(
            nc,
            param(config, "input.vertical.hq", nothing);
            sel = inds,
            defaults = 3.27,
            type = Float64,
        ) .* (Δt / basetimestep)
    alphanl = ncread(
        nc,
        param(config, "input.vertical.alphanl", nothing);
        sel = inds,
        defaults = 1.1,
        type = Float64,
    )
    perc =
        ncread(
            nc,
            param(config, "input.vertical.perc", nothing);
            sel = inds,
            defaults = 0.4,
            type = Float64,
        ) .* (Δt / basetimestep)
    cfr = ncread(
        nc,
        param(config, "input.vertical.cfr", nothing);
        sel = inds,
        defaults = 0.05,
        type = Float64,
    )
    pcorr = ncread(
        nc,
        param(config, "input.vertical.pcorr", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
    )
    rfcf = ncread(
        nc,
        param(config, "input.vertical.rfcf", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
    )
    sfcf = ncread(
        nc,
        param(config, "input.vertical.sfcf", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
    )
    cflux =
        ncread(
            nc,
            param(config, "input.vertical.cflux", nothing);
            sel = inds,
            defaults = 2.0,
            type = Float64,
        ) .* (Δt / basetimestep)
    icf = ncread(
        nc,
        param(config, "input.vertical.icf", nothing);
        sel = inds,
        defaults = 2.0,
        type = Float64,
    )
    cevpf = ncread(
        nc,
        param(config, "input.vertical.cevpf", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
    )
    epf = ncread(
        nc,
        param(config, "input.vertical.epf", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
    )
    ecorr = ncread(
        nc,
        param(config, "input.vertical.ecorr", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
    )

    # read x, y coordinates and calculate cell length [m]
    y_nc = "y" in keys(nc.dim) ? ncread(nc, "y") : ncread(nc, "lat")
    x_nc = "x" in keys(nc.dim) ? ncread(nc, "x") : ncread(nc, "lon")
    if dims_xy
        y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    else
        y = repeat(y_nc, outer = (1, length(x_nc)))[inds]
    end
    cellength = abs(mean(diff(x_nc)))

    xl = fill(mv, n)
    yl = fill(mv, n)
    for i = 1:n
        xl[i] = sizeinmetres ? cellength : lattometres(y[i])[1] * cellength
        yl[i] = sizeinmetres ? cellength : lattometres(y[i])[2] * cellength
    end

    threshold = fc .* lp

    altitude =
        ncread(nc, param(config, "input.vertical.altitude"); sel = inds, type = Float64)

    hbv = HBV{Float64}(
        Δt = tosecond(Δt),
        n = n,
        yl = yl,
        xl = xl,
        altitude = altitude,
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
        interceptionstorage = fill(0.0, n),
        snow = fill(0.0, n),
        snowwater = fill(0.0, n),
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
        ncread(nc, param(config, "input.lateral.river.width"); type = Float64, fill = 0)
    riverwidth = riverwidth_2d[inds]
    riverlength_2d =
        ncread(nc, param(config, "input.lateral.river.length"); type = Float64, fill = 0)
    riverlength = riverlength_2d[inds]

    inds_riv, rev_inds_riv = active_indices(river_2d, 0)
    nriv = length(inds_riv)

    # reservoirs
    pits = zeros(Bool, modelsize_2d)
    if do_reservoirs
        reservoirs, resindex, reservoir, pits =
            initialize_simple_reservoir(config, nc, inds_riv, nriv, pits)
    else
        reservoir = ()
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits =
            initialize_natural_lake(config, static_path, nc, inds_riv, nriv, pits)
    else
        lake = ()
    end

    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d = ncread(nc, param(config, "input.pits"); type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, inds)
    end

    βₗ = ncread(nc, param(config, "input.lateral.land.slope"); sel = inds, type = Float64)
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
        type = Float64,
    )

    alpha_pow = (2.0 / 3.0) * 0.6
    β = 0.6
    h = fill(0.0, n)
    alpha_term = pow.(n_land ./ sqrt.(βₗ), β)
    olf = SurfaceFlow(
        β = β,
        sl = βₗ,
        n = n_land,
        dl = dl,
        q = fill(0.0, n),
        qin = fill(0.0, n),
        q_av = fill(0.0, n),
        qlat = fill(0.0, n),
        h = h,
        h_av = fill(0.0, n),
        Δt = tosecond(Δt),
        its = kw_land_tstep > 0 ? Int(cld(tosecond(Δt), kw_land_tstep)) : kw_land_tstep,
        width = dw,
        wb_pit = pits[inds],
        alpha_term = alpha_term,
        alpha_pow = alpha_pow,
        α = alpha_term .* pow.(dw .+ 2.0 .* h, alpha_pow),
        eps = 1e-03,
        cel = fill(0.0, n),
        to_river = fill(0.0, n),
        rivercells = fill(false, n),
        reservoir_index = fill(0, n),
        lake_index = fill(0, n),
        reservoir = nothing,
        lake = nothing,
    )

    pcr_dir = dims_xy ? permute_indices(pcrdir) : pcrdir
    graph = flowgraph(ldd, inds, pcr_dir)

    riverslope = ncread(
        nc,
        param(config, "input.lateral.river.slope");
        sel = inds_riv,
        type = Float64,
    )
    clamp!(riverslope, 0.00001, Inf)
    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]
    n_river = ncread(
        nc,
        param(config, "input.lateral.river.n", nothing);
        sel = inds_riv,
        defaults = 0.036,
        type = Float64,
    )
    ldd_riv = ldd_2d[inds_riv]
    if do_pits
        ldd_riv = set_pit_ldd(pits_2d, ldd_riv, inds_riv)
    end
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, βₗ, n)

    h_river = fill(0.0, nriv)
    alpha_term = pow.(n_river ./ sqrt.(riverslope), β)
    rf = SurfaceFlow(
        β = β,
        sl = riverslope,
        n = n_river,
        dl = riverlength,
        q = fill(0.0, nriv),
        qin = fill(0.0, nriv),
        q_av = fill(0.0, nriv),
        qlat = fill(0.0, nriv),
        h = h_river,
        h_av = fill(0.0, nriv),
        Δt = tosecond(Δt),
        its = kw_river_tstep > 0 ? ceil(Int(tosecond(Δt) / kw_river_tstep)) :
              kw_river_tstep,
        width = riverwidth,
        wb_pit = pits[inds_riv],
        alpha_term = alpha_term,
        alpha_pow = alpha_pow,
        α = alpha_term .* pow.(riverwidth .+ 2.0 .* h_river, alpha_pow),
        eps = 1e-03,
        cel = fill(0.0, nriv),
        to_river = fill(0.0, nriv),
        reservoir_index = do_reservoirs ? resindex : fill(0, nriv),
        lake_index = do_lakes ? lakeindex : fill(0, nriv),
        reservoir = do_reservoirs ? reservoirs : nothing,
        lake = do_lakes ? lakes : nothing,
        rivercells = river,
    )

    state_ncnames = ncnames(config.state)

    reader = prepare_reader(dynamic_path, cyclic_path, config)

    modelmap = (vertical = hbv, lateral = (land = olf, river = rf))
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
    )
    writer = prepare_writer(
        config,
        reader,
        output_path,
        modelmap,
        state_ncnames,
        indices_reverse,
        x_nc,
        y_nc,
        dims_xy,
        nc,
    )

    # for each domain save the directed acyclic graph, the traversion order,
    # and the indices that map it back to the two dimensional grid
    land = (
        graph = graph,
        upstream_nodes = filter_upsteam_nodes(graph, olf.wb_pit),
        order = topological_sort_by_dfs(graph),
        indices = inds,
        reverse_indices = rev_inds,
    )
    river = (
        graph = graph_riv,
        upstream_nodes = filter_upsteam_nodes(graph_riv, rf.wb_pit),
        order = topological_sort_by_dfs(graph_riv),
        indices = inds_riv,
        reverse_indices = rev_inds_riv,
    )

    model = Model(
        config,
        (; land, river, reservoir, lake, index_river, frac_toriver),
        (land = olf, river = rf),
        hbv,
        Clock(config.starttime, 1, Δt),
        reader,
        writer,
    )

    # read and set states in model object if reinit=true
    if reinit == false
        set_states(instate_path, model, state_ncnames; type = Float64)
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
    kinwave_it = get(config.model, "kin_wave_iteration", false)

    update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end

    # vertical hbv concept is updated until snow state, after that (optional)
    # snow transport is possible
    update_until_snow(vertical, config)

    # lateral snow transport
    if get(config.model, "masswasting", false)
        snowflux_frac =
            min.(0.5, lateral.land.sl ./ 5.67) .* min.(1.0, vertical.snow ./ 10000.0)
        maxflux = snowflux_frac .* vertical.snow
        vertical.snow .= accucapacityflux(network.land, vertical.snow, maxflux)
        vertical.snowwater .= accucapacityflux(
            network.land,
            vertical.snowwater,
            vertical.snowwater .* snowflux_frac,
        )
    end
    
    # update vertical hbv concept
    update_after_snow(vertical, config)

    # determine lateral inflow for overland flow based on vertical runoff [mm] from vertical
    # hbv concept
    lateral.land.qlat .=
        (vertical.runoff .* vertical.xl .* vertical.yl .* 0.001) ./ lateral.land.Δt ./
        lateral.land.dl
    
    # run kinematic wave for overland flow
    update(
        lateral.land,
        network.land,
        frac_toriver = network.frac_toriver,
        do_iter = kinwave_it,
    )

    # determine lateral inflow (from overland flow) for river flow
    lateral.river.qlat .= (lateral.land.to_river[inds_riv]) ./ lateral.river.dl

    # run kinematic wave for river flow
    # check if reservoirs or lakes are defined, the inflow from overland flow is required
    if !isnothing(lateral.river.reservoir) || !isnothing(lateral.river.lake)
        update(
            lateral.river,
            network.river,
            inflow_wb = lateral.land.q_av[inds_riv],
            do_iter = kinwave_it,
            doy = dayofyear(clock.time),
        )
    else
        update(
            lateral.river,
            network.river,
            do_iter = kinwave_it,
            doy = dayofyear(clock.time),
        )
    end

    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
