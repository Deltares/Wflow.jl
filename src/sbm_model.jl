"""
    initialize_sbm_model(config::Config)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_sbm_model(config::Config)

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
    config_thicknesslayers = get(config.model, "thicknesslayers", 0.0)
    if length(config_thicknesslayers) > 0
        thicknesslayers = SVector(Tuple(push!(Float64.(config_thicknesslayers), mv)))
        sumlayers = pushfirst(cumsum(thicknesslayers), 0.0)
        maxlayers = length(thicknesslayers) # max number of soil layers
    else
        maxlayers = 1
    end
    do_reservoirs = get(config.model, "reservoirs", false)
    do_lakes = get(config.model, "lakes", false)
    do_snow = get(config.model, "snow", false)

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
    modelsize_2d = size(subcatch_2d)

    # alt = Wflow.symbols"input.vertical.altitude"
    # config
    # ncname = Wflow.param(config, alt, nothing)
    # nc[ncname]

    altitude =
        ncread(nc, param(config, "input.vertical.altitude"); sel = inds, type = Float64)
    river_2d = ncread(nc, param(config, "input.river_location"); type = Bool, fill = false)
    river = river_2d[inds]
    riverwidth_2d =
        ncread(nc, param(config, "input.lateral.river.width"); type = Float64, fill = 0)
    riverwidth = riverwidth_2d[inds]
    riverlength_2d =
        ncread(nc, param(config, "input.lateral.river.length"); type = Float64, fill = 0)
    riverlength = riverlength_2d[inds]

    # read x, y coordinates and calculate cell length [m]
    y_nc = "y" in keys(nc.dim) ? ncread(nc, "y") : ncread(nc, "lat")
    x_nc = "x" in keys(nc.dim) ? ncread(nc, "x") : ncread(nc, "lon")
    if dims_xy
        y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    else
        y = repeat(y_nc, outer = (1, length(x_nc)))[inds]
    end
    cellength = abs(mean(diff(x_nc)))


    cfmax = ncread(
        nc,
        param(config, "cfmax", nothing);
        sel = inds,
        defaults = 3.75653,
        type = Float64,
    ) .* (Δt / basetimestep)
    tt = ncread(
        nc,
        param(config, "input.vertical.tt", nothing);
        sel = inds,
        defaults = 0.0,
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
        defaults = 0.0,
        type = Float64,
    )
    whc = ncread(
        nc,
        param(config, "input.vertical.whc", nothing);
        sel = inds,
        defaults = 0.1,
        type = Float64,
    )
    w_soil = ncread(
        nc,
        param(config, "input.vertical.w_soil", nothing);
        sel = inds,
        defaults = 0.1125,
        type = Float64,
    ) .* (Δt / basetimestep)
    cf_soil = ncread(
        nc,
        param(config, "input.vertical.cf_soil", nothing);
        sel = inds,
        defaults = 0.038,
        type = Float64,
    )


    # soil parameters
    θₛ = ncread(
        nc,
        param(config, "input.vertical.θₛ", nothing);
        sel = inds,
        defaults = 0.6,
        type = Float64,
    )
    θᵣ = ncread(
        nc,
        param(config, "input.vertical.θᵣ", nothing);
        sel = inds,
        defaults = 0.01,
        type = Float64,
    )
    kv₀ = ncread(
        nc,
        param(config, "input.vertical.kv₀");
        sel = inds,
        defaults = 3000.0,
        type = Float64,
    ) .* (Δt / basetimestep)
    m = ncread(
        nc,
        param(config, "input.vertical.m", nothing);
        sel = inds,
        defaults = 300.0,
        type = Float64,
    )
    hb = ncread(
        nc,
        param(config, "input.vertical.hb", nothing);
        sel = inds,
        defaults = 10.0,
        type = Float64,
    )
    soilthickness = ncread(
        nc,
        param(config, "input.vertical.soilthickness", nothing);
        sel = inds,
        defaults = 2000.0,
        type = Float64,
    )
    infiltcappath = ncread(
        nc,
        param(config, "input.vertical.infiltcappath", nothing);
        sel = inds,
        defaults = 10.0,
        type = Float64,
    )
    infiltcapsoil = ncread(
        nc,
        param(config, "input.vertical.infiltcapsoil", nothing);
        sel = inds,
        defaults = 100.0,
        type = Float64,
    )
    maxleakage = ncread(
        nc,
        param(config, "input.vertical.maxleakage", nothing);
        sel = inds,
        defaults = 0.0,
        type = Float64,
    )

    c = ncread(
        nc,
        param(config, "input.vertical.c", nothing);
        sel = inds,
        defaults = 10.0,
        type = Float64,
        dimname = "layer",
    )
    kvfrac = ncread(
        nc,
        param(config, "input.vertical.khfrac", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
        dimname = "layer",
    )

    # fraction open water and compacted area (land cover)
    waterfrac = ncread(
        nc,
        param(config, "input.vertical.waterfrac", nothing);
        sel = inds,
        defaults = 0.0,
        type = Float64,
    )
    pathfrac = ncread(
        nc,
        param(config, "input.vertical.pathfrac", nothing);
        sel = inds,
        defaults = 0.01,
        type = Float64,
    )

    # vegetation parameters
    rootingdepth = ncread(
        nc,
        param(config, "input.vertical.rootingdepth", nothing);
        sel = inds,
        defaults = 750.0,
        type = Float64,
    )
    rootdistpar = ncread(
        nc,
        param(config, "input.vertical.rootdistpar", nothing);
        sel = inds,
        defaults = -500.0,
        type = Float64,
    )
    capscale = ncread(
        nc,
        param(config, "input.vertical.capscale", nothing);
        sel = inds,
        defaults = 100.0,
        type = Float64,
    )
    et_reftopot = ncread(
        nc,
        param(config, "input.vertical.et_reftopot", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
    )

    # if leaf area index climatology provided use sl, swood and kext to calculate cmax, e_r and canopygapfraction
    # TODO replace by something else
    if haskey(config.input.vertical, "leaf_area_index")
       # TODO confirm if leaf area index climatology is present in the NetCDF
       sl = ncread(
            nc,
            param(config, "input.vertical.specific_leaf");
            sel = inds,
            type = Float64,
        )
        swood = ncread(
            nc,
            param(config, "input.vertical.storage_wood");
            sel = inds,
            type = Float64,
        )
        kext = ncread(nc, param(config, "input.vertical.kext"); sel = inds, type = Float64)
        cmax = fill(mv, n)
        e_r = fill(mv, n)
        canopygapfraction = fill(mv, n)
    else
        # cmax, e_r, canopygapfraction only required when leaf area index climatology not provided
        cmax = ncread(
            nc,
            param(config, "input.vertical.cmax", nothing);
            sel = inds,
            defaults = 1.0,
            type = Float64,
        )
        e_r = ncread(
            nc,
            param(config, "input.vertical.eoverr", nothing);
            sel = inds,
            defaults = 0.1,
            type = Float64,
        )
        canopygapfraction = ncread(
            nc,
            param(config, "input.vertical.canopygapfraction", nothing);
            sel = inds,
            defaults = 0.1,
            type = Float64,
        )
        sl = fill(mv, n)
        swood = fill(mv, n)
        kext = fill(mv, n)
    end


    # these are filled in the loop below
    # TODO see if we can replace this approach
    nlayers = zeros(Int, n)
    act_thickl = zeros(Float64, maxlayers, n)
    s_layers = zeros(Float64, maxlayers + 1, n)
    xl = fill(mv, n)
    yl = fill(mv, n)
    riverfrac = fill(mv, n)

    for i = 1:n
        xl[i] = sizeinmetres ? cellength : lattometres(y[i])[1] * cellength
        yl[i] = sizeinmetres ? cellength : lattometres(y[i])[2] * cellength
        riverfrac[i] =
            river[i] ? min((riverlength[i] * riverwidth[i]) / (xl[i] * yl[i]), 1.0) : 0.0

        if length(config_thicknesslayers) > 0
            act_thickl_, nlayers_ =
                set_layerthickness(soilthickness[i], sumlayers, thicknesslayers)
            s_layers_ = pushfirst(cumsum(act_thickl_), 0.0)

            nlayers[i] = nlayers_
            act_thickl[:, i] = act_thickl_
            s_layers[:, i] = s_layers_
        else
            nlayers[i] = 1
            act_thickl[:, i] = SVector(soilthickness[i])
            s_layers[:, i] = pushfirst(cumsum(SVector(soilthickness[i])), 0.0)
        end
    end

    # needed for derived parameters below
    act_thickl = svectorscopy(act_thickl, Val{maxlayers}())
    θₑ = θₛ .- θᵣ
    soilwatercapacity = soilthickness .* θₑ
    satwaterdepth = 0.85 .* soilwatercapacity # cold state value for satwaterdepth
    zi = max.(0.0, soilthickness .- satwaterdepth ./ θₑ) # cold state value for zi

    # copied to array of sarray below
    vwc = fill(mv, maxlayers, n)
    vwc_perc = fill(mv, maxlayers, n)

    sbm = SBM{Float64,maxlayers,maxlayers + 1}(
        Δt = tosecond(Δt),
        maxlayers = maxlayers,
        n = n,
        nlayers = nlayers,
        n_unsatlayers = fill(0, n),
        yl = yl,
        xl = xl,
        riverfrac = riverfrac,
        θₛ = θₛ,
        θᵣ = θᵣ,
        θₑ = θₑ,
        kv₀ = kv₀,
        kvfrac = svectorscopy(kvfrac, Val{maxlayers}()),
        m = m,
        hb = hb,
        soilthickness = soilthickness,
        act_thickl = act_thickl,
        sumlayers = svectorscopy(s_layers, Val{maxlayers + 1}()),
        infiltcappath = infiltcappath,
        infiltcapsoil = infiltcapsoil,
        maxleakage = maxleakage,
        waterfrac = max.(waterfrac .- riverfrac, 0.0),
        pathfrac = pathfrac,
        altitude = altitude,
        rootingdepth = rootingdepth,
        rootdistpar = rootdistpar,
        capscale = capscale,
        et_reftopot = et_reftopot,
        c = svectorscopy(c, Val{maxlayers}()),
        stemflow = fill(mv, n),
        throughfall = fill(mv, n),
        f = θₑ ./ m,
        ustorelayerdepth = act_thickl .* 0.0,
        satwaterdepth = satwaterdepth,
        zi = zi,
        soilwatercapacity = soilwatercapacity,
        canopystorage = fill(0.0, n),
        cmax = cmax,
        canopygapfraction = canopygapfraction,
        e_r = e_r,
        precipitation = fill(mv, n),
        temperature = fill(mv, n),
        potential_evaporation = fill(mv, n),
        pottrans_soil = fill(mv, n),
        transpiration = fill(mv, n),
        ae_ustore = fill(mv, n),
        ae_sat = fill(mv, n),
        interception = fill(mv, n),
        soilevap = fill(mv, n),
        soilevapsat = fill(mv, n),
        actcapflux =  fill(mv, n),
        actevapsat = fill(mv, n),
        actevap = fill(mv, n),
        runoff_river = fill(mv, n),
        runoff_land = fill(mv, n),
        ae_openw_l = fill(mv, n),
        ae_openw_r = fill(mv, n),
        avail_forinfilt = fill(mv, n),
        actinfilt = fill(mv, n),
        actinfiltsoil = fill(mv, n),
        actinfiltpath = fill(mv, n),
        infiltsoilpath= fill(mv, n),
        infiltexcess = fill(mv, n),
        excesswater = fill(mv, n),
        exfiltsatwater = fill(mv, n),
        exfiltustore = fill(mv, n),
        excesswatersoil = fill(mv, n),
        excesswaterpath = fill(mv, n),
        runoff = fill(mv, n),
        net_runoff_river = fill(mv, n),
        vwc = svectorscopy(vwc, Val{maxlayers}()),
        vwc_perc = svectorscopy(vwc_perc, Val{maxlayers}()),
        rootstore = fill(mv, n),
        vwc_root = fill(mv, n),
        vwc_percroot = fill(mv, n),
        ustoredepth = fill(mv, n),
        transfer = fill(mv, n),
        recharge = fill(mv, n),
        actleakage = fill(mv, n),
        # snow parameters
        cfmax = cfmax,
        tt = tt,
        tti = tti,
        ttm = ttm,
        whc = whc,
        w_soil = w_soil,
        cf_soil = cf_soil,
        snow = fill(0.0, n),
        snowwater = fill(0.0, n),
        rainfallplusmelt = fill(mv, n),
        tsoil = fill(10.0, n),
        # Interception related to climatology (leaf_area_index)
        sl = sl,
        swood = swood,
        kext = kext,
        leaf_area_index = fill(mv, n),
        waterlevel_land = fill(mv, n),
        waterlevel_river = fill(0.0, n), #set to zero to account for cells outside river domain
    )

    # states sbm concept
    states = ()
    if do_snow
        for state in statevars(sbm, snow=true)
            states = (states...,(:vertical,state))
        end
    else
        for state in statevars(sbm)
            states = (states...,(:vertical,state))
        end
    end

    inds_riv, rev_inds_riv = active_indices(river_2d, 0)
    nriv = length(inds_riv)

    # reservoirs
    pits = zeros(Bool, modelsize_2d)
    if do_reservoirs
        reservoirs, resindex, reservoir, pits = initialize_simple_reservoir(config, nc, inds_riv, nriv, pits)
        for state in statevars(reservoirs)
            states = (states...,(:lateral,:river,:reservoir,state))
        end
    else
        reservoir = ()
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits = initialize_natural_lake(config, nc, inds_riv, nriv, pits)
        for state in statevars(lakes)
            states = (states...,(:lateral,:river,:lake,state))
        end
    else
        lake = ()
    end

    # lateral part sbm
    khfrac = ncread(
        nc,
        param(config, "input.lateral.subsurface.ksathorfrac", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
    )
    βₗ = ncread(nc, param(config, "input.lateral.land.slope"); sel = inds, type = Float64)
    clamp!(βₗ, 0.00001, Inf)
    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)
    ldd = ldd_2d[inds]
    kh₀ = khfrac .* kv₀
    dl = fill(mv, n)
    dw = fill(mv, n)
    sw = fill(mv, n)

    for i = 1:n
        dl[i] = detdrainlength(ldd[i], xl[i], yl[i])
        dw[i] = (xl[i] * yl[i])/dl[i]
        sw[i] = river[i] ? max(dw[i] - riverwidth[i], 0.0) : dw[i]
    end

    ssf = LateralSSF{Float64}(
        kh₀ = kh₀,
        f = sbm.f,
        zi = sbm.zi,
        soilthickness = soilthickness,
        θₑ = θₛ .- θᵣ,
        Δt = 1.0,
        βₗ = βₗ,
        dl = dl .* 1000.0,
        dw = dw .* 1000.0,
        ssf = ((kh₀ .* βₗ) ./ sbm.f) .* (exp.(-sbm.f .* sbm.zi) - exp.(-sbm.f .* soilthickness)) .* dw .*1000.0,
        wb_pit = pits[inds],
    )

    n_land = ncread(
        nc,
        param(config, "input.lateral.land.n", nothing);
        sel = inds,
        defaults = 0.072,
        type = Float64,
    )

    olf = SurfaceFlow(
        sl = βₗ,
        n = n_land,
        dl = dl,
        Δt = tosecond(Δt),
        width = sw,
        wb_pit = pits[inds],
    )

    pcr_dir = dims_xy ? permute_indices(Wflow.pcrdir) : Wflow.pcrdir
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
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = Wflow.fraction_runoff_toriver(graph, ldd, index_river, βₗ, n)

    rf = SurfaceFlow(
        sl = riverslope,
        n = n_river,
        dl = riverlength,
        Δt = tosecond(Δt),
        width = riverwidth,
        reservoir_index = do_reservoirs ? resindex : fill(0, nriv),
        lake_index = do_lakes ? lakeindex : fill(0, nriv),
        reservoir = do_reservoirs ? reservoirs : nothing,
        lake = do_lakes ? lakes : nothing,
        rivercells = river,
    )

    for state in statevars(ssf)
        states = (states...,(:lateral,:subsurface,state))
    end

    for state in statevars(rf)
        states = (states...,(:lateral,:river,state))
    end

    for state in statevars(olf)
        states = (states...,(:lateral,:land,state))
    end

    reader = prepare_reader(dynamic_path, cyclic_path, config)

    modelmap = (vertical = sbm, lateral = (subsurface = ssf, land = olf, river = rf))
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
        states,
        indices_reverse,
        x_nc,
        y_nc,
        dims_xy,
        nc,
        maxlayers = maxlayers,
    )

    # for each domain save the directed acyclic graph, the traversion order,
    # and the indices that map it back to the two dimensional grid
    # for reservoirs and lakes this is information is also available (except the graph) 
    # from the initialization functions
    land = (
        graph = graph,
        order = topological_sort_by_dfs(graph),
        indices = inds,
        reverse_indices = rev_inds,
    )
    river = (
        graph = graph_riv,
        order = topological_sort_by_dfs(graph_riv),
        indices = inds_riv,
        reverse_indices = rev_inds_riv,
    )

    model = Model(
        config,
        (; land, river, reservoir, lake, index_river, frac_toriver),
        (subsurface = ssf, land = olf, river = rf),
        sbm,
        Clock(config.starttime, 1, Δt),
        reader,
        writer,
    )

    # read and set states in model object if reinit=false
    if reinit  == false
        set_states(
            instate_path,
            model,
            states,
            inds,
            config;
            type = Float64
        )
    end

    # make sure the forcing is already loaded
    # it's fine to run twice, and may help catching errors earlier
    update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end
    return model
end

"update SBM model for a single timestep"
function update(model::Model{N,L,V,R,W}) where {N,L,V<:SBM,R,W}
    @unpack lateral, vertical, network, clock, config = model

    inds_riv = network.index_river
    kinwave_it = get(config.model, "kin_wave_iteration", false)

    update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end

    # extract water levels h_av [m] from the land and river domains
    # this is used to limit open water evaporation
    vertical.waterlevel_land .= lateral.land.h_av .* 1000.0
    vertical.waterlevel_river[inds_riv] .= lateral.river.h_av .* 1000.0

    # vertical sbm concept is updated until snow state, after that (optional)
    # snow transport is possible
    update_until_snow(vertical, config)

    # lateral snow transport 
    if get(config.model, "masswasting", false)
        snowflux_frac =
            min.(0.5, lateral.land.sl ./ 5.67) .* min.(1.0, vertical.snow ./ 10000.0)
        maxflux = snowflux_frac .* vertical.snow
        vertical.snow .= accucapacityflux(network.land, vertical.snow, maxflux)
        vertical.snowwater .= accucapacityflux(network.land, vertical.snowwater, vertical.snowwater .* snowflux_frac)
    end

    # update vertical sbm concept until recharge [mm] to the saturated store
    update_until_recharge(vertical, config)

    # exchange of recharge between vertical sbm concept and subsurface flow domain
    lateral.subsurface.recharge .= vertical.recharge
    lateral.subsurface.recharge .*= lateral.subsurface.dw
    lateral.subsurface.zi .= vertical.zi

    # update lateral subsurface flow domain
    update(lateral.subsurface, network.land, network.frac_toriver, lateral.river.rivercells)

    # update vertical sbm concept (runoff, ustorelayerdepth and satwaterdepth)
    update_after_lateralflow(
        vertical,
        lateral.subsurface.zi,
        lateral.subsurface.exfiltwater,
    )

    # update overland flow based on vertical runoff [mm] from vertical sbm concept
    lateral.land.qlat .=
        (vertical.runoff .* vertical.xl .* vertical.yl .* 0.001) ./ lateral.land.Δt ./
        lateral.land.dl
    update(
        lateral.land,
        network.land,
        frac_toriver = network.frac_toriver,
        river = lateral.river.rivercells,
        do_iter = kinwave_it,
    )

    # update river domain with net runoff from vertical sbm concept, overland flow
    # and lateral subsurface flow to the river cells
    net_runoff_river = 
        (vertical.net_runoff_river[inds_riv] .* vertical.xl[inds_riv] .* 
        vertical.yl[inds_riv] .* 0.001) ./ vertical.Δt
    lateral.river.qlat .= 
        (
            lateral.subsurface.to_river[inds_riv] ./ 1.0e9 ./ lateral.river.Δt .+
            lateral.land.to_river[inds_riv] .+ net_runoff_river
        ) ./ lateral.river.dl
    update(lateral.river, network.river, do_iter = kinwave_it, doy = dayofyear(clock.time))

    write_output(model, model.writer)

    # update the clock
    clock.iteration += 1
    clock.time += clock.Δt

    return model
end
