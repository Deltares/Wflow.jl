"""
    initialize_sbm_gwf_model(config::Config)

Initial part of the sbm_gwf model concept. The model contains:
    - the vertical SBM concept
    - the following lateral components:
        - 1-D kinematic wave for river flow
        - 1-D kinematic wave for overland flow
        - unconfined aquifer with groundwater flow in four directions (adjacent cells)

The unconfined aquifer contains a recharge, river and a drain (optional) boundary. 

The initial part reads the input settings and data as defined in the Config object. 
Will return a Model that is ready to run.
"""
function initialize_sbm_gwf_model(config::Config)

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
    do_drains = get(config.model, "drains", false)::Bool
    do_constanthead = get(config.model, "constanthead", false)::Bool

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool

    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, param(config, "input.subcatchment"); allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)
    modelsize_2d = size(subcatch_2d)

    river_2d = ncread(nc, param(config, "input.river_location"); type = Bool, fill = false)
    river = river_2d[inds]
    riverwidth_2d =
        ncread(nc, param(config, "input.lateral.river.width"); type = Float, fill = 0)
    riverwidth = riverwidth_2d[inds]
    riverlength_2d =
        ncread(nc, param(config, "input.lateral.river.length"); type = Float, fill = 0)
    riverlength = riverlength_2d[inds]

    altitude = ncread(nc, param(config, "input.altitude"); sel = inds, type = Float)
    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)
    riverfrac = river_fraction(river, riverlength, riverwidth, xl, yl)

    # initialize vertical SBM concept
    sbm = initialize_sbm(nc, config, riverfrac, inds)

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
            dirname(static_path),
            nc,
            inds_riv,
            nriv,
            pits,
            tosecond(Δt),
        )
    else
        lake = ()
    end

    # overland flow (kinematic wave)
    βₗ = ncread(nc, param(config, "input.lateral.land.slope"); sel = inds, type = Float)
    clamp!(βₗ, 0.00001, Inf)
    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)

    ldd = ldd_2d[inds]
    dl = fill(mv, n)
    dw = fill(mv, n)
    sw = fill(mv, n)

    for i = 1:n
        dl[i] = detdrainlength(ldd[i], xl[i], yl[i])
        dw[i] = (xl[i] * yl[i]) / dl[i]
        sw[i] = river[i] ? max(dw[i] - riverwidth[i], 0.0) : dw[i]
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
        width = sw,
        wb_pit = pits[inds],
        alpha_pow = alpha_pow,
        alpha_term = fill(mv, n),
        α = fill(mv, n),
        eps = Float(1e-03),
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

    # river flow (kinematic wave)
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
        q_av = zeros(Float, nriv),
        qin = zeros(Float, nriv),
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
        eps = Float(1e-03),
        cel = zeros(Float, nriv),
        to_river = zeros(Float, nriv),
        reservoir_index = do_reservoirs ? resindex : fill(0, nriv),
        lake_index = do_lakes ? lakeindex : fill(0, nriv),
        reservoir = do_reservoirs ? reservoirs : nothing,
        lake = do_lakes ? lakes : nothing,
        rivercells = river,
        kinwave_it = kinwave_it,
    )

    # unconfined aquifer
    if do_constanthead
        constanthead = ncread(
            nc,
            param(config, "input.lateral.subsurface.constant_head", nothing);
            sel = inds,
            type = Float,
            fill = mv,
        )
        index_constanthead = filter(i -> !isequal(constanthead[i], mv), 1:n)
        constant_head = ConstantHead(constanthead[index_constanthead], index_constanthead)
    else
        constant_head = ConstantHead{Float}(Float[], Int64[])
    end

    conductivity = ncread(
        nc,
        param(config, "input.lateral.subsurface.conductivity", nothing);
        sel = inds,
        type = Float,
    )
    specific_yield = ncread(
        nc,
        param(config, "input.lateral.subsurface.specific_yield", nothing);
        sel = inds,
        type = Float,
    )

    connectivity = Connectivity(inds, rev_inds, xl, yl)
    initial_head = altitude .- Float(0.10) # cold state for groundwater head
    initial_head[index_river] = altitude[index_river]

    if do_constanthead
        initial_head[constant_head.index] = constant_head.head
    end

    aquifer = UnconfinedAquifer(
        initial_head,
        conductivity,
        altitude,
        altitude .- sbm.soilthickness ./ Float(1000.0),
        xl .* yl,
        specific_yield,
        zeros(Float, connectivity.nconnection),  # conductance
    )

    # reset zi and satwaterdepth with groundwater head from unconfined aquifer 
    sbm.zi .= (altitude .- initial_head) .* 1000.0
    sbm.satwaterdepth .= (sbm.soilthickness .- sbm.zi) .* (sbm.θₛ .- sbm.θᵣ)

    # river boundary of unconfined aquifer
    infiltration_conductance = ncread(
        nc,
        param(config, "input.lateral.subsurface.infiltration_conductance", nothing);
        sel = inds_riv,
        type = Float,
    )
    exfiltration_conductance = ncread(
        nc,
        param(config, "input.lateral.subsurface.exfiltration_conductance", nothing);
        sel = inds_riv,
        type = Float,
    )
    river_bottom = ncread(
        nc,
        param(config, "input.lateral.subsurface.river_bottom", nothing);
        sel = inds_riv,
        type = Float,
    )

    river_flux = fill(mv, nriv)
    river_stage = fill(mv, nriv)
    river = River(
        river_stage,
        infiltration_conductance,
        exfiltration_conductance,
        river_bottom,
        river_flux,
        index_river,
    )

    # recharge boundary of unconfined aquifer
    r = fill(mv, n)
    recharge = Recharge(r, zeros(Float, n), collect(1:n))

    # drain boundary of unconfined aquifer (optional)
    if do_drains
        drain_2d = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain");
            type = Bool,
            fill = false,
        )

        drain = drain_2d[inds]
        # check if drain occurs where overland flow is not possible (sw = 0.0)
        # and correct if this is the case
        false_drain = filter(i -> !isequal(drain[i], 0) && sw[i] == Float(0), 1:n)
        n_false_drain = length(false_drain)
        if n_false_drain > 0
            drain_2d[inds[false_drain]] .= 0
            drain[false_drain] .= 0
            @info "$n_false_drain drain locations are removed that occur where overland flow
             is not possible (overland flow width is zero)"
        end
        inds_drain, rev_inds_drain = active_indices(drain_2d, 0)

        drain_elevation = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_elevation", nothing);
            sel = inds,
            type = Float,
            fill = mv,
        )
        drain_conductance = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_conductance", nothing);
            sel = inds,
            type = Float,
            fill = mv,
        )
        index_drain = filter(i -> !isequal(drain[i], 0), 1:n)
        drain_flux = fill(mv, length(index_drain))
        drains = Drainage(
            drain_elevation[index_drain],
            drain_conductance[index_drain],
            drain_flux,
            index_drain,
        )
        drain = (indices = inds_drain, reverse_indices = rev_inds_drain)
        aquifer_boundaries = AquiferBoundaryCondition[recharge, river, drains]
    else
        aquifer_boundaries = AquiferBoundaryCondition[recharge, river]
        drain = ()
    end

    gwf = GroundwaterFlow(aquifer, connectivity, constant_head, aquifer_boundaries)

    # map GroundwaterFlow and its boundaries
    if do_drains
        subsurface_map = (
            flow = gwf,
            recharge = gwf.boundaries[1],
            river = gwf.boundaries[2],
            drain = gwf.boundaries[3],
        )
    else
        subsurface_map =
            (flow = gwf, recharge = gwf.boundaries[1], river = gwf.boundaries[2])
    end

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

    modelmap =
        (vertical = sbm, lateral = (subsurface = subsurface_map, land = olf, river = rf))
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
        drain = isempty(drain) ? nothing : rev_inds_drain,
    )
    writer = prepare_writer(
        config,
        reader,
        modelmap,
        indices_reverse,
        x_nc,
        y_nc,
        nc,
        maxlayers = sbm.maxlayers,
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
        altitude = altitude,
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
        (; land, river, reservoir, lake, drain, index_river, frac_toriver),
        (subsurface = subsurface_map, land = olf, river = rf),
        sbm,
        clock,
        reader,
        writer,
    )

    # read and set states in model object if reinit=false
    if reinit == false
        instate_path = joinpath(tomldir, config.state.path_input)
        state_ncnames = ncnames(config.state)
        set_states(instate_path, model, state_ncnames, type = Float)
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


"update the sbm_gwf model for a single timestep"
function update_sbm_gwf(model)
    @unpack lateral, vertical, network, clock, config = model

    inds_riv = network.index_river
    do_drains = get(config.model, "drains", false)::Bool

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
    if get(config.model, "masswasting", false)::Bool
        lateral_snow_transport!(
            vertical.snow,
            vertical.snowwater,
            lateral.land.sl,
            network.land,
        )
    end

    # update vertical sbm concept until recharge [mm]
    update_until_recharge(vertical, config)

    # set river stage (groundwater) to average h from kinematic wave
    lateral.subsurface.river.stage .= lateral.river.h_av .+ lateral.subsurface.river.bottom

    # determine stable time step for groundwater flow
    Δt_gw = stable_timestep(lateral.subsurface.flow.aquifer) # time step in day (Float64)
    Δt_sbm = (vertical.Δt / 86400.0) # vertical.Δt is in seconds (Float64)
    if Δt_gw < Δt_sbm
        @warn("stable time step Δt $Δt for groundwater flow is smaller than sbm Δt $Δt_sbm")
    end

    Q = zeros(vertical.n)
    # exchange of recharge between vertical sbm concept and groundwater flow domain
    # recharge rate groundwater is required in units [m d⁻¹]
    lateral.subsurface.recharge.rate .= vertical.recharge ./ 1000.0 .* (1.0 / Δt_sbm)
    # update groundwater domain
    update(lateral.subsurface.flow, Q, Δt_sbm)

    # determine excess head [m] (exfiltwater) in groundwater domain (head > surface) and reset head
    exfiltwater =
        lateral.subsurface.flow.aquifer.head .-
        min.(lateral.subsurface.flow.aquifer.head, lateral.subsurface.flow.aquifer.top)
    lateral.subsurface.flow.aquifer.head .=
        min.(lateral.subsurface.flow.aquifer.head, lateral.subsurface.flow.aquifer.top)

    # update vertical sbm concept (runoff, ustorelayerdepth and satwaterdepth)
    update_after_subsurfaceflow(
        vertical,
        (network.land.altitude .- lateral.subsurface.flow.aquifer.head) .* 1000.0, # zi [mm] in vertical concept SBM
        exfiltwater .* 1000.0,
    )

    # determine lateral inflow for overland flow based on vertical runoff [mm] from vertical
    # sbm concept and drain flux (optional) from groundwater domain
    drainflux = zeros(vertical.n)
    if do_drains
        drainflux[lateral.subsurface.drain.index] = -lateral.subsurface.drain.flux
    end

    lateral.land.inwater .=
        (vertical.runoff .* network.land.xl .* network.land.yl .* 0.001 .+ drainflux) ./
        lateral.land.Δt
    lateral.land.qlat .= lateral.land.inwater ./ lateral.land.dl
    # run kinematic wave for overland flow
    update(lateral.land, network.land, frac_toriver = network.frac_toriver)

    # determine net runoff from vertical sbm concept in river cells
    net_runoff_river =
        (
            vertical.net_runoff_river[inds_riv] .* network.land.xl[inds_riv] .*
            network.land.yl[inds_riv] .* 0.001
        ) ./ vertical.Δt

    # determine lateral inflow to river from groundwater domain (river), overland flow
    # and net runoff from vertical sbm concept
    lateral.river.inwater .=
        -lateral.subsurface.river.flux ./ lateral.river.Δt .+
        lateral.land.to_river[inds_riv] .+ net_runoff_river
    lateral.river.qlat .= lateral.river.inwater ./ lateral.river.dl

    # run kinematic wave for river flow 
    # check if reservoirs or lakes are defined, then the inflow from overland flow is
    # required
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
