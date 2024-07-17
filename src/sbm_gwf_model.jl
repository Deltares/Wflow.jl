"""
    initialize_sbm_gwf_model(config::Config)

Initial part of the sbm_gwf model concept. The model contains:
    - the vertical SBM concept
    - unconfined aquifer with groundwater flow in four directions (adjacent cells)
    - the following surface routing options:
        - 1-D kinematic wave for river flow and 1-D kinematic wave for overland flow
        - 1-D local inertial model for river flow (optional floodplain) and 1-D kinematic wave for overland flow
        - 1-D local inertial model for river flow (optional floodplain) and 2-D local inertial model for overland flow

The unconfined aquifer contains a recharge, river and a drain (optional) boundary.

The initial part reads the input settings and data as defined in the Config object.
Will return a Model that is ready to run.
"""
function initialize_sbm_gwf_model(config::Config)

    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)

    reader = prepare_reader(config)
    clock = Clock(config, reader)
    dt = clock.dt

    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_drains = get(config.model, "drains", false)::Bool
    do_constanthead = get(config.model, "constanthead", false)::Bool

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool
    routing_options = ("kinematic-wave", "local-inertial")
    floodplain_1d = get(config.model, "floodplain_1d", false)::Bool
    river_routing = get_options(
        config.model,
        "river_routing",
        routing_options,
        "kinematic-wave",
    )::String
    land_routing =
        get_options(config.model, "land_routing", routing_options, "kinematic-wave")::String
    do_water_demand = haskey(config.model, "water_demand")

    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, config, "subcatchment"; optional = false, allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)
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

    altitude = ncread(nc, config, "altitude"; optional = false, sel = inds, type = Float)
    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)
    riverfrac = river_fraction(river, riverlength, riverwidth, xl, yl)

    inds_riv, rev_inds_riv = active_indices(river_2d, 0)
    nriv = length(inds_riv)

    # initialize vertical SBM concept
    sbm = initialize_sbm(nc, config, riverfrac, inds)

    # reservoirs
    pits = zeros(Bool, modelsize_2d)
    if do_reservoirs
        reservoirs, resindex, reservoir, pits =
            initialize_simple_reservoir(config, nc, inds_riv, nriv, pits, tosecond(dt))
    else
        reservoir = ()
        reservoirs = nothing
        resindex = fill(0, nriv)
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits =
            initialize_lake(config, nc, inds_riv, nriv, pits, tosecond(dt))
    else
        lake = ()
        lakes = nothing
        lakeindex = fill(0, nriv)
    end

    # overland flow (kinematic wave)
    landslope =
        ncread(nc, config, "lateral.land.slope"; optional = false, sel = inds, type = Float)
    clamp!(landslope, 0.00001, Inf)
    ldd_2d = ncread(nc, config, "ldd"; optional = false, allow_missing = true)

    ldd = ldd_2d[inds]

    dl = map(detdrainlength, ldd, xl, yl)
    dw = (xl .* yl) ./ dl
    sw = map(det_surfacewidth, dw, riverwidth, river)

    graph = flowgraph(ldd, inds, pcr_dir)
    ldd_riv = ldd_2d[inds_riv]
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, landslope, n)

    inds_allocation_areas = Vector{Int}[]
    inds_riv_allocation_areas = Vector{Int}[]
    if do_water_demand
        areas = unique(sbm.allocation.areas)
        for a in areas
            area_index = findall(x -> x == a, sbm.allocation.areas)
            push!(inds_allocation_areas, area_index)
            area_riv_index = findall(x -> x == a, sbm.allocation.areas[index_river])
            push!(inds_riv_allocation_areas, area_riv_index)
        end
    end

    if land_routing == "kinematic-wave"
        olf = initialize_surfaceflow_land(
            nc,
            config,
            inds;
            sl = landslope,
            dl,
            width = map(det_surfacewidth, dw, riverwidth, river),
            iterate = kinwave_it,
            tstep = kw_land_tstep,
            dt,
        )
    elseif land_routing == "local-inertial"
        index_river_nf = rev_inds_riv[inds] # not filtered (with zeros)
        olf, indices = initialize_shallowwater_land(
            nc,
            config,
            inds;
            modelsize_2d,
            indices_reverse = rev_inds,
            xlength = xl,
            ylength = yl,
            riverwidth = riverwidth_2d[inds_riv],
            graph_riv,
            ldd_riv,
            inds_riv,
            river,
            waterbody = !=(0).(resindex + lakeindex),
            dt,
        )
    end

    # river flow (kinematic wave)
    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]
    minimum(riverlength) > 0 || error("river length must be positive on river cells")
    minimum(riverwidth) > 0 || error("river width must be positive on river cells")

    if river_routing == "kinematic-wave"
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
            dt = dt,
        )
    elseif river_routing == "local-inertial"
        rf, nodes_at_link = initialize_shallowwater_river(
            nc,
            config,
            inds_riv;
            graph = graph_riv,
            ldd = ldd_riv,
            dl = riverlength,
            width = riverwidth,
            reservoir_index = resindex,
            reservoir = reservoirs,
            lake_index = lakeindex,
            lake = lakes,
            dt = dt,
            floodplain = floodplain_1d,
        )
    else
        error(
            """An unknown "river_routing" method is specified in the TOML file ($river_routing).
            This should be "kinematic-wave" or "local-inertial".
            """,
        )
    end

    # unconfined aquifer
    if do_constanthead
        constanthead = ncread(
            nc,
            config,
            "lateral.subsurface.constant_head";
            sel = inds,
            type = Float,
            fill = mv,
        )
        index_constanthead = filter(i -> !isequal(constanthead[i], mv), 1:n)
        constant_head = ConstantHead(constanthead[index_constanthead], index_constanthead)
    else
        constant_head = ConstantHead{Float}(Float[], Int64[])
    end

    conductivity =
        ncread(nc, config, "lateral.subsurface.conductivity"; sel = inds, type = Float)
    specific_yield =
        ncread(nc, config, "lateral.subsurface.specific_yield"; sel = inds, type = Float)
    gwf_f = ncread(
        nc,
        config,
        "lateral.subsurface.gwf_f";
        sel = inds,
        type = Float,
        defaults = 3.0,
    )
    conductivity_profile =
        get(config.input.lateral.subsurface, "conductivity_profile", "uniform")::String

    connectivity = Connectivity(inds, rev_inds, xl, yl)
    initial_head = altitude .- sbm.zi / 1000.0 # cold state for groundwater head based on SBM zi
    initial_head[index_river] = altitude[index_river]

    if do_constanthead
        initial_head[constant_head.index] = constant_head.head
    end

    bottom = altitude .- sbm.soilthickness ./ Float(1000.0)
    area = xl .* yl
    volume = @. (min(altitude, initial_head) - bottom) * area * specific_yield # total volume than can be released

    aquifer = UnconfinedAquifer(
        initial_head,
        conductivity,
        altitude,
        bottom,
        area,
        specific_yield,
        zeros(Float, connectivity.nconnection),  # conductance
        volume,
        gwf_f,
    )

    # river boundary of unconfined aquifer
    infiltration_conductance = ncread(
        nc,
        config,
        "lateral.subsurface.infiltration_conductance";
        sel = inds_riv,
        type = Float,
    )
    exfiltration_conductance = ncread(
        nc,
        config,
        "lateral.subsurface.exfiltration_conductance";
        sel = inds_riv,
        type = Float,
    )
    river_bottom =
        ncread(nc, config, "lateral.subsurface.river_bottom"; sel = inds_riv, type = Float)

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
        drain_2d = ncread(nc, config, "lateral.subsurface.drain"; type = Bool, fill = false)

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
            config,
            "lateral.subsurface.drain_elevation";
            sel = inds,
            type = Float,
            fill = mv,
        )
        drain_conductance = ncread(
            nc,
            config,
            "lateral.subsurface.drain_conductance";
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
    if land_routing == "kinematic-wave" || river_routing == "kinematic-wave"
        streamorder = stream_order(graph, toposort)
    end
    if land_routing == "kinematic-wave"
        toposort = topological_sort_by_dfs(graph)
        index_pit_land = findall(x -> x == 5, ldd)
        min_streamorder_land = get(config.model, "min_streamorder_land", 5)
        subbas_order, indices_subbas, topo_subbas = kinwave_set_subdomains(
            graph,
            toposort,
            index_pit_land,
            streamorder,
            min_streamorder_land,
        )
    end
    if river_routing == "kinematic-wave"
        min_streamorder_river = get(config.model, "min_streamorder_river", 6)
        toposort_riv = topological_sort_by_dfs(graph_riv)
        index_pit_river = findall(x -> x == 5, ldd_riv)
        subriv_order, indices_subriv, topo_subriv = kinwave_set_subdomains(
            graph_riv,
            toposort_riv,
            index_pit_river,
            streamorder[index_river],
            min_streamorder_river,
        )
    end

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
        modelmap,
        indices_reverse,
        x_nc,
        y_nc,
        nc,
        extra_dim = (name = "layer", value = Float64.(1:sbm.maxlayers)),
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
    if land_routing == "kinematic-wave"
        land = (
            graph = graph,
            upstream_nodes = filter_upsteam_nodes(graph, pits[inds]),
            subdomain_order = subbas_order,
            topo_subdomain = topo_subbas,
            indices_subdomain = indices_subbas,
            order = toposort,
            indices = inds,
            reverse_indices = rev_inds,
            area = xl .* yl,
            slope = landslope,
            altitude = altitude,
            indices_allocation_areas = inds_allocation_areas,
        )
    elseif land_routing == "local-inertial"
        land = (
            graph = graph,
            order = toposort,
            indices = inds,
            reverse_indices = rev_inds,
            area = xl .* yl,
            slope = landslope,
            altitude = altitude,
            index_river = index_river_nf,
            staggered_indices = indices,
            indices_allocation_areas = inds_allocation_areas,
        )
    end
    if do_water_demand
        # exclude waterbodies for local surface and ground water abstraction
        inds_riv_2d = copy(rev_inds_riv)
        inds_2d = ones(Bool, modelsize_2d)
        if !isempty(reservoir)
            inds_cov = collect(Iterators.flatten(reservoir.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 0
        end
        if !isempty(lake)
            inds_cov = collect(Iterators.flatten(lake.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 0
        end
        land = merge(land, (index_river_wb = inds_riv_2d[inds], index_wb = inds_2d[inds]))
    end
    if river_routing == "kinematic-wave"
        river = (
            graph = graph_riv,
            indices = inds_riv,
            reverse_indices = rev_inds_riv,
            # reservoir and lake index
            reservoir_index = resindex,
            lake_index = lakeindex,
            reservoir_index_f = filter(x -> x ≠ 0, resindex),
            lake_index_f = filter(x -> x ≠ 0, lakeindex),
            # specific for kinematic_wave
            upstream_nodes = filter_upsteam_nodes(graph_riv, pits[inds_riv]),
            subdomain_order = subriv_order,
            topo_subdomain = topo_subriv,
            indices_subdomain = indices_subriv,
            order = toposort_riv,
            # water allocation areas
            indices_allocation_areas = inds_riv_allocation_areas,
            area = xl[index_river] .* yl[index_river],
        )
    elseif river_routing == "local-inertial"
        river = (
            graph = graph_riv,
            indices = inds_riv,
            reverse_indices = rev_inds_riv,
            # reservoir and lake index
            reservoir_index = resindex,
            lake_index = lakeindex,
            reservoir_index_f = filter(x -> x ≠ 0, resindex),
            lake_index_f = filter(x -> x ≠ 0, lakeindex),
            # specific for local-inertial
            nodes_at_link = nodes_at_link,
            links_at_node = adjacent_links_at_node(graph_riv, nodes_at_link),
            # water allocation areas
            indices_allocation_areas = inds_riv_allocation_areas,
            area = xl[index_river] .* yl[index_river],
        )
    end

    model = Model(
        config,
        (; land, river, reservoir, lake, drain, index_river, frac_toriver),
        (subsurface = subsurface_map, land = olf, river = rf),
        sbm,
        clock,
        reader,
        writer,
        SbmGwfModel(),
    )

    model = set_states(model)

    return model
end


"update the sbm_gwf model for a single timestep"
function update(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmGwfModel}
    @unpack lateral, vertical, network, clock, config = model

    do_water_demand = haskey(config.model, "water_demand")
    inds_riv = network.index_river
    aquifer = lateral.subsurface.flow.aquifer
    constanthead = lateral.subsurface.flow.constanthead

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
            network.land.slope,
            network.land,
        )
    end

    # optional water demand and allocation
    if do_water_demand
        update_water_demand(vertical)
        update_water_allocation(model)
    end

    # update vertical sbm concept until recharge [mm]
    update_until_recharge(vertical, config)

    # set river stage (groundwater) to average h from kinematic wave
    lateral.subsurface.river.stage .= lateral.river.h_av .+ lateral.subsurface.river.bottom

    # determine stable time step for groundwater flow
    conductivity_profile =
        get(config.input.lateral.subsurface, "conductivity_profile", "uniform")
    dt_gw = stable_timestep(aquifer, conductivity_profile) # time step in day (Float64)
    dt_sbm = (vertical.dt / tosecond(basetimestep)) # vertical.dt is in seconds (Float64)
    if dt_gw < dt_sbm
        @warn(
            "stable time step dt $dt_gw for groundwater flow is smaller than sbm dt $dt_sbm"
        )
    end

    Q = zeros(vertical.n)
    # exchange of recharge between vertical sbm concept and groundwater flow domain
    # recharge rate groundwater is required in units [m d⁻¹]
    @. lateral.subsurface.recharge.rate = vertical.recharge / 1000.0 * (1.0 / dt_sbm)
    if do_water_demand
        @. lateral.subsurface.recharge.rate -=
            vertical.allocation.act_groundwater_abst / 1000.0 * (1.0 / dt_sbm)
    end
    # update groundwater domain
    update(lateral.subsurface.flow, Q, dt_sbm, conductivity_profile)

    # determine excess water depth [m] (exfiltwater) in groundwater domain (head > surface)
    # and reset head
    exfiltwater = (aquifer.head .- min.(aquifer.head, aquifer.top)) .* storativity(aquifer)
    aquifer.head .= min.(aquifer.head, aquifer.top)

    # Adjust for constant head boundary of groundwater domain
    exfiltwater[constanthead.index] .= 0
    aquifer.head[constanthead.index] .= constanthead.head

    # update vertical sbm concept (runoff, ustorelayerdepth and satwaterdepth)
    update_after_subsurfaceflow(
        vertical,
        (network.land.altitude .- aquifer.head) .* 1000.0, # zi [mm] in vertical concept SBM
        exfiltwater .* 1000.0,
    )

    ssf_toriver = zeros(vertical.n)
    ssf_toriver[inds_riv] = -lateral.subsurface.river.flux ./ lateral.river.dt
    surface_routing(model, ssf_toriver = ssf_toriver)

    return model
end
