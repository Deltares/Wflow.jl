"""
    initialize_sbm_gwf_model(config::Config)

Initial part of the sbm_gwf model concept. The model contains:
    - the land hydrology model with the SBM soil model
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

    # unconfined aquifer
    if do_constanthead
        constant_head = ConstantHead(dataset, config, indices)
    else
        variables = ConstantHeadVariables{Float}(; head = Float[])
        constant_head = ConstantHead{Float}(; variables, index = Int64[])
    end

    connectivity = Connectivity(indices, reverse_indices, x_length, y_length)

    initial_head = altitude .- land_hydrology.soil.variables.zi / 1000.0 # cold state for groundwater head based on SBM zi
    initial_head[inds_land_map2river] = altitude[inds_land_map2river]
    if do_constanthead
        initial_head[constant_head.index] = constant_head.variables.head
    end

    bottom = altitude .- land_hydrology.soil.parameters.soilthickness ./ Float(1000.0)
    area = x_length .* y_length
    conductance = zeros(Float, connectivity.nconnection)
    aquifer = UnconfinedAquifer(
        dataset,
        config,
        indices,
        altitude,
        bottom,
        area,
        conductance,
        initial_head,
    )

    # river boundary of unconfined aquifer
    river = River(dataset, config, inds_river, inds_land_map2river)

    # recharge boundary of unconfined aquifer
    recharge = Recharge(
        fill(mv, n_land_cells),
        zeros(Float, n_land_cells),
        collect(1:n_land_cells),
    )

    # drain boundary of unconfined aquifer (optional)
    if do_drains
        drain_2d = ncread(
            dataset,
            config,
            "routing.subsurface_flow.drain";
            type = Bool,
            fill = false,
        )

        drain = drain_2d[indices]
        # check if drain occurs where overland flow is not possible (surface_flow_width =
        # 0.0) and correct if this is the case
        false_drain = filter(
            i -> !isequal(drain[i], 0) && surface_flow_width[i] == Float(0),
            1:n_land_cells,
        )
        n_false_drain = length(false_drain)
        if n_false_drain > 0
            drain_2d[indices[false_drain]] .= 0
            drain[false_drain] .= 0
            @info "$n_false_drain drain locations are removed that occur where overland flow
             is not possible (overland flow width is zero)"
        end

        indices_drain, reverse_inds_drain = active_indices(drain_2d, 0)
        inds_land_map2drain = filter(i -> !isequal(drain[i], 0), 1:n_land_cells)

        drain = Drainage(dataset, config, indices, inds_land_map2drain)
        network_drain =
            NetworkDrain(; indices = indices_drain, reverse_indices = reverse_inds_drain)
        aquifer_boundaries = (; recharge, river, drain)
    else
        aquifer_boundaries = (; recharge, river)
        network_drain = NetworkDrain()
    end

    subsurface_flow = GroundwaterFlow{Float}(;
        aquifer,
        connectivity,
        constanthead = constant_head,
        boundaries = aquifer_boundaries,
    )

    # setup subdomains for the land and river kinematic wave domain, if nthreads = 1
    # subdomain is equal to the complete domain
    toposort = topological_sort_by_dfs(graph)
    if land_routing == "kinematic-wave" || river_routing == "kinematic-wave"
        streamorder = stream_order(graph, toposort)
    end
    if land_routing == "kinematic-wave"
        toposort = topological_sort_by_dfs(graph)
        land_pit_inds = findall(x -> x == 5, ldd)
        min_streamorder_land = get(config.model, "min_streamorder_land", 5)
        order_of_subdomains, subdomain_inds, toposort_subdomain = kinwave_set_subdomains(
            graph,
            toposort,
            land_pit_inds,
            streamorder,
            min_streamorder_land,
        )
    end
    if river_routing == "kinematic-wave"
        min_streamorder_river = get(config.model, "min_streamorder_river", 6)
        toposort_river = topological_sort_by_dfs(graph_river)
        river_pit_inds = findall(x -> x == 5, ldd_river)
        order_of_river_subdomains, river_subdomain_inds, toposort_river_subdomain =
            kinwave_set_subdomains(
                graph_river,
                toposort_river,
                river_pit_inds,
                streamorder[inds_land_map2river],
                min_streamorder_river,
            )
    end

    modelmap =
        (land = land_hydrology, routing = (; subsurface_flow, overland_flow, river_flow))
    indices_reverse = (
        land = reverse_indices,
        river = reverse_inds_river,
        reservoir = network_reservoir.reverse_indices,
        lake = isnothing(lake) ? nothing : lake.reverse_indices,
        drain = network_drain.reverse_indices,
    )
    writer = prepare_writer(
        config,
        modelmap,
        indices_reverse,
        x_coords,
        y_coords,
        dataset;
        extra_dim = (
            name = "layer",
            value = Float64.(1:(land_hydrology.soil.parameters.maxlayers)),
        ),
    )
    close(dataset)

    network_land = NetworkLand(;
        graph,
        order = toposort,
        indices,
        reverse_indices,
        area = x_length .* y_length,
        slope = land_slope,
        frac_to_river,
        altitude,
        allocation_area_indices = allocation_area_inds,
    )
    if land_routing == "kinematic-wave"
        @reset network_land.upstream_nodes = filter_upsteam_nodes(graph, pits[indices])
        @reset network_land.order_of_subdomains = order_of_subdomains
        @reset network_land.order_subdomain = toposort_subdomain
        @reset network_land.subdomain_indices = subdomain_inds
    elseif land_routing == "local-inertial"
        @reset network_land.river_indices = inds_river_map2land
        @reset network_land.edge_indices = edge_indices
    end
    if do_water_demand
        # exclude waterbodies for local surface and ground water abstraction
        inds_riv_2d = copy(reverse_inds_river)
        inds_2d = zeros(Bool, modelsize_2d)
        if do_reservoirs
            inds_cov = collect(Iterators.flatten(network_reservoir.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 1
        end
        if do_lakes
            inds_cov = collect(Iterators.flatten(network_lake.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 1
        end
        @reset network_land.river_inds_excl_waterbody = inds_riv_2d[indices]
        @reset network_land.waterbody = inds_2d[indices]
    end
    network_river = NetworkRiver(;
        graph = graph_river,
        indices = inds_river,
        reverse_indices = reverse_inds_river,
        reservoir_indices = inds_reservoir_map2river,
        lake_indices = inds_lake_map2river,
        land_indices = inds_land_map2river,
        allocation_area_indices = river_allocation_area_inds,
        cell_area = x_length[inds_land_map2river] .* y_length[inds_land_map2river],
    )
    if river_routing == "kinematic-wave"
        @reset network_river.upstream_nodes =
            filter_upsteam_nodes(graph_river, pits[inds_river])
        @reset network_river.order_of_subdomains = order_of_river_subdomains
        @reset network_river.order_subdomain = toposort_river_subdomain
        @reset network_river.subdomain_indices = river_subdomain_inds
        @reset network_river.order = toposort_river
    elseif river_routing == "local-inertial"
        @reset network_river.nodes_at_edge = NodesAtEdge(nodes_at_edge...)
        @reset network_river.edges_at_node =
            EdgesAtNode(adjacent_edges_at_node(graph_river, nodes_at_edge)...)
    end

    routing = Routing(; subsurface_flow, overland_flow, river_flow)
    return model
end

"update the sbm_gwf model for a single timestep"
function update!(model::AbstractModel{<:SbmGwfModel})
    (; routing, land, network, clock, config) = model
    (; soil, runoff, demand) = land

    do_water_demand = haskey(config.model, "water_demand")
    (; aquifer, boundaries) = routing.subsurface_flow
    dt = tosecond(clock.dt)

    update!(land, routing, network, config, dt)

    # set river stage (groundwater) to average h from kinematic wave
    boundaries.river.variables.stage .=
        routing.river_flow.variables.h_av .+ boundaries.river.parameters.bottom

    # determine stable time step for groundwater flow
    conductivity_profile =
        get(config.input.routing.subsurface_flow, "conductivity_profile", "uniform")
    dt_gw = stable_timestep(aquifer, conductivity_profile) # time step in day (Float64)
    dt_sbm = (dt / tosecond(basetimestep)) # dt is in seconds (Float64)
    if dt_gw < dt_sbm
        @warn(
            "stable time step dt $dt_gw for groundwater flow is smaller than `LandHydrologySBM` model dt $dt_sbm"
        )
    end

    Q = zeros(routing.subsurface_flow.connectivity.ncell)
    # exchange of recharge between SBM soil model and groundwater flow domain
    # recharge rate groundwater is required in units [m d⁻¹]
    @. boundaries.recharge.variables.rate =
        soil.variables.recharge / 1000.0 * (1.0 / dt_sbm)
    if do_water_demand
        @. boundaries.recharge.variables.rate -=
            land.allocation.variables.act_groundwater_abst / 1000.0 * (1.0 / dt_sbm)
    end
    # update groundwater domain
    update!(routing.subsurface_flow, Q, dt_sbm, conductivity_profile)

    # update SBM soil model (runoff, ustorelayerdepth and satwaterdepth)
    update!(soil, (; runoff, demand, subsurface_flow = routing.subsurface_flow))

    surface_routing!(model)

    return nothing
end
