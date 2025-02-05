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

    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)

    reader = prepare_reader(config)
    clock = Clock(config, reader)

    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_drains = get(config.model, "drains", false)::Bool
    do_constanthead = get(config.model, "constanthead", false)::Bool

    routing_options = ("kinematic-wave", "local-inertial")
    river_routing = get_options(
        config.model,
        "river_routing",
        routing_options,
        "kinematic-wave",
    )::String
    land_routing =
        get_options(config.model, "land_routing", routing_options, "kinematic-wave")::String
    do_water_demand = haskey(config.model, "water_demand")

    dataset = NCDataset(static_path)

    lens = lens_input(config, "subcatchment_location__count"; optional = false)
    subcatch_2d = ncread(dataset, config, lens; allow_missing = true)
    # indices based on catchment
    indices, reverse_indices = active_indices(subcatch_2d, missing)
    n_land_cells = length(indices)
    modelsize_2d = size(subcatch_2d)

    lens = lens_input(config, "river_location__mask"; optional = false)
    river_location_2d = ncread(dataset, config, lens; type = Bool, fill = false)
    river_location = river_location_2d[indices]

    lens = lens_input_parameter(config, "river__width"; optional = false)
    river_width_2d = ncread(dataset, config, lens; type = Float, fill = 0)
    river_width = river_width_2d[indices]

    lens = lens_input_parameter(config, "river__length"; optional = false)
    river_length_2d = ncread(dataset, config, lens; type = Float, fill = 0)
    river_length = river_length_2d[indices]

    lens = lens_input_parameter(config, "land_surface__elevation"; optional = false)
    altitude = ncread(dataset, config, lens; sel = indices, type = Float)

    # read x, y coordinates and calculate cell length [m]
    y_coords = read_y_axis(dataset)
    x_coords = read_x_axis(dataset)
    y = permutedims(repeat(y_coords; outer = (1, length(x_coords))))[indices]
    cell_length = abs(mean(diff(x_coords)))

    size_in_metres = get(config.model, "sizeinmetres", false)::Bool
    x_length, y_length = cell_lengths(y, cell_length, size_in_metres)
    river_fraction =
        get_river_fraction(river_location, river_length, river_width, x_length, y_length)

    inds_river, reverse_inds_river = active_indices(river_location_2d, 0)
    n_river_cells = length(inds_river)

    # initialize SBM concept
    land_hydrology = LandHydrologySBM(dataset, config, river_fraction, indices)

    # reservoirs
    pits = zeros(Bool, modelsize_2d)
    if do_reservoirs
        reservoir, reservoir_network, inds_reservoir_map2river, pits =
            SimpleReservoir(dataset, config, inds_river, n_river_cells, pits)
        network_reservoir = NetworkWaterBody(; reservoir_network...)
    else
        network_reservoir = NetworkWaterBody()
        inds_reservoir_map2river = fill(0, n_river_cells)
        reservoir = nothing
    end

    # lakes
    if do_lakes
        lake, lake_network, inds_lake_map2river, pits =
            Lake(dataset, config, inds_river, n_river_cells, pits)
        network_lake = NetworkWaterBody(; lake_network...)
    else
        network_lake = NetworkWaterBody()
        inds_lake_map2river = fill(0, n_river_cells)
        lake = nothing
    end

    # overland flow (kinematic wave)
    lens = lens_input_parameter(config, "land_surface__slope"; optional = false)
    land_slope = ncread(dataset, config, lens; sel = indices, type = Float)
    clamp!(land_slope, 0.00001, Inf)

    lens = lens_input(config, "local_drain_direction"; optional = false)
    ldd_2d = ncread(dataset, config, lens; allow_missing = true)
    ldd = ldd_2d[indices]

    flow_length = map(get_flow_length, ldd, x_length, y_length)
    flow_width = (x_length .* y_length) ./ flow_length
    surface_flow_width = map(det_surfacewidth, flow_width, river_width, river_location)

    graph = flowgraph(ldd, indices, pcr_dir)
    ldd_river = ldd_2d[inds_river]
    graph_river = flowgraph(ldd_river, inds_river, pcr_dir)

    # land indices where river is located
    inds_land_map2river = filter(i -> !isequal(river_location[i], 0), 1:n_land_cells)
    frac_to_river = fraction_runoff_to_river(graph, ldd, inds_land_map2river, land_slope)

    allocation_area_inds = Vector{Int}[]
    river_allocation_area_inds = Vector{Int}[]
    if do_water_demand
        areas = unique(land_hydrology.allocation.parameters.areas)
        for a in areas
            area_index = findall(x -> x == a, land_hydrology.allocation.parameters.areas)
            push!(allocation_area_inds, area_index)
            area_riv_index = findall(
                x -> x == a,
                land_hydrology.allocation.parameters.areas[inds_land_map2river],
            )
            push!(river_allocation_area_inds, area_riv_index)
        end
    end

    if land_routing == "kinematic-wave"
        overland_flow = KinWaveOverlandFlow(
            dataset,
            config,
            indices;
            slope = land_slope,
            flow_length,
            flow_width = surface_flow_width,
        )
    elseif land_routing == "local-inertial"
        inds_river_map2land = reverse_inds_river[indices] # not filtered (with zeros)
        overland_flow, edge_indices = LocalInertialOverlandFlow(
            dataset,
            config,
            indices;
            modelsize_2d,
            reverse_indices,
            x_length,
            y_length,
            river_width = river_width_2d[inds_river],
            graph_river,
            ldd_river,
            inds_river,
            river_location,
            waterbody = !=(0).(inds_reservoir_map2river + inds_lake_map2river),
        )
    end

    # river flow (kinematic wave)
    river_length = river_length_2d[inds_river]
    river_width = river_width_2d[inds_river]
    minimum(river_length) > 0 || error("river length must be positive on river cells")
    minimum(river_width) > 0 || error("river width must be positive on river cells")

    if river_routing == "kinematic-wave"
        river_flow = KinWaveRiverFlow(
            dataset,
            config,
            inds_river;
            river_length,
            river_width,
            reservoir,
            lake,
        )
    elseif river_routing == "local-inertial"
        river_flow, nodes_at_edge = LocalInertialRiverFlow(
            dataset,
            config,
            inds_river;
            graph_river,
            ldd_river,
            river_length,
            river_width,
            reservoir,
            lake,
            waterbody = !=(0).(inds_reservoir_map2river + inds_lake_map2river),
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
        lens = lens_input_parameter(config, "land_drain_location__flag")
        drain_2d = ncread(dataset, config, lens; type = Bool, fill = false)
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

    network = Network(;
        land = network_land,
        river = network_river,
        reservoir = network_reservoir,
        lake = network_lake,
        drain = network_drain,
    )

    routing = Routing(; subsurface_flow, overland_flow, river_flow)

    model = Model(
        config,
        network,
        routing,
        land_hydrology,
        clock,
        reader,
        writer,
        SbmGwfModel(),
    )

    set_states!(model)

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
    conductivity_profile = get(config.model, "conductivity_profile", "uniform")
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
