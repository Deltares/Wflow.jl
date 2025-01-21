const ROUTING_OPTIONS = (("kinematic-wave", "local-inertial"))

function Model(config::Config)::Model
    model_type = config.model.type

    if model_type ∉ ("sbm", "sbm_gwf", "sediment")
        error("Unknown model type $model_type.")
    end
    @info "Initialize model variables for model type `$model_type`."

    type = if model_type == "sbm"
        SbmModel()
    elseif model_type == "sbm_gwf"
        SbmGwfModel()
    elseif model_type == "sediment"
        SedimentModel()
    end

    return Model(config, type)
end

function Model(config::Config, type::AbstractSbmModelType)
    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)
    dataset = NCDataset(static_path)

    reader = prepare_reader(config)
    clock = Clock(config, reader)

    to_do = (;
        reservoirs = get(config.model, "reservoirs", false)::Bool,
        lakes = get(config.model, "lakes", false)::Bool,
        pits = get(config.model, "pits", false)::Bool,
        drains = get(config.model, "drains", false)::Bool,
        constanthead = get(config.model, "constanthead", false)::Bool,
        glacier = get(config.model, "glacier", false)::Bool,
        snow = get(config.model, "snow", false)::Bool,
        masswasting = get(config.model, "masswasting", false)::Bool,
        water_demand = haskey(config.model, "water_demand")::Bool,
        lateral_ssf = haskey(config.input.routing, "subsurface_flow")::Bool,
    )
    @info "General model settings" to_do...

    # Build nested NamedTuple with data that is required for multiple constructors
    sub_catchment_data = get_sub_catchment_data(config, dataset, type)
    cell_data = get_cell_data(config, dataset, sub_catchment_data)
    river_data = get_river_data(config, dataset, sub_catchment_data, cell_data)
    routing_types = get_routing_types(config)

    land_hydrology = LandHydrologySBM(
        dataset,
        config,
        river_data.river_fraction,
        sub_catchment_data.indices,
    )
    data = (; sub_catchment_data, cell_data, river_data, land_hydrology, routing_types)
    network, data = Network(config, dataset, data, to_do, type)
    routing, data = Routing(config, dataset, data, to_do, type)

    writer = get_writer(config, dataset, data)
    close(dataset)

    adjust_network_data!(data, to_do)

    model = Model(config, network, routing, land_hydrology, clock, reader, writer, type)
    @info "Initialized model"

    set_states!(model)
    return model
end

# Adjust for local-inertial and water demand
function adjust_network_data!(data::NamedTuple, to_do::NamedTuple)::Nothing
    (; network, routing_types) = data

    if routing_types.land_routing == "local-inertial"
        @reset network_land.river_indices = inds_river_map2land
        @reset network_land.edge_indices = edge_indices
    end
    if to_do.water_demand
        # exclude waterbodies for local surface and ground water abstraction
        inds_riv_2d = copy(reverse_inds_river)
        inds_2d = zeros(Bool, modelsize_2d)
        if !isnothing(reservoir)
            inds_cov = collect(Iterators.flatten(network_reservoir.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 1
        end
        if !isnothing(lake)
            inds_cov = collect(Iterators.flatten(lake_network.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 1
        end
        @reset network.land.river_inds_excl_waterbody = inds_riv_2d[indices]
        @reset network.land.waterbody = inds_2d[indices]
    end

    if routing_types.river_routing == "local-inertial"
        @reset network.river.nodes_at_edge = NodesAtEdge(; nodes_at_edge...)
        @reset network.river.edges_at_node =
            EdgesAtNode(; adjacent_edges_at_node(graph_river, nodes_at_edge)...)
    end
    return nothing
end

function get_routing_types(config::Config)
    land_routing =
        get_options(config.model, "land_routing", ROUTING_OPTIONS, "kinematic-wave")::String

    river_routing = get_options(
        config.model,
        "river_routing",
        ROUTING_OPTIONS,
        "kinematic-wave",
    )::String

    return (; land_routing, river_routing)
end

function get_writer(config::Config, dataset::NCDataset, data::NamedTuple)
    (;
        land_hydrology,
        routing,
        sub_catchment_data,
        cell_data,
        river_data,
        reservoir,
        lake,
        network,
    ) = data
    modelmap = (; land = land_hydrology, routing)
    indices_reverse = (
        land = sub_catchment_data.reverse_indices,
        river = river_data.reverse_inds_river,
        reservoir = isnothing(reservoir) ? nothing : network.reservoir.reverse_indices,
        lake = isnothing(lake) ? nothing : network.lake.reverse_indices,
    )
    (; maxlayers) = land_hydrology.soil.parameters
    writer = prepare_writer(
        config,
        modelmap,
        indices_reverse,
        cell_data.x_coords,
        cell_data.y_coords,
        dataset;
        extra_dim = (name = "layer", value = Float64.(1:(maxlayers))),
    )
    return writer
end

function get_subsurface_flow(
    config::Config,
    dataset::NCDataset,
    data::NamedTuple,
    to_do::NamedTuple,
    ::AbstractSbmModelType,
)
    (; sub_catchment_data, land_hydrology, cell_data, flow_data) = data
    (; indices, land_slope) = sub_catchment_data
    (; x_length, y_length) = cell_data
    (; flow_length, flow_width) = flow_data

    # check if lateral subsurface flow component is defined for the SBM model, when coupled
    # to another groundwater model, this component is not defined in the TOML file.
    if to_do.lateral_ssf
        subsurface_flow = LateralSSF(
            dataset,
            config,
            indices;
            soil = land_hydrology.soil,
            slope = land_slope,
            flow_length,
            flow_width,
            x_length,
            y_length,
        )
        # update variables `ssf`, `ssfmax` and `kh` (layered profile) based on ksat_profile
        kh_profile_type = get(config.input.land, "ksat_profile", "exponential")::String
        if kh_profile_type == "exponential" || kh_profile_type == "exponential_constant"
            initialize_lateralssf!(subsurface_flow, subsurface_flow.parameters.kh_profile)
        elseif kh_profile_type == "layered" || kh_profile_type == "layered_exponential"
            (; kv_profile) = land_hydrology.soil.parameters
            initialize_lateralssf!(
                subsurface_flow,
                land_hydrology.soil,
                kv_profile,
                tosecond(clock.dt),
            )
        end
    else
        # when the SBM model is coupled (BMI) to a groundwater model, the following
        # variables are expected to be exchanged from the groundwater model.
        subsurface_flow = GroundwaterExchange(n_land_cells)
    end
    return subsurface_flow
end

function get_overland_flow(
    config::Config,
    dataset::NCDataset,
    data::NamedTuple,
    ::AbstractSbmModelType,
)
    (; sub_catchment_data, flow_data, river_data, routing_types) = data
    (; land_routing) = routing_types
    (; indices, land_slope) = sub_catchment_data
    (; flow_length, flow_width) = flow_data
    (; river_width, river_location) = river_data

    if land_routing == "kinematic-wave"
        overland_flow = KinWaveOverlandFlow(
            dataset,
            config,
            indices;
            slope = land_slope,
            flow_length,
            flow_width = map(det_surfacewidth, flow_width, river_width, river_location),
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
            waterbody = !iszero.(inds_reservoir_map2river + inds_lake_map2river),
        )
    end
    return overland_flow
end

function get_river_flow(
    config::Config,
    dataset::NCDataset,
    data::NamedTuple,
    ::AbstractSbmModelType,
)
    (; river_data, reservoir, lake, routing_types) = data
    (; river_routing) = routing_types
    (; river_length, river_width, inds_river) = river_data

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
    return river_flow
end

function get_sub_catchment_data(
    config::Config,
    dataset::NCDataset,
    type::AbstractSbmModelType,
)
    subcatch_2d =
        ncread(dataset, config, "subcatchment"; optional = false, allow_missing = true)
    indices, reverse_indices = active_indices(subcatch_2d, missing)
    ldd_2d = ncread(dataset, config, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[indices]
    land_slope = get_land_slope(config, dataset, indices)
    altitude =
        type isa SbmModel ?
        ncread(dataset, config, "altitude"; optional = false, sel = indices, type = Float) :
        nothing
    return (; subcatch_2d, indices, reverse_indices, ldd_2d, ldd, land_slope, altitude)
end

function get_cell_data(
    config::Config,
    dataset::NCDataset,
    sub_catchment_data::NamedTuple,
)::NamedTuple
    (; indices) = sub_catchment_data
    y_coords = read_y_axis(dataset)
    x_coords = read_x_axis(dataset)
    y = permutedims(repeat(y_coords; outer = (1, length(x_coords))))[indices]
    cell_length = abs(mean(diff(x_coords)))

    size_in_metres = get(config.model, "sizeinmetres", false)::Bool
    x_length, y_length = cell_lengths(y, cell_length, size_in_metres)
    (; x_length, y_length, x_coords, y_coords)
end

function get_land_slope(
    config::Config,
    dataset::NCDataset,
    indices::Vector{CartesianIndex{2}},
)::Vector{Float64}
    land_slope = ncread(
        dataset,
        config,
        "routing.overland_flow.slope";
        optional = false,
        sel = indices,
        type = Float,
    )
    clamp!(land_slope, 0.00001, Inf)
    return land_slope
end

function get_kinwave_subdomains_land(
    config::Config,
    graph::DiGraph,
    data::NamedTuple,
    to_do::NamedTuple,
)::NamedTuple
    (; routing_types, sub_catchment_data) = data
    toposort = topological_sort_by_dfs(graph)
    streamorder = stream_order(graph, toposort)
    # setup subdomains for the land and river kinematic wave domain, if nthreads = 1
    # subdomain is equal to the complete domain
    if routing_types.land_routing == "kinematic-wave" || to_do.lateral_ssf
        land_pit_inds = findall(x -> x == 5, sub_catchment_data.ldd)
        min_streamorder_land = get(config.model, "min_streamorder_land", 5)
        order_of_subdomains, subdomain_inds, toposort_subdomain = kinwave_set_subdomains(
            graph,
            toposort,
            land_pit_inds,
            streamorder,
            min_streamorder_land,
        )
    else
        order_of_subdomains = Vector{Int}[]
        subdomain_inds = Vector{Int}[]
        toposort_subdomain = nothing
        min_streamorder_land = nothing
    end
    return (;
        toposort,
        order_of_subdomains,
        subdomain_inds,
        toposort_subdomain,
        min_streamorder_land,
        streamorder,
    )
end

function get_kinwave_subdomains_river(
    config::Config,
    graph_river::DiGraph,
    data::NamedTuple,
)
    (; routing_types, river_data, streamorder) = data
    (; ldd_river, inds_land_map2river) = river_data
    if routing_types.river_routing == "kinematic-wave"
        min_streamorder_river = get(config.model, "min_streamorder_river", 6)
        toposort_river = topological_sort_by_dfs(graph_river)
        river_pit_inds = findall(x -> x == 5, ldd_river)
        order_of_river_subdomains, river_subdomain_inds, toposort_river_subdomain =
            order_of_river_subdomains, river_subdomain_inds, toposort_river_subdomain =
                kinwave_set_subdomains(
                    graph_river,
                    toposort_river,
                    river_pit_inds,
                    streamorder[inds_land_map2river],
                    min_streamorder_river,
                )
    else
        order_of_river_subdomains = Vector{Int}[]
        river_subdomain_inds = Vector{Int}[]
        toposort_river = nothing
    end
    return (;
        min_streamorder_river,
        order_of_river_subdomains,
        river_subdomain_inds,
        toposort_river,
    )
end

function get_river_data(
    config::Config,
    dataset::NCDataset,
    sub_catchment_data::NamedTuple,
    cell_data::NamedTuple,
)::NamedTuple
    (; indices, ldd_2d) = sub_catchment_data
    river_location_2d = ncread(
        dataset,
        config,
        "river_location";
        optional = false,
        type = Bool,
        fill = false,
    )
    river_location = river_location_2d[indices]
    river_width_2d = ncread(
        dataset,
        config,
        "routing.river_flow.width";
        optional = false,
        type = Float,
        fill = 0,
    )
    river_length_2d = ncread(
        dataset,
        config,
        "routing.river_flow.length";
        optional = false,
        type = Float,
        fill = 0,
    )

    inds_land_map2river = findall(!iszero, river_location)
    inds_river, reverse_inds_river = active_indices(river_location_2d, 0)
    ldd_river = ldd_2d[inds_river]

    river_width = river_width_2d[indices]
    river_length = river_length_2d[indices]

    river_fraction = get_river_fraction(
        river_location,
        river_length,
        river_width,
        cell_data.x_length,
        cell_data.y_length,
    )

    @assert minimum(river_length) > 0 "river length must be positive on river cells"
    @assert minimum(river_width) > 0 "river width must be positive on river cells"

    (;
        river_location_2d,
        river_location,
        river_width_2d,
        river_width,
        river_length_2d,
        river_length,
        river_fraction,
        inds_land_map2river,
        inds_river,
        reverse_inds_river,
        ldd_river,
    )
end

function get_allocation_area_inds(
    land_hydrology::LandHydrologySBM,
    data::NamedTuple,
    to_do::NamedTuple,
)
    (; inds_land_map2river) = data.river_data
    allocation_area_inds = Vector{Int}[]
    river_allocation_area_inds = Vector{Int}[]
    if to_do.water_demand
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
    return allocation_area_inds, river_allocation_area_inds
end