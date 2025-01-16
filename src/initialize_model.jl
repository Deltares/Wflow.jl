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
    )
    @info "General model settings" to_do...

    # Build nested NamedTuple with data that is required for multiple constructors
    sub_catchment_data = get_sub_catchment_data(config, dataset)
    cell_data = get_cell_data(config, dataset, sub_catchment_data)
    river_data = get_river_data(config, dataset, sub_catchment_data, cell_data)
    data = (; sub_catchment_data, cell_data, river_data)

    network_land = NetworkLand(config, dataset, data, type)
    network_reservoir, inds_reservoir_map2river, reservoir =
        NetworkReservoir(config, dataset, to_do, data, type)
    network_lake, inds_lake_map2river, lake =
        NetworkLake(config, dataset, to_do, data, type)

    land_hydrology = LandHydrologySBM(
        dataset,
        config,
        river_data.river_fraction,
        sub_catchment_data.indices,
    )

    reservoir_data = (; inds_reservoir_map2river)
    lake_data = (; inds_lake_map2river)
    data = (;
        sub_catchment_data,
        cell_data,
        river_data,
        reservoir_data,
        lake_data,
        land_hydrology,
    )
    network_river = NetworkRiver(config, dataset, data, to_do, type)

    network_drain = NetworkDrain(config, dataset, data, to_do, type)

    network = Network(;
        land = network_land,
        reservoir = network_reservoir,
        lake = network_lake,
        river = network_river,
        drain = network_drain,
    )

    model = Model(config, network, routing, land, clock, reader, writer, type)
    @info "Initialized model"

    set_states!(model)
    return model
end

function NetworkLake(
    config::Config,
    dataset::NCDataset,
    to_do::NamedTuple,
    data::NamedTuple,
    ::SbmModel,
)
    (; inds_river) = data.river_data
    n_river_cells = length(inds_river)
    if to_do.lakes
        lake, lake_network, inds_lake_map2river, pits =
            Lake(dataset, config, inds_river, n_river_cells, pits)
        network_lake = NetworkWaterBody(; lake_network...)
    else
        network_lake = NetworkWaterBody()
        inds_lake_map2river = fill(0, n_river_cells)
        lake = nothing
    end
    network_lake, inds_lake_map2river, lake
end

function get_sub_catchment_data(config::Config, dataset::NCDataset)
    subcatch_2d =
        ncread(dataset, config, "subcatchment"; optional = false, allow_missing = true)
    indices, reverse_indices = active_indices(subcatch_2d, missing)
    ldd_2d = ncread(dataset, config, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[indices]
    return (; subcatch_2d, indices, reverse_indices, ldd_2d, ldd)
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
    (; x_length, y_length)
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

function NetworkLand(
    config::Config,
    dataset::NCDataset,
    data::NamedTuple,
    ::SbmModel,
)::NetworkLand
    (; sub_catchment_data, cell_data, river_data) = data
    (; ldd, indices, reverse_indices, subcatch_2d) = sub_catchment_data
    subcatch_2d =
        ncread(dataset, config, "subcatchment"; optional = false, allow_missing = true)

    # indices based on sub-catchments
    graph = flowgraph(ldd, indices, pcr_dir)
    pits = zeros(Bool, size(subcatch_2d))

    toposort = topological_sort_by_dfs(graph)
    land_pit_inds = findall(x -> x == 5, ldd)
    streamorder = stream_order(graph, toposort)
    min_streamorder_land = get(config.model, "min_streamorder_land", 5)
    order_of_subdomains, subdomain_indices, toposort_subdomain = kinwave_set_subdomains(
        graph,
        toposort,
        land_pit_inds,
        streamorder,
        min_streamorder_land,
    )

    land_slope = get_land_slope(config, dataset, indices)

    frac_to_river =
        fraction_runoff_to_river(graph, ldd, river_data.inds_land_map2river, land_slope)

    return NetworkLand(;
        graph,
        upstream_nodes = filter_upsteam_nodes(graph, pits[indices]),
        order_of_subdomains,
        order_subdomain = toposort_subdomain,
        subdomain_indices = subdomain_indices,
        order = toposort,
        indices,
        reverse_indices,
        area = cell_data.x_length .* cell_data.y_length,
        slope = land_slope,
        frac_to_river,
    )
end

function NetworkReservoir(
    config::Config,
    dataset::NCDataset,
    to_do::NamedTuple,
    data::NamedTuple,
    ::SbmModel,
)::Tuple{NetworkWaterBody, Vector{Int}, Union{Nothing}}
    n_river_cells = length(data.river_data.inds_river)
    if to_do.reservoirs
        reservoir, reservoir_network, inds_reservoir_map2river, pits =
            SimpleReservoir(dataset, config, inds_river, n_river_cells, pits)
        network_reservoir = NetworkWaterBody(; reservoir_network...)
    else
        network_reservoir = NetworkWaterBody()
        inds_reservoir_map2river = zeros(n_river_cells)
        reservoir = nothing
    end
    return network_reservoir, inds_reservoir_map2river, reservoir
end

function get_river_data(
    config::Config,
    dataset::NCDataset,
    sub_catchment_data::NamedTuple,
    cell_data::NamedTuple,
)::NamedTuple
    (; indices) = sub_catchment_data
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
    river_width = river_width_2d[indices]
    river_length_2d = ncread(
        dataset,
        config,
        "routing.river_flow.length";
        optional = false,
        type = Float,
        fill = 0,
    )
    river_length = river_length_2d[indices]
    river_fraction = get_river_fraction(
        river_location,
        river_length,
        river_width,
        cell_data.x_length,
        cell_data.y_length,
    )

    inds_land_map2river = findall(!iszero, river_location)
    inds_river, reverse_inds_river = active_indices(river_location_2d, 0)

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
    )
end

function get_allocation_area_inds(land_hydrology::LandHydrologySBM, to_do::NamedTuple)
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

function NetworkRiver(
    config::Config,
    dataset::NCDataset,
    data::NamedTuple,
    to_do::NamedTuple,
    ::SbmModel,
)::NetworkRiver
    (;
        sub_catchment_data,
        cell_data,
        reservoir_data,
        lake_data,
        river_data,
        land_hydrology,
    ) = data
    (; inds_land_map2river) = river_data
    (; ldd_2d) = sub_catchment_data
    inds_river, reverse_inds_river = active_indices(river_data.river_location_2d, 0)

    ldd_river = ldd_2d[inds_river]
    if to_do.pits
        ldd_river = set_pit_ldd(pits_2d, ldd_river, inds_river)
    end

    allocation_area_inds, river_allocation_area_inds =
        get_allocation_area_inds(land_hydrology, to_do)

    graph_river = flowgraph(ldd_river, inds_river, pcr_dir)
    return NetworkRiver(;
        graph = graph_river,
        indices = inds_river,
        reverse_indices = reverse_inds_river,
        reservoir_indices = reservoir_data.inds_reservoir_map2river,
        land_indices = inds_land_map2river,
        lake_indices = lake_data.inds_lake_map2river,
        allocation_area_indices = river_allocation_area_inds,
        cell_area = cell_data.x_length[inds_land_map2river] .*
                    cell_data.y_length[inds_land_map2river],
    )
end

NetworkReservoir(config, dataset, data, to_do, ::AbstractModelType) = NetworkReservoir()
NetworkLake(config, dataset, data, to_do, ::AbstractModelType) = NetworkLake()
NetworkDrain(config, dataset, data, to_do, ::AbstractModelType) = NetworkDrain()