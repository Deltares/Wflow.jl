"""
    initialize_sediment_model(config::Config)

Initial part of the sediment model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_sediment_model(config::Config)
    model_type = config.model.type::String
    @info "Initialize model variables for model type `$model_type`."

    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)

    reader = prepare_reader(config)
    clock = Clock(config, reader)
    dataset = NCDataset(static_path)

    lens = lens_input("subcatchment_location__count")
    subcatch_2d = ncread(dataset, config, lens; optional = false, allow_missing = true)
    # indices based on catchment
    indices, rev_indices = active_indices(subcatch_2d, missing)
    n = length(indices)

    lens = lens_input("river_location__mask")
    river_2d = ncread(dataset, config, lens; optional = false, type = Bool, fill = false)
    river = river_2d[indices]

    soilloss = SoilLoss(dataset, config, indices)

    # Get waterbodies mask
    do_reservoirs = get(config.model, "doreservoir", false)::Bool
    do_lakes = get(config.model, "dolake", false)::Bool
    waterbodies = fill(0.0, n)
    if do_reservoirs
        lens = lens_input("reservoir_area__count")
        reservoirs = ncread(
            dataset,
            config,
            lens;
            optional = false,
            sel = indices,
            type = Float,
            fill = 0,
        )
        waterbodies = waterbodies .+ reservoirs
    end
    if do_lakes
        lens = lens_input("lake_area__count")
        lakes = ncread(
            dataset,
            config,
            lens;
            optional = false,
            sel = indices,
            type = Float,
            fill = 0,
        )
        waterbodies = waterbodies .+ lakes
    end
    waterbodies = waterbodies .> 0

    lens = lens_input("local_drain_direction")
    ldd_2d = ncread(dataset, config, lens; optional = false, allow_missing = true)
    ldd = ldd_2d[indices]

    # # sediment in overland flow
    overland_flow_sediment =
        OverlandFlowSediment(dataset, soilloss, config, indices, waterbodies, river)

    graph = flowgraph(ldd, indices, pcr_dir)

    # River processes
    indices_riv, rev_indices_riv = active_indices(river_2d, 0)

    ldd_riv = ldd_2d[indices_riv]
    graph_riv = flowgraph(ldd_riv, indices_riv, pcr_dir)
    index_river = filter(i -> !isequal(river[i], 0), 1:n)

    river_sediment = RiverSediment(dataset, config, indices_riv, waterbodies)

    modelmap = (
        land = soilloss,
        routing = (overland_flow = overland_flow_sediment, river_flow = river_sediment),
    )
    indices_reverse = (land = rev_indices, river = rev_indices_riv)
    y_dataset = read_y_axis(dataset)
    x_dataset = read_x_axis(dataset)
    writer =
        prepare_writer(config, modelmap, indices_reverse, x_dataset, y_dataset, dataset)
    close(dataset)

    network_land = NetworkLand(;
        graph = graph,
        order = topological_sort_by_dfs(graph),
        indices = indices,
        reverse_indices = rev_indices,
    )
    network_river = NetworkRiver(;
        graph = graph_riv,
        order = topological_sort_by_dfs(graph_riv),
        indices = indices_riv,
        land_indices = index_river,
        reverse_indices = rev_indices_riv,
    )

    network = Network(; land = network_land, river = network_river)

    routing = Routing(; overland_flow = overland_flow_sediment, river_flow = river_sediment)

    model =
        Model(config, network, routing, soilloss, clock, reader, writer, SedimentModel())

    set_states!(model)
    @info "Initialized model"

    return model
end

"update sediment model for a single timestep"
function update!(model::AbstractModel{<:SedimentModel})
    (; routing, land, network, config, clock) = model
    dt = tosecond(clock.dt)

    # Soil erosion
    update!(land, dt)

    # Overland flow sediment transport
    update!(routing.overland_flow, land.soil_erosion, network.land, dt)

    # River sediment transport
    do_river = get(config.model, "run_river_model", false)::Bool
    if do_river
        update!(routing.river_flow, routing.overland_flow.to_river, network.river, dt)
    end

    return nothing
end

"set the initial states of the sediment model"
function set_states!(model::AbstractModel{<:SedimentModel})
    # read and set states in model object if reinit=false
    (; config) = model
    reinit = get(config.model, "reinit", true)::Bool
    if reinit == false
        instate_path = input_path(config, config.state.path_input)
        @info "Set initial conditions from state file `$instate_path`."
        set_states!(instate_path, model; type = Float)
    else
        @info "Set initial conditions from default values."
    end
    return nothing
end
