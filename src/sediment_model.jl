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

    lens = lens_input("subcatchment")
    subcatch_2d = ncread(dataset, config, lens; optional = false, allow_missing = true)
    # indices based on catchment
    indices, rev_indices = active_indices(subcatch_2d, missing)
    n = length(indices)

    lens = lens_input("river_location")
    river_2d = ncread(dataset, config, lens; optional = false, type = Bool, fill = false)
    river = river_2d[indices]

    # Needed to update the forcing
    reservoir = ()
    lake = ()

    soilloss = SoilLoss(dataset, config, indices)

    # Get waterbodies mask
    do_reservoirs = get(config.model, "doreservoir", false)::Bool
    do_lakes = get(config.model, "dolake", false)::Bool
    waterbodies = fill(0.0, n)
    if do_reservoirs
        lens = lens_input("reservoir_area__number")
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
        lens = lens_input("lake_area__number")
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

    lens = lens_input("ldd")
    ldd_2d = ncread(dataset, config, lens; optional = false, allow_missing = true)
    ldd = ldd_2d[indices]

    # # lateral part sediment in overland flow
    overland_flow_sediment =
        OverlandFlowSediment(dataset, soilloss, config, indices, waterbodies, river)

    graph = flowgraph(ldd, indices, pcr_dir)

    # River processes
    do_river = get(config.model, "run_river_model", false)::Bool
    # TODO: see if we can skip init if the river model is not needed
    # or if we leave it when we restructure the Wflow Model struct

    indices_riv, rev_indices_riv = active_indices(river_2d, 0)

    ldd_riv = ldd_2d[indices_riv]
    graph_riv = flowgraph(ldd_riv, indices_riv, pcr_dir)

    # Needed for frac_to_river?
    lens = lens_input_parameter("land_surface__slope")
    landslope = ncread(dataset, config, lens; optional = false, sel = indices, type = Float)
    clamp!(landslope, 0.00001, Inf)

    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_to_river(graph, ldd, index_river, landslope)

    river_sediment = RiverSediment(dataset, config, indices_riv, waterbodies)

    modelmap = (
        vertical = soilloss,
        lateral = (land = overland_flow_sediment, river = river_sediment),
    )
    indices_reverse = (
        land = rev_indices,
        river = rev_indices_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
    )
    y_dataset = read_y_axis(dataset)
    x_dataset = read_x_axis(dataset)
    writer =
        prepare_writer(config, modelmap, indices_reverse, x_dataset, y_dataset, dataset)
    close(dataset)

    # for each domain save the directed acyclic graph, the traversion order,
    # and the indices that map it back to the two dimensional grid
    land = (
        graph = graph,
        order = topological_sort_by_dfs(graph),
        indices = indices,
        reverse_indices = rev_indices,
    )
    river = (
        graph = graph_riv,
        order = topological_sort_by_dfs(graph_riv),
        indices = indices_riv,
        reverse_indices = rev_indices_riv,
    )

    model = Model(
        config,
        (; land, river, reservoir, lake, index_river, frac_toriver),
        (land = overland_flow_sediment, river = river_sediment),
        soilloss,
        clock,
        reader,
        writer,
        SedimentModel(),
    )

    set_states!(model)
    @info "Initialized model"

    return model
end

"update sediment model for a single timestep"
function update!(model::Model{N, L, V, R, W, T}) where {N, L, V, R, W, T <: SedimentModel}
    (; lateral, vertical, network, config, clock) = model
    dt = tosecond(clock.dt)

    # Soil erosion
    update!(vertical, dt)

    # Overland flow sediment transport
    update!(lateral.land, vertical.soil_erosion, network.land, dt)

    # River sediment transport
    do_river = get(config.model, "run_river_model", false)::Bool
    if do_river
        indices_riv = network.index_river
        update!(lateral.river, lateral.land.to_river, network.river, indices_riv, dt)
    end

    return nothing
end

"set the initial states of the sediment model"
function set_states!(
    model::Model{N, L, V, R, W, T},
) where {N, L, V, R, W, T <: SedimentModel}
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
