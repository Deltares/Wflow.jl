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
    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, config, "subcatchment"; optional = false, allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)

    river_2d =
        ncread(nc, config, "river_location"; optional = false, type = Bool, fill = false)
    river = river_2d[inds]

    # Needed to update the forcing
    reservoir = ()
    lake = ()

    soilloss = initialize_soil_loss(nc, config, inds)

    # Get waterbodies mask
    do_reservoirs = get(config.model, "doreservoir", false)::Bool
    do_lakes = get(config.model, "dolake", false)::Bool
    waterbodies = fill(0.0, n)
    if do_reservoirs
        reservoirs = ncread(
            nc,
            config,
            "reservoir_areas";
            optional = false,
            sel = inds,
            type = Float,
            fill = 0,
        )
        waterbodies = waterbodies .+ reservoirs
    end
    if do_lakes
        lakes = ncread(
            nc,
            config,
            "lake_areas";
            optional = false,
            sel = inds,
            type = Float,
            fill = 0,
        )
        waterbodies = waterbodies .+ lakes
    end
    waterbodies = waterbodies .> 0

    ldd_2d = ncread(nc, config, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[inds]

    # # lateral part sediment in overland flow
    overland_flow_sediment =
        initialize_overland_flow_sediment(nc, config, inds, waterbodies, river)

    graph = flowgraph(ldd, inds, pcr_dir)

    # River processes
    do_river = get(config.model, "runrivermodel", false)::Bool
    # TODO: see if we can skip init if the river model is not needed
    # or if we leave it when we restructure the Wflow Model struct

    inds_riv, rev_inds_riv = active_indices(river_2d, 0)

    ldd_riv = ldd_2d[inds_riv]
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # Needed for frac_to_river?
    landslope = ncread(
        nc,
        config,
        "vertical.geometry.slope";
        optional = false,
        sel = inds,
        type = Float,
    )
    clamp!(landslope, 0.00001, Inf)

    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, landslope, n)

    river_sediment = initialize_river_flow_sediment(nc, config, inds_riv, waterbodies)

    modelmap = (
        vertical = soilloss,
        lateral = (land = overland_flow_sediment, river = river_sediment),
    )
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
    )
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    writer = prepare_writer(config, modelmap, indices_reverse, x_nc, y_nc, nc)
    close(nc)

    # for each domain save the directed acyclic graph, the traversion order,
    # and the indices that map it back to the two dimensional grid
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

function update!(model::Model{N, L, V, R, W, T}) where {N, L, V, R, W, T <: SedimentModel}
    (; lateral, vertical, network, config, clock) = model
    dt = clock.dt

    # Soil erosion
    update!(vertical, dt)

    # Overland flow sediment transport
    update!(lateral.land, vertical.soil_erosion, network.land, dt)

    # River sediment transport
    do_river = get(config.model, "runrivermodel", false)::Bool
    if do_river
        inds_riv = network.index_river
        update!(lateral.river, lateral.land.to_river, network.river, inds_riv, dt)
    end

    return nothing
end

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
