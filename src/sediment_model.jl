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

    subcatch_2d =
        ncread(dataset, config, "subcatchment"; optional = false, allow_missing = true)
    # indices based on catchment
    indices_subcatch, reverse_indices_subcatch = active_indices(subcatch_2d, missing)
    number_of_cells = length(indices_subcatch)

    river_2d = ncread(
        dataset,
        config,
        "river_location";
        optional = false,
        type = Bool,
        fill = false,
    )
    river = river_2d[indices_subcatch]

    river_width_2d = ncread(
        dataset,
        config,
        "lateral.river.width";
        optional = false,
        type = Float,
        fill = 0,
    )
    river_width = river_width_2d[indices_subcatch]

    river_length_2d = ncread(
        dataset,
        config,
        "lateral.river.length";
        optional = false,
        type = Float,
        fill = 0,
    )
    river_length = river_length_2d[indices_subcatch]

    indices_river, reverse_indices_river = active_indices(river_2d, 0)

    # read x, y coordinates and calculate cell length [m]
    y_axis = read_y_axis(dataset)
    x_axis = read_x_axis(dataset)
    y = permutedims(repeat(y_axis; outer = (1, length(x_axis))))[indices_subcatch]
    cell_length = abs(mean(diff(x_axis)))

    size_in_metres = get(config.model, "sizeinmetres", false)::Bool
    x_length, y_length = cell_lengths(y, cell_length, size_in_metres)
    river_fraction_ = river_fraction(river, river_length, river_width, x_length, y_length)

    erosion = initialize_landsed(
        dataset,
        config,
        river,
        river_fraction_,
        x_length,
        y_length,
        indices_subcatch,
    )

    ldd_2d = ncread(dataset, config, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[indices_subcatch]

    # # lateral part sediment in overland flow
    overland_flow_sediment = init_overland_flow_sediment(river, number_of_cells)

    graph = flowgraph(ldd, indices_subcatch, pcr_dir)

    # River processes

    # the indices of the river cells in the land(+river) cell vector
    land_slope = ncread(
        dataset,
        config,
        "lateral.land.slope";
        optional = false,
        sel = indices_subcatch,
        type = Float,
    )
    clamp!(land_slope, 0.00001, Inf)

    river_length = river_length_2d[indices_river]
    minimum(river_length) > 0 || error("river length must be positive on river cells")

    river_width = river_width_2d[indices_river]
    minimum(river_width) > 0 || error("river width must be positive on river cells")

    ldd_river = ldd_2d[indices_river]
    graph_river = flowgraph(ldd_river, indices_river, pcr_dir)

    index_river = filter(i -> !isequal(river[i], 0), 1:number_of_cells)
    frac_to_river =
        fraction_runoff_toriver(graph, ldd, index_river, land_slope, number_of_cells)

    river_sediment =
        initialize_riversed(dataset, config, river_width, river_length, indices_river)

    model_map = (
        vertical = erosion,
        lateral = (land = overland_flow_sediment, river = river_sediment),
    )
    indices_reverse = (
        land = reverse_indices_subcatch,
        river = reverse_indices_river,
        reservoir = nothing,
        lake = nothing,
    )
    writer = prepare_writer(config, model_map, indices_reverse, x_axis, y_axis, dataset)
    close(dataset)

    # for each domain save the directed acyclic graph, the traversion order,
    # and the indices that map it back to the two dimensional grid
    land = (
        graph = graph,
        order = topological_sort_by_dfs(graph),
        indices = indices_subcatch,
        reverse_indices = reverse_indices_subcatch,
    )
    river = (
        graph = graph_river,
        order = topological_sort_by_dfs(graph_river),
        indices = indices_river,
        reverse_indices = reverse_indices_river,
    )

    sediment_model = SedimentModel()

    model = Model(
        config,
        (; land, river, reservoir = nothing, lake = nothing, index_river, frac_to_river),
        (land = overland_flow_sediment, river = river_sediment),
        erosion,
        clock,
        reader,
        writer,
        sediment_model,
    )

    model = set_states(model)
    @info "Initialized model"

    return model
end

function init_overland_flow_sediment(river, number_of_cells)
    river_cell = float(river)

    overland_flow_sediment = OverlandFlowSediment{Float}(;
        n = number_of_cells,
        rivcell = river_cell,
        soilloss = fill(mv, number_of_cells),
        erosclay = fill(mv, number_of_cells),
        erossilt = fill(mv, number_of_cells),
        erossand = fill(mv, number_of_cells),
        erossagg = fill(mv, number_of_cells),
        eroslagg = fill(mv, number_of_cells),
        TCsed = fill(mv, number_of_cells),
        TCclay = fill(mv, number_of_cells),
        TCsilt = fill(mv, number_of_cells),
        TCsand = fill(mv, number_of_cells),
        TCsagg = fill(mv, number_of_cells),
        TClagg = fill(mv, number_of_cells),
        olsed = fill(mv, number_of_cells),
        olclay = fill(mv, number_of_cells),
        olsilt = fill(mv, number_of_cells),
        olsand = fill(mv, number_of_cells),
        olsagg = fill(mv, number_of_cells),
        ollagg = fill(mv, number_of_cells),
        inlandsed = fill(mv, number_of_cells),
        inlandclay = fill(mv, number_of_cells),
        inlandsilt = fill(mv, number_of_cells),
        inlandsand = fill(mv, number_of_cells),
        inlandsagg = fill(mv, number_of_cells),
        inlandlagg = fill(mv, number_of_cells),
    )

    return overland_flow_sediment
end

function update(model::Model{N, L, V, R, W, T}) where {N, L, V, R, W, T <: SedimentModel}
    @unpack lateral, vertical, network, clock, config = model

    update_until_ols(vertical, config)
    update_until_oltransport(vertical, config)

    lateral.land.soilloss .= vertical.soilloss
    lateral.land.erosclay .= vertical.erosclay
    lateral.land.erossilt .= vertical.erossilt
    lateral.land.erossand .= vertical.erossand
    lateral.land.erossagg .= vertical.erossagg
    lateral.land.eroslagg .= vertical.eroslagg

    lateral.land.TCsed .= vertical.TCsed
    lateral.land.TCclay .= vertical.TCclay
    lateral.land.TCsilt .= vertical.TCsilt
    lateral.land.TCsand .= vertical.TCsand
    lateral.land.TCsagg .= vertical.TCsagg
    lateral.land.TClagg .= vertical.TClagg

    update(lateral.land, network.land, config)

    do_river = get(config.model, "runrivermodel", false)::Bool

    if do_river
        inds_riv = network.index_river
        lateral.river.inlandclay .= lateral.land.inlandclay[inds_riv]
        lateral.river.inlandsilt .= lateral.land.inlandsilt[inds_riv]
        lateral.river.inlandsand .= lateral.land.inlandsand[inds_riv]
        lateral.river.inlandsagg .= lateral.land.inlandsagg[inds_riv]
        lateral.river.inlandlagg .= lateral.land.inlandlagg[inds_riv]

        update(lateral.river, network.river, config)
    end

    return model
end

function set_states(
    model::Model{N, L, V, R, W, T},
) where {N, L, V, R, W, T <: SedimentModel}
    # read and set states in model object if reinit=false
    @unpack config = model
    reinit = get(config.model, "reinit", true)::Bool
    if reinit
        @info "Set initial conditions from default values."
    else
        instate_path = input_path(config, config.state.path_input)
        @info "Set initial conditions from state file `$instate_path`."
        set_states(instate_path, model; type = Float)
    end
    return model
end
