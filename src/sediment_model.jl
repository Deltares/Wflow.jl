"""
    initialize_sediment_model(config::Config)

Initial part of the sediment model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_sediment_model(config::Config)
    model_type = config.model.type::String
    @info "Initialize model variables for model type `$model_type`."

    # Unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)
    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, config, "subcatchment"; optional = false, allow_missing = true)
    # indices based on catchment
    indices_subcatch, reverse_indices_subcatch = active_indices(subcatch_2d, missing)
    n = length(indices_subcatch)

    river_2d =
        ncread(nc, config, "river_location"; optional = false, type = Bool, fill = false)
    river = river_2d[indices_subcatch]
    river_width_2d =
        ncread(nc, config, "lateral.river.width"; optional = false, type = Float, fill = 0)
    river_width = river_width_2d[indices_subcatch]
    river_length_2d =
        ncread(nc, config, "lateral.river.length"; optional = false, type = Float, fill = 0)
    river_length = river_length_2d[indices_subcatch]

    indices_river, reverse_indices_river = active_indices(river_2d, 0)

    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[indices_subcatch]
    cell_length = abs(mean(diff(x_nc)))

    size_in_metres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cell_length, size_in_metres)
    river_frac = river_fraction(river, river_length, river_width, xl, yl)

    ldd_2d = ncread(nc, config, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[indices_subcatch]

    # lateral part sediment in overland flow
    ols = overland_flow_sediment(river, mv, n)

    graph = flowgraph(ldd, indices_subcatch, pcr_dir)

    # River processes

    # the indices of the river cells in the land(+river) cell vector
    land_slope =
        ncread(nc, config, "lateral.land.slope"; optional = false, sel = indices_subcatch, type = Float)
    clamp!(land_slope, 0.00001, Inf)

    river_length = river_length_2d[indices_river]
    river_width = river_width_2d[indices_river]
    minimum(river_length) > 0 || error("river length must be positive on river cells")
    minimum(river_width) > 0 || error("river width must be positive on river cells")

    rs = initialize_riversed(nc, config, river_width, river_length, indices_river)

    ldd_riv = ldd_2d[indices_river]
    graph_riv = flowgraph(ldd_riv, indices_river, pcr_dir)

    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_to_river = fraction_runoff_toriver(graph, ldd, index_river, land_slope, n)

    river = (
        graph = graph_riv,
        order = topological_sort_by_dfs(graph_riv),
        indices = indices_river,
        reverse_indices = reverse_indices_river,
    )

    eros = initialize_landsed(nc, config, river, river_frac, xl, yl, indices_subcatch)

    modelmap = (vertical = eros, lateral = (land = ols, river = rs))
    indices_reverse = (
        land = reverse_indices_subcatch,
        river = reverse_indices_river,
        reservoir = nothing,
        lake = nothing,
    )

    # for each domain save the directed acyclic graph, the traversion order,
    # and the indices that map it back to the two dimensional grid
    land = (
        graph = graph,
        order = topological_sort_by_dfs(graph),
        indices = indices_subcatch,
        reverse_indices = reverse_indices_subcatch,
    )

    reader = prepare_reader(config)

    clock = Clock(config, reader)

    writer = prepare_writer(config, modelmap, indices_reverse, x_nc, y_nc, nc)

    close(nc)

    model = Model(
        config,
        (; land, river, reservoir=(), lake=(), index_river, frac_to_river),
        (land = ols, river = rs),
        eros,
        clock,
        reader,
        writer,
        SedimentModel(),
    )

    model = set_states(model)
    @info "Initialized model"

    return model
end

function overland_flow_sediment(T::Type{<:AbstractFloat}, river_cell, number_of_elements)
    return OverlandFlowSediment{T}(
        n = number_of_elements,
        rivcell = float(river_cell),
        soilloss = fill(T(NaN), number_of_elements),
        erosclay = fill(T(NaN), number_of_elements),
        erossilt = fill(T(NaN), number_of_elements),
        erossand = fill(T(NaN), number_of_elements),
        erossagg = fill(T(NaN), number_of_elements),
        eroslagg = fill(T(NaN), number_of_elements),
        TCsed = fill(T(NaN), number_of_elements),
        TCclay = fill(T(NaN), number_of_elements),
        TCsilt = fill(T(NaN), number_of_elements),
        TCsand = fill(T(NaN), number_of_elements),
        TCsagg = fill(T(NaN), number_of_elements),
        TClagg = fill(T(NaN), number_of_elements),
        olsed = fill(T(NaN), number_of_elements),
        olclay = fill(T(NaN), number_of_elements),
        olsilt = fill(T(NaN), number_of_elements),
        olsand = fill(T(NaN), number_of_elements),
        olsagg = fill(T(NaN), number_of_elements),
        ollagg = fill(T(NaN), number_of_elements),
        inlandsed = fill(T(NaN), number_of_elements),
        inlandclay = fill(T(NaN), number_of_elements),
        inlandsilt = fill(T(NaN), number_of_elements),
        inlandsand = fill(T(NaN), number_of_elements),
        inlandsagg = fill(T(NaN), number_of_elements),
        inlandlagg = fill(T(NaN), number_of_elements),
    )
end
overland_flow_sediment(river_cell, number_of_elements) = overland_flow_sediment(Float, river_cell, number_of_elements)

function update(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SedimentModel}
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

function set_states(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SedimentModel}
    # read and set states in model object if reinit=false
    @unpack config = model
    reinit = get(config.model, "reinit", true)::Bool
    if reinit == false
        instate_path = input_path(config, config.state.path_input)
        @info "Set initial conditions from state file `$instate_path`."
        set_states(instate_path, model; type = Float)
    else
        @info "Set initial conditions from default values."
    end
    return model
end
