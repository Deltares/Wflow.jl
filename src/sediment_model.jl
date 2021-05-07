"""
    initialize_sediment_model(config::Config)

Initial part of the sediment model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_sediment_model(config::Config)

    # unpack the paths to the NetCDF files
    tomldir = dirname(config)
    static_path = joinpath(tomldir, config.input.path_static)
    dynamic_path = joinpath(tomldir, config.input.path_forcing)

    reader = prepare_reader(dynamic_path, static_path, config)
    clock = Clock(config, reader)
    Δt = clock.Δt

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    reinit = get(config.model, "reinit", true)::Bool

    do_river = get(config.model, "runrivermodel", false)::Bool

    nc = NCDataset(static_path)
    dims = dimnames(nc[param(config, "input.subcatchment")])

    subcatch_2d = ncread(nc, param(config, "input.subcatchment"); allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)
    modelsize_2d = size(subcatch_2d)

    river_2d = ncread(nc, param(config, "input.river_location"); type = Bool, fill = false)
    river = river_2d[inds]
    riverwidth_2d =
        ncread(nc, param(config, "input.lateral.river.width"); type = Float, fill = 0)
    riverwidth = riverwidth_2d[inds]
    riverlength_2d =
        ncread(nc, param(config, "input.lateral.river.length"); type = Float, fill = 0)
    riverlength = riverlength_2d[inds]

    inds_riv, rev_inds_riv = active_indices(river_2d, 0)
    nriv = length(inds_riv)

    # Needed to update the forcing
    reservoir = ()
    lake = ()

    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    xl = fill(mv, n)
    yl = fill(mv, n)
    riverfrac = fill(mv, n)

    for i = 1:n
        xl[i] = sizeinmetres ? cellength : lattometres(y[i])[1] * cellength
        yl[i] = sizeinmetres ? cellength : lattometres(y[i])[2] * cellength
        riverfrac[i] =
            river[i] ? min((riverlength[i] * riverwidth[i]) / (xl[i] * yl[i]), 1.0) : 0.0
    end

    eros = initialize_landsed(nc, config, river, riverfrac, xl, yl, inds)

    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)
    ldd = ldd_2d[inds]

    # # lateral part sediment in overland flow
    rivcell = float(river)
    ols = OverlandFlowSediment{Float}(
        n = n,
        rivcell = rivcell,
        soilloss = fill(mv, n),
        erosclay = fill(mv, n),
        erossilt = fill(mv, n),
        erossand = fill(mv, n),
        erossagg = fill(mv, n),
        eroslagg = fill(mv, n),
        TCsed = fill(mv, n),
        TCclay = fill(mv, n),
        TCsilt = fill(mv, n),
        TCsand = fill(mv, n),
        TCsagg = fill(mv, n),
        TClagg = fill(mv, n),
        olsed = fill(mv, n),
        olclay = fill(mv, n),
        olsilt = fill(mv, n),
        olsand = fill(mv, n),
        olsagg = fill(mv, n),
        ollagg = fill(mv, n),
        inlandsed = fill(mv, n),
        inlandclay = fill(mv, n),
        inlandsilt = fill(mv, n),
        inlandsand = fill(mv, n),
        inlandsagg = fill(mv, n),
        inlandlagg = fill(mv, n),
    )

    graph = flowgraph(ldd, inds, pcr_dir)

    # River processes

    # the indices of the river cells in the land(+river) cell vector
    βₗ = ncread(nc, param(config, "input.lateral.land.slope"); sel = inds, type = Float)
    clamp!(βₗ, 0.00001, Inf)
    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]
    ldd_riv = ldd_2d[inds_riv]
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, βₗ, n)

    rs = initialize_riversed(nc, config, riverwidth, riverlength, inds_riv)

    modelmap = (vertical = eros, lateral = (land = ols, river = rs))
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
    )
    writer =
        prepare_writer(config, reader, modelmap, indices_reverse, x_nc, y_nc, nc)
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
        (land = ols, river = rs),
        eros,
        clock,
        reader,
        writer,
    )

    # read and set states in model object if reinit=false
    if reinit == false
        instate_path = joinpath(tomldir, config.state.path_input)
        state_ncnames = ncnames(config.state)
        set_states(instate_path, model, state_ncnames; type = Float)
    end

    # make sure the forcing is already loaded
    # it's fine to run twice, and may help catching errors earlier
    #update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end
    return model
end

function update(model::Model{N,L,V,R,W}) where {N,L,V<:LandSediment,R,W}
    @unpack lateral, vertical, network, clock, config = model

    update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end

    update_until_ols(vertical, config, network)
    update_until_oltransport(vertical, config, network)

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

    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
