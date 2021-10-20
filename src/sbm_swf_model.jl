"""
    initialize_sbm_swf_model(config::Config)

Initial part of the SBM + Shallow water flow model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""

# TODO (general): add docs/comments to functions and structs related to this new concept

@with_kw struct Indices
    xu::Vector{Int}
    xd::Vector{Int}
    yu::Vector{Int}
    yd::Vector{Int}
end

const dirs = ["yd", "xd", "xu", "yu"]

function initialize_sbm_swf_model(config::Config)

    # unpack the paths to the NetCDF files
    tomldir = dirname(config)
    static_path = joinpath(tomldir, config.input.path_static)
    dynamic_path = joinpath(tomldir, config.input.path_forcing)

    reader = prepare_reader(config)
    clock = Clock(config, reader)
    Δt = clock.Δt

    reinit = get(config.model, "reinit", true)::Bool
    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool

    nc = NCDataset(static_path)

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
    altitude_2d = ncread(nc, param(config, "input.altitude"); type = Float, fill = 0)
    elevation = altitude_2d[inds]
    river_elevation_2d = ncread(
        nc,
        param(config, "input.lateral.river.elevation");
        type = Wflow.Float,
        fill = 0,
    )

    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)
    riverfrac = river_fraction(river, riverlength, riverwidth, xl, yl)

    sbm = initialize_sbm(nc, config, riverfrac, inds)

    inds_riv, rev_inds_riv = active_indices(river_2d, 0)
    nriv = length(inds_riv)

    pits = zeros(Bool, modelsize_2d)

    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)
    ldd = ldd_2d[inds]

    βₗ = ncread(nc, param(config, "input.lateral.land.slope"); sel = inds, type = Float)
    clamp!(βₗ, 0.00001, Inf)
    dl = fill(mv, n)
    dw = fill(mv, n)

    for i = 1:n
        dl[i] = detdrainlength(ldd[i], xl[i], yl[i])
        dw[i] = (xl[i] * yl[i]) / dl[i]
    end

    # check if lateral subsurface flow component is defined for the SBM model, when coupled 
    # to another groundwater model, this component is not defined in the TOML file.
    if haskey(config.input.lateral, "subsurface")
        khfrac = ncread(
            nc,
            param(config, "input.lateral.subsurface.ksathorfrac", nothing);
            sel = inds,
            defaults = 1.0,
            type = Float,
        )

        # unit for lateral subsurface flow component is [m]
        kh₀ = khfrac .* sbm.kv₀ .* 0.001
        f = sbm.f .* 1000.0
        zi = sbm.zi .* 0.001
        soilthickness = sbm.soilthickness .* 0.001

        ssf = LateralSSF{Float}(
            kh₀ = kh₀,
            f = f,
            zi = zi,
            soilthickness = soilthickness,
            θₛ = sbm.θₛ,
            θᵣ = sbm.θᵣ,
            Δt = tosecond(Δt),
            t = 1.0,
            βₗ = βₗ,
            dl = dl,
            dw = dw,
            exfiltwater = fill(mv, n),
            recharge = fill(mv, n),
            ssf = ((kh₀ .* βₗ) ./ f) .* (exp.(-f .* zi) - exp.(-f .* soilthickness)) .* dw,
            ssfin = fill(mv, n),
            ssfmax = ((kh₀ .* βₗ) ./ f) .* (1.0 .- exp.(-f .* soilthickness)),
            to_river = zeros(n),
            wb_pit = pits[inds],
        )
    else
        # when the SBM model is coupled (BMI) to a groundwater model, the following
        # variables are expected to be exchanged from the groundwater model.
        ssf = GroundwaterExchange{Float}(
            Δt = tosecond(Δt),
            exfiltwater = fill(mv, n),
            zi = fill(mv, n),
            to_river = fill(mv, n),
            ssf = zeros(n),
        )
    end

    n_land = ncread(
        nc,
        param(config, "input.lateral.land.n", nothing);
        sel = inds,
        defaults = 0.072,
        type = Float,
    )

    graph = flowgraph(ldd, inds, pcr_dir)

    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]
    n_river = ncread(
        nc,
        param(config, "input.lateral.river.n", nothing);
        sel = inds_riv,
        defaults = 0.036,
        type = Float,
    )

    ldd_riv = ldd_2d[inds_riv]
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)
    river_elevation = river_elevation_2d[inds_riv]

    # the indices of the river cells in the vector with active cells (sub-catchment)
    index_river_f = filter(i -> !isequal(river[i], 0), 1:n) # filtered (without zeros)
    index_river = rev_inds_riv[inds] # not filtered (with zeros)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river_f, βₗ, n)

    # TODO: add more comments/ explanation for local inertial model

    # land part
    indices = Indices(xu = fill(0, n), xd = fill(0, n), yu = fill(0, n), yd = fill(0, n))

    nrow, ncol = modelsize_2d
    for (v, i) in enumerate(inds)
        for (m, neighbor) in enumerate(neighbors)
            j = i + neighbor
            if (1 <= j[1] <= nrow) && (1 <= j[2] <= ncol && rev_inds[j] != 0)
                param(indices, dirs[m])[v] = rev_inds[j]
            else
                param(indices, dirs[m])[v] = n + 1
            end
        end
    end

    zx_max = fill(Float(0), n)
    zy_max = fill(Float(0), n)
    for i = 1:n
        xu = indices.xu[i]
        if xu <= n
            zx_max[i] = max(elevation[i], elevation[xu])
        end
        yu = indices.yu[i]
        if yu <= n
            zy_max[i] = max(elevation[i], elevation[yu])
        end
    end

    # loop over cells containing subgrid river to set effective flow width for surface water
    # flow
    we_y = copy(yl)
    we_x = copy(xl)

    toposort_riv = topological_sort_by_dfs(graph_riv)

    for v in toposort_riv
        dst = outneighbors(graph_riv, v)
        isempty(dst) && continue
        w = min(riverwidth[v], riverwidth[only(dst)])
        dir = Wflow.pcr_dir[ldd_riv[v]]
        idx = rev_inds[inds_riv[v]]
        if dir == CartesianIndex(1, 1)
            we_x[idx] = we_x[idx] - 0.5 * w
            we_y[idx] = we_y[idx] - 0.5 * w
        elseif dir == CartesianIndex(-1, -1)
            we_x[indices.xd[idx]] = we_x[indices.xd[idx]] - 0.5 * w
            we_y[indices.yd[idx]] = we_y[indices.yd[idx]] - 0.5 * w
        elseif dir == CartesianIndex(1, 0)
            we_x[idx] = we_x[idx] - w
        elseif dir == CartesianIndex(0, 1)
            we_y[idx] = we_y[idx] - w
        elseif dir == CartesianIndex(-1, 0)
            we_x[indices.xd[idx]] = we_x[indices.xd[idx]] - w
        elseif dir == CartesianIndex(0, -1)
            we_y[indices.yd[idx]] = we_x[indices.yd[idx]] - w
        elseif dir == CartesianIndex(1, -1)
            we_x[idx] = we_x[idx] - 0.5 * w
            we_y[indices.yd[idx]] = we_x[indices.yd[idx]] - 0.5 * w
        elseif dir == CartesianIndex(-1, 1)
            we_x[indices.xd[idx]] = we_x[indices.xd[idx]] - 0.5 * w
            we_y[idx] = we_y[idx] - 0.5 * w
        end
    end

    bankheight = fill(Float(2), nriv)    # dummy value for bankheight
    bankvolume = bankheight .* riverwidth .* riverlength
    volume = fill(Float(0), n)

    # TODO: review field of struct
    sw_land = ShallowWaterLand{Float}(
        n = n,
        xl = xl, # also part of landnetwork...
        yl = yl, # also part of landnetwork...
        xwidth = we_x,
        ywidth = we_y,
        g = 9.80665,
        θ = 0.8,
        α = 0.7,
        qx0 = zeros(n + 1),
        qy0 = zeros(n + 1),
        qx = zeros(n + 1),
        qy = zeros(n + 1),
        zx_max = zx_max,
        zy_max = zy_max,
        mannings_n = n_land,
        volume = volume,
        runoff = zeros(n),
        h = zeros(n),
        h_av = zeros(n),
        slp_x = zeros(n),
        slp_y = zeros(n),
        z = elevation,
        froude = true,
        rivercells = river,
    )
    
    riverlength_bc = 1e05
    alpha = 0.7
    h_thresh = 1.0e-3
    froude_limit = true

    # set ghost points for boundary condition (downstream river outlet): river width and
    # manning n is copied from the upstream cell, river elevation and h are set at 0.0
    # (sea level). river length at boundary point is by default 1.0e5 m
    # (riverlength_bc).
    index_pit_river = findall(x -> x == 5, ldd_riv)
    n_ghost_points = length(index_pit_river)
    for (i, v) in enumerate(index_pit_river)
        add_vertex!(graph_riv)
        add_edge!(graph_riv, v, nriv + i)
        append!(river_elevation, 0.0)
        append!(riverwidth, riverwidth[v])
        append!(riverlength, riverlength_bc)
        append!(n_river, n_river[v])
    end

    # for each link the src and dst node is required
    nodes_at_link = adjacent_nodes_at_link(graph_riv)
    # for each node the src and dst link is required
    links_at_node = adjacent_links_at_node(graph_riv, nodes_at_link)

    _ne = ne(graph_riv)

    zmax = fill(Float(0), _ne)
    width_at_link = fill(Float(0), _ne)
    length_at_link = fill(Float(0), _ne)
    mannings_n = fill(Float(0), _ne)
    for i = 1:_ne
        zmax[i] = max(
            river_elevation[nodes_at_link.src[i]],
            river_elevation[nodes_at_link.dst[i]],
        )
        width_at_link[i] =
            min(riverwidth[nodes_at_link.dst[i]], riverwidth[nodes_at_link.src[i]])
        length_at_link[i] =
            0.5 * (riverlength[nodes_at_link.dst[i]] + riverlength[nodes_at_link.src[i]])
        mannings_n[i] =
            (
                n_river[nodes_at_link.dst[i]] * riverlength[nodes_at_link.dst[i]] +
                n_river[nodes_at_link.src[i]] * riverlength[nodes_at_link.src[i]]
            ) / (riverlength[nodes_at_link.dst[i]] + riverlength[nodes_at_link.src[i]])
    end

    sw_river = ShallowWaterRiver(
        n = nriv,
        ne = _ne,
        g = 9.80665,
        α = alpha,
        h_thresh = h_thresh,
        Δt = tosecond(Δt),
        q = zeros(_ne),
        q_av = zeros(_ne),
        zmax = zmax,
        mannings_n = mannings_n,
        h = fill(0.0, nriv + n_ghost_points),
        η_max = zeros(_ne),
        hf = zeros(_ne),
        h_av = zeros(nriv),
        width = riverwidth,
        width_at_link = width_at_link,
        a = zeros(_ne),
        r = zeros(_ne),
        volume = fill(0.0, nriv),
        error = zeros(Float, nriv),
        inwater = zeros(nriv),
        inwater0 = fill(mv, nriv),
        dl = riverlength,
        dl_at_link = length_at_link,
        bankvolume = bankvolume,
        bankheight = bankheight,
        z = river_elevation,
        froude_limit = froude_limit,
        reservoir_index = do_reservoirs ? resindex : fill(0, nriv),
        lake_index = do_lakes ? lakeindex : fill(0, nriv),
        reservoir = do_reservoirs ? reservoirs : nothing,
        lake = do_lakes ? lakes : nothing,
    )

    # setup subdomains for the land kinematic wave domain, if nthreads = 1
    # subdomain is equal to the complete domain
    toposort = topological_sort_by_dfs(graph)
    index_pit_land = findall(x -> x == 5, ldd)
    subbas_order, indices_subbas, topo_subbas =
        kinwave_set_subdomains(config, graph, toposort, index_pit_land)

    modelmap =
        (vertical = sbm, lateral = (subsurface = ssf, land = sw_land, river = sw_river))

    indices_reverse = (land = rev_inds, river = rev_inds_riv)

    writer = prepare_writer(
        config,
        reader,
        modelmap,
        indices_reverse,
        x_nc,
        y_nc,
        nc,
        maxlayers = sbm.maxlayers,
    )
    close(nc)

    # for each domain save:
    # - the directed acyclic graph (graph),
    # - the traversion order (order),
    # - upstream_nodes,
    # - subdomains for the kinematic wave domains for parallel execution (execution order of
    #   subbasins (subdomain_order), traversion order per subbasin (topo_subdomain) and
    #   Vector indices per subbasin matching the traversion order of the complete domain
    #   (indices_subdomain)) 
    # - the indices that map it back to the two dimensional grid (indices)

    # for the land domain the x and y length [m] of the grid cells are stored
    # for reservoirs and lakes indices information is available from the initialization
    # functions
    # TODO: review field of land and river network for local inertial model (1D and 2D) 
    land = (
        graph = graph,
        index_pit_land = index_pit_land,
        upstream_nodes = filter_upsteam_nodes(graph, ssf.wb_pit),
        subdomain_order = subbas_order,
        topo_subdomain = topo_subbas,
        indices_subdomain = indices_subbas,
        order = toposort,
        indices = inds,
        reverse_indices = rev_inds,
        staggered_indices = indices,
        xl = xl,
        yl = yl,
    )

    river = (
        graph = graph_riv,
        nodes_at_link = nodes_at_link,
        links_at_node = links_at_node,
        indices = inds_riv,
        reverse_indices = rev_inds_riv,
    )

    model = Model(
        config,
        (; land, river, index_river, index_river_f, frac_toriver),
        (subsurface = ssf, land = sw_land, river = sw_river),
        sbm,
        clock,
        reader,
        writer,
    )

    # read and set states in model object if reinit=false
    # TODO: add/check states for local inertial model 
    if reinit == false
        instate_path = joinpath(tomldir, config.state.path_input)
        state_ncnames = ncnames(config.state)
        set_states(instate_path, model, state_ncnames; type = Float)
        @unpack lateral, vertical = model
        # update zi for vertical sbm
        zi =
            max.(
                0.0,
                vertical.soilthickness .-
                vertical.satwaterdepth ./ (vertical.θₛ .- vertical.θᵣ),
            )
        vertical.zi .= zi
    end

    return model
end

"update the sbm_swf model for a single timestep"
function update_sbm_swf(model)

    @unpack lateral, vertical, network, clock, config = model
    inds_riv = network.index_river_f

    vertical.waterlevel_land .= lateral.land.h_av .* 1000.0
    vertical.waterlevel_river[inds_riv] .= lateral.river.h_av .* 1000.0

    # vertical sbm concept is updated until snow state, after that (optional)
    # snow transport is possible
    update_until_snow(vertical, config)

    # lateral snow transport 
    if get(config.model, "masswasting", false)::Bool
        lateral_snow_transport!(
            vertical.snow,
            vertical.snowwater,
            lateral.subsurface.βₗ,
            network.land,
        )
    end

    # update vertical sbm concept until recharge [mm] to the saturated store
    update_until_recharge(vertical, config)

    # exchange of recharge between vertical sbm concept and subsurface flow domain
    lateral.subsurface.recharge .= vertical.recharge ./ 1000.0
    lateral.subsurface.recharge .*= lateral.subsurface.dw
    lateral.subsurface.zi .= vertical.zi ./ 1000.0

    # update lateral subsurface flow domain (kinematic wave)
    update(lateral.subsurface, network.land, network.frac_toriver)

    # update vertical sbm concept (runoff, ustorelayerdepth and satwaterdepth)
    update_after_subsurfaceflow(
        vertical,
        lateral.subsurface.zi .* 1000.0,
        lateral.subsurface.exfiltwater .* 1000.0,
    )

    @. lateral.land.runoff = (
        (vertical.runoff / 1000.0) * (network.land.xl * network.land.yl) / vertical.Δt +
        lateral.subsurface.to_river / lateral.subsurface.Δt +
        # net_runoff_river
        (
            (vertical.net_runoff_river * network.land.xl * network.land.yl * 0.001) /
            vertical.Δt
        )
    )

    update(lateral.land, lateral.river, network)

    # update the clock
    advance!(clock)

    return model

end
