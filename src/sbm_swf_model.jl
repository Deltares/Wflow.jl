"""
    initialize_sbm_swf_model(config::Config)

Initial part of the SBM + Shallow water flow model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
@with_kw struct Indices
    xu :: Vector{Int}
    xd :: Vector{Int}
    yu :: Vector{Int}
    yd :: Vector{Int}
end

const dirs = ["yd", "xd", "xu", "yu"]

function initialize_sbm_swf_model(config::Config)

    # unpack the paths to the NetCDF files
    tomldir = dirname(config)
    static_path = joinpath(tomldir, config.input.path_static)
    dynamic_path = joinpath(tomldir, config.input.path_forcing)

    reader = prepare_reader(dynamic_path, static_path, config)
    clock = Clock(config, reader)
    Δt = clock.Δt

    reinit = get(config.model, "reinit", true)::Bool

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

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fill(Float(1), n)

    # the Almeida et al. (2012) algorithm makes use of a staggered grid (nodes and links).
    # the following is required:
    # for each node: - the source node (src_node) and destination node (dst_node)
    #                - the source link (src_link_node) and destination link (dst_link_node)
    # for each link: - the source link (src_link) and destination link (dst_link)
    links = collect(edges(graph_riv))
    src_node = src.(links)
    dst_node = dst.(links)
    
    src_link = Vector{Int}[]
    dst_link = copy(src_link)
    
    n_nodes = nv(graph_riv)
    n_links = ne(graph_riv)
    
    for i = 1:n_links
        push!(src_link, findall(isequal(src_node[i]), dst_node))
        push!(dst_link, findall(isequal(dst_node[i]), src_node))
    end
    
    nodes = vertices(graph_riv)
    src_link_node = Vector{Int}[]
    dst_link_node = copy(src_link_node)
    
    for i = 1:n_nodes
        push!(src_link_node, findall(isequal(nodes[i]), dst_node))
        push!(dst_link_node, findall(isequal(nodes[i]), src_node))
    end

    river_elevation = altitude_2d[inds_riv]
    # initialize subgrid river bed elevation at links
    zr_max = fill(Float(0), n_links)
    for i in 1:n_links
        zr_max[i] = max(river_elevation[src_node[i]], river_elevation[dst_node[i]])
    end

    # land part
    indices = Indices(
        xu = fill(0, n),
        xd = fill(0, n),
        yu = fill(0, n),
        yd = fill(0, n),
    )

    nrow, ncol  = modelsize_2d
    for (v,i) in enumerate(inds)
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
        yu = min(indices.yu[i],n)
        xu = min(indices.xu[i],n)
        zx_max[i] = max(elevation[i], elevation[xu])
        zy_max[i] = max(elevation[i], elevation[yu])
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
        if dir == CartesianIndex(1,1)
            we_x[idx] = we_x[idx] - 0.5 * w
            we_y[idx] = we_y[idx] - 0.5 * w
        elseif dir == CartesianIndex(-1,-1)
            we_x[indices.xd[idx]] = we_x[indices.xd[idx]] - 0.5 * w
            we_y[indices.yd[idx]] = we_y[indices.yd[idx]] - 0.5 * w
        elseif dir == CartesianIndex(1,0)
            we_x[idx] = we_x[idx] - w
        elseif dir == CartesianIndex(0,1)
            we_y[idx] = we_y[idx] - w
        elseif dir ==  CartesianIndex(-1,0)
            we_x[indices.xd[idx]] = we_x[indices.xd[idx]] - w
        elseif dir ==  CartesianIndex(0,-1)
            we_y[indices.yd[idx]] = we_x[indices.yd[idx]] - w
        elseif dir ==  CartesianIndex(1,-1)
            we_x[idx] = we_x[idx] - 0.5 * w
            we_y[indices.yd[idx]] = we_x[indices.yd[idx]] - 0.5 * w
        elseif dir ==  CartesianIndex(-1,1)
            we_x[indices.xd[idx]] = we_x[indices.xd[idx]] - 0.5 * w
            we_y[idx] = we_y[idx] - 0.5 * w
        end
    end

    bankheight = fill(Float(2), nriv)    # dummy value for bankheight
    bankvolume =  bankheight .* riverwidth .* riverlength
    volume = fill(Float(0), n)

    h_init = river .* 0.01  # cold state

    for i = 1:n
        if river[i]
            j = rev_inds_riv[inds][i]
            if h_init[i] <= bankheight[j]
                volume[i] = h_init[i] * riverwidth[j] * riverlength[j]
            else
                volume[i] = bankvolume[j] + (h_init[i] - bankheight[j]) * xl[i] * yl[i]
            end
        else
            volume[i] = h_init[i] * xl[i] * yl[i]
        end
    end
    
    sw_land = ShallowWaterLand{Float}(
        n = n,
        dx = xl,
        dy = yl,
        xwidth = we_x,
        ywidth = we_y,
        g = 9.80665,
        θ = 0.8,
        α = 0.2,
        qx0 = zeros(n+1),
        qy0 = zeros(n+1),
        qx = zeros(n+1),
        qy = zeros(n+1),
        zx_max = zx_max,
        zy_max = zy_max,
        mannings_n = n_land,
        volume = volume,
        runoff = zeros(n),
        h = h_init,
        slp_x = zeros(n),
        slp_y = zeros(n),
        z = elevation,
        froude = true,
        rivercells = river,
    )

    sw_river = ShallowWaterRiver{Float}(
        n = nriv,
        g = 9.80665,
        θ = 0.8,
        α = 0.2,
        qr0 = zeros(n_links),
        qr = zeros(n_links),
        zr_max = zr_max,
        mannings_n = n_river,
        h = h_init[index_river],
        riverwidth = riverwidth,
        riverlength = riverlength,
        bankvolume = bankvolume,
        bankheight = bankheight,
        slp = zeros(n_links),
        z_r = river_elevation,
        z_b = river_elevation .+ bankheight,
        froude = true,
    )

    # setup subdomains for the land kinematic wave domain, if nthreads = 1
    # subdomain is equal to the complete domain
    toposort = topological_sort_by_dfs(graph)
    index_pit_land = findall(x -> x == 5, ldd)
    subbas_order, indices_subbas, topo_subbas =
        kinwave_set_subdomains(config, graph, toposort, index_pit_land)

    modelmap = (vertical = sbm, lateral = (subsurface = ssf, land = sw_land, river = sw_river))
    
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
    )

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
    land = (
        graph = graph,
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
        links = links,
        src_node = src_node,
        dst_node = dst_node,
        src_link = src_link,
        dst_link = dst_link,
        src_link_node = src_link_node,
        dst_link_node = dst_link_node,
    )

    model = Model(
        config,
        (; land, river, index_river, frac_toriver),
        (subsurface = ssf, land = sw_land, river=sw_river),
        sbm,
        clock,
        reader,
        writer,
    )

    # read and set states in model object if reinit=false
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

    # make sure the forcing is already loaded
    # it's fine to run twice, and may help catching errors earlier
    update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end
    return model
end