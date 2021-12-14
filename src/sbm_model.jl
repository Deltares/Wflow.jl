"""
    initialize_sbm_model(config::Config)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_sbm_model(config::Config)

    model_type = config.model.type::String
    @info "Initializing of model variables for model type $model_type"

    # unpack the paths to the NetCDF files
    static_path = input_path(config, config.input.path_static)

    reader = prepare_reader(config)
    clock = Clock(config, reader)
    Δt = clock.Δt

    reinit = get(config.model, "reinit", true)::Bool
    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_pits = get(config.model, "pits", false)::Bool

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool
    river_routing = get(config.model, "river_routing", "kinematic-wave")

    snow = get(config.model, "snow", false)::Bool
    reservoirs = do_reservoirs
    lakes = do_lakes
    glacier = get(config.model, "glacier", false)::Bool
    masswasting = get(config.model, "masswasting", false)::Bool
    @info "General model settings" reservoirs lakes snow masswasting glacier

    nc = NCDataset(static_path)

    subcatch_2d =
        ncread(nc, config.input, "subcatchment"; optional = false, allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)
    modelsize_2d = size(subcatch_2d)

    river_2d = ncread(
        nc,
        config.input,
        "river_location";
        optional = false,
        type = Bool,
        fill = false,
    )
    river = river_2d[inds]
    riverwidth_2d = ncread(
        nc,
        config.input,
        "lateral.river.width";
        optional = false,
        type = Float,
        fill = 0,
    )
    riverwidth = riverwidth_2d[inds]
    riverlength_2d = ncread(
        nc,
        config.input,
        "lateral.river.length";
        optional = false,
        type = Float,
        fill = 0,
    )
    riverlength = riverlength_2d[inds]

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

    # reservoirs
    pits = zeros(Bool, modelsize_2d)
    if do_reservoirs
        reservoirs, resindex, reservoir, pits =
            initialize_simple_reservoir(config, nc, inds_riv, nriv, pits, tosecond(Δt))
    else
        reservoir = ()
        reservoirs = nothing
        resindex = fill(0, nriv)
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits =
            initialize_natural_lake(config, nc, inds_riv, nriv, pits, tosecond(Δt))
    else
        lake = ()
        lakes = nothing
        lakeindex = fill(0, nriv)
    end

    ldd_2d = ncread(nc, config.input, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d =
            ncread(nc, config.input, "pits"; optional = false, type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, inds)
    end

    βₗ = ncread(
        nc,
        config.input,
        "lateral.land.slope";
        optional = false,
        sel = inds,
        type = Float,
    )
    clamp!(βₗ, 0.00001, Inf)

    dl = map(detdrainlength, ldd, xl, yl)
    dw = (xl .* yl) ./ dl

    # check if lateral subsurface flow component is defined for the SBM model, when coupled 
    # to another groundwater model, this component is not defined in the TOML file.
    if haskey(config.input.lateral, "subsurface")
        khfrac = ncread(
            nc,
            config.input,
            "lateral.subsurface.ksathorfrac";
            optional = false,
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

    olf = initialize_surfaceflow_land(
        nc,
        config,
        inds;
        sl = βₗ,
        dl = dl,
        width = map(det_surfacewidth, dw, riverwidth, river),
        wb_pit = pits[inds],
        iterate = kinwave_it,
        tstep = kw_land_tstep,
        Δt = Δt,
    )

    graph = flowgraph(ldd, inds, pcr_dir)

    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]

    ldd_riv = ldd_2d[inds_riv]
    if do_pits
        ldd_riv = set_pit_ldd(pits_2d, ldd_riv, inds_riv)
    end
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, βₗ, n)

    if river_routing == "kinematic-wave"
        rf = initialize_surfaceflow_river(
            nc,
            config,
            inds_riv;
            dl = riverlength,
            width = riverwidth,
            wb_pit = pits[inds_riv],
            reservoir_index = resindex,
            reservoir = reservoirs,
            lake_index = lakeindex,
            lake = lakes,
            river = river,
            iterate = kinwave_it,
            tstep = kw_river_tstep,
            Δt = Δt,
        )
    elseif river_routing == "local-inertial"
        rf, nodes_at_link = initialize_shallowwater_river(
            nc,
            config,
            inds_riv;
            graph = graph_riv,
            ldd = ldd_riv,
            dl = riverlength,
            width = riverwidth,
            reservoir_index = resindex,
            reservoir = reservoirs,
            lake_index = lakeindex,
            lake = lakes,
            Δt = Δt,
        )
    else
        error(
            """An unknown "river_routing" method is specified in the TOML file ($river_routing). 
            This should be "kinematic_wave" or "local-inertial".
            """,
        )
    end

    # setup subdomains for the land and river kinematic wave domain, if nthreads = 1
    # subdomain is equal to the complete domain
    toposort = topological_sort_by_dfs(graph)
    index_pit_land = findall(x -> x == 5, ldd)
    subbas_order, indices_subbas, topo_subbas =
        kinwave_set_subdomains(config, graph, toposort, index_pit_land)
    if river_routing == "kinematic-wave"
        toposort_riv = topological_sort_by_dfs(graph_riv)
        index_pit_river = findall(x -> x == 5, ldd_riv)
        subriv_order, indices_subriv, topo_subriv =
            kinwave_set_subdomains(config, graph_riv, toposort_riv, index_pit_river)
    end

    if nthreads() > 1
        @info "Parallel execution of kinematic wave, minimum stream order = $(get(config.model, "min_streamorder", 4))"
    end

    modelmap = (vertical = sbm, lateral = (subsurface = ssf, land = olf, river = rf))
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
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
        upstream_nodes = filter_upsteam_nodes(graph, olf.wb_pit),
        subdomain_order = subbas_order,
        topo_subdomain = topo_subbas,
        indices_subdomain = indices_subbas,
        order = toposort,
        indices = inds,
        reverse_indices = rev_inds,
        xl = xl,
        yl = yl,
    )
    if river_routing == "kinematic-wave"
        river = (
            graph = graph_riv,
            upstream_nodes = filter_upsteam_nodes(graph_riv, rf.wb_pit),
            subdomain_order = subriv_order,
            topo_subdomain = topo_subriv,
            indices_subdomain = indices_subriv,
            order = toposort_riv,
            indices = inds_riv,
            reverse_indices = rev_inds_riv,
        )
    elseif river_routing == "local-inertial"
        river = (
            graph = graph_riv,
            nodes_at_link = nodes_at_link,
            links_at_node = adjacent_links_at_node(graph, nodes_at_link),
            indices = inds_riv,
            reverse_indices = rev_inds_riv,
        )
    end

    model = Model(
        config,
        (; land, river, reservoir, lake, index_river, frac_toriver),
        (subsurface = ssf, land = olf, river = rf),
        sbm,
        clock,
        reader,
        writer,
        SbmModel(),
    )

    # read and set states in model object if reinit=false
    if reinit == false
        instate_path = input_path(config, config.state.path_input)
        @info "Read and set initial conditions from state file $instate_path"
        state_ncnames = ncnames(config.state)
        set_states(instate_path, model, state_ncnames; type = Float)
        @unpack lateral, vertical = model
        # update zi for vertical sbm and kinematic wave volume for river and land domain
        zi =
            max.(
                0.0,
                vertical.soilthickness .-
                vertical.satwaterdepth ./ (vertical.θₛ .- vertical.θᵣ),
            )
        vertical.zi .= zi
        # makes sure land cells with zero flow width are set to zero q and h
        for i in eachindex(lateral.land.width)
            if lateral.land.width[i] <= 0.0
                lateral.land.q[i] = 0.0
                lateral.land.h[i] = 0.0
            end
        end
        lateral.land.volume .= lateral.land.h .* lateral.land.width .* lateral.land.dl
        # only set active cells for river (ignore boundary conditions/ghost points)
        lateral.river.volume[1:nriv] .=
            lateral.river.h[1:nriv] .* lateral.river.width[1:nriv] .*
            lateral.river.dl[1:nriv]

        if do_lakes
            # storage must be re-initialized after loading the state with the current
            # waterlevel otherwise the storage will be based on the initial water level
            lakes.storage .=
                initialize_storage(lakes.storfunc, lakes.area, lakes.waterlevel, lakes.sh)
        end
    else
        @info "Setting initial conditions to default"
    end

    @info "End of model initialization"

    return model
end

"update SBM model for a single timestep"
function update(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}

    @unpack lateral, vertical, network, clock, config = model
    model = update_until_recharge(model)
    # exchange of recharge between vertical sbm concept and subsurface flow domain
    lateral.subsurface.recharge .= vertical.recharge ./ 1000.0
    lateral.subsurface.recharge .*= lateral.subsurface.dw
    lateral.subsurface.zi .= vertical.zi ./ 1000.0
    # update lateral subsurface flow domain (kinematic wave)
    update(lateral.subsurface, network.land, network.frac_toriver)
    model = update_after_subsurfaceflow(model)
end

"""
    update_until_recharge(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}

Update SBM model until recharge for a single timestep. This function is also accessible
through BMI, to couple the SBM model to an external groundwater model.
"""
function update_until_recharge(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}
    @unpack lateral, vertical, network, clock, config = model

    inds_riv = network.index_river

    # extract water levels h_av [m] from the land and river domains
    # this is used to limit open water evaporation
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
            lateral.land.sl,
            network.land,
        )
    end

    # update vertical sbm concept until recharge [mm] to the saturated store
    update_until_recharge(vertical, config)

    return model
end

"""
    update_after_subsurfaceflow(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}

Update SBM model after subsurface flow for a single timestep. This function is also 
accessible through BMI, to couple the SBM model to an external groundwater model.
"""
function update_after_subsurfaceflow(
    model::Model{N,L,V,R,W,T},
) where {N,L,V,R,W,T<:SbmModel}
    @unpack lateral, vertical, network, clock, config = model

    inds_riv = network.index_river

    # update vertical sbm concept (runoff, ustorelayerdepth and satwaterdepth)
    update_after_subsurfaceflow(
        vertical,
        lateral.subsurface.zi * 1000.0,
        lateral.subsurface.exfiltwater * 1000.0,
    )

    ssf_toriver = lateral.subsurface.to_river[inds_riv] ./ lateral.subsurface.Δt
    surface_routing(model, ssf_toriver = ssf_toriver)

    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
