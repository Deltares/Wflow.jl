"""
    initialize_sbm_model(config::Config)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_sbm_model(config::Config)

    # unpack the paths to the NetCDF files
    tomldir = dirname(config)
    static_path = joinpath(tomldir, config.input.path_static)
    dynamic_path = joinpath(tomldir, config.input.path_forcing)

    reader = prepare_reader(dynamic_path, static_path, config)
    clock = Clock(config, reader)
    Δt = clock.Δt

    reinit = get(config.model, "reinit", true)::Bool
    do_snow = get(config.model, "snow", false)::Bool
    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_pits = get(config.model, "pits", false)::Bool

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)

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
            initialize_simple_reservoir(config, nc, inds_riv, nriv, pits)
    else
        reservoir = ()
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits =
            initialize_natural_lake(config, dirname(static_path), nc, inds_riv, nriv, pits)
    else
        lake = ()
    end

    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d = ncread(nc, param(config, "input.pits"); type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, inds)
    end

    βₗ = ncread(nc, param(config, "input.lateral.land.slope"); sel = inds, type = Float)
    clamp!(βₗ, 0.00001, Inf)
    dl = fill(mv, n)
    dw = fill(mv, n)
    sw = fill(mv, n)

    for i = 1:n
        dl[i] = detdrainlength(ldd[i], xl[i], yl[i])
        dw[i] = (xl[i] * yl[i]) / dl[i]
        sw[i] = river[i] ? max(dw[i] - riverwidth[i], 0.0) : dw[i]
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

        kh₀ = khfrac .* sbm.kv₀

        ssf = LateralSSF{Float}(
            kh₀ = kh₀,
            f = sbm.f,
            zi = sbm.zi,
            soilthickness = sbm.soilthickness,
            θₛ = sbm.θₛ,
            θᵣ = sbm.θᵣ,
            Δt = tosecond(Δt),
            t = 1.0,
            βₗ = βₗ,
            dl = dl .* 1000.0,
            dw = dw .* 1000.0,
            exfiltwater = fill(mv, n),
            recharge = fill(mv, n),
            ssf = ((kh₀ .* βₗ) ./ sbm.f) .*
                  (exp.(-sbm.f .* sbm.zi) - exp.(-sbm.f .* sbm.soilthickness)) .* dw .*
                  1000.0,
            ssfin = fill(mv, n),
            ssfmax = ((kh₀ .* βₗ) ./ sbm.f) .* (1.0 .- exp.(-sbm.f .* sbm.soilthickness)),
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

    alpha_pow = Float((2.0 / 3.0) * 0.6)
    β = Float(0.6)
    olf = SurfaceFlow(
        β = β,
        sl = βₗ,
        n = n_land,
        dl = dl,
        q = zeros(Float, n),
        qin = zeros(Float, n),
        q_av = zeros(Float, n),
        qlat = zeros(Float, n),
        inwater = zeros(Float, n),
        volume = zeros(Float, n),
        h = zeros(Float, n),
        h_av = zeros(Float, n),
        Δt = Float(tosecond(Δt)),
        its = kw_land_tstep > 0 ? Int(cld(tosecond(Δt), kw_land_tstep)) : kw_land_tstep,
        width = sw,
        wb_pit = pits[inds],
        alpha_pow = alpha_pow,
        α = fill(mv, n),
        eps = Float(1e-3),
        cel = zeros(Float, n),
        to_river = zeros(Float, n),
        rivercells = fill(false, n),
        reservoir_index = fill(0, n),
        lake_index = fill(0, n),
        reservoir = nothing,
        lake = nothing,
    )

    graph = flowgraph(ldd, inds, pcr_dir)


    riverslope =
        ncread(nc, param(config, "input.lateral.river.slope"); sel = inds_riv, type = Float)
    clamp!(riverslope, 0.00001, Inf)
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
    if do_pits
        ldd_riv = set_pit_ldd(pits_2d, ldd_riv, inds_riv)
    end
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, βₗ, n)

    rf = SurfaceFlow(
        β = β,
        sl = riverslope,
        n = n_river,
        dl = riverlength,
        q = zeros(Float, nriv),
        qin = zeros(Float, nriv),
        q_av = zeros(Float, nriv),
        qlat = zeros(Float, nriv),
        inwater = zeros(Float, nriv),
        volume = zeros(Float, nriv),
        h = zeros(Float, nriv),
        h_av = zeros(Float, nriv),
        Δt = Float(tosecond(Δt)),
        its = kw_river_tstep > 0 ? ceil(Int(tosecond(Δt) / kw_river_tstep)) :
              kw_river_tstep,
        width = riverwidth,
        wb_pit = pits[inds_riv],
        alpha_pow = alpha_pow,
        α = fill(mv, nriv),
        eps = Float(1e-3),
        cel = zeros(Float, nriv),
        to_river = zeros(Float, nriv),
        reservoir_index = do_reservoirs ? resindex : fill(0, nriv),
        lake_index = do_lakes ? lakeindex : fill(0, nriv),
        reservoir = do_reservoirs ? reservoirs : nothing,
        lake = do_lakes ? lakes : nothing,
        rivercells = river,
    )

    # setup subdomains for the land and river kinematic wave domain, if nthreads = 1
    # subdomain is equal to the complete domain
    toposort = topological_sort_by_dfs(graph)
    toposort_riv = topological_sort_by_dfs(graph_riv)
    index_pit_land = findall(x -> x == 5, ldd)
    index_pit_river = findall(x -> x == 5, ldd_riv)
    subbas_order, indices_subbas, topo_subbas =
        kinwave_set_subdomains(config, graph, toposort, index_pit_land)
    subriv_order, indices_subriv, topo_subriv =
        kinwave_set_subdomains(config, graph_riv, toposort_riv, index_pit_river)

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

    model = Model(
        config,
        (; land, river, reservoir, lake, index_river, frac_toriver),
        (subsurface = ssf, land = olf, river = rf),
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
        # update zi for vertical sbm and kinematic wave volume for river and land domain
        zi =
            max.(
                0.0,
                vertical.soilthickness .-
                vertical.satwaterdepth ./ (vertical.θₛ .- vertical.θᵣ),
            )
        vertical.zi .= zi
        lateral.land.volume .= lateral.land.h .* lateral.land.width .* lateral.land.dl
        lateral.river.volume .= lateral.river.h .* lateral.river.width .* lateral.river.dl
    end

    # make sure the forcing is already loaded
    # it's fine to run twice, and may help catching errors earlier
    update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end
    return model
end

"update SBM model for a single timestep"
function update(model::Model{N,L,V,R,W}) where {N,L,V<:SBM,R,W}
    @unpack lateral, vertical, network, clock, config = model
    model = update_until_recharge(model)
    # exchange of recharge between vertical sbm concept and subsurface flow domain
    lateral.subsurface.recharge .= vertical.recharge
    lateral.subsurface.recharge .*= lateral.subsurface.dw
    lateral.subsurface.zi .= vertical.zi
    # update lateral subsurface flow domain (kinematic wave)
    update(lateral.subsurface, network.land, network.frac_toriver)
    model = update_after_subsurfaceflow(model)
    return model
end

"""
    update_until_recharge(model::Model{N,L,V,R,W}) where {N,L,V<:SBM,R,W}

Update SBM model until recharge for a single timestep. This function is also accessible
through BMI, to couple the SBM model to an external groundwater model.
"""
function update_until_recharge(model::Model{N,L,V,R,W}) where {N,L,V<:SBM,R,W}
    @unpack lateral, vertical, network, clock, config = model

    inds_riv = network.index_river

    update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end

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
    update_after_subsurfaceflow(model::Model{N,L,V,R,W}) where {N,L,V<:SBM,R,W}

Update SBM model after subsurface flow for a single timestep. This function is also 
accessible through BMI, to couple the SBM model to an external groundwater model.
"""
function update_after_subsurfaceflow(model::Model{N,L,V,R,W}) where {N,L,V<:SBM,R,W}
    @unpack lateral, vertical, network, clock, config = model

    inds_riv = network.index_river
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool

    # update vertical sbm concept (runoff, ustorelayerdepth and satwaterdepth)
    update_after_subsurfaceflow(
        vertical,
        lateral.subsurface.zi,
        lateral.subsurface.exfiltwater,
    )

    # determine lateral inflow for overland flow based on vertical runoff [mm] from vertical
    # sbm concept
    @. lateral.land.inwater =
        (vertical.runoff * network.land.xl * network.land.yl * 0.001) / lateral.land.Δt
    lateral.land.qlat .= lateral.land.inwater ./ lateral.land.dl

    # run kinematic wave for overland flow
    update(
        lateral.land,
        network.land,
        frac_toriver = network.frac_toriver,
        do_iter = kinwave_it,
    )

    # determine net runoff from vertical sbm concept in river cells, and lateral inflow from
    # overland flow lateral subsurface flow and net runoff to the river cells
    @. lateral.river.inwater = (
        lateral.subsurface.to_river[inds_riv] / 1.0e9 / lateral.river.Δt +
        lateral.land.to_river[inds_riv] +
        # net_runoff_river
        (
            (
                vertical.net_runoff_river[inds_riv] *
                network.land.xl[inds_riv] *
                network.land.yl[inds_riv] *
                0.001
            ) / vertical.Δt
        )
    )
    lateral.river.qlat .= lateral.river.inwater ./ lateral.river.dl

    # run kinematic wave for river flow
    # check if reservoirs or lakes are defined, the inflow from lateral subsurface and
    # overland flow is required
    if !isnothing(lateral.river.reservoir) || !isnothing(lateral.river.lake)
        inflow_wb =
            lateral.subsurface.ssf[inds_riv] ./ 1.0e9 ./ lateral.river.Δt .+
            lateral.land.q_av[inds_riv]
        update(
            lateral.river,
            network.river,
            inflow_wb = inflow_wb,
            do_iter = kinwave_it,
            doy = dayofyear(clock.time),
        )
    else
        update(
            lateral.river,
            network.river,
            do_iter = kinwave_it,
            doy = dayofyear(clock.time),
        )
    end
    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
