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

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    reinit = get(config.model, "reinit", true)::Bool
    do_snow = get(config.model, "snow", false)::Bool
    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_pits = get(config.model, "pits", false)::Bool

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)

    nc = NCDataset(static_path)
    dims = dimnames(nc[param(config, "input.subcatchment")])

    # There is no need to permute the dimensions of the data, since the active indices are
    # correctly calculated in both ways.
    # The dimension order only needs to be known for interpreting the LDD directions
    # and creating the coordinate maps.
    dims_xy = dims[2] in ("y", "lat")

    subcatch_2d = ncread(nc, param(config, "input.subcatchment"); allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)
    modelsize_2d = size(subcatch_2d)

    river_2d = ncread(nc, param(config, "input.river_location"); type = Bool, fill = false)
    river = river_2d[inds]
    riverwidth_2d =
        ncread(nc, param(config, "input.lateral.river.width"); type = Float64, fill = 0)
    riverwidth = riverwidth_2d[inds]
    riverlength_2d =
        ncread(nc, param(config, "input.lateral.river.length"); type = Float64, fill = 0)
    riverlength = riverlength_2d[inds]

    # read x, y coordinates and calculate cell length [m]
    y_nc = "y" in keys(nc.dim) ? ncread(nc, "y") : ncread(nc, "lat")
    x_nc = "x" in keys(nc.dim) ? ncread(nc, "x") : ncread(nc, "lon")
    if dims_xy
        y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    else
        y = repeat(y_nc, outer = (1, length(x_nc)))[inds]
    end
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

    sbm = initialize_sbm(nc, config, riverfrac, xl, yl, inds)

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
            initialize_natural_lake(config, static_path, nc, inds_riv, nriv, pits)
    else
        lake = ()
    end

    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d = ncread(nc, param(config, "input.pits"); type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, inds)
    end


    # lateral part sbm
    khfrac = ncread(
        nc,
        param(config, "input.lateral.subsurface.ksathorfrac", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float64,
    )
    βₗ = ncread(nc, param(config, "input.lateral.land.slope"); sel = inds, type = Float64)
    clamp!(βₗ, 0.00001, Inf)
    kh₀ = khfrac .* sbm.kv₀
    dl = fill(mv, n)
    dw = fill(mv, n)
    sw = fill(mv, n)

    for i = 1:n
        dl[i] = detdrainlength(ldd[i], xl[i], yl[i])
        dw[i] = (xl[i] * yl[i]) / dl[i]
        sw[i] = river[i] ? max(dw[i] - riverwidth[i], 0.0) : dw[i]
    end

    ssf = LateralSSF{Float64}(
        kh₀ = kh₀,
        f = sbm.f,
        zi = sbm.zi,
        soilthickness = sbm.soilthickness,
        θₑ = sbm.θₛ .- sbm.θᵣ,
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

    n_land = ncread(
        nc,
        param(config, "input.lateral.land.n", nothing);
        sel = inds,
        defaults = 0.072,
        type = Float64,
    )

    alpha_pow = (2.0 / 3.0) * 0.6
    β = 0.6
    h = fill(0.0, n)
    alpha_term = pow.(n_land ./ sqrt.(βₗ), β)
    olf = SurfaceFlow(
        β = β,
        sl = βₗ,
        n = n_land,
        dl = dl,
        q = fill(0.0, n),
        qin = fill(0.0, n),
        q_av = fill(0.0, n),
        qlat = fill(0.0, n),
        h = h,
        h_av = fill(0.0, n),
        Δt = tosecond(Δt),
        its = kw_land_tstep > 0 ? Int(cld(tosecond(Δt), kw_land_tstep)) : kw_land_tstep,
        width = sw,
        wb_pit = pits[inds],
        alpha_term = alpha_term,
        alpha_pow = alpha_pow,
        α = alpha_term .* pow.(sw .+ 2.0 .* h, alpha_pow),
        eps = 1e-03,
        cel = fill(0.0, n),
        to_river = fill(0.0, n),
        rivercells = fill(false, n),
        reservoir_index = fill(0, n),
        lake_index = fill(0, n),
        reservoir = nothing,
        lake = nothing,
    )

    pcr_dir = dims_xy ? permute_indices(pcrdir) : pcrdir
    graph = flowgraph(ldd, inds, pcr_dir)


    riverslope = ncread(
        nc,
        param(config, "input.lateral.river.slope");
        sel = inds_riv,
        type = Float64,
    )
    clamp!(riverslope, 0.00001, Inf)
    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]
    n_river = ncread(
        nc,
        param(config, "input.lateral.river.n", nothing);
        sel = inds_riv,
        defaults = 0.036,
        type = Float64,
    )
    ldd_riv = ldd_2d[inds_riv]
    if do_pits
        ldd_riv = set_pit_ldd(pits_2d, ldd_riv, inds_riv)
    end
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, βₗ, n)

    h_river = fill(0.0, nriv)
    alpha_term = pow.(n_river ./ sqrt.(riverslope), β)
    rf = SurfaceFlow(
        β = β,
        sl = riverslope,
        n = n_river,
        dl = riverlength,
        q = fill(0.0, nriv),
        qin = fill(0.0, nriv),
        q_av = fill(0.0, nriv),
        qlat = fill(0.0, nriv),
        h = h_river,
        h_av = fill(0.0, nriv),
        Δt = tosecond(Δt),
        its = kw_river_tstep > 0 ? ceil(Int(tosecond(Δt) / kw_river_tstep)) :
              kw_river_tstep,
        width = riverwidth,
        wb_pit = pits[inds_riv],
        alpha_term = alpha_term,
        alpha_pow = alpha_pow,
        α = alpha_term .* pow.(riverwidth .+ 2.0 .* h_river, alpha_pow),
        eps = 1e-03,
        cel = fill(0.0, nriv),
        to_river = fill(0.0, nriv),
        reservoir_index = do_reservoirs ? resindex : fill(0, nriv),
        lake_index = do_lakes ? lakeindex : fill(0, nriv),
        reservoir = do_reservoirs ? reservoirs : nothing,
        lake = do_lakes ? lakes : nothing,
        rivercells = river,
    )

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
        dims_xy,
        nc,
        maxlayers = sbm.maxlayers,
    )
    close(nc)

    # for each domain save the directed acyclic graph, the traversion order,
    # and the indices that map it back to the two dimensional grid
    # for reservoirs and lakes this is information is also available (except the graph) 
    # from the initialization functions
    land = (
        graph = graph,
        upstream_nodes = filter_upsteam_nodes(graph, olf.wb_pit),
        order = topological_sort_by_dfs(graph),
        indices = inds,
        reverse_indices = rev_inds,
    )
    river = (
        graph = graph_riv,
        upstream_nodes = filter_upsteam_nodes(graph_riv, rf.wb_pit),
        order = topological_sort_by_dfs(graph_riv),
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
        set_states(instate_path, model, state_ncnames; type = Float64)
        @unpack lateral, vertical = model
        # update zi for vertical sbm and α for river and overland flow
        zi = max.(0.0, vertical.soilthickness .- vertical.satwaterdepth ./ vertical.θₑ)
        vertical.zi .= zi
        lateral.river.α .=
            lateral.river.alpha_term .*
            pow.(lateral.river.width .+ 2.0 .* lateral.river.h, lateral.river.alpha_pow)
        lateral.land.α .=
            lateral.land.alpha_term .*
            pow.(lateral.land.width .+ 2.0 .* lateral.land.h, lateral.land.alpha_pow)
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

    inds_riv = network.index_river
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool

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

    # exchange of recharge between vertical sbm concept and subsurface flow domain
    lateral.subsurface.recharge .= vertical.recharge
    lateral.subsurface.recharge .*= lateral.subsurface.dw
    lateral.subsurface.zi .= vertical.zi

    # update lateral subsurface flow domain (kinematic wave)
    update(lateral.subsurface, network.land, network.frac_toriver)

    # update vertical sbm concept (runoff, ustorelayerdepth and satwaterdepth)
    update_after_lateralflow(
        vertical,
        lateral.subsurface.zi,
        lateral.subsurface.exfiltwater,
    )

    # determine lateral inflow for overland flow based on vertical runoff [mm] from vertical
    # sbm concept
    lateral.land.qlat .=
        (vertical.runoff .* vertical.xl .* vertical.yl .* 0.001) ./ lateral.land.Δt ./
        lateral.land.dl
    # run kinematic wave for overland flow
    update(
        lateral.land,
        network.land,
        frac_toriver = network.frac_toriver,
        do_iter = kinwave_it,
    )

    # determine net runoff from vertical sbm concept in river cells, and lateral inflow from
    # overland flow lateral subsurface flow and net runoff to the river cells
    net_runoff_river =
        (
            vertical.net_runoff_river[inds_riv] .* vertical.xl[inds_riv] .*
            vertical.yl[inds_riv] .* 0.001
        ) ./ vertical.Δt
    lateral.river.qlat .=
        (
            lateral.subsurface.to_river[inds_riv] ./ 1.0e9 ./ lateral.river.Δt .+
            lateral.land.to_river[inds_riv] .+ net_runoff_river
        ) ./ lateral.river.dl

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
