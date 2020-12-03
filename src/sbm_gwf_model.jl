"""
    initialize_sbm_gwf_model(config::Config)

Initial part of the sbm_gwf model concept. The model contains:
    - the vertical SBM concept
    - the following lateral components:
        - 1-D kinematic wave for river flow
        - 1-D kinematic wave for overland flow
        - unconfined aquifer with groundwater flow in four directions (adjacent cells)

The unconfined aquifer contains a recharge, river and a drain (optional) boundary. 

The initial part reads the input settings and data as defined in the Config object. 
Will return a Model that is ready to run.
"""
function initialize_sbm_gwf_model(config::Config)

    # unpack the paths to the NetCDF files
    tomldir = dirname(config)
    static_path = joinpath(tomldir, config.input.path_static)
    cyclic_path = joinpath(tomldir, config.input.path_static)
    dynamic_path = joinpath(tomldir, config.input.path_forcing)
    instate_path = joinpath(tomldir, config.state.path_input)
    output_path = joinpath(tomldir, config.output.path)

    Δt = Second(config.timestepsecs)
    sizeinmetres = get(config.model, "sizeinmetres", false)
    reinit = get(config.model, "reinit", true)
    do_snow = get(config.model, "snow", false)
    do_reservoirs = get(config.model, "reservoirs", false)
    do_lakes = get(config.model, "lakes", false)
    do_drains = get(config.model, "drains", false)
    do_constanthead = get(config.model, "constanthead", false)

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

    # initialize vertical SBM concept
    sbm = initialize_sbm(nc, config, riverfrac, xl, yl, inds)

    inds_riv, rev_inds_riv = active_indices(river_2d, 0)
    nriv = length(inds_riv)

    # reservoirs
    pits = zeros(Bool, modelsize_2d)
    if do_reservoirs
        reservoirs, resindex, reservoir, pits = initialize_simple_reservoir(config, nc, inds_riv, nriv, pits)
    else
        reservoir = ()
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits = initialize_natural_lake(config, nc, inds_riv, nriv, pits)
    else
        lake = ()
    end

    # overland flow (kinematic wave)
    βₗ = ncread(nc, param(config, "input.lateral.land.slope"); sel = inds, type = Float64)
    clamp!(βₗ, 0.00001, Inf)
    ldd_2d = ncread(nc, param(config, "input.ldd"); allow_missing = true)
    ldd = ldd_2d[inds]
    dl = fill(mv, n)
    dw = fill(mv, n)
    sw = fill(mv, n)

    for i = 1:n
        dl[i] = detdrainlength(ldd[i], xl[i], yl[i])
        dw[i] = (xl[i] * yl[i])/dl[i]
        sw[i] = river[i] ? max(dw[i] - riverwidth[i], 0.0) : dw[i]
    end

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

    # river flow (kinematic wave)
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
        q_av = fill(0.0, nriv),
        qlat = fill(0.0, nriv),
        h = h_river,
        h_av = fill(0.0, nriv),
        Δt = tosecond(Δt),
        its = kw_river_tstep > 0 ? ceil(Int(tosecond(Δt) / kw_river_tstep)) : kw_river_tstep,
        width = riverwidth,
        wb_pit = fill(false, nriv),
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

    # unconfined aquifer
    if do_constanthead
        constanthead = ncread(
            nc,
            param(config, "input.lateral.subsurface.constant_head", nothing);
            sel = inds,
            type = Float64,
            fill = mv,
        )
        index_constanthead = filter(i -> !isequal(constanthead[i], mv), 1:n)
        constant_head = ConstantHead(constanthead[index_constanthead], 
                                    index_constanthead,)

    else
        constant_head = ConstantHead[]
    end

    conductivity = ncread(
        nc,
        param(config, "input.lateral.subsurface.conductivity", nothing);
        sel = inds,
        type = Float64,
    )
    specific_yield = ncread(
        nc,
        param(config, "input.lateral.subsurface.specific_yield", nothing);
        sel = inds,
        type = Float64,
    )

    connectivity = Connectivity(inds, rev_inds, sbm.xl, sbm.yl)
    initial_head = sbm.altitude .- 0.10 # cold state for groundwater head
    initial_head[index_river] = sbm.altitude[index_river]

    if do_constanthead
        initial_head[constant_head.index] = constant_head.head
    end

    aquifer = UnconfinedAquifer(
        initial_head,
        conductivity,
        sbm.altitude,
        sbm.altitude .- sbm.soilthickness ./ 1000.0,
        xl .* yl,
        specific_yield,
        fill(0.0, connectivity.nconnection),  # conductance
    )

    # reset zi and satwaterdepth with groundwater head from unconfined aquifer 
    sbm.zi .= (sbm.altitude .- initial_head) .* 1000.0
    sbm.satwaterdepth .= (sbm.soilthickness .- sbm.zi) .* sbm.θₑ

    # river boundary of unconfined aquifer
    infiltration_conductance = ncread(
        nc,
        param(config, "input.lateral.subsurface.infiltration_conductance", nothing);
        sel = inds_riv,
        type = Float64,
    )
    exfiltration_conductance = ncread(
        nc,
        param(config, "input.lateral.subsurface.exfiltration_conductance", nothing);
        sel = inds_riv,
        type = Float64,
    )
    river_bottom = ncread(
        nc,
        param(config, "input.lateral.subsurface.river_bottom", nothing);
        sel = inds_riv,
        type = Float64,
    )

    river_flux = fill(mv,nriv)
    river_stage = fill(mv,nriv)
    river = River(
        river_stage,
        infiltration_conductance,
        exfiltration_conductance,
        river_bottom,
        river_flux,    
        index_river,
    )

    # drain boundary of unconfined aquifer (optional)
    if do_drains
        drain_2d = ncread(nc, param(config, "input.lateral.subsurface.drain"); type = Bool, fill=false)
        inds_drain, rev_inds_drain = active_indices(drain_2d, 0)
        drain = drain_2d[inds]
        drain_elevation = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_elevation", nothing);
            sel = inds,
            type = Float64,
            fill = mv,
        )
        drain_conductance = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_conductance", nothing);
            sel = inds,
            type = Float64,
            fill = mv,
        )
        index_drain = filter(i -> !isequal(drain[i], 0), 1:n)
        drain_flux = fill(mv,length(index_drain))
        drains = Drainage(
            drain_elevation[index_drain],
            drain_conductance[index_drain],
            drain_flux,
            index_drain,
        )
        drain = (indices = inds_drain, reverse_indices = rev_inds_drain)
    else
        drains = Drainage[]
        drain = ()
    end

    # recharge boundary of unconfined aquifer
    r = fill(mv,n)
    recharge = Recharge(r, fill(0.0, n), collect(1:n))

    gwf = GroundwaterFlow(
        aquifer,
        connectivity,
        constant_head,
        AquiferBoundaryCondition[recharge, river, drains],
    )

    state_ncnames = ncnames(config.state)

    reader = prepare_reader(dynamic_path, cyclic_path, config)

    modelmap = (vertical = sbm,
                lateral = (subsurface = (flow = gwf, recharge = gwf.boundaries[1], river = gwf.boundaries[2], drain = gwf.boundaries[3]),
                land = olf,
                river = rf))
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
        drain = isempty(drain) ? nothing : rev_inds_drain,
    )
    writer = prepare_writer(
        config,
        reader,
        output_path,
        modelmap,
        state_ncnames,
        indices_reverse,
        x_nc,
        y_nc,
        dims_xy,
        nc,
        maxlayers = sbm.maxlayers,
    )

    # for each domain save the directed acyclic graph, the traversion order,
    # and the indices that map it back to the two dimensional grid
    # for reservoirs and lakes this is information is also available (except the graph) 
    # from the initialization functions
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
        (; land, river, reservoir, lake, drain, index_river, frac_toriver),
        (subsurface = (flow = gwf, recharge = gwf.boundaries[1], river = gwf.boundaries[2], drain = gwf.boundaries[3]), 
        land = olf, 
        river = rf),
        sbm,
        Clock(config.starttime, 1, Δt),
        reader,
        writer,
    )

    # read and set states in model object if reinit=false
    if reinit  == false
        set_states(
            instate_path,
            model,
            state_ncnames,
            type = Float64
        )
    end

    # make sure the forcing is already loaded
    # it's fine to run twice, and may help catching errors earlier
    update_forcing!(model)
    if haskey(config.input, "cyclic")
        update_cyclic!(model)
    end

    return model  
end


"update the sbm_gwf model for a single timestep"
function update_sbm_gwf(model)
    @unpack lateral, vertical, network, clock, config = model

    inds_riv = network.index_river
    kinwave_it = get(config.model, "kin_wave_iteration", false)

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
    if get(config.model, "masswasting", false)
        snowflux_frac =
            min.(0.5, lateral.land.sl ./ 5.67) .* min.(1.0, vertical.snow ./ 10000.0)
        maxflux = snowflux_frac .* vertical.snow
        vertical.snow .= accucapacityflux(network.land, vertical.snow, maxflux)
        vertical.snowwater .= accucapacityflux(network.land, vertical.snowwater, vertical.snowwater .* snowflux_frac)
    end

    # update vertical sbm concept until recharge [mm]
    update_until_recharge(vertical, config)

    # set river stage (groundwater) to average h from kinematic wave
    lateral.subsurface.river.stage .= lateral.river.h_av .+ lateral.subsurface.river.bottom

    # determine stable time step for groundwater flow
    Δt_gw = stable_timestep(lateral.subsurface.flow.aquifer) # time step in day (Float64)
    Δt_sbm = (vertical.Δt / 86400.0) # vertical.Δt is in seconds (Float64)
    if Δt_gw < Δt_sbm
        @warn("stable time step Δt $Δt for groundwater flow is smaller than sbm Δt $Δt_sbm")
    end

    Q = zeros(vertical.n)
    # exchange of recharge between vertical sbm concept and groundwater flow domain
    # recharge rate groundwater [m d⁻¹]
    lateral.subsurface.recharge.rate .= vertical.recharge ./ 1000.0 .* (1.0/Δt_sbm)
    # update groundwater domain
    update(lateral.subsurface.flow, Q, Δt_sbm)
    
    # determine excess head [m] (exfiltwater) in groundwater domain (head > surface) and reset head
    exfiltwater = lateral.subsurface.flow.aquifer.head .- min.(lateral.subsurface.flow.aquifer.head, lateral.subsurface.flow.aquifer.top)
    lateral.subsurface.flow.aquifer.head .= min.(lateral.subsurface.flow.aquifer.head, lateral.subsurface.flow.aquifer.top)

    # update vertical sbm concept (runoff, ustorelayerdepth and satwaterdepth)
    update_after_lateralflow(
        vertical,
        (vertical.altitude .- lateral.subsurface.flow.aquifer.head) .* 1000.0, # zi [mm] in vertical concept SBM
        exfiltwater .* 1000.0,
    )

    # update overland flow based on vertical runoff [mm] from vertical sbm concept
    lateral.land.qlat .=
        (vertical.runoff .* vertical.xl .* vertical.yl .* 0.001) ./ lateral.land.Δt ./
        lateral.land.dl
    update(
        lateral.land,
        network.land,
        frac_toriver = network.frac_toriver,
        river = lateral.river.rivercells,
        do_iter = kinwave_it,
    )

    # update river domain with net runoff from vertical sbm concept, overland flow
    # and groundwater flow to the river cells
    net_runoff_river = 
        (vertical.net_runoff_river[inds_riv] .* vertical.xl[inds_riv] .* 
        vertical.yl[inds_riv] .* 0.001) ./ vertical.Δt
    
    # flux from groundwater domain to river (Q to river from drains and groundwater)
    flux_gw = zeros(vertical.n)
    flux_gw[lateral.subsurface.river.index] = -lateral.subsurface.river.flux
    flux_gw[lateral.subsurface.drain.index] = flux_gw[lateral.subsurface.drain.index] - lateral.subsurface.drain.flux
    
    lateral.river.qlat .= 
        (
            flux_gw[inds_riv] ./ lateral.river.Δt .+
            lateral.land.to_river[inds_riv] .+ net_runoff_river
        ) ./ lateral.river.dl
    update(lateral.river, network.river, do_iter = kinwave_it, doy = dayofyear(clock.time))

    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
