"""
    initialize_sbm_gwf_meta_model(config::Config)

Initial part of the sbm_gwf model concept. The model contains:
    - the vertical SBM concept
    - the following lateral components:
        - unconfined aquifer with groundwater flow in four directions (adjacent cells)

The unconfined aquifer contains a recharge, river and a drain (optional) boundary. 

The initial part reads the input settings and data as defined in the Config object. 
Will return a Model that is ready to run.
"""
function initialize_sbm_gwf_meta_model(config::Config)

    # unpack the paths to the NetCDF files
    tomldir = dirname(config)
    static_path = joinpath(tomldir, config.input.path_static)

    reader = prepare_reader(config)
    clock = Clock(config, reader)

    reinit = get(config.model, "reinit", true)::Bool
    do_drains = get(config.model, "drains", false)::Bool
    do_constanthead = get(config.model, "constanthead", false)::Bool
    do_underwaterdrains = get(config.model, "underwaterdrains", false)::Bool

    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, param(config, "input.subcatchment"); allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)

    altitude = ncread(nc, param(config, "input.altitude"); sel = inds, type = Float)
    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)
    riverfrac = zeros(n) # set riverfrac to zero for now, all open water handled by waterfrac

    # initialize vertical SBM concept
    sbm = initialize_sbm(nc, config, riverfrac, inds)

    # unconfined aquifer
    if do_constanthead
        constanthead = ncread(
            nc,
            param(config, "input.lateral.subsurface.constant_head", nothing);
            sel = inds,
            type = Float,
            fill = mv,
        )
        index_constanthead = filter(i -> !isequal(constanthead[i], mv), 1:n)
        constant_head = ConstantHead(constanthead[index_constanthead], index_constanthead)
    else
        constant_head = ConstantHead{Float}(Float[], Int64[])
    end

    conductivity = ncread(
        nc,
        param(config, "input.lateral.subsurface.conductivity", nothing);
        sel = inds,
        type = Float,
    )
    specific_yield = ncread(
        nc,
        param(config, "input.lateral.subsurface.specific_yield", nothing);
        sel = inds,
        type = Float,
    )

    connectivity = Connectivity(inds, rev_inds, xl, yl)
    initial_head = altitude .- Float(0.10) # cold state for groundwater head

    if do_constanthead
        initial_head[constant_head.index] = constant_head.head
    end

    aquifer = UnconfinedAquifer(
        initial_head,
        conductivity,
        altitude,
        altitude .- sbm.soilthickness ./ Float(1000.0),
        xl .* yl,
        specific_yield,
        zeros(Float, connectivity.nconnection),  # conductance
    )

    # reset zi and satwaterdepth with groundwater head from unconfined aquifer 
    sbm.zi .= (altitude .- initial_head) .* 1000.0
    sbm.satwaterdepth .= (sbm.soilthickness .- sbm.zi) .* (sbm.θₛ .- sbm.θᵣ)

    # recharge boundary of unconfined aquifer
    r = fill(mv, n)
    recharge = Recharge(r, zeros(Float, n), collect(1:n))

    # surface water boundaries of unconfined aquifer
    # main water system (h)
    river_h_2d = ncread(nc, param(config, "input.river_h_location"); type = Bool, fill = false)
    inds_riv_h, rev_inds_riv_h = active_indices(river_h_2d, false)
    infiltcond_h = ncread(
        nc,
        param(config, "input.lateral.subsurface.infiltration_conductance_h", nothing);
        sel = inds_riv_h,
        type = Float,
    )
    exfiltcond_h = ncread(
        nc,
        param(config, "input.lateral.subsurface.exfiltration_conductance_h", nothing);
        sel = inds_riv_h,
        type = Float,
    )
    bottom_h = ncread(
        nc,
        param(config, "input.lateral.subsurface.river_bottom_h", nothing);
        sel = inds_riv_h,
        type = Float,
    )
    stage_h = ncread(
        nc,
        param(config, "input.lateral.subsurface.river_stage_h", nothing);
        sel = inds_riv_h,
        type = Float,
    )
    n_h = length(inds_riv_h)
    index_h = findall(x-> x==true, river_h_2d[inds])
    river_h = River(
        stage_h,
        infiltcond_h,
        exfiltcond_h,
        bottom_h,
        fill(mv, n_h),
        index_h,
    )
    # primary water system (p)
    river_p_2d = ncread(nc, param(config, "input.river_p_location"); type = Bool, fill = false)
    inds_riv_p, rev_inds_riv_p = active_indices(river_p_2d, false)
    infiltcond_p = ncread(
        nc,
        param(config, "input.lateral.subsurface.infiltration_conductance_p", nothing);
        sel = inds_riv_p,
        type = Float,
    )
    exfiltcond_p = ncread(
        nc,
        param(config, "input.lateral.subsurface.exfiltration_conductance_p", nothing);
        sel = inds_riv_p,
        type = Float,
    )
    n_p = length(inds_riv_p)
    index_p = findall(x-> x==true, river_p_2d[inds])
    river_p = River(
        fill(mv, n_p),
        infiltcond_p,
        exfiltcond_p,
        fill(mv, n_p),
        fill(mv, n_p),
        index_p,
    )
    # secondary water system (s)
    river_s_2d = ncread(nc, param(config, "input.river_s_location"); type = Bool, fill = false)
    inds_riv_s, rev_inds_riv_s = active_indices(river_s_2d, false)
    infiltcond_s = ncread(
        nc,
        param(config, "input.lateral.subsurface.infiltration_conductance_s", nothing);
        sel = inds_riv_s,
        type = Float,
    )
    exfiltcond_s = ncread(
        nc,
        param(config, "input.lateral.subsurface.exfiltration_conductance_s", nothing);
        sel = inds_riv_s,
        type = Float,
    )
    n_s = length(inds_riv_s)
    index_s = findall(x-> x==true, river_s_2d[inds])
    river_s = River(
        fill(mv, n_s),
        infiltcond_s,
        exfiltcond_s,
        fill(mv, n_s),
        fill(mv, n_s),
        index_s,
    )
    # tertiary water system (t)
    river_t_2d = ncread(nc, param(config, "input.river_t_location"); type = Bool, fill = false)
    inds_riv_t, rev_inds_riv_t = active_indices(river_t_2d, false)
    infiltcond_t = ncread(
        nc,
        param(config, "input.lateral.subsurface.infiltration_conductance_t", nothing);
        sel = inds_riv_t,
        type = Float,
    )
    exfiltcond_t = ncread(
        nc,
        param(config, "input.lateral.subsurface.exfiltration_conductance_t", nothing);
        sel = inds_riv_t,
        type = Float,
    )
    n_t = length(inds_riv_t)
    index_t = findall(x-> x==true, river_t_2d[inds])
    river_t = River(
        fill(mv, n_t),
        infiltcond_t,
        exfiltcond_t,
        fill(mv, n_t),
        fill(mv, n_t),
        index_t,
    )
    index_river = vcat(index_h, index_p, index_s, index_t)
    index_river = unique(index_river)
    initial_head[index_river] = altitude[index_river]

    # under water drains (uwd)
    if do_underwaterdrains
        uwd_2d = ncread(
            nc,
            param(config, "input.lateral.subsurface.uwd_location");
            type = Bool,
            fill = false,
        )
        inds_uwd, rev_inds_uwd = active_indices(uwd_2d, false)
        cond_uwd = ncread(
            nc,
            param(config, "input.lateral.subsurface.uwd_conductance", nothing);
            sel = inds_uwd,
            type = Float,
        )
        elevation_uwd = ncread(
            nc,     
            param(config, "input.lateral.subsurface.uwd_elevation", nothing);
            sel = inds_uwd,
            type = Float,
        )

        n_uwd = length(inds_uwd)
        index_uwd = findall(x-> x==true, uwd_2d[inds])
        uwds = River(
            fill(mv, n_uwd),
            cond_uwd,
            cond_uwd,
            elevation_uwd,
            fill(mv, n_uwd),
            index_uwd,
        )

        uwd = (indices = inds_uwd, reverse_indices = rev_inds_uwd)

    end

    # drain boundaries of unconfined aquifer (optional)
    if do_drains
        # ditch drains (d)
        drain_d_2d = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_d");
            type = Bool,
            fill = false,
        )
        inds_drain_d, rev_inds_drain_d = active_indices(drain_d_2d, false)
        drain_elevation_d = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_elevation_d", nothing);
            sel = inds_drain_d,
            type = Float,
        )
        drain_conductance_d = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_conductance_d", nothing);
            sel = inds_drain_d,
            type = Float,
        )
        n_d = length(inds_drain_d)
        index_d = findall(x-> x==true, drain_d_2d[inds])
        drains_d = Drainage(
            drain_elevation_d,
            drain_conductance_d,
            fill(mv, n_d),
            index_d,
        )
        # pipe drains (p)
        drain_p_2d = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_p");
            type = Bool,
            fill = false,
        )
        inds_drain_p, rev_inds_drain_p = active_indices(drain_p_2d, false)
        drain_elevation_p = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_elevation_p", nothing);
            sel = inds_drain_p,
            type = Float,
        )
        drain_conductance_p = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_conductance_p", nothing);
            sel = inds_drain_p,
            type = Float,
        )
        n_p = length(inds_drain_p)
        index_p = findall(x-> x==true, drain_p_2d[inds])        
        drains_p = Drainage(
            drain_elevation_p,
            drain_conductance_p,
            fill(mv, n_p),
            index_p,
        )
        # surface drainage (s)
        drain_s_2d = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_s");
            type = Bool,
            fill = false,
        )
        inds_drain_s, rev_inds_drain_s = active_indices(drain_s_2d, false)
        drain_elevation_s = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_elevation_s", nothing);
            sel = inds_drain_s,
            type = Float,
        )
        drain_conductance_s = ncread(
            nc,
            param(config, "input.lateral.subsurface.drain_conductance_s", nothing);
            sel = inds_drain_s,
            type = Float,
        )
        n_s = length(inds_drain_s)
        index_s = findall(x-> x==true, drain_s_2d[inds])
        drains_s = Drainage(
            drain_elevation_s,
            drain_conductance_s,
            fill(mv, n_s),
            index_s,
        )
        drain_d = (indices = inds_drain_d, reverse_indices = rev_inds_drain_d)
        drain_p = (indices = inds_drain_p, reverse_indices = rev_inds_drain_p)
        drain_s = (indices = inds_drain_s, reverse_indices = rev_inds_drain_s)

        if do_underwaterdrains
            aquifer_boundaries = AquiferBoundaryCondition[recharge, river_h, river_p, river_s, river_t, drains_d, drains_p, drains_s, uwds]
        else
            aquifer_boundaries = AquiferBoundaryCondition[recharge, river_h, river_p, river_s, river_t, drains_d, drains_p, drains_s]
            uwd = ()
        end
    else
        aquifer_boundaries = AquiferBoundaryCondition[recharge, river_h, river_p, river_s, river_t]
        drain_d = ()
        drain_p = ()
        drain_s = ()
        uwd = ()
    end

    gwf = GroundwaterFlow(aquifer, connectivity, constant_head, aquifer_boundaries)

    # map GroundwaterFlow and its boundaries
    if do_drains && do_underwaterdrains
        subsurface_map = (
            flow = gwf,
            recharge = gwf.boundaries[1],
            river_h = gwf.boundaries[2],
            river_p = gwf.boundaries[3],
            river_s = gwf.boundaries[4],
            river_t = gwf.boundaries[5],
            drain_d = gwf.boundaries[6],
            drain_p = gwf.boundaries[7],
            drain_s = gwf.boundaries[8],
            uwd = gwf.boundaries[9]
        )
    elseif do_drains
        subsurface_map = (
            flow = gwf,
            recharge = gwf.boundaries[1],
            river_h = gwf.boundaries[2],
            river_p = gwf.boundaries[3],
            river_s = gwf.boundaries[4],
            river_t = gwf.boundaries[5],
            drain_d = gwf.boundaries[6],
            drain_p = gwf.boundaries[7],
            drain_s = gwf.boundaries[8],
        )
    else
        subsurface_map =
            (flow = gwf, 
            recharge = gwf.boundaries[1], 
            river_h = gwf.boundaries[2],
            river_p = gwf.boundaries[3],
            river_s = gwf.boundaries[4],
            river_t = gwf.boundaries[5],    
        )
    end

    modelmap =
        (vertical = sbm, lateral = (; subsurface = subsurface_map))
    indices_reverse = (
        land = rev_inds,
        river_h = rev_inds_riv_h,
        river_p = rev_inds_riv_p,
        river_s = rev_inds_riv_s,
        river_t = rev_inds_riv_t,
        drain_d = isempty(drain_d) ? nothing : rev_inds_drain_d,
        drain_p = isempty(drain_p) ? nothing : rev_inds_drain_p,
        drain_s = isempty(drain_s) ? nothing : rev_inds_drain_s,
        uwd = isempty(uwd) ? nothing : rev_inds_uwd,

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

    land = (
        indices = inds,
        reverse_indices = rev_inds,
        xl = xl,
        yl = yl,
        altitude = altitude,
    )
    river_h = (
        indices = inds_riv_h,
        reverse_indices = rev_inds_riv_h,
    )
    river_p = (
        indices = inds_riv_p,
        reverse_indices = rev_inds_riv_p,
    )
    river_s = (
        indices = inds_riv_s,
        reverse_indices = rev_inds_riv_s,
    )
    river_t = (
        indices = inds_riv_t,
        reverse_indices = rev_inds_riv_t,
    )

    model = Model(
        config,
        (; land, river_h, river_p, river_s, river_t, drain_d, drain_p, drain_s, uwd),
        (; subsurface = subsurface_map),
        sbm,
        clock,
        reader,
        writer,
    )

    # read and set states in model object if reinit=false
    if reinit == false
        instate_path = joinpath(tomldir, config.state.path_input)
        state_ncnames = ncnames(config.state)
        set_states(instate_path, model, state_ncnames, type = Float)
    end

    return model
end

"update the sbm_gwf_meta model for a single timestep"
function update_sbm_gwf_meta(model)
    @unpack lateral, vertical, network, clock, config = model

    # extract water levels h_av [m] from the land domain
    # this is used to limit open water evaporation
    vertical.waterlevel_land .= 1000.0

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

    # update vertical sbm concept until recharge [mm]
    update_until_recharge(vertical, config)

    # determine stable time step for groundwater flow
    Δt_gw = stable_timestep(lateral.subsurface.flow.aquifer) # time step in day (Float64)
    Δt_sbm = (vertical.Δt / tosecond(basetimestep)) # vertical.Δt is in seconds (Float64)
    if Δt_gw < Δt_sbm
        @warn(
            "stable time step Δt $Δt_gw for groundwater flow is smaller than sbm Δt $Δt_sbm"
        )
    end

    Q = zeros(vertical.n)
    # exchange of recharge between vertical sbm concept and groundwater flow domain
    # recharge rate groundwater is required in units [m d⁻¹]
    lateral.subsurface.recharge.rate .= vertical.recharge ./ 1000.0 .* (1.0 / Δt_sbm)
    # update groundwater domain
    update(lateral.subsurface.flow, Q, Δt_sbm)

    # determine excess water depth [m] (exfiltwater) in groundwater domain (head > surface)
    # and reset head
    exfiltwater =
        (
            lateral.subsurface.flow.aquifer.head .-
            min.(lateral.subsurface.flow.aquifer.head, lateral.subsurface.flow.aquifer.top)
        ) .* lateral.subsurface.flow.aquifer.specific_yield
    lateral.subsurface.flow.aquifer.head .=
        min.(lateral.subsurface.flow.aquifer.head, lateral.subsurface.flow.aquifer.top)

    # update vertical sbm concept (runoff, ustorelayerdepth and satwaterdepth)
    update_after_subsurfaceflow(
        vertical,
        (network.land.altitude .- lateral.subsurface.flow.aquifer.head) .* 1000.0, # zi [mm] in vertical concept SBM
        exfiltwater .* 1000.0,
    )

    write_output(model)

    # update the clock
    advance!(clock)

    return model
end
