"""
    initialize_sbm_model(config::Config)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_sbm_model(config::Config)
    model_type = config.model.type::String
    @info "Initialize model variables for model type `$model_type`."

    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)

    reader = prepare_reader(config)
    clock = Clock(config, reader)
    dt = clock.dt

    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_pits = get(config.model, "pits", false)::Bool

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool
    routing_options = ("kinematic-wave", "local-inertial")
    floodplain_1d = get(config.model, "floodplain_1d", false)::Bool
    river_routing = get_options(
        config.model,
        "river_routing",
        routing_options,
        "kinematic-wave",
    )::String
    land_routing =
        get_options(config.model, "land_routing", routing_options, "kinematic-wave")::String
    do_water_demand = haskey(config.model, "water_demand")

    snow = get(config.model, "snow", false)::Bool
    reservoirs = do_reservoirs
    lakes = do_lakes
    glacier = get(config.model, "glacier", false)::Bool
    masswasting = get(config.model, "masswasting", false)::Bool
    @info "General model settings" reservoirs lakes snow masswasting glacier

    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, config, "subcatchment"; optional = false, allow_missing = true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)
    modelsize_2d = size(subcatch_2d)

    river_2d =
        ncread(nc, config, "river_location"; optional = false, type = Bool, fill = false)
    river = river_2d[inds]
    riverwidth_2d =
        ncread(nc, config, "lateral.river.width"; optional = false, type = Float, fill = 0)
    riverwidth = riverwidth_2d[inds]
    riverlength_2d =
        ncread(nc, config, "lateral.river.length"; optional = false, type = Float, fill = 0)
    riverlength = riverlength_2d[inds]

    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc; outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)
    riverfrac = river_fraction(river, riverlength, riverwidth, xl, yl)

    lhm = LandHydrologySBM(nc, config, riverfrac, inds)

    inds_riv, rev_inds_riv = active_indices(river_2d, 0)
    nriv = length(inds_riv)

    # reservoirs
    pits = zeros(Bool, modelsize_2d)
    if do_reservoirs
        reservoirs, resindex, reservoir, pits =
            initialize_simple_reservoir(config, nc, inds_riv, nriv, pits, tosecond(dt))
    else
        reservoir = ()
        reservoirs = nothing
        resindex = fill(0, nriv)
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits =
            initialize_lake(config, nc, inds_riv, nriv, pits, tosecond(dt))
    else
        lake = ()
        lakes = nothing
        lakeindex = fill(0, nriv)
    end

    ldd_2d = ncread(nc, config, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d = ncread(nc, config, "pits"; optional = false, type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, inds)
    end

    landslope =
        ncread(nc, config, "lateral.land.slope"; optional = false, sel = inds, type = Float)
    clamp!(landslope, 0.00001, Inf)

    dl = map(detdrainlength, ldd, xl, yl)
    dw = (xl .* yl) ./ dl

    # check if lateral subsurface flow component is defined for the SBM model, when coupled
    # to another groundwater model, this component is not defined in the TOML file.
    subsurface_flow = haskey(config.input.lateral, "subsurface")
    if subsurface_flow
        khfrac = ncread(
            nc,
            config,
            "lateral.subsurface.ksathorfrac";
            sel = inds,
            defaults = 1.0,
            type = Float,
        )

        (; theta_s, theta_r, soilthickness) = lhm.soil.parameters
        (; zi) = lhm.soil.variables
        ssf_soilthickness = soilthickness .* 0.001
        ssf_zi = zi .* 0.001

        kh_profile_type = get(config.input.vertical, "ksat_profile", "exponential")::String
        if kh_profile_type == "exponential"
            (; kv_0, f) = lhm.soil.parameters.kv_profile
            kh_0 = khfrac .* kv_0 .* 0.001 .* (basetimestep / dt)
            kh_profile = KhExponential(kh_0, f .* 1000.0)
        elseif kh_profile_type == "exponential_constant"
            (; z_exp) = lhm.soil.parameters.kv_profile
            (; kv_0, f) = lhm.soil.parameters.kv_profile.exponential
            kh_0 = khfrac .* kv_0 .* 0.001 .* (basetimestep / dt)
            exp_profile = KhExponential(kh_0, f .* 1000.0)
            kh_profile = KhExponentialConstant(exp_profile, z_exp .* 0.001)
        elseif kh_profile_type == "layered" || kh_profile_type == "layered_exponential"
            kh_profile = KhLayered(fill(mv, n))
        end

        # unit for lateral subsurface flow component is [m³ d⁻¹], kv_0 [mm Δt⁻¹]
        ssf = LateralSSF{Float, typeof(kh_profile)}(;
            kh_profile = kh_profile,
            khfrac = khfrac,
            zi = ssf_zi,
            soilthickness = ssf_soilthickness,
            theta_s,
            theta_r,
            dt = dt / basetimestep,
            slope = landslope,
            dl = dl,
            dw = dw,
            exfiltwater = fill(mv, n),
            recharge = fill(mv, n),
            ssf = fill(mv, n),
            ssfin = fill(mv, n),
            ssfmax = fill(mv, n),
            to_river = zeros(n),
            volume = (theta_s .- theta_r) .* (ssf_soilthickness .- ssf_zi) .* (xl .* yl),
        )
        # update variables `ssf`, `ssfmax` and `kh` (layered profile) based on ksat_profile
        if kh_profile_type == "exponential" || kh_profile_type == "exponential_constant"
            initialize_lateralssf!(ssf, kh_profile)
        elseif kh_profile_type == "layered" || kh_profile_type == "layered_exponential"
            (; kv_profile) = lhm.soil.parameters
            initialize_lateralssf!(ssf, lhm.soil, kv_profile, tosecond(dt))
        end
    else
        # when the SBM model is coupled (BMI) to a groundwater model, the following
        # variables are expected to be exchanged from the groundwater model.
        ssf = GroundwaterExchange{Float}(;
            dt = dt / basetimestep,
            exfiltwater = fill(mv, n),
            zi = fill(mv, n),
            to_river = fill(mv, n),
            ssf = zeros(n),
        )
    end

    graph = flowgraph(ldd, inds, pcr_dir)
    ldd_riv = ldd_2d[inds_riv]
    if do_pits
        ldd_riv = set_pit_ldd(pits_2d, ldd_riv, inds_riv)
    end
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, landslope, n)

    inds_allocation_areas = Vector{Int}[]
    inds_riv_allocation_areas = Vector{Int}[]
    if do_water_demand
        areas = unique(lhm.allocation.parameters.areas)
        for a in areas
            area_index = findall(x -> x == a, lhm.allocation.parameters.areas)
            push!(inds_allocation_areas, area_index)
            area_riv_index =
                findall(x -> x == a, lhm.allocation.parameters.areas[index_river])
            push!(inds_riv_allocation_areas, area_riv_index)
        end
    end

    if land_routing == "kinematic-wave"
        olf = initialize_surfaceflow_land(
            nc,
            config,
            inds;
            sl = landslope,
            dl,
            width = map(det_surfacewidth, dw, riverwidth, river),
            iterate = kinwave_it,
            tstep = kw_land_tstep,
            dt,
        )
    elseif land_routing == "local-inertial"
        index_river_nf = rev_inds_riv[inds] # not filtered (with zeros)
        olf, indices = initialize_shallowwater_land(
            nc,
            config,
            inds;
            modelsize_2d,
            indices_reverse = rev_inds,
            xlength = xl,
            ylength = yl,
            riverwidth = riverwidth_2d[inds_riv],
            graph_riv,
            ldd_riv,
            inds_riv,
            river,
            waterbody = !=(0).(resindex + lakeindex),
            dt,
        )
    end

    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]
    minimum(riverlength) > 0 || error("river length must be positive on river cells")
    minimum(riverwidth) > 0 || error("river width must be positive on river cells")
    if river_routing == "kinematic-wave"
        rf = initialize_surfaceflow_river(
            nc,
            config,
            inds_riv;
            dl = riverlength,
            width = riverwidth,
            reservoir_index = resindex,
            reservoir = reservoirs,
            lake_index = lakeindex,
            lake = lakes,
            iterate = kinwave_it,
            tstep = kw_river_tstep,
            dt = dt,
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
            dt = dt,
            floodplain = floodplain_1d,
        )
    else
        error(
            """An unknown "river_routing" method is specified in the TOML file ($river_routing).
            This should be "kinematic-wave" or "local-inertial".
            """,
        )
    end

    # setup subdomains for the land and river kinematic wave domain, if nthreads = 1
    # subdomain is equal to the complete domain
    toposort = topological_sort_by_dfs(graph)
    if land_routing == "kinematic-wave" ||
       river_routing == "kinematic-wave" ||
       subsurface_flow
        streamorder = stream_order(graph, toposort)
    end
    if land_routing == "kinematic-wave" || subsurface_flow
        toposort = topological_sort_by_dfs(graph)
        index_pit_land = findall(x -> x == 5, ldd)
        min_streamorder_land = get(config.model, "min_streamorder_land", 5)
        subbas_order, indices_subbas, topo_subbas = kinwave_set_subdomains(
            graph,
            toposort,
            index_pit_land,
            streamorder,
            min_streamorder_land,
        )
    end
    if river_routing == "kinematic-wave"
        min_streamorder_river = get(config.model, "min_streamorder_river", 6)
        toposort_riv = topological_sort_by_dfs(graph_riv)
        index_pit_river = findall(x -> x == 5, ldd_riv)
        subriv_order, indices_subriv, topo_subriv = kinwave_set_subdomains(
            graph_riv,
            toposort_riv,
            index_pit_river,
            streamorder[index_river],
            min_streamorder_river,
        )
    end

    if nthreads() > 1
        if river_routing == "kinematic-wave"
            @info "Parallel execution of kinematic wave" min_streamorder_land min_streamorder_river
        elseif land_routing == "kinematic-wave" || subsurface_flow
            @info "Parallel execution of kinematic wave" min_streamorder_land
        end
    end

    modelmap = (vertical = lhm, lateral = (subsurface = ssf, land = olf, river = rf))
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
    )
    (; maxlayers) = lhm.soil.parameters
    writer = prepare_writer(
        config,
        modelmap,
        indices_reverse,
        x_nc,
        y_nc,
        nc;
        extra_dim = (name = "layer", value = Float64.(1:(maxlayers))),
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
        upstream_nodes = filter_upsteam_nodes(graph, pits[inds]),
        subdomain_order = subbas_order,
        topo_subdomain = topo_subbas,
        indices_subdomain = indices_subbas,
        order = toposort,
        indices = inds,
        reverse_indices = rev_inds,
        area = xl .* yl,
        slope = landslope,
        indices_allocation_areas = inds_allocation_areas,
    )
    if land_routing == "local-inertial"
        land = merge(land, (index_river = index_river_nf, staggered_indices = indices))
    end
    if do_water_demand
        # exclude waterbodies for local surface and ground water abstraction
        inds_riv_2d = copy(rev_inds_riv)
        inds_2d = ones(Bool, modelsize_2d)
        if !isempty(reservoir)
            inds_cov = collect(Iterators.flatten(reservoir.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 0
        end
        if !isempty(lake)
            inds_cov = collect(Iterators.flatten(lake.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 0
        end
        land = merge(land, (index_river_wb = inds_riv_2d[inds], index_wb = inds_2d[inds]))
    end
    if river_routing == "kinematic-wave"
        river = (
            graph = graph_riv,
            indices = inds_riv,
            reverse_indices = rev_inds_riv,
            # reservoir and lake index
            reservoir_index = resindex,
            lake_index = lakeindex,
            reservoir_index_f = filter(x -> x ≠ 0, resindex),
            lake_index_f = filter(x -> x ≠ 0, lakeindex),
            # specific for kinematic_wave
            upstream_nodes = filter_upsteam_nodes(graph_riv, pits[inds_riv]),
            subdomain_order = subriv_order,
            topo_subdomain = topo_subriv,
            indices_subdomain = indices_subriv,
            order = toposort_riv,
            # water allocation areas
            indices_allocation_areas = inds_riv_allocation_areas,
            area = xl[index_river] .* yl[index_river],
        )
    elseif river_routing == "local-inertial"
        river = (
            graph = graph_riv,
            indices = inds_riv,
            reverse_indices = rev_inds_riv,
            # reservoir and lake index
            reservoir_index = resindex,
            lake_index = lakeindex,
            reservoir_index_f = filter(x -> x ≠ 0, resindex),
            lake_index_f = filter(x -> x ≠ 0, lakeindex),
            # specific for local-inertial
            nodes_at_link = nodes_at_link,
            links_at_node = adjacent_links_at_node(graph_riv, nodes_at_link),
            # water allocation areas
            indices_allocation_areas = inds_riv_allocation_areas,
            area = xl[index_river] .* yl[index_river],
        )
    end

    model = Model(
        config,
        (; land, river, reservoir, lake, index_river, frac_toriver),
        (subsurface = ssf, land = olf, river = rf),
        lhm,
        clock,
        reader,
        writer,
        SbmModel(),
    )

    set_states!(model)

    @info "Initialized model"
    return model
end

"update SBM model for a single timestep"
function update!(model::Model{N, L, V, R, W, T}) where {N, L, V, R, W, T <: SbmModel}
    (; lateral, vertical, network, config) = model
    do_water_demand = haskey(config.model, "water_demand")
    (; kv_profile) = vertical.soil.parameters

    update_until_recharge!(model)
    # exchange of recharge between SBM soil model and subsurface flow domain
    lateral.subsurface.recharge .= vertical.soil.variables.recharge ./ 1000.0
    if do_water_demand
        @. lateral.subsurface.recharge -=
            vertical.allocation.variables.act_groundwater_abst / 1000.0
    end
    lateral.subsurface.recharge .*= lateral.subsurface.dw
    lateral.subsurface.zi .= vertical.soil.variables.zi ./ 1000.0
    # update lateral subsurface flow domain (kinematic wave)
    kh_layered_profile!(vertical.soil, lateral.subsurface, kv_profile, vertical.dt)
    update!(lateral.subsurface, network.land, network.frac_toriver)
    update_after_subsurfaceflow!(model)
    update_total_water_storage!(model)
    return nothing
end

"""
    update_until_recharge!(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}

Update SBM model until recharge for a single timestep. This function is also accessible
through BMI, to couple the SBM model to an external groundwater model.
"""
function update_until_recharge!(
    model::Model{N, L, V, R, W, T},
) where {N, L, V, R, W, T <: SbmModel}
    (; lateral, vertical, network, config) = model
    update!(vertical, lateral, network, config)
    return nothing
end

"""
    update_after_subsurfaceflow!(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}

Update SBM model after subsurface flow for a single timestep. This function is also
accessible through BMI, to couple the SBM model to an external groundwater model.
"""
function update_after_subsurfaceflow!(
    model::Model{N, L, V, R, W, T},
) where {N, L, V, R, W, T <: SbmModel}
    (; lateral, vertical) = model
    (; soil, runoff, demand) = vertical
    (; subsurface) = lateral

    # update SBM soil model (runoff, ustorelayerdepth and satwaterdepth)
    update!(soil, (; runoff, demand, subsurface))

    ssf_toriver = lateral.subsurface.to_river ./ tosecond(basetimestep)
    surface_routing!(model; ssf_toriver = ssf_toriver)

    return nothing
end

"""
Update of the total water storage at the end of each timestep per model cell.

This is done here at model level.
"""
function update_total_water_storage!(
    model::Model{N, L, V, R, W, T},
) where {N, L, V, R, W, T <: SbmModel}
    (; lateral, vertical, network) = model

    # Update the total water storage based on vertical states
    # TODO Maybe look at routing in the near future
    update_total_water_storage!(
        vertical,
        network.index_river,
        network.land.area,
        lateral.river,
        lateral.land,
    )
    return nothing
end

function set_states!(
    model::Model{N, L, V, R, W, T},
) where {N, L, V, R, W, T <: Union{SbmModel, SbmGwfModel}}
    (; lateral, vertical, network, config) = model

    reinit = get(config.model, "reinit", true)::Bool
    routing_options = ("kinematic-wave", "local-inertial")
    land_routing =
        get_options(config.model, "land_routing", routing_options, "kinematic-wave")::String
    do_lakes = get(config.model, "lakes", false)::Bool
    floodplain_1d = get(config.model, "floodplain_1d", false)::Bool

    # read and set states in model object if reinit=false
    if reinit == false
        nriv = length(network.river.indices)
        instate_path = input_path(config, config.state.path_input)
        @info "Set initial conditions from state file `$instate_path`."
        set_states!(instate_path, model; type = Float, dimname = :layer)
        # update zi for SBM soil model
        zi =
            max.(
                0.0,
                vertical.soil.parameters.soilthickness .-
                vertical.soil.variables.satwaterdepth ./
                (vertical.soil.parameters.theta_s .- vertical.soil.parameters.theta_r),
            )
        vertical.soil.variables.zi .= zi
        if land_routing == "kinematic-wave"
            # make sure land cells with zero flow width are set to zero q and h
            for i in eachindex(lateral.land.width)
                if lateral.land.width[i] <= 0.0
                    lateral.land.q[i] = 0.0
                    lateral.land.h[i] = 0.0
                end
            end
            lateral.land.volume .= lateral.land.h .* lateral.land.width .* lateral.land.dl
        elseif land_routing == "local-inertial"
            for i in eachindex(lateral.land.volume)
                if lateral.land.rivercells[i]
                    j = network.land.index_river[i]
                    if lateral.land.h[i] > 0.0
                        lateral.land.volume[i] =
                            lateral.land.h[i] * lateral.land.xl[i] * lateral.land.yl[i] +
                            lateral.river.bankfull_volume[j]
                    else
                        lateral.land.volume[i] =
                            lateral.river.h[j] *
                            lateral.river.width[j] *
                            lateral.river.dl[j]
                    end
                else
                    lateral.land.volume[i] =
                        lateral.land.h[i] * lateral.land.xl[i] * lateral.land.yl[i]
                end
            end
        end
        # only set active cells for river (ignore boundary conditions/ghost points)
        lateral.river.volume[1:nriv] .=
            lateral.river.h[1:nriv] .* lateral.river.width[1:nriv] .*
            lateral.river.dl[1:nriv]

        if floodplain_1d
            initialize_volume!(lateral.river, nriv)
        end

        if do_lakes
            # storage must be re-initialized after loading the state with the current
            # waterlevel otherwise the storage will be based on the initial water level
            lakes = lateral.river.lake
            lakes.storage .=
                initialize_storage(lakes.storfunc, lakes.area, lakes.waterlevel, lakes.sh)
        end
    else
        @info "Set initial conditions from default values."
    end
    return nothing
end
