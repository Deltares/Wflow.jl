"""
    initialize_sbm_model(config::Config)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_sbm_model(config::Config)
    routing_options = ("kinematic-wave", "local-inertial")
    river_routing = get_options(
        config.model,
        "river_routing",
        routing_options,
        "kinematic-wave",
    )::String
    land_routing =
        get_options(config.model, "land_routing", routing_options, "kinematic-wave")::String

    reservoirs = do_reservoirs
    lakes = do_lakes
    glacier = get(config.model, "glacier", false)::Bool
    masswasting = get(config.model, "masswasting", false)::Bool
    @info "General model settings" reservoirs lakes snow masswasting glacier

    n_land_cells = length(indices)

    river_location_2d = ncread(
        dataset,
        config,
        "river_location";
        optional = false,
        type = Bool,
        fill = false,
    )
    river_location = river_location_2d[indices]
    river_width_2d = ncread(
        dataset,
        config,
        "routing.river_flow.width";
        optional = false,
        type = Float,
        fill = 0,
    )
    river_width = river_width_2d[indices]
    river_length_2d = ncread(
        dataset,
        config,
        "routing.river_flow.length";
        optional = false,
        type = Float,
        fill = 0,
    )
    river_length = river_length_2d[indices]
    river_fraction =
        get_river_fraction(river_location, river_length, river_width, x_length, y_length)

    land_hydrology = LandHydrologySBM(dataset, config, river_fraction, indices)

    inds_river, reverse_inds_river = active_indices(river_location_2d, 0)
    n_river_cells = length(inds_river)

    # reservoirs
    if do_reservoirs
        reservoir, reservoir_network, inds_reservoir_map2river, pits =
            SimpleReservoir(dataset, config, inds_river, n_river_cells, pits)
        network_reservoir = NetworkWaterBody(; reservoir_network...)
    else
        network_reservoir = NetworkWaterBody()
        inds_reservoir_map2river = fill(0, n_river_cells)
        reservoir = nothing
    end

    ldd_2d = ncread(dataset, config, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[indices]
    if do_pits
        pits_2d =
            ncread(dataset, config, "pits"; optional = false, type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, indices)
    end
    flow_length = map(get_flow_length, ldd, x_length, y_length)
    flow_width = (x_length .* y_length) ./ flow_length

    # check if lateral subsurface flow component is defined for the SBM model, when coupled
    # to another groundwater model, this component is not defined in the TOML file.
    do_lateral_ssf = haskey(config.input.routing, "subsurface_flow")
    if do_lateral_ssf
        subsurface_flow = LateralSSF(
            dataset,
            config,
            indices;
            soil = land_hydrology.soil,
            slope = land_slope,
            flow_length,
            flow_width,
            x_length,
            y_length,
        )
        # update variables `ssf`, `ssfmax` and `kh` (layered profile) based on ksat_profile
        kh_profile_type = get(config.input.land, "ksat_profile", "exponential")::String
        if kh_profile_type == "exponential" || kh_profile_type == "exponential_constant"
            initialize_lateralssf!(subsurface_flow, subsurface_flow.parameters.kh_profile)
        elseif kh_profile_type == "layered" || kh_profile_type == "layered_exponential"
            (; kv_profile) = land_hydrology.soil.parameters
            initialize_lateralssf!(
                subsurface_flow,
                land_hydrology.soil,
                kv_profile,
                tosecond(clock.dt),
            )
        end
    else
        # when the SBM model is coupled (BMI) to a groundwater model, the following
        # variables are expected to be exchanged from the groundwater model.
        subsurface_flow = GroundwaterExchange(n_land_cells)
    end

    graph = flowgraph(ldd, indices, pcr_dir)
    ldd_river = ldd_2d[inds_river]
    if do_pits
        ldd_river = set_pit_ldd(pits_2d, ldd_river, inds_river)
    end
    graph_river = flowgraph(ldd_river, inds_river, pcr_dir)

    # land indices where river is located
    inds_land_map2river = filter(i -> !isequal(river_location[i], 0), 1:n_land_cells)

    if land_routing == "kinematic-wave"
        overland_flow = KinWaveOverlandFlow(
            dataset,
            config,
            indices;
            slope = land_slope,
            flow_length,
            flow_width = map(det_surfacewidth, flow_width, river_width, river_location),
        )
    elseif land_routing == "local-inertial"
        inds_river_map2land = reverse_inds_river[indices] # not filtered (with zeros)
        overland_flow, edge_indices = LocalInertialOverlandFlow(
            dataset,
            config,
            indices;
            modelsize_2d,
            reverse_indices,
            x_length,
            y_length,
            river_width = river_width_2d[inds_river],
            graph_river,
            ldd_river,
            inds_river,
            river_location,
            waterbody = !=(0).(inds_reservoir_map2river + inds_lake_map2river),
        )
    end

    river_length = river_length_2d[inds_river]
    river_width = river_width_2d[inds_river]
    minimum(river_length) > 0 || error("river length must be positive on river cells")
    minimum(river_width) > 0 || error("river width must be positive on river cells")
    if river_routing == "kinematic-wave"
        river_flow = KinWaveRiverFlow(
            dataset,
            config,
            inds_river;
            river_length,
            river_width,
            reservoir = reservoir,
            lake = lake,
        )
    elseif river_routing == "local-inertial"
        river_flow, nodes_at_edge = LocalInertialRiverFlow(
            dataset,
            config,
            inds_river;
            graph_river,
            ldd_river,
            river_length,
            river_width,
            reservoir,
            lake,
            waterbody = !=(0).(inds_reservoir_map2river + inds_lake_map2river),
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
    if land_routing == "kinematic-wave" ||
       river_routing == "kinematic-wave" ||
       do_lateral_ssf
    end
    if land_routing == "kinematic-wave" || do_lateral_ssf
    end
    if river_routing == "kinematic-wave"
        min_streamorder_river = get(config.model, "min_streamorder_river", 6)
        toposort_river = topological_sort_by_dfs(graph_river)
        river_pit_inds = findall(x -> x == 5, ldd_river)
        order_of_river_subdomains, river_subdomain_inds, toposort_river_subdomain =
            kinwave_set_subdomains(
                graph_river,
                toposort_river,
                river_pit_inds,
                streamorder[inds_land_map2river],
                min_streamorder_river,
            )
    end

    if nthreads() > 1
        if river_routing == "kinematic-wave"
            @info "Parallel execution of kinematic wave" min_streamorder_land min_streamorder_river
        elseif land_routing == "kinematic-wave" || do_lateral_ssf
            @info "Parallel execution of kinematic wave" min_streamorder_land
        end
    end

    modelmap =
        (land = land_hydrology, routing = (; subsurface_flow, overland_flow, river_flow))
    indices_reverse = (
        land = reverse_indices,
        river = reverse_inds_river,
        reservoir = isnothing(reservoir) ? nothing : network_reservoir.reverse_indices,
        lake = isnothing(lake) ? nothing : lake_network.reverse_indices,
    )
    (; maxlayers) = land_hydrology.soil.parameters
    writer = prepare_writer(
        config,
        modelmap,
        indices_reverse,
        x_coords,
        y_coords,
        dataset;
        extra_dim = (name = "layer", value = Float64.(1:(maxlayers))),
    )
    close(dataset)

    if land_routing == "local-inertial"
        @reset network_land.river_indices = inds_river_map2land
        @reset network_land.edge_indices = edge_indices
    end
    if do_water_demand
        # exclude waterbodies for local surface and ground water abstraction
        inds_riv_2d = copy(reverse_inds_river)
        inds_2d = zeros(Bool, modelsize_2d)
        if !isnothing(reservoir)
            inds_cov = collect(Iterators.flatten(network_reservoir.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 1
        end
        if !isnothing(lake)
            inds_cov = collect(Iterators.flatten(lake_network.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 1
        end
        @reset network_land.river_inds_excl_waterbody = inds_riv_2d[indices]
        @reset network_land.waterbody = inds_2d[indices]
    end
    if river_routing == "kinematic-wave"
        @reset network_river.upstream_nodes =
            filter_upsteam_nodes(graph_river, pits[inds_river])
        @reset network_river.order_of_subdomains = order_of_river_subdomains
        @reset network_river.order_subdomain = toposort_river_subdomain
        @reset network_river.subdomain_indices = river_subdomain_inds
        @reset network_river.order = toposort_river
    elseif river_routing == "local-inertial"
        @reset network_river.nodes_at_edge = NodesAtEdge(; nodes_at_edge...)
        @reset network_river.edges_at_node =
            EdgesAtNode(; adjacent_edges_at_node(graph_river, nodes_at_edge)...)
    end

    routing = Routing(; subsurface_flow, overland_flow, river_flow)
end

"update SBM model for a single timestep"
function update!(model::AbstractModel{<:SbmModel})
    (; routing, land, network, clock, config) = model
    dt = tosecond(clock.dt)
    do_water_demand = haskey(config.model, "water_demand")
    (; kv_profile) = land.soil.parameters

    update_until_recharge!(model)
    # exchange of recharge between SBM soil model and subsurface flow domain
    routing.subsurface_flow.boundary_conditions.recharge .=
        land.soil.variables.recharge ./ 1000.0
    if do_water_demand
        @. routing.subsurface_flow.boundary_conditions.recharge -=
            land.allocation.variables.act_groundwater_abst / 1000.0
    end
    routing.subsurface_flow.boundary_conditions.recharge .*=
        routing.subsurface_flow.parameters.flow_width
    routing.subsurface_flow.variables.zi .= land.soil.variables.zi ./ 1000.0
    # update lateral subsurface flow domain (kinematic wave)
    kh_layered_profile!(land.soil, routing.subsurface_flow, kv_profile, dt)
    update!(routing.subsurface_flow, network.land, clock.dt / basetimestep)
    update_after_subsurfaceflow!(model)
    update_total_water_storage!(model)
    return nothing
end

"""
    update_until_recharge!model::AbstractModel{<:SbmModel})

Update SBM model until recharge for a single timestep. This function is also accessible
through BMI, to couple the SBM model to an external groundwater model.
"""
function update_until_recharge!(model::AbstractModel{<:SbmModel})
    (; routing, land, network, clock, config) = model
    dt = tosecond(clock.dt)
    update!(land, routing, network, config, dt)
    return nothing
end

"""
    update_after_subsurfaceflow!(model::AbstractModel{<:SbmModel})

Update SBM model after subsurface flow for a single timestep. This function is also
accessible through BMI, to couple the SBM model to an external groundwater model.
"""
function update_after_subsurfaceflow!(model::AbstractModel{<:SbmModel})
    (; routing, land) = model
    (; soil, runoff, demand) = land
    (; subsurface_flow) = routing

    # update SBM soil model (runoff, ustorelayerdepth and satwaterdepth)
    update!(soil, (; runoff, demand, subsurface_flow))

    surface_routing!(model)

    return nothing
end

"""
Update of the total water storage at the end of each timestep per model cell.

This is done here at model level.
"""
function update_total_water_storage!(model::AbstractModel{<:SbmModel})
    (; routing, land, network) = model

    # Update the total water storage based on land states
    # TODO Maybe look at routing in the near future
    update_total_water_storage!(
        land,
        network.river.land_indices,
        network.land.area,
        routing.river_flow,
        routing.overland_flow,
    )
    return nothing
end

function set_states!(model::AbstractModel{<:Union{SbmModel, SbmGwfModel}})
    (; routing, land, network, config) = model
    land_v = routing.overland_flow.variables
    land_p = routing.overland_flow.parameters
    river_v = routing.river_flow.variables
    river_p = routing.river_flow.parameters

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
                land.soil.parameters.soilthickness .-
                land.soil.variables.satwaterdepth ./
                (land.soil.parameters.theta_s .- land.soil.parameters.theta_r),
            )
        land.soil.variables.zi .= zi
        if land_routing == "kinematic-wave"
            # make sure land cells with zero flow width are set to zero q and h
            for i in eachindex(land_p.flow_width)
                if land_p.flow_width[i] <= 0.0
                    land_v.q[i] = 0.0
                    land_v.h[i] = 0.0
                end
            end
            land_v.volume .= land_v.h .* land_p.flow_width .* land_p.flow_length
        elseif land_routing == "local-inertial"
            for i in eachindex(routing.overland_flow.volume)
                if land_p.rivercells[i]
                    j = network.land.index_river[i]
                    if land_v.h[i] > 0.0
                        land_v.volume[i] =
                            land_v.h[i] * land_p.xl[i] * land_p.yl[i] +
                            land_p.bankfull_volume[j]
                    else
                        land_v.volume[i] =
                            river_v.h[j] * river_p.flow_width[j] * river_p.flow_length[j]
                    end
                else
                    routing.overland_flow.volume[i] =
                        routing.overland_flow.h[i] *
                        routing.overland_flow.xl[i] *
                        routing.overland_flow.yl[i]
                end
            end
        end
        # only set active cells for river (ignore boundary conditions/ghost points)
        river_v.volume[1:nriv] .=
            river_v.h[1:nriv] .* river_p.flow_width[1:nriv] .* river_p.flow_length[1:nriv]

        if floodplain_1d
            initialize_volume!(routing.river_flow, nriv)
        end

        if do_lakes
            # storage must be re-initialized after loading the state with the current
            # waterlevel otherwise the storage will be based on the initial water level
            lakes = routing.river_flow.boundary_conditions.lake
            lakes.storage .=
                initialize_storage(lakes.storfunc, lakes.area, lakes.waterlevel, lakes.sh)
        end
    else
        @info "Set initial conditions from default values."
    end
    return nothing
end
