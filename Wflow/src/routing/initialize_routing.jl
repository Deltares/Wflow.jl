
"Initialize subsurface flow routing for model type `sbm`"
function initialize_subsurface_flow(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil::SbmSoilModel,
    ::SbmModel,
)
    (; parameters) = domain.land
    subsurface_flow = LateralSSF(dataset, config, domain.land, soil)

    kh_profile_type = config.model.saturated_hydraulic_conductivity_profile

    if kh_profile_type == VerticalConductivityProfile.exponential ||
       kh_profile_type == VerticalConductivityProfile.exponential_constant
        initialize_lateral_ssf!(
            subsurface_flow,
            parameters,
            subsurface_flow.parameters.kh_profile,
        )
    elseif kh_profile_type == VerticalConductivityProfile.layered ||
           kh_profile_type == VerticalConductivityProfile.layered_exponential
        (; kv_profile) = soil.parameters
        dt = Second(config.time.timestepsecs)
        initialize_lateral_ssf!(subsurface_flow, soil, parameters, kv_profile, tosecond(dt))
    end
    return subsurface_flow
end

"Initialize subsurface flow routing for model type `sbm_gwf`"
function initialize_subsurface_flow(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil::SbmSoilModel,
    ::SbmGwfModel,
)
    (; land, river, drain) = domain

    (; indices, reverse_indices) = land.network
    (; x_length, y_length, area) = land.parameters

    n_cells = length(indices)

    elevation = ncread(
        dataset,
        config,
        "land_surface__elevation";
        optional = false,
        sel = indices,
        type = Float64,
    )

    # unconfined aquifer
    if config.model.constanthead__flag
        constant_head = ConstantHead(dataset, config, indices)
    else
        variables = ConstantHeadVariables(; head = Float64[])
        constant_head = ConstantHead(; variables, index = Int64[])
    end

    connectivity = Connectivity(indices, reverse_indices, x_length, y_length)

    # cold state for groundwater head based on water table depth zi
    initial_head = elevation .- soil.variables.zi / 1000.0
    initial_head[river.network.land_indices] = elevation[river.network.land_indices]
    if config.model.constanthead__flag
        initial_head[constant_head.index] = constant_head.variables.head
    end
    # reset soil (cold) state and related variables based on initial_head (river cells and constanthead)
    if config.model.cold_start__flag
        (;
            zi,
            satwaterdepth,
            ustorecapacity,
            ustorelayerthickness,
            n_unsatlayers,
            total_soilwater_storage,
        ) = soil.variables
        (; theta_s, theta_r, soilthickness, soilwatercapacity, sumlayers, act_thickl) =
            soil.parameters

        @. zi = (elevation - min(elevation, initial_head)) * 1000.0
        @. satwaterdepth = (soilthickness - zi) * (theta_s - theta_r)
        @. ustorecapacity = soilwatercapacity - satwaterdepth
        @. ustorelayerthickness = set_layerthickness(zi, sumlayers, act_thickl)
        @. n_unsatlayers = number_of_active_layers.(ustorelayerthickness)
        @. total_soilwater_storage = satwaterdepth
    end

    bottom = elevation .- soil.parameters.soilthickness ./ 1000.0
    conductance = zeros(connectivity.nconnection)
    aquifer = UnconfinedAquifer(
        dataset,
        config,
        indices,
        elevation,
        bottom,
        area,
        conductance,
        initial_head,
    )

    # river boundary of unconfined aquifer
    gwf_river = GwfRiver(dataset, config, river.network.indices, river.network.land_indices)

    # recharge boundary of unconfined aquifer
    gwf_recharge = Recharge(; n = n_cells)

    # drain boundary of unconfined aquifer (optional)
    if config.model.drain__flag
        gwf_drain = Drainage(dataset, config, indices, drain.network.land_indices)
        aquifer_boundaries = AquiferBoundaries(;
            recharge = gwf_recharge,
            river = gwf_river,
            drain = gwf_drain,
        )
    else
        aquifer_boundaries = AquiferBoundaries(; recharge = gwf_recharge, river = gwf_river)
    end

    cfl = config.model.subsurface_water_flow__alpha_coefficient
    @info "Numerical stability coefficient for groundwater flow `alpha`: `$cfl`."
    timestepping = TimeStepping(; cfl)

    subsurface_flow = GroundwaterFlow(;
        timestepping,
        aquifer,
        connectivity,
        constanthead = constant_head,
        boundaries = aquifer_boundaries,
    )
    return subsurface_flow
end

"Initialize kinematic wave or local inertial overland flow routing"
function initialize_overland_flow(dataset::NCDataset, config::Config, domain::Domain)
    (; land_routing) = config.model

    if land_routing == RoutingType.kinematic_wave
        overland_flow = KinWaveOverlandFlow(dataset, config, domain.land)
    elseif land_routing == RoutingType.local_inertial
        overland_flow = LocalInertialOverlandFlow(dataset, config, domain)
    end
    return overland_flow
end

"""
Initialize kinematic wave or local inertial overland flow routing, including optional
reservoirs.
"""
function initialize_river_flow(dataset::NCDataset, config::Config, domain::Domain)
    (; river_routing) = config.model

    reservoir =
        config.model.reservoir__flag ?
        Reservoir(dataset, config, domain.reservoir.network) : nothing

    if river_routing == RoutingType.kinematic_wave
        river_flow = KinWaveRiverFlow(dataset, config, domain.river, reservoir)
    elseif river_routing == RoutingType.local_inertial
        river_flow = LocalInertialRiverFlow(dataset, config, domain.river, reservoir)
    end
end

"Initialize `Routing` for model types `sbm` and `sbm_gwf`"
function Routing(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil::SbmSoilModel,
    type,
)
    subsurface_flow = initialize_subsurface_flow(dataset, config, domain, soil, type)
    overland_flow = initialize_overland_flow(dataset, config, domain)
    river_flow = initialize_river_flow(dataset, config, domain)

    routing = Routing(; subsurface_flow, overland_flow, river_flow)
    return routing
end

"Initialize `Routing` for model type `sediment`"
function Routing(dataset::NCDataset, config::Config, domain::Domain, soil::SoilLoss)
    overland_flow = OverlandFlowSediment(dataset, config, domain.land, soil)
    river_flow = RiverSediment(dataset, config, domain.river)

    routing = Routing(; overland_flow, river_flow)
    return routing
end
