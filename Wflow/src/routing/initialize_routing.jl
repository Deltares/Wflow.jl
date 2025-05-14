
"Initialize subsurface flow routing for model type `sbm`"
function initialize_subsurface_flow(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil::SbmSoilModel,
    type::SbmModel,
)
    (; parameters) = domain.land
    subsurface_flow = LateralSSF(dataset, config, domain.land, soil)

    kh_profile_type =
        get(config.model, "saturated_hydraulic_conductivity_profile", "exponential")::String

    if kh_profile_type == "exponential" || kh_profile_type == "exponential_constant"
        initialize_lateral_ssf!(
            subsurface_flow,
            parameters,
            subsurface_flow.parameters.kh_profile,
        )
    elseif kh_profile_type == "layered" || kh_profile_type == "layered_exponential"
        (; kv_profile) = soil.parameters
        dt = Second(config.time.timestepsecs)
        initialize_lateral_ssf!(
            subsurface_flow,
            soil,
            parameters,
            kv_profile,
            Float(tosecond(dt)),
        )
    end
    return subsurface_flow
end

"Initialize subsurface flow routing for model type `sbm_gwf`"
function initialize_subsurface_flow(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil::SbmSoilModel,
    type::SbmGwfModel,
)
    do_drains = get(config.model, "drain__flag", false)::Bool
    do_constanthead = get(config.model, "constanthead__flag", false)::Bool

    (; land, river, drain) = domain

    (; indices, reverse_indices) = land.network
    (; x_length, y_length, area) = land.parameters

    n_cells = length(indices)

    lens = lens_input_parameter(config, "land_surface__elevation"; optional = false)
    elevation = ncread(dataset, config, lens; sel = indices, type = Float)

    # unconfined aquifer
    if do_constanthead
        constant_head = ConstantHead(dataset, config, indices)
    else
        variables = ConstantHeadVariables(; head = Float[])
        constant_head = ConstantHead(; variables, index = Int[])
    end

    connectivity = Connectivity(indices, reverse_indices, x_length, y_length)

    # cold state for groundwater head based on water table depth zi
    initial_head = elevation .- soil.variables.zi / 1000.0
    initial_head[river.network.land_indices] = elevation[river.network.land_indices]
    if do_constanthead
        initial_head[constant_head.index] = constant_head.variables.head
    end

    bottom = elevation .- soil.parameters.soilthickness ./ 1000.0
    conductance = zeros(Float, connectivity.nconnection)
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
    gwf_recharge =
        Recharge(fill(MISSING_VALUE, n_cells), zeros(n_cells), collect(1:n_cells))

    # drain boundary of unconfined aquifer (optional)
    if do_drains
        gwf_drain = Drainage(dataset, config, indices, drain.network.land_indices)
        aquifer_boundaries =
            (; recharge = gwf_recharge, river = gwf_river, drain = gwf_drain)
    else
        aquifer_boundaries = (; recharge = gwf_recharge, river = gwf_river)
    end

    subsurface_flow = GroundwaterFlow(;
        aquifer,
        connectivity,
        constanthead = constant_head,
        boundaries = aquifer_boundaries,
    )
    return subsurface_flow
end

"Initialize kinematic wave or local inertial overland flow routing"
function initialize_overland_flow(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    routing_types::NamedTuple,
)
    if routing_types.land == "kinematic-wave"
        overland_flow = KinWaveOverlandFlow(dataset, config, domain.land)
    elseif routing_types.land == "local-inertial"
        overland_flow = LocalInertialOverlandFlow(dataset, config, domain)
    else
        error(
            """An unknown "land_routing" method is specified in the TOML file
            ($(routing_types.land)). This should be "kinematic-wave" or "local-inertial".
            """,
        )
    end
    return overland_flow
end

"""
Initialize kinematic wave or local inertial overland flow routing, including optional
reservoirs and lakes.
"""
function initialize_river_flow(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    routing_types::NamedTuple,
)
    do_reservoirs = get(config.model, "reservoir__flag", false)::Bool
    if do_reservoirs
        reservoir = SimpleReservoir(dataset, config, domain.reservoir.network)
    else
        reservoir = nothing
    end

    do_lakes = get(config.model, "lake__flag", false)::Bool
    if do_lakes
        lake = Lake(dataset, config, domain.lake.network)
    else
        lake = nothing
    end

    if routing_types.river == "kinematic-wave"
        river_flow = KinWaveRiverFlow(dataset, config, domain.river, reservoir, lake)
    elseif routing_types.river == "local-inertial"
        river_flow = LocalInertialRiverFlow(dataset, config, domain.river, reservoir, lake)
    else
        error(
            """An unknown "river_routing" method is specified in the TOML file
             ($(routing_types.river)). This should be "kinematic-wave" or "local-inertial".
            """,
        )
    end
end

"Initialize `Routing` for model types `sbm` and `sbm_gwf`"
function Routing(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil::SbmSoilModel,
    routing_types::NamedTuple,
    type,
)
    subsurface_flow = initialize_subsurface_flow(dataset, config, domain, soil, type)
    overland_flow = initialize_overland_flow(dataset, config, domain, routing_types)
    river_flow = initialize_river_flow(dataset, config, domain, routing_types)

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