
"Initialize subsurface flow routing for model type `sbm`"
function initialize_subsurface_flow(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil::SbmSoilModel,
    ::SbmModel,
)
    (; parameters) = domain.land
    subsurface_flow = LateralSSF(dataset, config, domain, soil)

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
    subsurface_flow = GroundwaterFlow(dataset, config, domain, soil)
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
