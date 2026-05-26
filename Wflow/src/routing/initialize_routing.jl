
"Initialize subsurface flow routing for model type `sbm`"
function initialize_subsurface_flow_model(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil_model::SbmSoilModel,
    ::SbmModel;
    data_lookup::DataLookup = DataLookup(),
)
    (; parameters) = domain.land
    subsurface_flow_model =
        LateralSSFModel(dataset, config, domain, soil_model; data_lookup)

    kh_profile_type = config.model.saturated_hydraulic_conductivity_profile

    if kh_profile_type == VerticalConductivityProfile.exponential ||
       kh_profile_type == VerticalConductivityProfile.exponential_constant
        initialize_lateral_ssf_model!(
            subsurface_flow_model,
            parameters,
            subsurface_flow_model.parameters.kh_profile,
        )
    elseif kh_profile_type == VerticalConductivityProfile.layered ||
           kh_profile_type == VerticalConductivityProfile.layered_exponential
        (; kv_profile) = soil_model.parameters
        dt = Second(config.time.timestepsecs)
        initialize_lateral_ssf_model!(
            subsurface_flow_model,
            soil_model,
            parameters,
            kv_profile,
            tosecond(dt),
        )
    end
    return subsurface_flow_model
end

"Initialize subsurface flow routing for model type `sbm_gwf`"
function initialize_subsurface_flow_model(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil_model::SbmSoilModel,
    ::SbmGwfModel;
    data_lookup::DataLookup = DataLookup(),
)
    subsurface_flow_model =
        GroundwaterFlowModel(dataset, config, domain, soil_model; data_lookup)
    return subsurface_flow_model
end

"Initialize kinematic wave or local inertial overland flow routing"
function initialize_overland_flow(
    dataset::NCDataset,
    config::Config,
    domain::Domain;
    data_lookup::DataLookup = DataLookup(),
)
    (; land_routing) = config.model

    if land_routing == RoutingType.kinematic_wave
        overland_flow = KinWaveOverlandFlowModel(dataset, config, domain.land; data_lookup)
    elseif land_routing == RoutingType.local_inertial
        overland_flow = LocalInertialOverlandFlowModel(dataset, config, domain; data_lookup)
    end
    return overland_flow
end

"""
Initialize kinematic wave or local inertial overland flow routing, including optional
reservoirs.
"""
function initialize_river_flow(
    dataset::NCDataset,
    config::Config,
    domain::Domain;
    data_lookup::DataLookup = DataLookup(),
)
    (; river_routing) = config.model

    reservoir =
        config.model.reservoir__flag ?
        ReservoirModel(dataset, config, domain.reservoir.network; data_lookup) : nothing

    if river_routing == RoutingType.kinematic_wave
        river_flow =
            KinWaveRiverFlowModel(dataset, config, domain.river, reservoir; data_lookup)
    elseif river_routing == RoutingType.local_inertial
        river_flow = LocalInertialRiverFlowModel(
            dataset,
            config,
            domain.river,
            reservoir;
            data_lookup,
        )
    end
end

"Initialize `Routing` for model types `sbm` and `sbm_gwf`"
function Routing(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil_model::SbmSoilModel,
    type;
    data_lookup::DataLookup = DataLookup(),
)
    subsurface_flow = initialize_subsurface_flow_model(
        dataset,
        config,
        domain,
        soil_model,
        type;
        data_lookup,
    )
    overland_flow = initialize_overland_flow(dataset, config, domain; data_lookup)
    river_flow = initialize_river_flow(dataset, config, domain; data_lookup)

    routing = Routing(; subsurface_flow, overland_flow, river_flow)
    return routing
end

"Initialize `Routing` for model type `sediment`"
function Routing(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
    soil_loss_model::SoilLossModel;
    data_lookup::DataLookup = DataLookup(),
)
    overland_flow = OverlandFlowSedimentModel(
        dataset,
        config,
        domain.land,
        soil_loss_model;
        data_lookup,
    )
    river_flow = RiverSedimentModel(dataset, config, domain.river; data_lookup)

    routing = Routing(; overland_flow, river_flow)
    return routing
end
