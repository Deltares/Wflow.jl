
function initialize_lateral_ssf(
    dataset::NCDataset,
    config::Config,
    soil::SbmSoilModel,
    network::NetworkLand,
    parameters::LandParameters,
)
    do_subsurface_flow = get(config.model, "kinematic-wave_subsurface", true)::Bool
    if do_subsurface_flow
        subsurface_flow = LateralSSF(dataset, config, network, soil, parameters)

        kh_profile_type = get(
            config.model,
            "saturated_hydraulic_conductivity_profile",
            "exponential",
        )::String

        if kh_profile_type == "exponential" || kh_profile_type == "exponential_constant"
            initialize_lateral_ssf!(subsurface_flow, subsurface_flow.parameters.kh_profile)
        elseif kh_profile_type == "layered" || kh_profile_type == "layered_exponential"
            (; kv_profile) = soil.parameters
            dt = Second(config.time.timestepsecs)
            initialize_lateral_ssf!(subsurface_flow, soil, kv_profile, tosecond(dt))
        end
    else
        # when the SBM model is coupled (BMI) to a groundwater model, the following
        # variables are expected to be exchanged from the groundwater model.
        n = length(network.indices)
        subsurface_flow = GroundwaterExchange(n)
    end
    return subsurface_flow
end

function initialize_overland_flow(
    dataset::NCDataset,
    config::Config,
    network::Network,
    parameters::SharedParameters,
    routing_types::NamedTuple,
)
    if routing_types.land == "kinematic-wave"
        overland_flow = KinWaveOverlandFlow(dataset, config, network.land, parameters.land)
    elseif routing_types.land == "local-inertial"
        overland_flow = LocalInertialOverlandFlow(dataset, config, network, parameters)
    else
        error(
            """An unknown "land_routing" method is specified in the TOML file
            ($(routing_types.land)). This should be "kinematic-wave" or "local-inertial".
            """,
        )
    end
    return overland_flow
end

function initialize_river_flow(
    dataset::NCDataset,
    config::Config,
    network::Network,
    parameters::RiverParameters,
    routing_types::NamedTuple,
)
    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    if do_reservoirs
        reservoir = SimpleReservoir(dataset, config, network.reservoir)
    else
        reservoir = nothing
    end

    do_lakes = get(config.model, "lakes", false)::Bool
    if do_lakes
        lake = Lake(dataset, config, network.lake)
    else
        lake = nothing
    end

    if routing_types.river == "kinematic-wave"
        river_flow =
            KinWaveRiverFlow(dataset, config, network.river, parameters, reservoir, lake)
    elseif routing_types.river == "local-inertial"
        river_flow = LocalInertialRiverFlow(
            dataset,
            config,
            network.river,
            parameters,
            reservoir,
            lake,
        )
    else
        error("""An unknown "river_routing" method is specified in the TOML file
              ($river_routing). This should be "kinematic-wave" or "local-inertial".
              """)
    end
end

function Routing(
    dataset::NCDataset,
    config::Config,
    soil::SbmSoilModel,
    network::Network,
    parameters::SharedParameters,
    routing_types::NamedTuple,
)
    subsurface_flow =
        initialize_lateral_ssf(dataset, config, soil, network.land, parameters.land)
    overland_flow =
        initialize_overland_flow(dataset, config, network, parameters, routing_types)
    river_flow =
        initialize_river_flow(dataset, config, network, parameters.river, routing_types)

    routing = Routing(; subsurface_flow, overland_flow, river_flow)
    return routing
end