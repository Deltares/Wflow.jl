
function initialize_subsurface_flow(
    dataset::NCDataset,
    config::Config,
    soil::SbmSoilModel,
    network::Network,
    parameters::LandParameters,
    type::SbmModel,
)
    do_subsurface_flow = get(config.model, "kinematic-wave_subsurface", true)::Bool
    if do_subsurface_flow
        subsurface_flow = LateralSSF(dataset, config, network.land, soil, parameters)

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
        n = length(network.land.indices)
        subsurface_flow = GroundwaterExchange(n)
    end
    return subsurface_flow
end

function initialize_subsurface_flow(
    dataset::NCDataset,
    config::Config,
    soil::SbmSoilModel,
    network::Network,
    parameters::LandParameters,
    type::SbmGwfModel,
)
    do_drains = get(config.model, "drains", false)::Bool
    do_constanthead = get(config.model, "constanthead", false)::Bool

    (; indices, reverse_indices) = network.land
    (; x_length, y_length, area) = parameters

    n_cells = length(indices)

    lens = lens_input_parameter(config, "land_surface__elevation"; optional = false)
    elevation = ncread(dataset, config, lens; sel = indices, type = Float64)

    # unconfined aquifer
    if do_constanthead
        constant_head = ConstantHead(dataset, config, indices)
    else
        variables = ConstantHeadVariables(; head = Float64[])
        constant_head = ConstantHead(; variables, index = Int64[])
    end

    connectivity = Connectivity(indices, reverse_indices, x_length, y_length)

    # cold state for groundwater head based on water table depth zi
    initial_head = elevation .- soil.variables.zi / 1000.0
    initial_head[network.river.land_indices] = elevation[network.river.land_indices]
    if do_constanthead
        initial_head[constant_head.index] = constant_head.variables.head
    end

    bottom = elevation .- soil.parameters.soilthickness ./ 1000.0
    conductance = zeros(Float64, connectivity.nconnection)
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
    river = GwfRiver(dataset, config, network.river.indices, network.river.land_indices)

    # recharge boundary of unconfined aquifer
    recharge = Recharge(fill(MISSING_VALUE, n_cells), zeros(n_cells), collect(1:n_cells))

    # drain boundary of unconfined aquifer (optional)
    if do_drains
        drain = Drainage(dataset, config, indices, network.drain.land_indices)
        aquifer_boundaries = (; recharge, river, drain)
    else
        aquifer_boundaries = (; recharge, river)
    end

    subsurface_flow = GroundwaterFlow(;
        aquifer,
        connectivity,
        constanthead = constant_head,
        boundaries = aquifer_boundaries,
    )
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
    type,
)
    subsurface_flow =
        initialize_subsurface_flow(dataset, config, soil, network, parameters.land, type)
    overland_flow =
        initialize_overland_flow(dataset, config, network, parameters, routing_types)
    river_flow =
        initialize_river_flow(dataset, config, network, parameters.river, routing_types)

    routing = Routing(; subsurface_flow, overland_flow, river_flow)
    return routing
end