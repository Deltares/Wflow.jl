"""
    initialize_sbm_gwf_model(config::Config)

Initial part of the sbm_gwf model concept. The model contains:
    - the land hydrology model with the SBM soil model
    - unconfined aquifer with groundwater flow in four directions (adjacent cells)
    - the following surface routing options:
        - 1-D kinematic wave for river flow and 1-D kinematic wave for overland flow
        - 1-D local inertial model for river flow (optional floodplain) and 1-D kinematic wave for overland flow
        - 1-D local inertial model for river flow (optional floodplain) and 2-D local inertial model for overland flow

The unconfined aquifer contains a recharge, river and a drain (optional) boundary.

The initial part reads the input settings and data as defined in the Config object.
Will return a Model that is ready to run.
"""
function Model(config::Config, type::SbmGwfModel)

    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)
    dataset = NCDataset(static_path)

    reader = prepare_reader(config)
    clock = Clock(config, reader)

    modelsettings = (;
        snow = get(config.model, "snow", false)::Bool,
        gravitational_snow_transport = get(
            config.model,
            "gravitational_snow_transport",
            false,
        )::Bool,
        glacier = get(config.model, "glacier", false)::Bool,
        reservoirs = get(config.model, "reservoirs", false)::Bool,
        lakes = get(config.model, "lakes", false)::Bool,
        drains = get(config.model, "drains", false)::Bool,
        constanthead = get(config.model, "constanthead", false)::Bool,
        water_demand = haskey(config.model, "water_demand"),
        pits = get(config.model, "pits", false)::Bool,
        min_streamorder_river = get(config.model, "min_streamorder_river", 6),
        min_streamorder_land = get(config.model, "min_streamorder_land", 5),
    )

    @info "General model settings" modelsettings[keys(modelsettings)[1:8]]...

    routing_types = get_routing_types(config)
    domain = Domain(dataset, config, modelsettings, routing_types)

    land_hydrology = LandHydrologySBM(dataset, config, domain.land)
    routing = Routing(dataset, config, domain, land_hydrology.soil, routing_types, type)

    modelmap = (land = land_hydrology, routing)
    (; maxlayers) = land_hydrology.soil.parameters
    writer = prepare_writer(
        config,
        modelmap,
        domain,
        dataset;
        extra_dim = (name = "layer", value = Float64.(1:(maxlayers))),
    )
    close(dataset)

    model =
        Model(config, domain, routing, land_hydrology, clock, reader, writer, SbmGwfModel())

    set_states!(model)

    return model
end

"update the sbm_gwf model for a single timestep"
function update!(model::AbstractModel{<:SbmGwfModel})
    (; routing, land, domain, clock, config) = model
    (; soil, runoff, demand) = land

    do_water_demand = haskey(config.model, "water_demand")
    (; aquifer, boundaries) = routing.subsurface_flow
    dt = tosecond(clock.dt)

    update!(land, routing, domain, config, dt)

    # set river stage and storage (groundwater boundary) based on river flow routing
    # variables
    for i in eachindex(boundaries.river.variables.stage)
        boundaries.river.variables.stage[i] =
            routing.river_flow.variables.h[i] + boundaries.river.parameters.bottom[i]
        boundaries.river.variables.storage[i] = routing.river_flow.variables.storage[i]
    end

    # determine stable time step for groundwater flow
    conductivity_profile = get(config.model, "conductivity_profile", "uniform")
    dt_gw = stable_timestep(aquifer, conductivity_profile) # time step in day (Float64)
    dt_sbm = (dt / tosecond(BASETIMESTEP)) # dt is in seconds (Float64)
    if dt_gw < dt_sbm
        @warn(
            "stable time step dt $dt_gw for groundwater flow is smaller than `LandHydrologySBM` model dt $dt_sbm"
        )
    end

    Q = zeros(routing.subsurface_flow.connectivity.ncell)
    # exchange of recharge between SBM soil model and groundwater flow domain
    # recharge rate groundwater is required in units [m d⁻¹]
    @. boundaries.recharge.variables.rate =
        soil.variables.recharge / 1000.0 * (1.0 / dt_sbm)
    if do_water_demand
        @. boundaries.recharge.variables.rate -=
            land.allocation.variables.act_groundwater_abst / 1000.0 * (1.0 / dt_sbm)
    end
    # update groundwater domain
    update!(routing.subsurface_flow, Q, dt_sbm, conductivity_profile)

    # update SBM soil model (runoff, ustorelayerdepth and satwaterdepth)
    update!(soil, (; runoff, demand, subsurface_flow = routing.subsurface_flow))

    surface_routing!(model)

    return nothing
end
