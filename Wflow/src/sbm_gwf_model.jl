"""
    Model(config::Config, type::SbmGwfModel)

Initial part of the `sbm_gwf` model concept. The model contains the land hydrology model
with the `SBM` soil model, an unconfined aquifer with groundwater flow in four directions
(adjacent cells). The unconfined aquifer contains a recharge, river and a drain (optional)
boundary. The initial part reads the input settings and data as defined in the Config
object. Will return a Model that is ready to run.
"""
function Model(config::Config, type::SbmGwfModel)

    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)
    dataset = NCDataset(static_path)

    reader = NCReader(config)
    clock = Clock(config, reader)

    @info "General model settings." (;
        snow = config.model.snow__flag,
        gravitational_snow_transport = config.model.snow_gravitational_transport__flag,
        glacier = config.model.glacier__flag,
        reservoirs = config.model.reservoir__flag,
        drains = config.model.drain__flag,
        constanthead = config.model.constanthead__flag,
        water_demand = do_water_demand(config),
    )...

    domain = Domain(dataset, config, type)

    land_hydrology = LandHydrologySBM(dataset, config, domain.land)
    routing = Routing(dataset, config, domain, land_hydrology.soil, type)
    mass_balance = HydrologicalMassBalance(domain, config)

    modelmap = (land = land_hydrology, routing, mass_balance)
    (; maxlayers) = land_hydrology.soil.parameters
    writer = Writer(
        config,
        modelmap,
        domain,
        dataset;
        extra_dim = (name = "layer", value = Float64.(1:(maxlayers))),
    )
    close(dataset)

    model = Model(
        config,
        domain,
        routing,
        land_hydrology,
        mass_balance,
        clock,
        reader,
        writer,
        SbmGwfModel(),
    )

    set_states!(model)

    return model
end

"update the `sbm_gwf` model type for a single timestep"
function update!(model::AbstractModel{<:SbmGwfModel})
    (; routing, land, domain, clock, config) = model
    (; soil, runoff, demand) = land
    (; boundaries) = routing.subsurface_flow

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
    dt_gwf = (dt / tosecond(BASETIMESTEP)) # dt is in seconds (Float64)

    # exchange of recharge between SBM soil model and groundwater flow domain
    # recharge rate groundwater is required in units [m d⁻¹]
    @. boundaries.recharge.variables.rate =
        soil.variables.recharge / 1000.0 * (1.0 / dt_gwf)
    if do_water_demand(config)
        @. boundaries.recharge.variables.rate -=
            land.allocation.variables.act_groundwater_abst / 1000.0 * (1.0 / dt_gwf)
    end

    # update groundwater domain
    update!(routing.subsurface_flow, dt_gwf, config.model.conductivity_profile)

    # update SBM soil model (runoff, ustorelayerdepth and satwaterdepth)
    update!(soil, (; runoff, demand, subsurface_flow = routing.subsurface_flow))

    surface_routing!(model)

    return nothing
end
