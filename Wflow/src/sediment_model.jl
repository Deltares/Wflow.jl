"""
    Model(config::Config, type::SedimentModel)

Initial part of the sediment model concept. Reads the input settings and data as defined in
the Config object. Will return a Model that is ready to run.
"""
function Model(config::Config, type::SedimentModel)
    model_type = config.model.type
    @info "Initialize model variables for model type `$model_type`."

    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)
    dataset = NCDataset(static_path)

    reader = NCReader(config)
    clock = Clock(config, reader)

    @info "General model settings." reservoirs = config.model.reservoir__flag

    domain = Domain(dataset, config, type)
    soilloss = SoilLoss(dataset, config, domain.land.network.indices)
    routing = Routing(dataset, config, domain, soilloss)
    mass_balance = NoMassBalance()

    modelmap = (land = soilloss, routing, mass_balance)
    writer = Writer(config, modelmap, domain, dataset)
    close(dataset)

    model = Model(
        config,
        domain,
        routing,
        soilloss,
        mass_balance,
        clock,
        reader,
        writer,
        SedimentModel(),
    )

    set_states!(model)
    @info "Initialized model"

    return model
end

"update `sediment` model for a single timestep"
function update!(model::AbstractModel{<:SedimentModel})
    (; routing, land, domain, config, clock) = model
    dt = tosecond(clock.dt)

    # Soil erosion
    update!(land, domain.land.parameters, dt)

    # Overland flow sediment transport
    update!(routing.overland_flow, land.soil_erosion, domain.land, dt)

    # River sediment transport
    if config.model.run_river_model__flag
        update!(routing.river_flow, routing.overland_flow.to_river, domain.river, dt)
    end

    return nothing
end

"set the initial states of the `sediment` model"
function set_states!(model::AbstractModel{<:SedimentModel})
    # read and set states in model object if cold_start=false
    (; config) = model
    if !config.model.cold_start__flag
        instate_path = input_path(config, config.state.path_input)
        @info "Set initial conditions from state file `$instate_path`."
        set_states!(instate_path, model; type = Float64)
    else
        @info "Set initial conditions from default values."
    end
    return nothing
end
