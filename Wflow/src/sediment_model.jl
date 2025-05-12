"""
    Model(config::Config, type::SedimentModel)

Initial part of the sediment model concept. Reads the input settings and data as defined in
the Config object. Will return a Model that is ready to run.
"""
function Model(config::Config, type::SedimentModel)
    model_type = config.model.type::String
    @info "Initialize model variables for model type `$model_type`."

    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)
    dataset = NCDataset(static_path)

    reader = prepare_reader(config)
    clock = Clock(config, reader)

    modelsettings = (;
        reservoirs = get(config.model, "reservoir__flag", false)::Bool,
        lakes = get(config.model, "lake__flag", false)::Bool,
        pits = get(config.model, "pit__flag", false)::Bool,
    )

    @info "General model settings" modelsettings[keys(modelsettings)[1:2]]...

    domain = Domain(dataset, config, modelsettings)
    soilloss = SoilLoss(dataset, config, domain.land.network.indices)
    routing = Routing(dataset, config, domain, soilloss)

    modelmap = (land = soilloss, routing)
    writer = prepare_writer(config, modelmap, domain, dataset)
    close(dataset)

    model = Model(config, domain, routing, soilloss, clock, reader, writer, SedimentModel())

    set_states!(model)
    @info "Initialized model"

    return model
end

"update `sediment` model for a single timestep"
function update!(model::AbstractModel{<:SedimentModel})
    (; routing, land, domain, config, clock) = model
    dt = Float(tosecond(clock.dt))

    # Soil erosion
    update!(land, domain.land.parameters, dt)

    # Overland flow sediment transport
    update!(routing.overland_flow, land.soil_erosion, domain.land, dt)

    # River sediment transport
    do_river = get(config.model, "run_river_model__flag", false)::Bool
    if do_river
        update!(routing.river_flow, routing.overland_flow.to_river, domain.river, dt)
    end

    return nothing
end

"set the initial states of the `sediment` model"
function set_states!(model::AbstractModel{<:SedimentModel})
    # read and set states in model object if cold_start=false
    (; config) = model
    cold_start = get(config.model, "cold_start__flag", true)::Bool
    if cold_start == false
        instate_path = input_path(config, config.state.path_input)
        @info "Set initial conditions from state file `$instate_path`."
        set_states!(instate_path, model; type = Float)
    else
        @info "Set initial conditions from default values."
    end
    return nothing
end
