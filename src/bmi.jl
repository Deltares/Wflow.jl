function BMI.initialize(config_file)
    config = Config(config_file)
    model = if modeltype == "sbm"
        initialize_sbm_model(config)
    elseif modeltype == "sbm_gwf"
        initialize_sbm_gwf_model(config)
    elseif modeltype == "hbv"
        initialize_hbv_model(config)
    else
        error("unknown model type")
    end
    return model
end

function BMI.update(model::Model)
    @unpack network, config = model
    update_func = config.model.type == "sbm_gwf" ? update_sbm_gwf : update
    return update_func(model)
end

function BMI.update_until(model::Model, time::Float64)
    @unpack network, config = model
    update_func = config.model.type == "sbm_gwf" ? update_sbm_gwf : update
    n_iter = Int(time)
    for i in 1:n_iter
        update_func(model)
    end
    return model
end

function BMI.finalize(model::Model)
    @unpack config, writer, clock = model
    rewind!(clock)
    write_netcdf_timestep(model, writer.state_dataset, writer.state_parameters)
    reset_clock!(model.clock, config)
    close_files(model, delete_output = false)
end

function BMI.get_component_name(model::Model)
    @unpack config = model
    return config.model.type
end

function BMI.get_input_item_count(model::Model)
    length(BMI.get_input_var_names(model))
end

function BMI.get_output_item_count(model::Model)
    length(BMI.get_output_var_names(model))
end

function BMI.get_input_var_names(model::Model)
    @unpack config = model
    if haskey(config,"API")
        var_names = Vector{String}()
        for c in config.API.components
            append!(var_names, collect(string.(c,".", fieldnames(typeof(param(model, c))))))
        end
        return var_names
    else
        @warn("TOML file does not contain section [API] to extract model var names")
        return []
    end
end

function BMI.get_output_var_names(model::Model)
    BMI.get_input_var_names(model)
end

function BMI.get_var_grid(model::Model, name::String)
    if name in BMI.get_input_var_names(model)
        return 1
    end
    return -1
end

function BMI.get_var_type(model::Model, name::String)
    string(typeof(param(model,name)))
end

function BMI.get_var_units(model::Model, name::String)
    #TODO implement properly
    "mm"
end

function BMI.get_var_itemsize(model::Model, name::String)
    sizeof(Wflow.param(model, name)[1])
end

function BMI.get_var_nbytes(model::Model, name::String)
    sizeof(Wflow.param(model, name))
end

function BMI.get_var_location(model::Model, name::String)
    if name in BMI.get_input_var_names(model)
        return "node"
    end
end

function BMI.get_current_time(model::Model)
    datetime2unix(model.clock.time)
end

function BMI.get_start_time(model::Model)
    datetime2unix(model.clock.time)
end

function BMI.get_start_time(model::Model)
    datetime2unix(model.config.starttime)
end

function BMI.get_end_time(model::Model)
    datetime2unix(model.config.endtime)
end

function BMI.get_time_units()
    string("seconds since ", unix2datetime(0))
end

function BMI.get_time_step(model::Model)
    Float64(model.config.timestepsecs)
end