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
    end
end

function BMI.get_output_var_names(model::Model)
    BMI.get_input_var_names(model)
end
