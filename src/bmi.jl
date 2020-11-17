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
