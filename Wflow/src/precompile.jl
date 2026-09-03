using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    configured_dir = get(ENV, "WFLOW_PRECOMPILE_CONFIG_DIR", nothing)
    config_dir = if isnothing(configured_dir)
        normpath(@__DIR__, "..", "test")
    else
        configured_dir
    end
    config_paths = map(
        config_name -> joinpath(config_dir, config_name),
        ("sbm_config.toml", "sbm_gwf_config.toml", "sediment_config.toml"),
    )
    missing_paths = filter(!isfile, config_paths)
    data_dir = joinpath(config_dir, "data", "input")

    if !isnothing(configured_dir)
        isempty(missing_paths) ||
            error("Wflow precompile workload not found at $(join(missing_paths, ", "))")
        isdir(data_dir) || error("Wflow precompile workload data not found at $data_dir")
    elseif !isempty(missing_paths) || !isdir(data_dir)
        return
    end

    @compile_workload begin
        for config_path in config_paths
            run(config_path)
        end
    end
end
