using PrecompileTools: @compile_workload, @setup_workload

# The CLI build overrides Wflow's default-off preference for this workload.
@setup_workload begin
    # Keep the workload next to the package so JuliaC can resolve the local checkout.
    config_dir = normpath(@__DIR__, "..", "test")
    config_paths = map(
        config_name -> joinpath(config_dir, config_name),
        ("sbm_config.toml", "sbm_gwf_config.toml", "sediment_config.toml"),
    )
    for config_path in config_paths
        isfile(config_path) || error("Wflow precompile workload not found at $config_path")
    end
    data_dir = joinpath(config_dir, "data", "input")
    isdir(data_dir) || error("Wflow precompile workload data not found at $data_dir")

    @compile_workload begin
        for config_path in config_paths
            run(config_path)
        end
    end
end
