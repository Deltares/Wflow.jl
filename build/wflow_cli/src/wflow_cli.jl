module wflow_cli

using PrecompileTools: @compile_workload, @setup_workload
using Wflow

@setup_workload begin
    # JuliaC copies this package, while Wflow still resolves to the local checkout.
    config_dir = joinpath(pkgdir(Wflow), "test")
    config_paths = map(
        config_name -> joinpath(config_dir, config_name),
        ("sbm_config.toml", "sbm_gwf_config.toml", "sediment_config.toml"),
    )
    for config_path in config_paths
        isfile(config_path) || error("Wflow precompile workload not found at $config_path")
    end
    data_dir = joinpath(config_dir, "data", "input")
    isdir(data_dir) || error("Wflow precompile workload data not found at $data_dir")

    # Use a workload, Wflow.run on its own will not cache the precompilation of the package.
    @compile_workload begin
        for config_path in config_paths
            Wflow.run(config_path)
        end
    end
end

function help(x)::Cint
    println(x)
    println("Usage: wflow_cli path/to/config.toml")
    return 1
end

function (@main)(_)::Cint
    n = length(ARGS)
    if n != 1
        return help("Only 1 argument expected, got $n")
    end
    toml_path = only(ARGS)
    if !isfile(toml_path)
        return help("File not found: $toml_path")
    end

    try
        Wflow.run(toml_path)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end

    return 0
end

end # module
