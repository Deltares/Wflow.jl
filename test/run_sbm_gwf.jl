# simple example of running a model simulation, not part of the tests
using Wflow

tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)
#config["model"]["reinit"] = true
model = Wflow.initialize_sbm_gwf_model(config)
Wflow.run_simulation(model; close_files = true)