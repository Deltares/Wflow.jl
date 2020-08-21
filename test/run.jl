# simple example of running a model simulation, not part of the tests
using Wflow

tomlpath = joinpath(@__DIR__, "config.toml")
config = Wflow.Config(tomlpath)
config["model"]["reinit"] = true

model = Wflow.initialize_sbm_model(config)
@time Wflow.run_simulation(model; close_files = true)
# Wflow.close_files(model)
