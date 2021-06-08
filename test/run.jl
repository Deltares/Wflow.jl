using Wflow
using Dates

tomlpath = joinpath(@__DIR__, "sbm_simple.toml")
Wflow.run(tomlpath)

# test whether restarted runs get the same results as continuus ones, i.e. state is captured
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

# January at once, cold start
config.starttime = DateTime("2000-01-01T00:00:00")
config.endtime   = DateTime("2000-02-01T00:00:00")
config.model.reinit = true  # cold start
# note that this needs to be relative to the tomlpath
config.state.path_output = joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january.nc")
config.output.path = joinpath(dirname(tomlpath), "data/state-test/output-moselle-january.nc")
model = Wflow.initialize_sbm_model(config)
Wflow.run(model)

# first half of January, cold start
config.starttime = DateTime("2000-01-01T00:00:00")
config.endtime   = DateTime("2000-01-15T00:00:00")
config.model.reinit = true  # cold start
config.state.path_output = joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-1of2.nc")
config.output.path = joinpath(dirname(tomlpath), "data/state-test/output-moselle-january-1of2.nc")
model = Wflow.initialize_sbm_model(config)
Wflow.run(model)

# second half of January, warm start
config.starttime = DateTime("2000-01-15T00:00:00")
config.endtime   = DateTime("2000-02-01T00:00:00")
config.model.reinit = false  # cold start
config.state.path_input = joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-1of2.nc")
config.state.path_output = joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-2of2.nc")
config.output.path = joinpath(dirname(tomlpath), "data/state-test/output-moselle-january-2of2.nc")
model = Wflow.initialize_sbm_model(config)
Wflow.run(model)
