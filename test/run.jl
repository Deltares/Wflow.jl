using Test
using Wflow
using Dates
using NCDatasets

tomlpath = joinpath(@__DIR__, "sbm_simple.toml")
Wflow.run(tomlpath; silent = true)

# test whether restarted runs get the same results as continuous ones, i.e. state is captured
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

# January at once, cold start
config.starttime = DateTime("2000-01-01T00:00:00")
config.endtime = DateTime("2000-02-01T00:00:00")
config.model.reinit = true  # cold start
# note that this needs to be relative to the tomlpath
config.state.path_output =
    joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january.nc")
config.output.path =
    joinpath(dirname(tomlpath), "data/state-test/output-moselle-january.nc")
model = Wflow.initialize_sbm_model(config)
Wflow.run(model)

# first half of January, cold start
config.starttime = DateTime("2000-01-01T00:00:00")
config.endtime = DateTime("2000-01-15T00:00:00")
config.model.reinit = true  # cold start
config.state.path_output =
    joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-1of2.nc")
config.output.path =
    joinpath(dirname(tomlpath), "data/state-test/output-moselle-january-1of2.nc")
model = Wflow.initialize_sbm_model(config)
Wflow.run(model)

# second half of January, warm start
config.starttime = DateTime("2000-01-15T00:00:00")
config.endtime = DateTime("2000-02-01T00:00:00")
config.model.reinit = false  # warm start
config.state.path_input =
    joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-1of2.nc")
config.state.path_output =
    joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-2of2.nc")
config.output.path =
    joinpath(dirname(tomlpath), "data/state-test/output-moselle-january-2of2.nc")
model = Wflow.initialize_sbm_model(config)
Wflow.run(model)

# verify that there are minimal differences in the end state of the two runs
endstate_one_run_path =
    joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january.nc")
endstate_restart_path =
    joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-2of2.nc")
endstate_one_run = NCDataset(endstate_one_run_path)
endstate_restart = NCDataset(endstate_restart_path)

varnames = setdiff(keys(endstate_restart), keys(endstate_restart.dim))
@testset "restart" begin
    @test length(varnames) > 10
    for varname in varnames
        a = endstate_one_run[varname][:]
        b = endstate_restart[varname][:]
        maxdiff = maximum(abs.(skipmissing(b - a)))
        @test maxdiff < 1e-9
    end
end
