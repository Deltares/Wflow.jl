@testitem "Run" begin
    using Test
    using Wflow
    using Dates
    using NCDatasets

    # test whether restarted runs get the same results as continuous ones, i.e. state is captured
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)

    # January at once, cold start
    config.time.starttime = DateTime("2000-01-01T00:00:00")
    config.time.endtime = DateTime("2000-02-01T00:00:00")
    config.model.cold_start__flag = true  # cold start
    # note that this needs to be relative to the tomlpath
    config.state.path_output =
        joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january.nc")
    model = Wflow.Model(config)
    Wflow.run!(model)

    # first half of January, cold start
    config.time.starttime = DateTime("2000-01-01T00:00:00")
    config.time.endtime = DateTime("2000-01-15T00:00:00")
    config.model.cold_start__flag = true  # cold start
    config.state.path_output =
        joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-1of2.nc")
    model = Wflow.Model(config)
    Wflow.run!(model)

    # second half of January, warm start
    config.time.starttime = DateTime("2000-01-15T00:00:00")
    config.time.endtime = DateTime("2000-02-01T00:00:00")
    config.model.cold_start__flag = false  # warm start
    config.state.path_input =
        joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-1of2.nc")
    config.state.path_output =
        joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-2of2.nc")
    model = Wflow.Model(config)
    Wflow.run!(model)

    # second half of January, warm start, fews_run set to true, and starttime set one day earlier
    # to match endtime of part 1
    config.time.starttime = DateTime("2000-01-14T00:00:00")
    config.time.endtime = DateTime("2000-02-01T00:00:00")
    config.model.cold_start__flag = false  # warm start
    config.fews_run__flag = true
    config.state.path_input =
        joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-1of2.nc")
    config.state.path_output = joinpath(
        dirname(tomlpath),
        "data/state-test/outstates-moselle-january-2of2-fews_run.nc",
    )
    model = Wflow.Model(config)
    Wflow.run!(model)

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
            tol = varname == "volume_reservoir" ? 1e-8 : 1e-9
            @test maxdiff < tol
        end
    end

    # the fews_run restart should match the other restart exactly
    endstate_fewsrun_path = joinpath(
        dirname(tomlpath),
        "data/state-test/outstates-moselle-january-2of2-fews_run.nc",
    )
    endstate_restart_path =
        joinpath(dirname(tomlpath), "data/state-test/outstates-moselle-january-2of2.nc")
    endstate_fewsrun = NCDataset(endstate_fewsrun_path)
    endstate_restart = NCDataset(endstate_restart_path)

    varnames = setdiff(keys(endstate_restart), keys(endstate_restart.dim))
    @testset "fews_run" begin
        for varname in varnames
            a = endstate_fewsrun[varname][:]
            b = endstate_restart[varname][:]
            @test all(skipmissing(a) .== skipmissing(b))
        end
    end
end
