using NCDatasets
using Dates
using CFTime
using Random
using UnPack

tomlpath = joinpath(@__DIR__, "config.toml")
parsed_toml = Wflow.parsefile(tomlpath)
config = Wflow.Config(tomlpath)

@testset "configuration file" begin
    @test parsed_toml isa Dict{String,Any}
    @test config isa Wflow.Config
    @test Dict(config) == parsed_toml
    @test pathof(config) == tomlpath
    @test dirname(config) == dirname(tomlpath)

    # test if the values are parsed as expected
    @test config.casename == "testcase"
    @test config.starttime === DateTime(2000)
    @test config.endtime === DateTime(2000, 2)
    @test config.output.path == "data/output_moselle.nc"
    @test config.output isa Wflow.Config
    @test collect(keys(config.output)) == ["lateral", "vertical", "path"]
end

@testset "checkdims" begin
    @test_throws AssertionError Wflow.checkdims(("z", "lat", "time"))
    @test_throws AssertionError Wflow.checkdims(("z", "lat"))
    @test Wflow.checkdims(("lon", "lat", "time")) == ("lon", "lat", "time")
    @test Wflow.checkdims(("time", "lon", "lat")) == ("time", "lon", "lat")
    @test Wflow.checkdims(("lat", "lon", "time")) == ("lat", "lon", "time")
    @test Wflow.checkdims(("time", "lat", "lon")) == ("time", "lat", "lon")
end

@testset "timecycles" begin
    @test Wflow.timecycles([Date(2020, 4, 21), Date(2020, 10, 21)]) == [(4, 21), (10, 21)]
    @test_throws ErrorException Wflow.timecycles([Date(2020, 4, 21), Date(2021, 10, 21)])
    @test_throws ErrorException Wflow.timecycles(collect(1:400))
    @test Wflow.timecycles(collect(1:12)) == collect(zip(1:12, fill(1, 12)))
    @test Wflow.timecycles(collect(1:366)) ==
          monthday.(Date(2000, 1, 1):Day(1):Date(2000, 12, 31))
end

# test reading and setting of warm states (reinit=true)
# modify existing config an re-initialize the model
@test !config.model.reinit
config["model"]["reinit"] = true
@test config.model.reinit
model = Wflow.initialize_sbm_model(config)

@unpack vertical, clock, reader, writer = model

@testset "output and state names" begin
    ncdims = ("lon", "lat", "layer", "time")
    @test dimnames(writer.dataset["ustorelayerdepth"]) == ncdims
    ncvars = [k for k in keys(writer.dataset) if !in(k, ncdims)]
    @test "snow" in ncvars
    @test "q" in ncvars
    @test writer.states == (
        Wflow.symbols"vertical.satwaterdepth",
        Wflow.symbols"vertical.snow",
        Wflow.symbols"vertical.tsoil",
        Wflow.symbols"vertical.ustorelayerdepth",
        Wflow.symbols"vertical.snowwater",
        Wflow.symbols"vertical.canopystorage",
        Wflow.symbols"lateral.river.reservoir.volume",
        Wflow.symbols"lateral.subsurface.ssf",
        Wflow.symbols"lateral.river.q",
        Wflow.symbols"lateral.river.h",
        Wflow.symbols"lateral.land.q",
        Wflow.symbols"lateral.land.h",
    )
end

@testset "warm states" begin
    @test Wflow.param(model, "lateral.river.reservoir.volume")[2] ≈ 2.7801042e7
    @test Wflow.param(model, "vertical.satwaterdepth")[26625] ≈ 168.40777587890625
    @test Wflow.param(model, "vertical.snow")[26625] ≈ 0.9690762758255005
    @test Wflow.param(model, "vertical.tsoil")[26625] ≈ 2.9367642402648926
    @test Wflow.param(model, "vertical.ustorelayerdepth")[1][1] ≈ 3.31813645362854
    @test Wflow.param(model, "vertical.snowwater")[26625] ≈ 0.011469922959804535
    @test Wflow.param(model, "vertical.canopystorage")[1] ≈ 0.0
    @test Wflow.param(model, "lateral.subsurface.ssf")[39308] ≈ 2.49775357952e11
    @test Wflow.param(model, "lateral.river.q")[4061] ≈ 1.3987680673599243
    @test Wflow.param(model, "lateral.river.h")[4061] ≈ 0.628973126411438
    @test Wflow.param(model, "lateral.land.q")[39308] ≈ 0.08347763121128082
    @test Wflow.param(model, "lateral.land.h")[39308] ≈ 0.008668862283229828
end

@testset "reducer" begin
    V = [6, 5, 4, 1]
    @test Wflow.reducerfunction("maximum")(V) == 6
    @test Wflow.reducerfunction("mean")(V) == 4
    @test Wflow.reducerfunction("median")(V) == 4.5
    @test Wflow.reducerfunction("first")(V) == 6
    @test Wflow.reducerfunction("last")(V) == 1
    @test_throws ErrorException Wflow.reducerfunction("other")
end

@testset "network" begin
    @unpack network = model
    @unpack indices, reverse_indices = model.network.land
    # test if the reverse index reverses the index
    linear_index = 100
    cartesian_index = indices[linear_index]
    @test cartesian_index === CartesianIndex(115, 6)
    @test reverse_indices[cartesian_index] === linear_index
end

Wflow.close_files(model, delete_output = true)
