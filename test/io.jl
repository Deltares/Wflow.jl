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
    @test config.output.parameters isa Wflow.Config
    @test collect(keys(config.output.parameters)) == [
        "snow",
        "soilthickness",
        "snowwater",
        "satwaterdepth",
        "q",
        "ustorelayerdepth",
        "canopystorage",
        "tsoil",
    ]
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
    @test writer.statenames == (
        "satwaterdepth",
        "snow",
        "tsoil",
        "ustorelayerdepth",
        "snowwater",
        "canopystorage",
        "volume_reservoir",
        "ssf",
        "q_river",
        "h_river",
        "q_land",
        "h_land",
    )
end

@testset "warm states" begin
    @test get(model, Wflow.paramap["volume_reservoir"])[2] ≈ 2.7801042e7
    @test get(model, Wflow.paramap["satwaterdepth"])[26625] ≈ 168.40777587890625
    @test get(model, Wflow.paramap["snow"])[26625] ≈ 0.9690762758255005
    @test get(model, Wflow.paramap["tsoil"])[26625] ≈ 2.9367642402648926
    @test get(model, Wflow.paramap["ustorelayerdepth"])[1][1] ≈ 3.31813645362854
    @test get(model, Wflow.paramap["snowwater"])[26625] ≈ 0.011469922959804535
    @test get(model, Wflow.paramap["canopystorage"])[1] ≈ 0.0
    @test get(model, Wflow.paramap["ssf"])[39308] ≈ 2.49775357952e11
    @test get(model, Wflow.paramap["q_river"])[4061] ≈ 1.3987680673599243
    @test get(model, Wflow.paramap["h_river"])[4061] ≈ 0.628973126411438
    @test get(model, Wflow.paramap["q_land"])[39308] ≈ 0.08347763121128082
    @test get(model, Wflow.paramap["h_land"])[39308] ≈ 0.008668862283229828
end

@testset "reducer" begin
    @test Wflow.reducer(Dict("reducer" => "maximum"))([6, 5, 2]) == 6
    @test Wflow.reducer(Dict("reducer" => "mean"))([6, 5, 4, 1]) == 4
    @test Wflow.reducer(Dict("reducer" => "median"))([6, 5, 4, 1]) == 4.5
    @test Wflow.reducer(Dict("reducer" => "first"))([6, 5, 4, 1]) == 6
    @test Wflow.reducer(Dict("reducer" => "last"))([6, 5, 4, 1]) == 1
    @test Wflow.reducer(Dict("index" => 2))([6, 5, 4, 1]) == 5
end

Wflow.close_files(model, delete_output = true)
