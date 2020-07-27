using NCDatasets
using Dates
using CFTime
using Random
using UnPack

@testset "configuration file" begin
    tomlpath = joinpath(@__DIR__, "config.toml")
    parsed_toml = Wflow.parsefile(tomlpath)
    config = Wflow.Config(tomlpath)
    @test parsed_toml isa Dict{String,Any}
    @test config isa Wflow.Config
    @test Dict(config) == parsed_toml
    @test pathof(config) == tomlpath
    @test dirname(config) == dirname(tomlpath)

    # test if the values are parsed as expected
    @test config.casename == "testcase"
    @test config.λ == 1.2
    @test config.starttime === DateTime(2000)
    @test config.endtime === DateTime(2000, 2)
    @test config.output.path == "data/specified_output.nc"
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

tomlpath = joinpath(@__DIR__, "config.toml")
tomldir = dirname(tomlpath)
config = Wflow.Config(tomlpath)

# initialize a vector of SBM structs
model = Wflow.initialize_sbm_model(config)

@unpack vertical, clock, reader, writer = model
@unpack dataset, buffer, inds = reader

@testset "output" begin
    ncdims = ("lon", "lat", "layer", "time")
    @test dimnames(writer.dataset["ustorelayerdepth"]) == ncdims
    ncvars = [k for k in keys(writer.dataset) if !in(k, ncdims)]
    @test "snow" in ncvars
    @test "q" in ncvars
end

# get the output path before it's closed, and remove up the file
output_path = path(model.writer.dataset)
Wflow.close_files(model)
rm(output_path)
