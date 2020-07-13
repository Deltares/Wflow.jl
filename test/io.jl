using NCDatasets
using Dates
using CFTime
using Random
using UnPack

@testset "configuration file" begin
    tomlpath = joinpath(@__DIR__, "config.toml")
    parsed_toml = Wflow.parsefile(tomlpath)
    config = Wflow.Config(parsed_toml)
    @test parsed_toml isa Dict{String,Any}
    @test config isa Wflow.Config
    @test getfield(config, :dict) === parsed_toml

    # test if the values are parsed as expected
    @test config.casename == "testcase"
    @test config.λ == 1.2
    @test config.input.starttime === DateTime(2000)
    @test config.input.endtime === DateTime(2000, 2)
    @test config.output.path == "data/specified_output.nc"
    @test config.output.parameters isa Vector
    @test config.output.parameters == [
        "satwaterdepth",
        "snow",
        "tsoil",
        "ustorelayerdepth",
        "snowwater",
        "canopystorage",
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
config = Wflow.Config(Wflow.parsefile(tomlpath))

# initialize a vector of SBM structs
model = Wflow.initialize_sbm_model(
    config,
    joinpath(tomldir, config.input.staticmaps),
    joinpath(tomldir, config.input.cyclic),
    joinpath(tomldir, config.input.forcing),
    joinpath(tomldir, config.output.path),
)

@unpack vertical, clock, reader, writer = model
@unpack dataset, buffer, inds = reader

# get the output path before it's closed, and remove up the file
output_path = path(model.writer.dataset)
Wflow.close_files(model)
rm(output_path)
