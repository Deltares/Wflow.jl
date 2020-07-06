using NCDatasets
using Dates
using CFTime
using Random
using UnPack

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
config = Wflow.Config(TOML.parsefile(tomlpath))

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

# close both input and output datasets
close(reader.dataset)
close(reader.cyclic_dataset)
output_path = path(writer.dataset)
close(writer.dataset)
rm(output_path)
