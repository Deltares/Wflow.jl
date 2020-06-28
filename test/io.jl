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

tomlpath = joinpath(@__DIR__, "config.toml")
tomldir = dirname(tomlpath)
config = Wflow.Config(TOML.parsefile(tomlpath))

# initialize a vector of SBM structs
model = Wflow.initialize_sbm_model(
    config,
    joinpath(tomldir, config.input.staticmaps),
    joinpath(tomldir, config.input.leafarea),
    joinpath(tomldir, config.input.forcing),
    joinpath(tomldir, config.output.path),
)

@unpack vertical, clock, reader, writer = model
@unpack dataset, buffer, inds = reader

Wflow.update_forcing!(model)

for parameter in writer.parameters
    A = rand(Float32, size(buffer))
    Wflow.grow_netcdf!(writer.dataset, parameter, clock.time, A)
end

# close both input and output datasets
close(reader.dataset)
output_path = path(writer.dataset)
close(writer.dataset)
rm(output_path)
