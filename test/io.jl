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

# TODO remove random string from the filename
# this makes it easier to develop for now, since we don't run into issues with open files
base, ext = splitext(config.output.path)
randomized_path = string(base, '_', randstring('a':'z', 4), ext)
output_path = joinpath(tomldir, randomized_path)

# initialize a vector of SBM structs
model = Wflow.initialize_sbm_model(
    joinpath(tomldir, config.input.staticmaps),
    joinpath(tomldir, config.input.leafarea),
    joinpath(tomldir, config.input.forcing),
    output_path,
)

@unpack vertical, clock, reader = model
@unpack dataset, buffer, inds = reader

Wflow.update_forcing!(model)

## write output function

nclon = Float64.(nomissing(reader.dataset["lon"][:]))
nclat = Float64.(nomissing(reader.dataset["lat"][:]))

# writer = model.writer
# writer = Wflow.setup_netcdf(output_path, nclon, nclat)

# q = rand(291, 313)
# grow_netcdf!(writer, "q", 1, q)

# close both input and output datasets
close(reader.dataset)
# close(writer)
# rm(output_path)
