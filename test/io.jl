using NCDatasets
using Dates
using CFTime
using Random

tomlpath = joinpath(@__DIR__, "config.toml")
tomldir = dirname(tomlpath)
config = Wflow.Config(TOML.parsefile(tomlpath))

reader = NCDataset(joinpath(tomldir, config.input.forcing))

nclon = Float64.(nomissing(reader["lon"][:]))
nclat = Float64.(nomissing(reader["lat"][:]))

# TODO remove random string from the filename
# this makes it easier to develop for now, since we don't run into issues with open files
base, ext = splitext(config.output.path)
randomized_path = string(base, '_', randstring('a':'z', 4), ext)
output_path = joinpath(tomldir, randomized_path)

writer = Wflow.setup_netcdf(output_path, nclon, nclat)

# initialize a vector of SBM structs
model = Wflow.initialize_sbm_model(
    joinpath(tomldir, config.input.staticmaps),
    joinpath(tomldir, config.input.leafarea),
    reader,
    writer,
)

model.inds

Wflow.update_forcing!(model)


## write output function

# q = rand(291, 313)
# grow_netcdf!(writer, "q", 1, q)

# close both input and output datasets
close(reader)
close(writer)
rm(output_path)
