## load test dependencies and set paths to testing data
using Wflow
using Test
using NCDatasets
using Dates
using LightGraphs
using Pkg.TOML

# ensure test data is present
datadir = joinpath(@__DIR__, "data")
isdir(datadir) || mkdir(datadir)
staticmaps_rhine_path = joinpath(datadir, "staticmaps-rhine.nc")
isfile(staticmaps_rhine_path) || download("https://github.com/visr/wflow-artifacts/releases/download/v0.1.0/staticmaps.nc", staticmaps_rhine_path)
staticmaps_moselle_path = joinpath(datadir, "staticmaps-moselle.nc")
isfile(staticmaps_moselle_path) || download("https://github.com/visr/wflow-artifacts/releases/download/v0.2.0/staticmaps.nc", staticmaps_moselle_path)
leafarea_moselle_path = joinpath(datadir, "lai_clim-moselle.nc")
isfile(leafarea_moselle_path) || download("https://github.com/visr/wflow-artifacts/releases/download/v0.2.0/lai_clim.nc", leafarea_moselle_path)

## run all tests

@testset "Wflow.jl" begin
    include("config.jl")
    include("kinematic_wave.jl")
    include("sbm.jl")
end
