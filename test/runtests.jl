using Wflow
using Test
using NCDatasets
using LightGraphs

# ensure test data is present
datadir = joinpath(@__DIR__, "data")
isdir(datadir) || mkdir(datadir)
staticmaps_path = joinpath(datadir, "staticmaps.nc")
isfile(staticmaps_path) || download("https://github.com/visr/wflow-artifacts/releases/download/v0.1.0/staticmaps.nc", staticmaps_path)

@testset "Wflow.jl" begin
    include("kinematic_wave.jl")
end
