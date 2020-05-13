## load test dependencies and set paths to testing data
using Dates
using LightGraphs
using NCDatasets
using Pkg.TOML
using Setfield: setproperties
using StaticArrays
using Test
using UnPack
using Wflow

# ensure test data is present
const datadir = joinpath(@__DIR__, "data")
isdir(datadir) || mkdir(datadir)

"Download a test data file if it does not already exist"
function testdata(version, source_filename, target_filename)
    target_path = joinpath(datadir, target_filename)
    base_url = "https://github.com/visr/wflow-artifacts/releases/download"
    url = joinpath(base_url, string('v', version), source_filename)
    isfile(target_path) || download(url, target_path)
    return target_path
end

staticmaps_rhine_path = testdata(v"0.1", "staticmaps.nc", "staticmaps-rhine.nc")
staticmaps_moselle_path = testdata(v"0.2", "staticmaps.nc", "staticmaps-moselle.nc")
leafarea_moselle_path = testdata(v"0.2", "lai_clim.nc", "lai_clim-moselle.nc")

## run all tests

@testset "Wflow.jl" begin
    include("config.jl")
    include("kinematic_wave.jl")
    include("sbm.jl")
    include("vertical_process.jl")
    include("run_sbm.jl")
end
