## load test dependencies and set paths to testing data
using Dates
using LightGraphs
using NCDatasets
using StaticArrays
using Test
using UnPack
using Wflow
using CSV
using BenchmarkTools
using Base.MathConstants: eulergamma
import Polynomials

# ensure test data is present
testdir = @__DIR__
datadir = joinpath(testdir, "data")
isdir(datadir) || mkdir(datadir)

"Download a test data file if it does not already exist"
function testdata(version, source_filename, target_filename)
    target_path = joinpath(datadir, target_filename)
    base_url = "https://github.com/visr/wflow-artifacts/releases/download"
    url = string(base_url, '/', string('v', version), '/', source_filename)
    isfile(target_path) || download(url, target_path)
    return target_path
end

staticmaps_rhine_path = testdata(v"0.1", "staticmaps.nc", "staticmaps-rhine.nc")
staticmaps_moselle_path = testdata(v"0.2.1", "staticmaps.nc", "staticmaps-moselle.nc")
staticmaps_lahn_path = testdata(v"0.2.1", "staticmaps-lahn.nc", "staticmaps-lahn.nc")
forcing_moselle_path = testdata(v"0.2", "forcing-2000.nc", "forcing-moselle.nc")
forcing_lahn_path = testdata(v"0.2", "forcing-lahn.nc", "forcing-lahn.nc")
instates_moselle_path = testdata(v"0.2.1", "instates-moselle.nc", "instates-moselle.nc")
forcing_sbm_gw_path = testdata(v"0.2.1", "forcing-sbm-groundwater.nc", "forcing-sbm-groundwater.nc")
staticmaps_sbm_gw_path = testdata(v"0.2.1", "staticmaps-sbm-groundwater.nc", "staticmaps-sbm-groundwater.nc")
lake_sh_1_path = testdata(v"0.2.1", "lake_sh_1.csv", "lake_sh_1.csv")
lake_sh_2_path = testdata(v"0.2.1", "lake_sh_2.csv", "lake_sh_2.csv")
lake_hq_2_path = testdata(v"0.2.1", "lake_hq_2.csv", "lake_hq_2.csv")

include("testing_utils.jl")

## run all tests
@testset "Wflow.jl" begin
    include("horizontal_process.jl")
    include("io.jl")
    include("vertical_process.jl")
    include("reservoir_lake.jl")
    include("run_sbm.jl")
    include("run_hbv.jl")
    include("run_sbm_gwf.jl")
    include("run.jl")
    include("groundwater.jl")
    include("utils.jl")
end
