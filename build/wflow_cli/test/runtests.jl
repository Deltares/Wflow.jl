using Test
using Wflow

extension = Sys.iswindows() ? ".exe" : ""
wflow_exe = normpath(@__DIR__, "../../create_binaries/wflow_bundle/bin/wflow_cli" * extension)

# this assumes that the Wflow tests have already been run, so the data has been downloaded
testdir = abspath(dirname(pathof(Wflow)), "..", "test")

# ensure test data is present
# this code is copied from runtests.jl, and is a temporary solution to get the data in place
datadir = joinpath(testdir, "data")
outputdir = joinpath(datadir, "output")
rm(outputdir; force=true, recursive=true)

@testset "SBM" begin
    toml = normpath(testdir, "sbm_config.toml")
    run(`$wflow_exe $toml`)
end
