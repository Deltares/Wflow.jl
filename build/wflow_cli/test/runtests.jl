using Test
using Wflow

extension = Sys.iswindows() ? ".exe" : ""
wflow_exe =
    normpath(@__DIR__, "../../create_binaries/wflow_bundle/bin/wflow_cli" * extension)

# this assumes that the Wflow tests have already been run, so the data has been downloaded
testdir = abspath(dirname(pathof(Wflow)), "..", "test")

# ensure test data is present
# this code is copied from runtests.jl, and is a temporary solution to get the data in place
datadir = joinpath(testdir, "data")
outputdir = joinpath(datadir, "output")

@testset "no_config" begin
    # Clean output directory
    rm(outputdir; force = true, recursive = true)
    @test_throws ProcessFailedException run(`$wflow_exe`)
    # Check if no files are being created
    @test !(isdir(outputdir))
end

@testset "wflow_sbm" begin
    # Clean directory
    rm(outputdir; force = true, recursive = true)
    # Run cli with the toml
    toml = normpath(testdir, "sbm_config.toml")
    run(`$wflow_exe $toml`)

    # List of files to check if they exists
    out_files = [
        "output_moselle.csv",
        "output_moselle.nc",
        "output_scalar_moselle.nc",
        "outstates-moselle.nc",
        "log.txt",
        "sbm_config.toml",
    ]
    for file in out_files
        @test filesize(joinpath(outputdir, file)) > 0
    end
end

@testset "wflow_sbm-gwf" begin
    # Clean directory
    rm(outputdir; force = true, recursive = true)
    # Run cli with the toml
    toml = normpath(testdir, "sbm_gwf_config.toml")
    run(`$wflow_exe $toml`)

    # List of files to check if they exists
    out_files = [
        "outstates-example-sbm-gwf.nc",
        "output_example-sbm-gwf.nc",
        "output_example-sbm-gwf.csv",
        "log.txt",
        "sbm_gwf_config.toml",
    ]
    for file in out_files
        @test filesize(joinpath(outputdir, file)) > 0
    end
end

@testset "wflow_sediment" begin
    # Clean directory
    rm(outputdir; force = true, recursive = true)
    # Run cli with the toml
    toml = normpath(testdir, "sediment_config.toml")
    run(`$wflow_exe $toml`)

    # List of files to check if they exists
    out_files = [
        "outstates-moselle-sed.nc",
        "output-moselle-sed.nc",
        "output-moselle-sediment.csv",
        "sediment_config.toml",
        "log.txt",
    ]
    for file in out_files
        @test filesize(joinpath(outputdir, file)) > 0
    end
end

@testset "wflow_sbm_timing" begin
    toml = normpath(testdir, "sbm_config.toml")
    time_sbm_1thread =
        @elapsed run(Cmd(`$wflow_exe $toml`; env = ("JULIA_NUM_THREADS" => "1",)))

    time_sbm_4thread =
        @elapsed run(Cmd(`$wflow_exe $toml`; env = ("JULIA_NUM_THREADS" => "4",)))

    # Test if run with more threads is indeed faster
    @test time_sbm_4thread < time_sbm_1thread
    # Test timings of different runs (very depending on machine, be careful with numbers!)
    @test time_sbm_1thread < 35
    @test time_sbm_4thread < 25
end
