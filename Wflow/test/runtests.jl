## load test dependencies and set paths to testing data
using Accessors
using Dates
using Graphs
using NCDatasets
using StaticArrays
using Statistics
using Test
using Wflow
using Base.MathConstants: eulergamma
using Base.Threads
using BasicModelInterface
using Polynomials: Polynomials
using DelimitedFiles
using LoggingExtras
using QuadGK
using Aqua: Aqua

const BMI = BasicModelInterface

include("testing_utils.jl")

@info "testing Wflow with" nthreads() VERSION

# disable logging output during testing
with_logger(NullLogger()) do
    ## run all tests
    @testset "Wflow.jl" begin
        include("routing_process.jl")
        include("io.jl")
        include("land_process.jl")
        include("reservoir.jl")
        include("run_sbm.jl")
        include("run_sbm_piave.jl")
        include("run_sbm_gwf_piave.jl")
        include("run_sbm_gwf.jl")
        include("run.jl")
        include("groundwater.jl")
        include("utils.jl")
        include("bmi.jl")
        include("run_sediment.jl")
        include("subdomains.jl")
        Aqua.test_all(Wflow; ambiguities = false, persistent_tasks = false)
    end
end
