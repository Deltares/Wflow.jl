using ZMQ
using JSON3
using StructTypes
using Wflow
using WflowServer
using Statistics
using Logging
using Test

# The tests depend on wflow_sbm model files (Moselle), which are downloaded as part of the
# Wflow.jl tests, so these need to run first.
with_logger(NullLogger()) do
    @testset "Test client server Wflow ZMQ Server" begin
        include("client.jl")
    end
end