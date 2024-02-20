import ZMQ
import JSON3
import StructTypes
import Wflow
import WflowServer
import Statistics: mean
import Logging: with_logger, NullLogger
import Test: @testset, @test

# The tests depend on wflow_sbm model files (Moselle), which are downloaded as part of the
# Wflow.jl tests, so these need to run first.
with_logger(NullLogger()) do
    @testset "Test client server Wflow ZMQ Server" begin
        include("client.jl")
    end
end
