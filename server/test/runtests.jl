using ZMQ: ZMQ
using JSON3: JSON3
using StructTypes: StructTypes
using Wflow: Wflow
using WflowServer: WflowServer
import Statistics: mean
import Logging: with_logger, NullLogger
import Test: @testset, @test

with_logger(NullLogger()) do
    @testset "Test client server Wflow ZMQ Server" begin
        include("client.jl")
    end
end
