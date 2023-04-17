using ZMQ
using JSON3
using StructTypes
using Wflow
using Statistics
using Logging
using Test

with_logger(NullLogger()) do
    @testset "Test Wflow client server Wflow (BMI and ZMQ)" begin
        include("client.jl")
    end
end