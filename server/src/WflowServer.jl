module WflowServer
using ZMQ
using JSON3
using StructTypes
using Wflow

include("bmi_service.jl")
include("server.jl")
end
