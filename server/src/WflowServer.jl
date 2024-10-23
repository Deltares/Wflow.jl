module WflowServer
using ZMQ: ZMQ
using JSON3: JSON3
using StructTypes: StructTypes
using Wflow: Wflow

include("bmi_service.jl")
include("server.jl")
end
