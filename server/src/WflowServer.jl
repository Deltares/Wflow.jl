module WflowServer
using ZMQ: ZMQ
using JSON: JSON
using Wflow: Wflow

include("bmi_service.jl")
include("server.jl")
end
