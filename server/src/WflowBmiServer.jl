module WflowBmiServer
    using ZMQ
    using JSON3
    using StructTypes
    using Wflow

    include("bmi_api.jl")
    include("server.jl")
end