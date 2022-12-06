using ZMQ
using JSON3

@info("start ZMQ server for Wflow (@async)")
@async begin
    include("../src/WflowBmiServer.jl")
end
sleep(2)
@info("Started ZMQ server for Wflow")

context = Context()

@info("Connecting to the ZMQ server for Wflow...")
socket = Socket(context, REQ)
ZMQ.connect(socket, "tcp://localhost:5555")

@info("Request Wflow initialization...")
message = Dict{String, Any}("fn" => "initialize", "config_file" => joinpath(@__DIR__, "../../test/sbm_config.toml"))
ZMQ.send(socket, JSON3.write(message))
ret_value = String(ZMQ.recv(socket))
println("Received reply [ $ret_value]")

@info("Request start time...")
message = Dict{String, Any}("fn" => "get_start_time")
ZMQ.send(socket, JSON3.write(message))
ret_value = String(ZMQ.recv(socket))
println("Received reply [ $ret_value]")

ZMQ.close(socket)
ZMQ.close(context)