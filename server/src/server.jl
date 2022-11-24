
const map_structs = Dict(
    "initialize" => Initialize
)

context = Context()
socket = Socket(context, REP)
ZMQ.bind(socket, "tcp://*:5555")

mutable struct ModelHandler
    model::Wflow.Model
    ModelHandler() = new()
end

handler = ModelHandler()

try
    while true
        # Wait for next request from client
        req = ZMQ.recv(socket)
        json = JSON3.read(req)
        f = StructTypes.constructfrom(map_structs[json.fn], json)
        if f.fn === "terminate"
            ZMQ.send(socket, "terminate")
            break
        else
            try
                handler.model = wflow_bmi(f)
                ZMQ.send(socket, "done")
            catch e
                ZMQ.send(socket, "error")
                @error "Wflow BMI failed" exception = (e, catch_backtrace())
            end
            
        end
    end
catch e
    @warn "ZMQ Wflow Server: exception in process"
finally
    ZMQ.close(socket)
    ZMQ.close(context)
end