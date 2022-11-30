const map_structs = Dict("initialize" => Initialize)

context = Context()
socket = Socket(context, REP)
ZMQ.bind(socket, "tcp://*:5555")

mutable struct ModelHandler
    model::Wflow.Model
    ModelHandler() = new()
end


function shutdown(s::Socket, ctx::Context)
    ZMQ.close(s)
    ZMQ.close(ctx)
end

function response(err::AbstractString, s::Socket)
    resp = Dict{String,String}("status" => "ERROR", "error" => err)
    ZMQ.send(s, JSON3.write(resp))
end

function response(s::Socket)
    resp = Dict{String,String}("status" => "OK")
    ZMQ.send(s, JSON3.write(resp))
end

function valid_request(json)
    for f in fieldnames(map_structs[json.fn])
        if f ∉ keys(json)
            return f
            break
        end
    end
end

function wflow_bmi(s::Socket, handler::ModelHandler, f)
    try
        handler.model = wflow_bmi(f)
        response(s)
    catch e
        @error "Wflow BMI $(f.fn) failed" exception = (e, catch_backtrace())
        err = string(
            "Wflow BMI $(f.fn) failed\n exception =\n ",
            sprint(showerror, e, catch_backtrace()),
        )
        response(err, s)
    end
end


handler = ModelHandler()
try
    while true
        # Wait for next request from client
        req = ZMQ.recv(socket)
        json = JSON3.read(req)

        if haskey(map_structs, json.fn)
            v = valid_request(json)
            if isnothing(v)
                f = StructTypes.constructfrom(map_structs[json.fn], json)
                wflow_bmi(socket, handler, f)
            else
                err = ("""At least one required argument name 
                        ($v) not available for function: $(json.fn)""")
                @error err
                response(err, socket)
            end
        elseif json.fn === "shutdown"
            response(socket)
            shutdown(socket, context)
        else
            err = "Received invalid function: $(json.fn)"
            @error err
            response(err, socket)
        end
    end
catch e
    err = "ZMQ Wflow Server: exception in process"
    @error err
    response(err, socket)
finally
    shutdown(socket, context)
end
