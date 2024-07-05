# map JSON function name to Struct (bmi_service.jl)
const map_structs = Dict(
    "initialize" => Initialize,
    "get_component_name" => GetComponentName,
    "get_input_item_count" => GetInputItemCount,
    "get_output_item_count" => GetOutputItemCount,
    "get_start_time" => GetStartTime,
    "get_start_unix_time" => GetStartUnixTime,
    "get_end_time" => GetEndTime,
    "get_time_step" => GetTimeStep,
    "get_time_units" => GetTimeUnits,
    "get_current_time" => GetCurrentTime,
    "update_until" => UpdateUntil,
    "update" => Update,
    "get_input_var_names" => GetInputVarNames,
    "get_output_var_names" => GetOutputVarNames,
    "get_var_itemsize" => GetVarItemSize,
    "get_var_type" => GetVarType,
    "get_var_units" => GetVarUnits,
    "get_var_location" => GetVarLocation,
    "get_var_nbytes" => GetVarNbytes,
    "get_value" => GetValue,
    "get_value_ptr" => GetValuePtr,
    "get_value_at_indices" => GetValueAtIndices,
    "set_value" => SetValue,
    "set_value_at_indices" => SetValueAtIndices,
    "get_grid_type" => GetGridType,
    "get_var_grid" => GetVarGrid,
    "get_grid_rank" => GetGridRank,
    "get_grid_size" => GetGridSize,
    "get_grid_x" => GetGridX,
    "get_grid_y" => GetGridY,
    "get_grid_node_count" => GetGridNodeCount,
    "get_grid_edge_count" => GetGridEdgeCount,
    "get_grid_edge_nodes" => GetGridEdgeNodes,
    "finalize" => Finalize,
    "load_state" => LoadState,
    "save_state" => SaveState,
)

mutable struct ModelHandler
    model::Union{Wflow.Model,Nothing}
end

"Shutdown ZMQ server"
function shutdown(s::ZMQ.Socket, ctx::ZMQ.Context)
    @info "Shutting down Wflow ZMQ server on request..."
    ZMQ.close(s)
    return ZMQ.close(ctx)
end

"Error response ZMQ server"
function response(err::AbstractString, s::ZMQ.Socket)
    @info "Send error response"
    resp = Dict{String,String}("status" => "ERROR", "error" => err)
    return ZMQ.send(s, JSON3.write(resp))
end

"Status response ZMQ server"
function response(s::ZMQ.Socket)
    @info "Send status response"
    resp = Dict{String,String}("status" => "OK")
    return ZMQ.send(s, JSON3.write(resp))
end

"Validate JSON request against mapped Struct"
function valid_request(json)
    for f in fieldnames(map_structs[json.fn])
        if f ∉ keys(json)
            return f
            break
        end
    end
end

"""
    wflow_bmi(s::Socket, handler::ModelHandler, f)

Run a Wflow function through `wflow.bmi(f, handler.model)` and update Wflow Model `handler`
if required, depending on return type of `wflow.bmi(f, handler.model)`.
"""
function wflow_bmi(s::ZMQ.Socket, handler::ModelHandler, f)
    try
        ret = wflow_bmi(f, handler.model)
        if typeof(ret) <: Wflow.Model # update of Wflow model
            handler.model = ret
            response(s)
        elseif isnothing(ret) # for SetValue and SetValueAtIndices
            response(s)
        else
            @info "Send response including output from Wflow function `$(f.fn)`"
            ZMQ.send(s, JSON3.write(ret))
        end
    catch e
        @error "Wflow function `$(f.fn)` failed" exception = (e, catch_backtrace())
        err = string(
            "Wflow function `$(f.fn)` failed\n exception =\n ",
            sprint(showerror, e, catch_backtrace()),
        )
        response(err, s)
    end
end

main() = main(ARGS)

"""
    main(ARGS::Vector{String})
    main()

This is the main entry point of the Wflow ZMQ Server. Performs argument parsing and starts
the Wflow ZMQ Server, with `WflowServer.start(port::Int)`.
"""
function main(ARGS::Vector{String})
    n = length(ARGS)
    if n == 0
        port = 5555
    elseif startswith(ARGS[1], "--port")
        if occursin("=", ARGS[1])
            port = parse(Int, split(ARGS[1], "=")[2])
        else
            port = parse(Int, ARGS[2])
        end
    else
        throw(
            ArgumentError(
                "One argument is allowed to specify the port number: `--port=<port>` or `--port <port>`, 
                where `<port>` refers to the port number.",
            ),
        )
    end
    return start(port)
end

"""
    start(port::Int)

Start the Wflow ZMQ Server using port number `port`.
"""
function start(port::Int)
    @info "Start Wflow ZMQ Server..."

    # initialize Wflow model handler
    handler = ModelHandler(nothing)

    # set up a ZMQ context, with optional port number (default = 5555) argument
    context = ZMQ.Context()
    socket = ZMQ.Socket(context, ZMQ.REP)
    ZMQ.bind(socket, "tcp://*:$port")

    try
        while true
            # Wait for next request from client
            req = ZMQ.recv(socket)
            json = JSON3.read(req)
            @info "Received request to run function `$(json.fn)`..."

            if haskey(map_structs, json.fn)
                v = valid_request(json)
                if isnothing(v)
                    f = StructTypes.constructfrom(map_structs[json.fn], json)
                    wflow_bmi(socket, handler, f)
                else
                    err = (
                        """At least one required argument name (`$v`) not available for function: `$(json.fn)`"""
                    )
                    @error err
                    response(err, socket)
                end
            elseif json.fn === "shutdown"
                response(socket)
                break
            else
                err = "Received invalid Wflow function: `$(json.fn)`"
                @error err
                response(err, socket)
            end
        end
    catch e
        err = "Wflow ZMQ Server: exception in process"
        @error err
        response(err, socket)
    finally
        shutdown(socket, context)
    end
end
