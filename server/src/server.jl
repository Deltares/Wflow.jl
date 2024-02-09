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
function shutdown(s::Socket, ctx::Context)
    ZMQ.close(s)
    ZMQ.close(ctx)
end

"Error response ZMQ server"
function response(err::AbstractString, s::Socket)
    resp = Dict{String,String}("status" => "ERROR", "error" => err)
    ZMQ.send(s, JSON3.write(resp))
end

"Status response ZMQ server"
function response(s::Socket)
    resp = Dict{String,String}("status" => "OK")
    ZMQ.send(s, JSON3.write(resp))
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

Run a Wflow.BMI (or Wflow) function through `wflow.bmi(f, handler.model)` and update Wflow
Model `handler` if required, depending on return type of `wflow.bmi(f, handler.model)`.
"""
function wflow_bmi(s::Socket, handler::ModelHandler, f)
    try
        ret = wflow_bmi(f, handler.model)
        if typeof(ret) <: Wflow.Model
            handler.model = ret
            response(s)
        elseif isnothing(ret)
            response(s)
        else
            ZMQ.send(s, JSON3.write(ret))
        end
    catch e
        @error "Wflow BMI $(f.fn) failed" exception = (e, catch_backtrace())
        err = string(
            "Wflow BMI $(f.fn) failed\n exception =\n ",
            sprint(showerror, e, catch_backtrace()),
        )
        response(err, s)
    end
end

# initialize Wflow model handler
handler = ModelHandler(nothing)

# set up a ZMQ context, with optional port number (default = 5555) argument
context = Context()
socket = Socket(context, REP)
n = length(ARGS)
if n == 0
    port = 5555
elseif n == 1
    port = parse(Int, ARGS[1])
else
    throw(
        ArgumentError(
            "More than one argument provided, while only one port number is allowed.",
        ),
    )
end
ZMQ.bind(socket, "tcp://*:$port")

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
            break
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
