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

function request(message)
    ZMQ.send(socket, JSON3.write(message))
    ret_value = JSON3.read(ZMQ.recv(socket), Dict)
    return ret_value
end

@testset "Client requests" begin
    msg = Dict(
        "fn" => "initialize",
        "config_file" => joinpath(@__DIR__, "../../test/sbm_config.toml"),
    )
    @test request(msg) == Dict("status" => "OK")

    msg = Dict("fn" => "get_end_time")
    @test request(msg) == Dict("end_time" => 949363200)

    msg = Dict("fn" => "get_start_time")
    @test request(msg) == Dict("start_time" => 946771200)

    msg = Dict("fn" => "update_until", "time" => 9.468576e8)
    @test request(msg) == Dict("status" => "OK")

    msg = Dict("fn" => "get_current_time")
    @test request(msg) == Dict("current_time" => 946857600)

    msg = Dict("fn" => "update")
    @test request(msg) == Dict("status" => "OK")

    msg = Dict("fn" => "get_component_name")
    @test request(msg) == Dict("component_name" => "sbm")

    msg = Dict("fn" => "get_input_item_count")
    @test request(msg) == Dict("input_item_count" => 186)

    msg = Dict("fn" => "get_output_item_count")
    @test request(msg) == Dict("output_item_count" => 186)

    msg = Dict("fn" => "get_time_step")
    @test request(msg) == Dict("time_step" => 86400)

    msg = Dict("fn" => "get_time_units")
    @test request(msg) == Dict("time_units" => "seconds since 1970-01-01T00:00:00")

    msg = Dict("fn" => "get_input_var_names")
    @test request(msg)["input_var_names"][[1, 5, 120, 186]] == [
        "vertical.Δt",
        "vertical.n_unsatlayers",
        "lateral.land.q_av",
        "lateral.river.reservoir.evaporation",
    ]
    msg = Dict("fn" => "get_output_var_names")
    @test request(msg)["output_var_names"][[1, 5, 120, 186]] == [
        "vertical.Δt",
        "vertical.n_unsatlayers",
        "lateral.land.q_av",
        "lateral.river.reservoir.evaporation",
    ]
    msg = Dict("fn" => "get_var_itemsize", "name" => "lateral.subsurface.ssf")
    @test request(msg) == Dict("var_itemsize" => sizeof(Wflow.Float))

    msg = Dict("fn" => "get_var_type", "name" => "vertical.n")
    @test request(msg) == Dict("var_type" => "Int64")

    msg = Dict("fn" => "get_var_units", "name" => "vertical.θₛ")
    @test request(msg) == Dict("var_units" => "-")

    msg = Dict("fn" => "get_var_location", "name" => "lateral.river.q")
    @test request(msg) == Dict("var_location" => "node")

    msg = Dict("fn" => "get_var_nbytes", "name" => "vertical.n")
    @test request(msg) == Dict("var_nbytes" => 8)

    msg = Dict("fn" => "get_value", "name" => "vertical.zi")
    @test mean(request(msg)["value"]) ≈ 278.1510965581235

    msg = Dict("fn" => "get_value_ptr", "name" => "vertical.θₛ")
    @test mean(request(msg)["value_ptr"]) ≈ 0.4409211971535584

    msg = Dict(
        "fn" => "get_value_at_indices",
        "name" => "lateral.river.q",
        "inds" => [1, 5, 10],
    )
    @test request(msg)["value_at_indices"] ≈
          [2.1901434445889123, 2.6778265820894545, 3.470059871798756]

    msg = Dict("fn" => "get_grid_size", "grid" => 4)
    n = request(msg)["grid_size"]

    msg = Dict("fn" => "set_value", "name" => "vertical.zi", "src" => fill(300.0, n))
    @test request(msg) == Dict("status" => "OK")

    msg = Dict("fn" => "get_value", "name" => "vertical.zi")
    @test mean(request(msg)["value"]) == 300.0

    msg = Dict(
        "fn" => "set_value_at_indices",
        "name" => "vertical.zi",
        "src" => [250.0, 350.0],
        "inds" => [1, 2],
    )
    @test request(msg) == Dict("status" => "OK")

    msg = Dict("fn" => "get_value_at_indices", "name" => "vertical.zi", "inds" => [1, 2, 3])
    @test request(msg)["value_at_indices"] == [250.0, 350.0, 300.0]

    msg = Dict("fn" => "get_grid_type", "grid" => 0)
    @test request(msg) == Dict("grid_type" => "unstructured")

    msg = Dict("fn" => "get_var_grid", "name" => "lateral.river.h")
    @test request(msg) == Dict("var_grid" => 2)

    msg = Dict("fn" => "get_grid_shape", "grid" => 4)
    @test request(msg) == Dict("grid_shape" => 50063)

    msg = Dict("fn" => "get_grid_rank", "grid" => 0)
    @test request(msg) == Dict("grid_rank" => 2)

    msg = Dict("fn" => "get_grid_x", "grid" => 0)
    @test request(msg)["grid_x"] ≈ [5.760000000000002f0, 5.918333333333336f0]

    msg = Dict("fn" => "get_grid_y", "grid" => 0)
    @test request(msg)["grid_y"] ≈ [48.92583333333333f0, 49.909166666666664f0]

    msg = Dict("fn" => "get_grid_node_count", "grid" => 0)
    @test request(msg) == Dict("grid_node_count" => 2)

    msg = Dict("fn" => "load_state")
    @test request(msg) == Dict("status" => "OK")

    msg = Dict("fn" => "save_state")
    @test request(msg) == Dict("status" => "OK")

    msg = Dict("fn" => "finalize")
    @test request(msg) == Dict("status" => "OK")
end

@testset "Error handling and shutdown" begin
    msg = Dict("fn" => "not_existing_function")
    @test request(msg) == Dict(
        "status" => "ERROR",
        "error" => "Received invalid function: not_existing_function",
    )

    msg = Dict("fn" => "initialize")
    @test request(msg) == Dict(
        "status" => "ERROR",
        "error" => "At least one required argument name \n(config_file) not available for function: initialize",
    )

    msg = Dict("fn" => "shutdown")
    @test request(msg) == Dict("status" => "OK")
end
