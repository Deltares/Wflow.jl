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

    @testset "initialization and time functions" begin
        msg = Dict(
            "fn" => "initialize",
            "config_file" => joinpath(@__DIR__, "../../test/sbm_config.toml"),
        )
        @test request(msg) == Dict("status" => "OK")
        msg = Dict("fn" => "get_end_time")
        @test request(msg) == Dict("end_time" => 2678400)
        msg = Dict("fn" => "get_start_time")
        @test request(msg) == Dict("start_time" => 0)
        msg = Dict("fn" => "get_start_unix_time")
        @test request(msg) == Dict("start_unix_time" => 946684800)
        msg = Dict("fn" => "get_time_step")
        @test request(msg) == Dict("time_step" => 86400)
        msg = Dict("fn" => "get_time_units")
        @test request(msg) == Dict("time_units" => "s")
    end

    @testset "update functions" begin
        msg = Dict("fn" => "update_until", "time" => 86400.0)
        @test request(msg) == Dict("status" => "OK")
        msg = Dict("fn" => "get_current_time")
        @test request(msg) == Dict("current_time" => 86400)
        msg = Dict("fn" => "update")
        @test request(msg) == Dict("status" => "OK")
    end

    @testset "model information functions" begin
        msg = Dict("fn" => "get_component_name")
        @test request(msg) == Dict("component_name" => "sbm")
        msg = Dict("fn" => "get_input_item_count")
        @test request(msg) == Dict("input_item_count" => 180)
        msg = Dict("fn" => "get_output_item_count")
        @test request(msg) == Dict("output_item_count" => 180)
        msg = Dict("fn" => "get_input_var_names")
        @test request(msg)["input_var_names"][[1, 5, 151, 175]] == [
            "vertical.nlayers",
            "vertical.θᵣ",
            "lateral.river.q",
            "lateral.river.reservoir.outflow",
        ]
        msg = Dict("fn" => "get_output_var_names")
        @test request(msg)["output_var_names"][[1, 5, 151, 175]] == [
            "vertical.nlayers",
            "vertical.θᵣ",
            "lateral.river.q",
            "lateral.river.reservoir.outflow",
        ]
    end

    zi_size = 0
    vwc_1_size = 0
    @testset "variable information functions" begin
        msg = Dict("fn" => "get_var_itemsize", "name" => "lateral.subsurface.ssf")
        @test request(msg) == Dict("var_itemsize" => sizeof(Wflow.Float))
        msg = Dict("fn" => "get_var_type", "name" => "vertical.n")
        @test request(msg)["status"] == "ERROR"
        msg = Dict("fn" => "get_var_units", "name" => "vertical.θₛ")
        @test request(msg) == Dict("var_units" => "-")
        msg = Dict("fn" => "get_var_location", "name" => "lateral.river.q")
        @test request(msg) == Dict("var_location" => "node")
        msg = Dict("fn" => "get_var_nbytes", "name" => "vertical.zi")
        zi_nbytes = request(msg)["var_nbytes"]
        @test zi_nbytes == 400504
        msg = Dict("fn" => "get_var_itemsize", "name" => "vertical.zi")
        zi_itemsize = request(msg)["var_itemsize"]
        zi_size = Int(zi_nbytes / zi_itemsize)
        msg = Dict("fn" => "get_var_nbytes", "name" => "vertical.vwc[1]")
        vwc_1_nbytes = request(msg)["var_nbytes"]
        @test vwc_1_nbytes == 400504
        msg = Dict("fn" => "get_var_itemsize", "name" => "vertical.vwc[1]")
        vwc_1_itemsize = request(msg)["var_itemsize"]
        vwc_1_size = Int(vwc_1_nbytes / vwc_1_itemsize)
        msg = Dict("fn" => "get_var_grid", "name" => "lateral.river.h")
        @test request(msg) == Dict("var_grid" => 3)
    end

    @testset "get and set functions" begin
        msg =
            Dict("fn" => "get_value", "name" => "vertical.zi", "dest" => fill(0.0, zi_size))
        @test mean(request(msg)["value"]) ≈ 278.1510965581235
        msg = Dict("fn" => "get_value_ptr", "name" => "vertical.θₛ")
        @test mean(request(msg)["value_ptr"]) ≈ 0.4409211971535584
        msg = Dict(
            "fn" => "get_value_at_indices",
            "name" => "lateral.river.q",
            "dest" => [0.0, 0.0, 0.0],
            "inds" => [1, 5, 10],
        )
        @test request(msg)["value_at_indices"] ≈
              [2.1901434445889123, 2.6778265820894545, 3.470059871798756]
        msg = Dict(
            "fn" => "set_value",
            "name" => "vertical.zi",
            "src" => fill(300.0, zi_size),
        )
        @test request(msg) == Dict("status" => "OK")
        msg =
            Dict("fn" => "get_value", "name" => "vertical.zi", "dest" => fill(0.0, zi_size))
        @test mean(request(msg)["value"]) == 300.0
        msg = Dict(
            "fn" => "set_value_at_indices",
            "name" => "vertical.zi",
            "src" => [250.0, 350.0],
            "inds" => [1, 2],
        )
        @test request(msg) == Dict("status" => "OK")
        msg = Dict(
            "fn" => "get_value_at_indices",
            "name" => "vertical.zi",
            "dest" => [0.0, 0.0, 0.0],
            "inds" => [1, 2, 3],
        )
        @test request(msg)["value_at_indices"] == [250.0, 350.0, 300.0]
        msg = Dict(
            "fn" => "get_value",
            "name" => "vertical.vwc[1]",
            "dest" => fill(0.0, vwc_1_size),
        )
        @test mean(request(msg)["value"]) ≈ 0.1845159140308566f0
        msg = Dict(
            "fn" => "get_value_at_indices",
            "name" => "vertical.vwc[1]",
            "dest" => [0.0, 0.0, 0.0],
            "inds" => [1, 2, 3],
        )
        @test request(msg)["value_at_indices"] ≈
              [0.12089607119560242f0, 0.11967185322648062f0, 0.14503555864288548f0]
        msg = Dict(
            "fn" => "set_value",
            "name" => "vertical.vwc[1]",
            "src" => fill(0.3, vwc_1_size),
        )
        @test request(msg) == Dict("status" => "OK")
        msg = Dict(
            "fn" => "get_value",
            "name" => "vertical.vwc[1]",
            "dest" => fill(0.0, vwc_1_size),
        )
        @test mean(request(msg)["value"]) ≈ 0.3f0
        msg = Dict(
            "fn" => "get_value_at_indices",
            "name" => "vertical.vwc[1]",
            "dest" => [0.0, 0.0, 0.0],
            "inds" => [1, 2, 3],
        )
        @test request(msg)["value_at_indices"] == [0.3, 0.3, 0.3]
        msg = Dict(
            "fn" => "set_value_at_indices",
            "name" => "vertical.vwc[1]",
            "src" => [0.1, 0.25],
            "inds" => [1, 2],
        )
        @test request(msg) == Dict("status" => "OK")
        msg = Dict(
            "fn" => "get_value_at_indices",
            "name" => "vertical.vwc[1]",
            "dest" => [0.0, 0.0],
            "inds" => [1, 2],
        )
        @test request(msg)["value_at_indices"] == [0.1, 0.25]
    end

    @testset "model grid functions" begin
        msg = Dict("fn" => "get_grid_type", "grid" => 0)
        @test request(msg) == Dict("grid_type" => "points")
        msg = Dict("fn" => "get_grid_rank", "grid" => 0)
        @test request(msg) == Dict("grid_rank" => 2)
        msg = Dict("fn" => "get_grid_size", "grid" => 4)
        grid_size = request(msg)["grid_size"]
        @test grid_size == 50063
        msg = Dict("fn" => "get_grid_x", "grid" => 4, "x" => fill(0.0, grid_size))
        @test request(msg)["grid_x"][1:3] ≈
              [6.826666666666673, 6.810000000000006, 6.81833333333334]
        msg = Dict("fn" => "get_grid_y", "grid" => 4, "y" => fill(0.0, grid_size))
        @test request(msg)["grid_y"][1:3] ≈
              [47.8175, 47.825833333333335, 47.825833333333335]
        msg = Dict("fn" => "get_grid_node_count", "grid" => 0)
        @test request(msg) == Dict("grid_node_count" => 2)
        msg = Dict("fn" => "get_grid_edge_count", "grid" => 3)
        @test request(msg)["grid_edge_count"] == 5808
        msg = Dict(
            "fn" => "get_grid_edge_nodes",
            "grid" => 3,
            "edge_nodes" => fill(0, 2 * 5808),
        )
        @test request(msg)["grid_edge_nodes"][1:6] == [1, 5, 2, 1, 3, 2]
    end

    @testset "model states and finalize functions" begin
        msg = Dict("fn" => "load_state")
        @test request(msg) == Dict("status" => "OK")
        msg = Dict("fn" => "save_state")
        @test request(msg) == Dict("status" => "OK")
        msg = Dict("fn" => "finalize")
        @test request(msg) == Dict("status" => "OK")
    end
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
