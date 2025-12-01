@testitem "Client server Wflow ZMQ Server" begin
    import ZMQ, JSON3, StructTypes, Wflow
    using Statistics: mean
    using Logging: with_logger, NullLogger

    # start Wflow ZMQ server (@async)
    @async begin
        WflowServer.main()
    end

    # Connecting to the Wflow ZMQ Server
    context = ZMQ.Context()
    socket = ZMQ.Socket(context, ZMQ.REQ)
    ZMQ.connect(socket, "tcp://localhost:5555")

    function request(message)
        ZMQ.send(socket, JSON3.write(message))
        ret_value = JSON3.read(ZMQ.recv(socket), Dict; allow_inf = true)
        return ret_value
    end

    @testset "initialization and time functions" begin
        msg = (
            fn = "initialize",
            config_file = joinpath(@__DIR__, "../../Wflow/test/sbm_config.toml"),
        )
        @test request(msg) == Dict("status" => "OK")
        @test request((fn = "get_end_time",)) == Dict("end_time" => 2678400)
        @test request((fn = "get_start_time",)) == Dict("start_time" => 0)
        @test request((fn = "get_start_unix_time",)) == Dict("start_unix_time" => 946684800)
        @test request((fn = "get_time_step",)) == Dict("time_step" => 86400)
        @test request((fn = "get_time_units",)) == Dict("time_units" => "s")
    end

    @testset "Reading and writing NaN values allowed" begin
        msg = (
            fn = "get_value",
            name = "soil_layer_1_water__volume_fraction",
            dest = fill(0.0, 50063),
        )
        @test isnan(mean(request(msg)["value"]))
    end

    @testset "update functions" begin
        @test request((fn = "update_until", time = 86400.0)) == Dict("status" => "OK")
        @test request((fn = "get_current_time",)) == Dict("current_time" => 86400)
        @test request((fn = "update",)) == Dict("status" => "OK")
    end

    @testset "model information functions" begin
        @test request((fn = "get_component_name",)) == Dict("component_name" => "sbm")
        @test request((fn = "get_input_item_count",)) == Dict("input_item_count" => 7)
        @test request((fn = "get_output_item_count",)) == Dict("output_item_count" => 7)
        to_check = [
            "river_water__volume_flow_rate",
            "soil_water_unsaturated_zone__depth",
            "soil_water__transpiration_volume_flux",
            "soil_layer_2_water_unsaturated_zone__depth",
        ]
        retrieved_vars = request((fn = "get_input_var_names",))["input_var_names"]
        @test all(x -> x in retrieved_vars, to_check)
        retrieved_vars = request((fn = "get_output_var_names",))["output_var_names"]
        @test all(x -> x in retrieved_vars, to_check)
    end

    zi_size = 0
    vwc_1_size = 0
    @testset "variable information and get and set functions" begin
        @test request((
            fn = "get_var_itemsize",
            name = "subsurface_water__volume_flow_rate",
        )) == Dict("var_itemsize" => sizeof(Float64))
        @test request((fn = "get_var_units", name = "river_water__volume_flow_rate")) ==
              Dict("var_units" => "m3 s-1")
        @test request((fn = "get_var_location", name = "river_water__volume_flow_rate")) ==
              Dict("var_location" => "node")
        zi_nbytes =
            request((fn = "get_var_nbytes", name = "soil_water_saturated_zone_top__depth"))["var_nbytes"]
        @test zi_nbytes == 400504
        zi_itemsize = request((
            fn = "get_var_itemsize",
            name = "soil_water_saturated_zone_top__depth",
        ))["var_itemsize"]
        zi_size = Int(zi_nbytes / zi_itemsize)
        vwc_1_nbytes =
            request((fn = "get_var_nbytes", name = "soil_layer_1_water__volume_fraction"))["var_nbytes"]
        @test vwc_1_nbytes == 400504
        vwc_1_itemsize = request((
            fn = "get_var_itemsize",
            name = "soil_layer_1_water__volume_fraction",
        ))["var_itemsize"]
        vwc_1_size = Int(vwc_1_nbytes / vwc_1_itemsize)
        @test request((fn = "get_var_grid", name = "river_water__depth")) ==
              Dict("var_grid" => 2)
        msg = (
            fn = "get_value",
            name = "soil_water_saturated_zone_top__depth",
            dest = fill(0.0, zi_size),
        )
        @test mean(request(msg)["value"]) ≈ 277.88881025600784
        msg = (fn = "get_value_ptr", name = "soil_water_root_zone__depth")
        @test mean(request(msg)["value_ptr"]) ≈ 28.89577387100501
        msg = (
            fn = "get_value_at_indices",
            name = "river_water__instantaneous_volume_flow_rate",
            dest = [0.0, 0.0, 0.0],
            inds = [1, 5, 10],
        )
        @test request(msg)["value_at_indices"] ≈
              [2.2552403727984323, 2.7512381415155462, 3.5566676103297965]
        msg = (
            fn = "set_value",
            name = "soil_water_saturated_zone_top__depth",
            src = fill(300.0, zi_size),
        )
        @test request(msg) == Dict("status" => "OK")
        msg = (
            fn = "get_value",
            name = "soil_water_saturated_zone_top__depth",
            dest = fill(0.0, zi_size),
        )
        @test mean(request(msg)["value"]) == 300.0
        msg = (
            fn = "set_value_at_indices",
            name = "soil_water_saturated_zone_top__depth",
            src = [250.0, 350.0],
            inds = [1, 2],
        )
        @test request(msg) == Dict("status" => "OK")
        msg = (
            fn = "get_value_at_indices",
            name = "soil_water_saturated_zone_top__depth",
            dest = [0.0, 0.0, 0.0],
            inds = [1, 2, 3],
        )
        @test request(msg)["value_at_indices"] == [250.0, 350.0, 300.0]
        msg = (
            fn = "get_value",
            name = "soil_layer_1_water__volume_fraction",
            dest = fill(0.0, vwc_1_size),
        )
        @test mean(request(msg)["value"]) ≈ 0.18607159416148555
        msg = (
            fn = "get_value_at_indices",
            name = "soil_layer_1_water__volume_fraction",
            dest = [0.0, 0.0, 0.0],
            inds = [1, 2, 3],
        )
        @test request(msg)["value_at_indices"] ≈
              [0.12089607119560242, 0.11968416924304527, 0.14602328618707333]
        msg = (
            fn = "set_value",
            name = "soil_layer_1_water__volume_fraction",
            src = fill(0.3, vwc_1_size),
        )
        @test request(msg) == Dict("status" => "OK")
        msg = (
            fn = "get_value",
            name = "soil_layer_1_water__volume_fraction",
            dest = fill(0.0, vwc_1_size),
        )
        @test mean(request(msg)["value"]) ≈ 0.3
        msg = (
            fn = "get_value_at_indices",
            name = "soil_layer_1_water__volume_fraction",
            dest = [0.0, 0.0, 0.0],
            inds = [1, 2, 3],
        )
        @test request(msg)["value_at_indices"] == [0.3, 0.3, 0.3]
        msg = (
            fn = "set_value_at_indices",
            name = "soil_layer_1_water__volume_fraction",
            src = [0.1, 0.25],
            inds = [1, 2],
        )
        @test request(msg) == Dict("status" => "OK")
        msg = (
            fn = "get_value_at_indices",
            name = "soil_layer_1_water__volume_fraction",
            dest = [0.0, 0.0],
            inds = [1, 2],
        )
        @test request(msg)["value_at_indices"] == [0.1, 0.25]
    end

    @testset "model grid functions" begin
        @test request((fn = "get_grid_type", grid = 0)) == Dict("grid_type" => "points")
        @test request((fn = "get_grid_rank", grid = 0)) == Dict("grid_rank" => 2)
        grid_size = request((fn = "get_grid_size", grid = 4))["grid_size"]
        @test grid_size == 50063
        msg = (fn = "get_grid_x", grid = 4, x = fill(0.0, grid_size))
        @test request(msg)["grid_x"][1:3] ≈
              [6.826666666666673, 6.810000000000006, 6.81833333333334]
        msg = (fn = "get_grid_y", grid = 4, y = fill(0.0, grid_size))
        @test request(msg)["grid_y"][1:3] ≈
              [47.8175, 47.825833333333335, 47.825833333333335]
        @test request((fn = "get_grid_node_count", grid = 0)) ==
              Dict("grid_node_count" => 2)
        @test request((fn = "get_grid_edge_count", grid = 3))["grid_edge_count"] == 5808
        msg = (fn = "get_grid_edge_nodes", grid = 3, edge_nodes = fill(0, 2 * 5808))
        @test request(msg)["grid_edge_nodes"][1:6] == [1, 5, 2, 1, 3, 2]
    end

    @testset "model states and finalize functions" begin
        @test request((fn = "load_state",)) == Dict("status" => "OK")
        @test request((fn = "save_state",)) == Dict("status" => "OK")
        @test request((fn = "finalize",)) == Dict("status" => "OK")
    end

    @testset "Error handling and shutdown" begin
        msg = (fn = "initialize", config_file = joinpath(@__DIR__, "not_existing.toml"))
        @test request(msg)["status"] == "ERROR"
        @test split(request(msg)["error"], "\n")[1] == "Wflow function `initialize` failed"
        @test request((fn = "not_existing_function",)) == Dict(
            "status" => "ERROR",
            "error" => "Received invalid Wflow function: `not_existing_function`",
        )
        @test request((fn = "initialize",)) == Dict(
            "status" => "ERROR",
            "error" => "At least one required argument name (`config_file`) not available for function: `initialize`",
        )
        @test request((fn = "shutdown",)) == Dict("status" => "OK")
    end
end
