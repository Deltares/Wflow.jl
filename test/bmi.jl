
tomlpath = joinpath(@__DIR__, "sbm_config.toml")


@testset "BMI" begin

    @testset "BMI functions" begin

        model = BMI.initialize(Wflow.Model, tomlpath)

        @testset "initialization and time functions" begin
            @test BMI.get_time_units(model) == "s"
            @test BMI.get_time_step(model) == 86400.0
            @test BMI.get_start_time(model) == 0.0
            @test BMI.get_current_time(model) == 0.0
            @test BMI.get_end_time(model) == 31 * 86400.0
            model.config.endtime = "2000-02-01T00:00:00"
            @test BMI.get_end_time(model) == 31 * 86400.0
        end

        @testset "model information functions" begin
            @test BMI.get_component_name(model) == "sbm"
            @test BMI.get_input_item_count(model) == 193
            @test BMI.get_output_item_count(model) == 193
            to_check = [
                "vertical.nlayers",
                "vertical.θᵣ",
                "lateral.river.q",
                "lateral.river.reservoir.outflow",
            ]
            retrieved_vars = BMI.get_input_var_names(model)
            @test all(x -> x in retrieved_vars, to_check)
            retrieved_vars = BMI.get_output_var_names(model)
            @test all(x -> x in retrieved_vars, to_check)
        end

        @testset "variable information functions" begin
            @test BMI.get_var_grid(model, "vertical.θₛ") == 6
            @test BMI.get_var_grid(model, "lateral.river.h") == 3
            @test BMI.get_var_grid(model, "lateral.river.reservoir.inflow") == 0
            @test_throws ErrorException BMI.get_var_grid(model, "lateral.river.lake.volume")
            @test BMI.get_var_type(model, "lateral.river.reservoir.inflow") == "$Float"
            @test BMI.get_var_units(model, "vertical.θₛ") == "-"
            @test BMI.get_var_itemsize(model, "lateral.subsurface.ssf") == sizeof(Float)
            @test BMI.get_var_nbytes(model, "lateral.river.q") ==
                  length(model.lateral.river.q) * sizeof(Float)
            @test BMI.get_var_location(model, "lateral.river.q") == "node"
            @test_throws ErrorException(
                "lateral.land.alpha_pow not listed as variable for BMI exchange",
            ) BMI.get_var_itemsize(model, "lateral.land.alpha_pow")
        end

        model = BMI.update(model)

        @testset "update and get and set functions" begin
            @test BMI.get_current_time(model) == 86400.0
            @test_throws ErrorException BMI.get_value_ptr(model, "vertical.")
            dest = zeros(Float, size(model.vertical.zi))
            BMI.get_value(model, "vertical.zi", dest)
            @test mean(dest) ≈ 276.16325589542333
            @test BMI.get_value_at_indices(
                model,
                "vertical.vwc[1]",
                zeros(Float, 3),
                [1, 2, 3],
            ) ≈ getindex.(model.vertical.vwc, 1)[1:3]
            BMI.set_value_at_indices(
                model,
                "vertical.vwc[2]",
                [1, 2, 3],
                [0.10, 0.15, 0.20],
            ) ≈ getindex.(model.vertical.vwc, 2)[1:3]
            @test BMI.get_value_at_indices(
                model,
                "lateral.river.q",
                zeros(Float, 3),
                [1, 100, 5617],
            ) ≈ [0.623325399343309, 5.227139951657074, 0.02794287432778194]
            BMI.set_value(model, "vertical.zi", fill(300.0, length(model.vertical.zi)))
            @test mean(
                BMI.get_value(model, "vertical.zi", zeros(Float, size(model.vertical.zi))),
            ) == 300.0
            BMI.set_value_at_indices(model, "vertical.zi", [1], [250.0])
            @test BMI.get_value_at_indices(model, "vertical.zi", zeros(Float, 2), [1, 2]) ==
                  [250.0, 300.0]
        end

        @testset "model grid functions" begin
            @test BMI.get_grid_type(model, 0) == "points"
            @test BMI.get_grid_type(model, 2) == "points"
            @test BMI.get_grid_type(model, 6) == "unstructured"
            @test_throws ErrorException BMI.get_grid_rank(model, 8)
            @test BMI.get_grid_rank(model, 0) == 2
            @test BMI.get_grid_rank(model, 6) == 2
            @test_throws ErrorException BMI.get_grid_rank(model, 8)
            @test BMI.get_grid_node_count(model, 0) == 2
            @test BMI.get_grid_node_count(model, 3) == 5809
            @test BMI.get_grid_node_count(model, 4) == 50063
            @test BMI.get_grid_node_count(model, 5) == 50063
            @test BMI.get_grid_size(model, 0) == 2
            @test BMI.get_grid_size(model, 3) == 5809
            @test BMI.get_grid_size(model, 4) == 50063
            @test BMI.get_grid_size(model, 5) == 50063
            @test minimum(BMI.get_grid_x(model, 5, zeros(Float, 50063))) ≈
                  5.426666666666667f0
            @test maximum(BMI.get_grid_x(model, 5, zeros(Float, 50063))) ≈
                  7.843333333333344f0
            @test BMI.get_grid_x(model, 0, zeros(Float, 2)) ≈
                  [5.760000000000002f0, 5.918333333333336f0]
            @test BMI.get_grid_y(model, 0, zeros(Float, 2)) ≈
                  [48.92583333333333f0, 49.909166666666664f0]
            @test BMI.get_grid_node_count(model, 0) == 2
            @test BMI.get_grid_edge_count(model, 3) == 5808
            @test BMI.get_grid_edge_nodes(model, 3, fill(0, 2 * 5808))[1:6] ==
                  [1, 5, 2, 1, 3, 2]
        end

        @testset "update until and finalize" begin
            time = BMI.get_current_time(model) + 2 * BMI.get_time_step(model)
            model = BMI.update_until(model, time)
            @test model.clock.iteration == 3
            time_off = BMI.get_current_time(model) + 1 * BMI.get_time_step(model) + 1e-06
            @test_throws ErrorException model = BMI.update_until(model, time_off)
            @test_throws ErrorException model =
                BMI.update_until(model, time - BMI.get_time_step(model))
            BMI.finalize(model)
        end

    end

    @testset "BMI grid edges" begin
        tomlpath = joinpath(@__DIR__, "sbm_swf_config.toml")
        model = BMI.initialize(Wflow.Model, tomlpath)
        @test BMI.get_var_grid(model, "lateral.land.qx") == 4
        @test BMI.get_var_grid(model, "lateral.land.qy") == 5
        @test BMI.get_grid_edge_count(model, 4) == 50063
        @test BMI.get_grid_edge_count(model, 5) == 50063
        @test BMI.get_grid_edge_nodes(model, 4, fill(0, 2 * 50063))[1:4] == [1, -999, 2, 3]
        @test BMI.get_grid_edge_nodes(model, 5, fill(0, 2 * 50063))[1:4] == [1, 4, 2, 10]
        @test_logs (
            :warn,
            "edges are not provided for grid type 2 (variables are located at nodes)",
        ) BMI.get_grid_edge_nodes(model, 2, fill(0, 2 * 50063))
        @test_throws ErrorException BMI.get_grid_edge_nodes(model, 7, fill(0, 2 * 50063))
        BMI.finalize(model)
    end

    @testset "BMI run SBM in parts" begin
        tomlpath = joinpath(@__DIR__, "sbm_gw.toml")
        model = BMI.initialize(Wflow.Model, tomlpath)

        # update the recharge part of the SBM model
        model = BMI.update(model, run = "sbm_until_recharge")

        @testset "recharge part of SBM" begin
            sbm = model.vertical
            @test sbm.interception[1] ≈ 0.32734913737568716f0
            @test sbm.ustorelayerdepth[1][1] ≈ 0.0f0
            @test sbm.snow[1] ≈ 3.484789961176288f0
            @test sbm.recharge[5] ≈ -0.0f0
            @test sbm.zi[5] ≈ 300.0f0
        end

        # set zi and exfiltwater from external source (e.g. a groundwater model)
        BMI.set_value(
            model,
            "lateral.subsurface.zi",
            fill(0.25, BMI.get_grid_node_count(model, 6)),
        )
        BMI.set_value(
            model,
            "lateral.subsurface.exfiltwater",
            fill(1.0e-5, BMI.get_grid_node_count(model, 6)),
        )
        # update SBM after subsurface flow
        model = BMI.update(model, run = "sbm_after_subsurfaceflow")

        @testset "SBM after subsurface flow" begin
            sbm = model.vertical
            sub = model.lateral.subsurface
            @test sbm.interception[1] ≈ 0.32734913737568716f0
            @test sbm.ustorelayerdepth[1][1] ≈ 0.0f0
            @test sbm.snow[1] ≈ 3.484789961176288f0
            @test sbm.recharge[5] ≈ 0.0f0
            @test sbm.zi[5] ≈ 250.0f0
            @test sub.zi[5] ≈ 0.25f0
            @test sub.exfiltwater[1] ≈ 1.0f-5
            @test sub.ssf[1] ≈ 0.0f0
        end

        BMI.finalize(model)
    end

end

@testset "BMI extension functions" begin

    model = BMI.initialize(Wflow.Model, tomlpath)
    @test Wflow.get_start_unix_time(model) == 9.466848e8
    satwaterdepth = mean(model.vertical.satwaterdepth)
    model.config.model.reinit = false
    model = Wflow.load_state(model)
    @test satwaterdepth ≠ mean(model.vertical.satwaterdepth)
    @test_logs (
        :info,
        "Write output states to netCDF file `$(model.writer.state_nc_path)`.",
    ) Wflow.save_state(model)
    @test !isopen(model.writer.state_dataset)
end
