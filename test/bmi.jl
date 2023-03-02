
tomlpath = joinpath(@__DIR__, "sbm_config.toml")


@testset "BMI" begin

    @testset "BMI functions" begin

        model = BMI.initialize(Wflow.Model, tomlpath)

        @testset "initialization and time functions" begin
            @test BMI.get_time_units(model) == "s"
            @test BMI.get_time_step(model) == 86400.0
            @test BMI.get_start_time(model) == 9.467712e8
            @test BMI.get_current_time(model) == 9.467712e8
            @test BMI.get_end_time(model) == 9.493632e8
        end

        @testset "model information functions" begin
            @test BMI.get_component_name(model) == "sbm"
            @test BMI.get_input_item_count(model) == 186
            @test BMI.get_output_item_count(model) == 186
            @test BMI.get_input_var_names(model)[[1, 5, 120, 186]] == [
                "vertical.Δt",
                "vertical.n_unsatlayers",
                "lateral.land.q_av",
                "lateral.river.reservoir.evaporation",
            ]
        end

        @testset "variable information functions" begin
            @test BMI.get_var_grid(model, "vertical.θₛ") == 4
            @test BMI.get_var_grid(model, "lateral.river.h") == 2
            @test BMI.get_var_grid(model, "lateral.river.reservoir.inflow") == 0
            @test_throws ErrorException BMI.get_var_grid(model, "lateral.river.lake.volume")
            # Vector{Float64} printing on Julia 1.6+
            @test BMI.get_var_type(model, "lateral.river.reservoir.inflow") in
                  ["Array{$Float,1}", "Vector{$Float}"]
            @test BMI.get_var_type(model, "vertical.n") == "Int64"
            @test BMI.get_var_units(model, "vertical.θₛ") == "-"
            @test BMI.get_var_itemsize(model, "lateral.subsurface.ssf") == sizeof(Float)
            @test BMI.get_var_nbytes(model, "vertical.n") == 8
            @test BMI.get_var_nbytes(model, "lateral.river.q") ==
                  length(model.lateral.river.q) * sizeof(Float)
            @test BMI.get_var_location(model, "lateral.river.q") == "node"
        end

        model = BMI.update(model)

        @testset "update and get and set functions" begin
            @test BMI.get_current_time(model) == 9.468576e8
            @test mean(BMI.get_value(model, "vertical.zi")) ≈ 276.3767651555451
            @test BMI.get_value_at_indices(model, "lateral.river.q", [1, 100, 5617]) ≈
                  [0.6211503865184697, 5.219305686635002, 0.026163746306482282]
            BMI.set_value(model, "vertical.zi", fill(300.0, length(model.vertical.zi)))
            @test mean(BMI.get_value(model, "vertical.zi")) == 300.0
            BMI.set_value_at_indices(model, "vertical.zi", [1], [250.0])
            @test BMI.get_value_at_indices(model, "vertical.zi", [1, 2]) == [250.0, 300.0]
        end

        @testset "model grid functions" begin
            @test BMI.get_grid_type(model, 0) == "unstructured"
            @test BMI.get_grid_rank(model, 0) == 2
            @test BMI.get_grid_size(model, 0) == 2
            @test BMI.get_grid_size(model, 2) == 5809
            @test BMI.get_grid_size(model, 4) == 50063
            @test BMI.get_grid_shape(model, 4) == 50063
            @test minimum(BMI.get_grid_x(model, 4)) ≈ 5.426666666666667f0
            @test maximum(BMI.get_grid_x(model, 4)) ≈ 7.843333333333344f0
            @test BMI.get_grid_x(model, 0) ≈ [5.760000000000002f0, 5.918333333333336f0]
            @test BMI.get_grid_y(model, 0) ≈ [48.92583333333333f0, 49.909166666666664f0]
            @test BMI.get_grid_node_count(model, 0) == 2
        end

        @testset "update until and finalize" begin
            time = BMI.get_current_time(model) + 2 * BMI.get_time_step(model)
            @test_logs (:info, "update model until 9.470304e8")
            model = BMI.update_until(model, time)
            BMI.finalize(model)
        end

    end

    @testset "BMI run SBM in parts" begin
        tomlpath = joinpath(@__DIR__, "sbm_gw.toml")
        model = BMI.initialize(Wflow.Model, tomlpath)

        # update the recharge part of the SBM model
        model = BMI.update(model, run = "sbm_until_recharge")

        @testset "recharge part of SBM" begin
            sbm = model.vertical
            @test sbm.interception[1] ≈ 0.6299999952316284f0
            @test sbm.ustorelayerdepth[1][1] ≈ 0.0f0
            @test sbm.snow[1] ≈ 3.1912317735997524f0
            @test sbm.recharge[5] ≈ -0.0727941579914808f0
            @test sbm.zi[5] ≈ 300.0f0
        end

        # set zi and exfiltwater from external source (e.g. a groundwater model)
        BMI.set_value(
            model,
            "lateral.subsurface.zi",
            fill(0.25, BMI.get_grid_size(model, 4)),
        )
        BMI.set_value(
            model,
            "lateral.subsurface.exfiltwater",
            fill(1.0e-5, BMI.get_grid_size(model, 4)),
        )
        # update SBM after subsurface flow
        model = BMI.update(model, run = "sbm_after_subsurfaceflow")

        @testset "SBM after subsurface flow" begin
            sbm = model.vertical
            sub = model.lateral.subsurface
            @test sbm.interception[1] ≈ 0.6299999952316284
            @test sbm.ustorelayerdepth[1][1] ≈ 0.0f0
            @test sbm.snow[1] ≈ 3.1912317735997524f0
            @test sbm.recharge[5] ≈ -0.0727941579914808f0
            @test sbm.zi[5] ≈ 250.0f0
            @test sub.zi[5] ≈ 0.25f0
            @test sub.exfiltwater[1] ≈ 1.0f-5
            @test sub.ssf[1] ≈ 0.0f0
        end

        BMI.finalize(model)
    end
end
