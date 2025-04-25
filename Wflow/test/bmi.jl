
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
            model.config.time.endtime = "2000-02-01T00:00:00"
            @test BMI.get_end_time(model) == 31 * 86400.0
        end

        @testset "model information functions" begin
            @test BMI.get_component_name(model) == "sbm"
            @test BMI.get_input_item_count(model) == 7
            @test BMI.get_output_item_count(model) == 7
            to_check = [
                "river_water__volume_flow_rate",
                "soil_water_unsat-zone__depth",
                "soil_water__transpiration_volume_flux",
                "soil_layer~2_water_unsat-zone__depth",
            ]
            retrieved_vars = BMI.get_input_var_names(model)
            @test all(x -> x in retrieved_vars, to_check)
            retrieved_vars = BMI.get_output_var_names(model)
            @test all(x -> x in retrieved_vars, to_check)
        end

        @testset "variable information functions" begin
            @test BMI.get_var_grid(model, "soil_water__infiltration_volume_flux") == 6
            @test BMI.get_var_grid(model, "river_water__volume_flow_rate") == 3
            @test BMI.get_var_grid(model, "reservoir_water~outgoing__volume_flow_rate") == 0
            @test BMI.get_var_type(model, "reservoir_water~incoming__volume_flow_rate") ==
                  "Float64"
            @test BMI.get_var_units(model, "river_water__volume_flow_rate") == "m3 s-1"
            @test BMI.get_var_itemsize(model, "subsurface_water__volume_flow_rate") ==
                  sizeof(Float64)
            @test BMI.get_var_nbytes(
                model,
                "river_water__instantaneous_volume_flow_rate",
            ) == length(model.routing.river_flow.variables.q) * sizeof(Float64)
            @test BMI.get_var_location(model, "river_water__volume_flow_rate") == "node"
        end

        BMI.update(model)

        @testset "update and get and set functions" begin
            @test BMI.get_current_time(model) == 86400.0
            dest = zeros(Float64, size(model.land.soil.variables.zi))
            BMI.get_value(model, "soil_water_sat-zone_top__depth", dest)
            @test mean(dest) ≈ 276.1625022866973
            @test BMI.get_value_at_indices(
                model,
                "soil_layer~1_water__volume_fraction",
                zeros(Float64, 3),
                [1, 2, 3],
            ) ≈ getindex.(model.land.soil.variables.vwc, 1)[1:3]
            BMI.set_value_at_indices(
                model,
                "soil_layer~2_water__volume_fraction",
                [1, 2, 3],
                [0.10, 0.15, 0.20],
            ) ≈ getindex.(model.land.soil.variables.vwc, 2)[1:3]
            @test BMI.get_value_at_indices(
                model,
                "river_water__instantaneous_volume_flow_rate",
                zeros(Float64, 3),
                [1, 100, 5617],
            ) ≈ [0.6525634030110335, 7.493560511070567, 0.023197145979653312]
            BMI.set_value(
                model,
                "soil_water_sat-zone_top__depth",
                fill(300.0, length(model.land.soil.variables.zi)),
            )
            @test mean(
                BMI.get_value(
                    model,
                    "soil_water_sat-zone_top__depth",
                    zeros(Float64, size(model.land.soil.variables.zi)),
                ),
            ) == 300.0
            BMI.set_value_at_indices(model, "soil_water_sat-zone_top__depth", [1], [250.0])
            @test BMI.get_value_at_indices(
                model,
                "soil_water_sat-zone_top__depth",
                zeros(Float64, 2),
                [1, 2],
            ) == [250.0, 300.0]
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
            @test minimum(BMI.get_grid_x(model, 5, zeros(Float64, 50063))) ≈
                  5.426666666666667
            @test maximum(BMI.get_grid_x(model, 5, zeros(Float64, 50063))) ≈
                  7.843333333333344
            @test BMI.get_grid_x(model, 0, zeros(Float64, 2)) ≈
                  [5.760000000000002, 5.918333333333336]
            @test BMI.get_grid_y(model, 0, zeros(Float64, 2)) ≈
                  [48.92583333333333, 49.909166666666664]
            @test BMI.get_grid_node_count(model, 0) == 2
            @test BMI.get_grid_edge_count(model, 3) == 5808
            @test BMI.get_grid_edge_nodes(model, 3, fill(0, 2 * 5808))[1:6] ==
                  [1, 5, 2, 1, 3, 2]
        end

        @testset "update until and finalize" begin
            time = BMI.get_current_time(model) + 2 * BMI.get_time_step(model)
            BMI.update_until(model, time)
            @test model.clock.iteration == 3
            time_off = BMI.get_current_time(model) + 1 * BMI.get_time_step(model) + 1e-06
            @test_throws ErrorException BMI.update_until(model, time_off)
            @test_throws ErrorException BMI.update_until(
                model,
                time - BMI.get_time_step(model),
            )
            BMI.finalize(model)
        end
    end

    @testset "BMI grid edges and element type" begin
        tomlpath = joinpath(@__DIR__, "sbm_river-land-local-inertial_config.toml")
        model = BMI.initialize(Wflow.Model, tomlpath)
        @test BMI.get_var_grid(
            model,
            "land_surface_water__x_component_of_instantaneous_volume_flow_rate",
        ) == 4
        @test BMI.get_var_grid(
            model,
            "land_surface_water__y_component_of_instantaneous_volume_flow_rate",
        ) == 5
        @test BMI.get_grid_edge_count(model, 4) == 50063
        @test BMI.get_grid_edge_count(model, 5) == 50063
        @test_logs (
            :warn,
            "edges are not provided for grid type 2 (variables are located at nodes)",
        ) BMI.get_grid_edge_count(model, 2)
        @test_throws ErrorException BMI.get_grid_edge_count(model, 7)
        @test BMI.get_grid_edge_nodes(model, 4, fill(0, 2 * 50063))[1:4] == [1, -999, 2, 3]
        @test BMI.get_grid_edge_nodes(model, 5, fill(0, 2 * 50063))[1:4] == [1, 4, 2, 10]
        @test_logs (
            :warn,
            "edges are not provided for grid type 2 (variables are located at nodes)",
        ) BMI.get_grid_edge_nodes(model, 2, fill(0, 2 * 50063))
        @test_throws ErrorException BMI.get_grid_edge_nodes(model, 7, fill(0, 2 * 50063))
        @test BMI.get_var_location(model, "river_water__volume_flow_rate") == "edge"
        @test BMI.get_var_location(
            model,
            "land_surface_water__y_component_of_instantaneous_volume_flow_rate",
        ) == "edge"
        @test BMI.get_var_location(model, "river_water__depth") == "node"
        BMI.finalize(model)
    end

    @testset "BMI run SBM in parts" begin
        tomlpath = joinpath(@__DIR__, "sbm_gw.toml")
        model = BMI.initialize(Wflow.Model, tomlpath)

        # update the recharge part of the SBM model
        BMI.update(model; run = "sbm_until_recharge")

        @testset "recharge part of SBM" begin
            sbm = model.land
            @test sbm.interception.variables.interception_rate[1] ≈ 0.32734913737568716
            @test sbm.soil.variables.ustorelayerdepth[1][1] ≈ 0.0
            @test sbm.snow.variables.snow_storage[1] ≈ 3.4847899611762876
            @test sbm.soil.variables.recharge[5] ≈ 0.0
            @test sbm.soil.variables.zi[5] ≈ 300.0
        end

        # set zi and exfiltwater from external source (e.g. a groundwater model)
        BMI.set_value(
            model,
            "subsurface_water_sat-zone_top__depth",
            fill(0.25, BMI.get_grid_node_count(model, 6)),
        )
        BMI.set_value(
            model,
            "subsurface_water__exfiltration_volume_flux",
            fill(1.0e-5, BMI.get_grid_node_count(model, 6)),
        )
        # update SBM after subsurface flow
        BMI.update(model; run = "sbm_after_subsurfaceflow")

        @testset "SBM after subsurface flow" begin
            sbm = model.land
            sub = model.routing.subsurface_flow
            @test sbm.interception.variables.interception_rate[1] ≈ 0.32734913737568716
            @test sbm.soil.variables.ustorelayerdepth[1][1] ≈ 0.0
            @test sbm.snow.variables.snow_storage[1] ≈ 3.4847899611762876
            @test sbm.soil.variables.recharge[5] ≈ 0.0
            @test sbm.soil.variables.zi[5] ≈ 250.0
            @test sub.variables.zi[5] ≈ 0.25
            @test sub.variables.exfiltwater[1] ≈ 1.0e-5
            @test sub.variables.ssf[1] ≈ 0.0
        end

        BMI.finalize(model)
    end
end

@testset "BMI extension functions" begin
    model = BMI.initialize(Wflow.Model, tomlpath)
    @test Wflow.get_start_unix_time(model) == 9.466848e8
    satwaterdepth = mean(model.land.soil.variables.satwaterdepth)
    model.config.model.cold_start = false
    Wflow.load_state(model)
    @test satwaterdepth ≠ mean(model.land.soil.variables.satwaterdepth)
    @test_logs (
        :info,
        "Write output states to netCDF file `$(model.writer.state_nc_path)`.",
    ) Wflow.save_state(model)
    @test !isopen(model.writer.state_dataset)
end
