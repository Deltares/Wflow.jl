
tomlpath = joinpath(@__DIR__, "sbm_config.toml")

@testset "bmi" begin

    model = Wflow.BMI.initialize(tomlpath)

    @testset "initialization and time functions" begin 
        @test Wflow.BMI.get_time_units() == "seconds since 1970-01-01T00:00:00"
        @test Wflow.BMI.get_time_step(model) == 86400.0
        @test Wflow.BMI.get_start_time(model) == 9.466848e8
        @test Wflow.BMI.get_current_time(model) == 9.466848e8
        @test Wflow.BMI.get_end_time(model) == 9.493632e8
    end

    @testset "model information functions" begin
        @test Wflow.BMI.get_component_name(model) == "sbm"
        @test Wflow.BMI.get_input_item_count(model) == 179
        @test Wflow.BMI.get_output_item_count(model) == 179
        @test Wflow.BMI.get_input_var_names(model)[[1,5,120,179]] ==
            ["vertical.Δt",
            "vertical.n_unsatlayers",
            "lateral.land.sl",
            "lateral.river.reservoir.evaporation",
            ]
    end

    @testset "variable information functions" begin
        @test Wflow.BMI.get_var_grid(model, "vertical.θₛ") == 4
        @test Wflow.BMI.get_var_grid(model, "lateral.river.h") == 2
        @test Wflow.BMI.get_var_grid(model, "lateral.river.reservoir.inflow") == 0
        @test_logs (:warn, "Model variable lateral.river.lake.volume not listed as input of " *
        "output variable") Wflow.BMI.get_var_grid(model, "lateral.river.lake.volume")
        @test Wflow.BMI.get_var_type(model, "lateral.river.reservoir.inflow") == 
        "Array{Float64,1}"
        @test Wflow.BMI.get_var_type(model, "vertical.n") == "Int64"
        @test Wflow.BMI.get_var_units(model, "vertical.θₛ") == "mm mm-1"
        @test Wflow.BMI.get_var_itemsize(model,"lateral.subsurface.ssf") == 8
        @test Wflow.BMI.get_var_nbytes(model,"vertical.n") == 8
        @test Wflow.BMI.get_var_nbytes(model,"vertical.xl") == 400560
        @test Wflow.BMI.get_var_nbytes(model,"lateral.river.q") == 45240
        @test Wflow.BMI.get_var_location(model,"lateral.river.q") == "node"
    end

    model = Wflow.BMI.update(model)

    @testset "update and get and set functions" begin
        @test Wflow.BMI.get_current_time(model) == 9.467712e8
        @test mean(Wflow.BMI.get_value(model, "vertical.zi")) ≈ 276.2528389914268
        @test Wflow.BMI.get_value_at_indices(model, "lateral.river.q", [1,100,5617]) ≈ 
        [0.0019494398840044726,0.0552464309192915,2.662842827356121]
        Wflow.BMI.set_value(model, "vertical.zi" , fill(300.0,50070))
        @test mean(Wflow.BMI.get_value(model, "vertical.zi")) == 300.0
        Wflow.BMI.set_value_at_indices(model, "vertical.zi", [1], [250.0])
        @test Wflow.BMI.get_value_at_indices(model, "vertical.zi", [1,2]) == [250.0, 300.0]
    end

    @testset "model grid functions" begin
        @test Wflow.BMI.get_grid_type(0) == "unstructured"
        @test Wflow.BMI.get_grid_rank(0) == 2
        @test Wflow.BMI.get_grid_size(model,0) == 2
        @test Wflow.BMI.get_grid_size(model,2) == 5655
        @test Wflow.BMI.get_grid_size(model,4) == 50070
        @test Wflow.BMI.get_grid_shape(model,4) == 50070
        @test minimum(Wflow.BMI.get_grid_x(model,4)) ≈ 5.4291666666666725
        @test maximum(Wflow.BMI.get_grid_x(model,4)) ≈ 7.845833333333333
        @test Wflow.BMI.get_grid_x(model,0) ≈ [5.920833333333338, 5.7625000000000055]
        @test Wflow.BMI.get_grid_y(model,0) ≈ [49.9125, 48.920834] 
        @test Wflow.BMI.get_grid_node_count(model, 0) == 2
    end

    @testset "update until and finalize" begin
        time = Wflow.BMI.get_current_time(model) + 2 * Wflow.BMI.get_time_step(model)
        @test_logs (:info, "update model until 9.46944e8" )
        model = Wflow.BMI.update_until(model, time)
        Wflow.BMI.finalize(model)
    end

end