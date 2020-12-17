
tomlpath = joinpath(@__DIR__, "sbm_config.toml")

@testset "BMI" begin

    model = BMI.initialize(Wflow.Model, tomlpath)

    @testset "initialization and time functions" begin 
        @test BMI.get_time_units(Wflow.Model) == "seconds since 1970-01-01T00:00:00"
        @test BMI.get_time_step(model) == 86400.0
        @test BMI.get_start_time(model) == 9.466848e8
        @test BMI.get_current_time(model) == 9.466848e8
        @test BMI.get_end_time(model) == 9.493632e8
    end

    @testset "model information functions" begin
        @test BMI.get_component_name(model) == "sbm"
        @test BMI.get_input_item_count(model) == 179
        @test BMI.get_output_item_count(model) == 179
        @test BMI.get_input_var_names(model)[[1,5,120,179]] ==
            ["vertical.Δt",
            "vertical.n_unsatlayers",
            "lateral.land.sl",
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
            ["Array{Float64,1}", "Vector{Float64}"]
        @test BMI.get_var_type(model, "vertical.n") == "Int64"
        @test BMI.get_var_units(model, "vertical.θₛ") == "mm mm-1"
        @test BMI.get_var_itemsize(model,"lateral.subsurface.ssf") == 8
        @test BMI.get_var_nbytes(model,"vertical.n") == 8
        @test BMI.get_var_nbytes(model,"vertical.xl") == 400560
        @test BMI.get_var_nbytes(model,"lateral.river.q") == 45240
        @test BMI.get_var_location(model,"lateral.river.q") == "node"
    end

    model = BMI.update(model)

    @testset "update and get and set functions" begin
        @test BMI.get_current_time(model) == 9.467712e8
        @test mean(BMI.get_value(model, "vertical.zi")) ≈ 276.2528389914268
        @test BMI.get_value_at_indices(model, "lateral.river.q", [1,100,5617]) ≈ 
        [0.0019494398840044726,0.0552464309192915,2.662842827356121]
        BMI.set_value(model, "vertical.zi" , fill(300.0,50070))
        @test mean(BMI.get_value(model, "vertical.zi")) == 300.0
        BMI.set_value_at_indices(model, "vertical.zi", [1], [250.0])
        @test BMI.get_value_at_indices(model, "vertical.zi", [1,2]) == [250.0, 300.0]
    end

    @testset "model grid functions" begin
        @test BMI.get_grid_type(model, 0) == "unstructured"
        @test BMI.get_grid_rank(model, 0) == 2
        @test BMI.get_grid_size(model, 0) == 2
        @test BMI.get_grid_size(model, 2) == 5655
        @test BMI.get_grid_size(model, 4) == 50070
        @test BMI.get_grid_shape(model, 4) == 50070
        @test minimum(BMI.get_grid_x(model, 4)) ≈ 5.4291666666666725
        @test maximum(BMI.get_grid_x(model, 4)) ≈ 7.845833333333333
        @test BMI.get_grid_x(model, 0) ≈ [5.920833333333338, 5.7625000000000055]
        @test BMI.get_grid_y(model, 0) ≈ [49.9125, 48.920834] 
        @test BMI.get_grid_node_count(model, 0) == 2
    end

    @testset "update until and finalize" begin
        time = BMI.get_current_time(model) + 2 * BMI.get_time_step(model)
        @test_logs (:info, "update model until 9.46944e8" )
        model = BMI.update_until(model, time)
        BMI.finalize(model)
    end

end