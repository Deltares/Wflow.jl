tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.initialize_sbm_gwf_model(config)
Wflow.run_timestep!(model)

(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    model.land.allocation.variables
@testset "piave water demand and allocation first timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1016.4268126167575
    @test sum(act_groundwater_abst) ≈ 182.2057678312209
    @test paddy.variables.h[[45, 76, 296]] ≈
          [33.55659436283413, 44.11663357735189, 35.232731550849486]
    @test paddy.parameters.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[45, 76, 296]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[10, 33, 1293]] ≈ [3.3014913197447964, 0.0, 0.0]
    @test industry.demand.demand_gross[[1, end]] ≈ [0.2105557769536972, 0.0485190823674202]
    @test industry.demand.demand_net[[1, end]] ≈ [0.05265098437666893, 0.012132546864449978]
    @test industry.variables.returnflow[[1, end]] ≈
          [0.15790479257702827, 0.03638653550297022]
    @test livestock.demand.demand_gross[[1, end]] ≈
          [9.896758274408057e-5, 6.352497439365834e-5]
    @test livestock.demand.demand_net[[1, end]] ≈
          [9.896758274408057e-5, 6.352497439365834e-5]
    @test livestock.variables.returnflow[[1, end]] ≈ [0.0, 0.0]
    @test domestic.demand.demand_gross[[1, end]] ≈ [0.5389957427978516, 0.0]
    @test domestic.demand.demand_net[[1, end]] ≈ [0.33949509263038635, 0.0]
    @test domestic.variables.returnflow[[1, end]] ≈ [0.1995006501674652, 0.0]
end

Wflow.run_timestep!(model)
(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    model.land.allocation.variables

@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1058.8052175705648
    @test sum(act_groundwater_abst) ≈ 189.5040823394481
    @test paddy.variables.h[[45, 76, 296]] ≈
          [28.197082339552537, 24.547545670400083, 30.066801639786547]
    @test paddy.parameters.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[45, 76, 296]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[10, 33, 1293]] ≈
          [4.059144161330735, 0.0, 1.9399078662788196]
end

Wflow.close_files(model; delete_output = false)
