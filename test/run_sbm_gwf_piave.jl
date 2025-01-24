tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
Wflow.run_timestep!(model)

(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    model.land.allocation.variables
@testset "piave water demand and allocation first timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1016.4268126167575f0
    @test sum(act_groundwater_abst) ≈ 182.2057678312209f0
    @test paddy.variables.h[[45, 76, 296]] ≈
          [33.55659436283413f0, 44.11663357735189f0, 35.232731550849486f0]
    @test paddy.parameters.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[45, 76, 296]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[10, 33, 1293]] ≈
          [3.3014913197447964f0, 0.0f0, 0.0f0]
    @test industry.demand.demand_gross[[1, end]] ≈
          [0.2105557769536972f0, 0.0485190823674202f0]
    @test industry.demand.demand_net[[1, end]] ≈
          [0.05265098437666893f0, 0.012132546864449978f0]
    @test industry.variables.returnflow[[1, end]] ≈
          [0.15790479257702827f0, 0.03638653550297022f0]
    @test livestock.demand.demand_gross[[1, end]] ≈
          [9.896758274408057f-5, 6.352497439365834f-5]
    @test livestock.demand.demand_net[[1, end]] ≈
          [9.896758274408057f-5, 6.352497439365834f-5]
    @test livestock.variables.returnflow[[1, end]] ≈ [0.0f0, 0.0f0]
    @test domestic.demand.demand_gross[[1, end]] ≈ [0.5389957427978516f0, 0.0f0]
    @test domestic.demand.demand_net[[1, end]] ≈ [0.33949509263038635f0, 0.0f0]
    @test domestic.variables.returnflow[[1, end]] ≈ [0.1995004952035704f0, 0.0f0]
end

Wflow.run_timestep!(model)
(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    model.land.allocation.variables

@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1057.2573971583527f0
    @test sum(act_groundwater_abst) ≈ 189.12436534660375f0
    @test paddy.variables.h[[45, 76, 296]] ≈
          [28.197082339552537f0, 25.873022895247782f0, 30.066801639786547f0]
    @test paddy.parameters.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[45, 76, 296]] ≈ [0.0f0, 0.0f0, 0.0f0]
    @test nonpaddy.parameters.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[10, 33, 1293]] ≈
          [4.059144161330735f0, 0.0f0, 1.9399078662788196f0]
end

Wflow.close_files(model; delete_output = false)
