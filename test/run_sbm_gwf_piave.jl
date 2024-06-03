tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.initialize_sbm_gwf_model(config)
model = Wflow.run_timestep(model)
sbm = model.vertical

@testset "piave water demand and allocation first timestep" begin
    sum_total_alloc = sum(sbm.waterallocation.total_alloc)
    @test sum(sbm.waterallocation.irri_alloc) + sum(sbm.waterallocation.nonirri_alloc) ≈
          sum_total_alloc
    @test sum(sbm.waterallocation.surfacewater_alloc) ≈ 1016.4268126167575f0
    @test sum(sbm.waterallocation.act_groundwater_abst) ≈ 182.2057678312209f0
    @test sbm.paddy.h[[45, 76, 296]] ≈
          [33.55659436283413f0, 44.11663357735189f0, 35.232731550849486f0]
    @test sbm.paddy.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test sbm.paddy.demand_gross[[45, 76, 296]] ≈ [0.0, 0.0, 0.0]
    @test sbm.nonpaddy.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test sbm.nonpaddy.demand_gross[[10, 33, 1293]] ≈ [3.3014913197447964f0, 0.0f0, 0.0f0]
    @test sbm.industry.demand_gross[[1, end]] ≈ [0.2105557769536972f0, 0.0485190823674202f0]
    @test sbm.industry.demand_net[[1, end]] ≈
          [0.05265098437666893f0, 0.012132546864449978f0]
    @test sbm.industry.returnflow[[1, end]] ≈ [0.15790479257702827f0, 0.03638653550297022f0]
    @test sbm.livestock.demand_gross[[1, end]] ≈
          [9.896758274408057f-5, 6.352497439365834f-5]
    @test sbm.livestock.demand_net[[1, end]] ≈ [9.896758274408057f-5, 6.352497439365834f-5]
    @test sbm.livestock.returnflow[[1, end]] ≈ [0.0f0, 0.0f0]
    @test sbm.domestic.demand_gross[[1, end]] ≈ [0.5389957427978516f0, 0.0f0]
    @test sbm.domestic.demand_net[[1, end]] ≈ [0.33949509263038635f0, 0.0f0]
    @test sbm.domestic.returnflow[[1, end]] ≈ [0.1995004952035704f0, 0.0f0]
end

model = Wflow.run_timestep(model)
sbm = model.vertical

@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(sbm.waterallocation.total_alloc)
    @test sum(sbm.waterallocation.irri_alloc) + sum(sbm.waterallocation.nonirri_alloc) ≈
          sum_total_alloc
    @test sum(sbm.waterallocation.surfacewater_alloc) ≈ 1040.0177693609326f0
    @test sum(sbm.waterallocation.act_groundwater_abst) ≈ 183.93121472382694f0
    @test sbm.paddy.h[[45, 76, 296]] ≈
          [28.002849362721747f0, 25.880301981405225f0, 29.880048223290544f0]
    @test sbm.paddy.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test sbm.paddy.demand_gross[[45, 76, 296]] ≈ [0.0f0, 0.0f0, 0.0f0]
    @test sbm.nonpaddy.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test sbm.nonpaddy.demand_gross[[10, 33, 1293]] ≈ [0.0f0, 0.0f0, 2.1321618457666403f0]
end

Wflow.close_files(model, delete_output = false)
