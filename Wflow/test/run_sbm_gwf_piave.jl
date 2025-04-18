tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
Wflow.run_timestep!(model)

(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    model.land.allocation.variables
(; lake, reservoir) = model.routing.river_flow.boundary_conditions

@testset "piave water demand and allocation first timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1206.2599733337572
    @test sum(act_groundwater_abst) ≈ 222.28168654037432
    @test paddy.variables.h[[25, 42, 45]] ≈
          [43.221297434296346, 1.6725004267537358, 34.33337318442007]
    @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 25.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈ [0.0, 0.0, 0.0074476682445751684]
    @test industry.demand.demand_gross[[1, end]] ≈ [0.2105557769536972, 0.0485190823674202]
    @test industry.demand.demand_net[[1, end]] ≈ [0.05265098437666893, 0.012132546864449978]
    @test industry.variables.returnflow[[1, end]] ≈
          [0.15790479257702827, 0.03638653550297022]
    @test livestock.demand.demand_gross[[1, end]] ≈
          [9.896758274408057e-5, 6.352497439365834e-5]
    @test livestock.demand.demand_net[[1, end]] ≈
          [9.896758274408057e-5, 6.352497439365834e-5]
    @test livestock.variables.returnflow[[1, end]] ≈ [0.0, 0.0]
    @test domestic.demand.demand_gross[[1, end]] ≈ [0.6012673377990723, 0.0]
    @test domestic.demand.demand_net[[1, end]] ≈ [0.3802947998046875, 0.0]
    @test domestic.variables.returnflow[[1, end]] ≈ [0.2209725379943848, 0.0]
    @test lake.variables.waterlevel_av ≈ [24.143207490889875]
    @test lake.variables.storage_av ≈ [1.5644798454096633e8]
    @test lake.variables.outflow_av ≈ [3.2987797843728366]
    @test reservoir.variables.storage_av ≈ [4.2798909526192345e7, 7.159752155889235e7]
    @test reservoir.variables.outflow_av ≈ [5.006045233380978, 14.032226692922634]
    @test reservoir.variables.demandrelease ≈ [1.182999968528626, 7.902500152587073]
end

Wflow.run_timestep!(model)

@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1252.8380255930115
    @test sum(act_groundwater_abst) ≈ 226.94288689148232
    @test paddy.variables.h[[25, 42, 45]] ≈ [39.29489908470131, 0.0, 28.928631484478558]
    @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 25.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
          [0.0, 3.9938486455693436, 4.6316316016074115]
    @test lake.variables.waterlevel_av ≈ [24.132227749090518]
    @test lake.variables.storage_av ≈ [1.5637683581410643e8]
    @test lake.variables.outflow_av ≈ [3.295780210974653]
    @test reservoir.variables.storage_av ≈ [4.279864320599316e7, 7.159769075666577e7]
    @test reservoir.variables.outflow_av ≈ [4.670430399032648, 13.21257202829881]
    @test reservoir.variables.demandrelease ≈ [1.182999968528626, 7.902500152587073]
end

Wflow.close_files(model; delete_output = false)
