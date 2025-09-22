tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
Wflow.run_timestep!(model)

(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    model.land.allocation.variables
(; reservoir) = model.routing.river_flow.boundary_conditions

@testset "piave water demand and allocation first timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1747.401299613164
    @test sum(act_groundwater_abst) ≈ 392.3593761085565
    @test paddy.variables.h[[25, 42, 45]] ≈ [42.968584878704704, 0.0, 33.2318007065323]
    @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 25.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈ [0.0, 0.0, 4.487607239115634]
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
    @test reservoir.variables.waterlevel_av ≈
          [23.95914400555985, 32.68535935572136, 39.969135230891744]
    @test reservoir.variables.storage_av ≈
          [1.5525525315602788e8, 4.27990593597248e7, 7.15981208511133e7]
    @test reservoir.variables.outflow_av ≈
          [3.248673046140208, 8.352196766583088, 29.02990124474297]
end

Wflow.run_timestep!(model)

@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1884.1479041244102
    @test sum(act_groundwater_abst) ≈ 426.57564117049156
    @test paddy.variables.h[[25, 42, 45]] ≈ [38.99648725170036, 0.0, 27.60970255170497]
    @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 25.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
          [0.0, 3.9965040974684207, 5.44810857188258]
    @test reservoir.variables.waterlevel_av ≈
          [23.95845092987531, 32.68554507549814, 39.96988968054929]
    @test reservoir.variables.storage_av ≈
          [1.552507620255922e8, 4.2799302546029136e7, 7.159947232337902e7]
    @test reservoir.variables.outflow_av ≈
          [3.2484850081729024, 9.328049956914716, 38.06870720301024]
end

Wflow.close_files(model; delete_output = false)
