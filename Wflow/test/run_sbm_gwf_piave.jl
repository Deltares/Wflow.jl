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
          [23.956777464091655, 32.68528833371174, 39.96881815487859]
    @test reservoir.variables.storage_av ≈
          [1.5523991796731403e8, 4.279896636165852e7, 7.159755286167501e7]
    @test reservoir.variables.outflow_av ≈
          [3.248031210948757, 5.04778640418958, 14.259912531748482]
    @test reservoir.variables.demandrelease[[2, 3]] ≈ [1.182999968528626, 7.902500152587073]
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
          [23.946156606576015, 32.68503511233224, 39.968913991701285]
    @test reservoir.variables.storage_av ≈
          [1.5517109481061247e8, 4.2798634787000515e7, 7.159772453755285e7]
    @test reservoir.variables.outflow_av ≈
          [3.2451518657160476, 4.654493902558098, 13.362835228238662]
    @test reservoir.variables.demandrelease[[2, 3]] ≈ [1.182999968528626, 7.902500152587073]
end

Wflow.close_files(model; delete_output = false)
