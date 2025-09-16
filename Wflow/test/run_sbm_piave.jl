
function run_piave(model, steps)
    q = zeros(steps)
    ssf_storage = zeros(steps)
    riv_storage = zeros(steps)
    for i in 1:steps
        Wflow.run_timestep!(model)
        ssf_storage[i] = mean(model.routing.subsurface_flow.variables.storage)
        riv_storage[i] = mean(model.routing.river_flow.variables.storage)
        q[i] = model.routing.river_flow.variables.q_av[1]
    end
    return q, riv_storage, ssf_storage
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
q_demand, riv_storage_demand, ssf_storage_demand = run_piave(model, 30)
Wflow.close_files(model; delete_output = false)

tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
q_, riv_storage, ssf_storage = run_piave(model, 30)
Wflow.close_files(model; delete_output = false)

@testset "piave with and without water demand" begin
    idx = 1:3:28
    @test q_demand[idx] ≈ [
        227.48037896764353,
        201.78742424759918,
        277.79682561233307,
        167.55786362554602,
        173.52980435020592,
        159.07569567125717,
        134.13850871542112,
        119.48343702580726,
        190.2621015758869,
        116.56593315061009,
    ]
    @test q_[idx] ≈ [
        230.50893814063687,
        208.82584824313827,
        284.81955178174405,
        175.73953774993535,
        186.33517285165013,
        166.40698312714562,
        145.0657059985612,
        127.31965453809961,
        192.85958279090423,
        123.05180814279916,
    ]
    @test riv_storage_demand[idx] ≈ [
        61400.322078051824,
        57191.663377309276,
        63162.60070846872,
        51134.248213288876,
        54101.732744489476,
        46317.84990590213,
        44750.600812209865,
        42700.577549691414,
        46707.77542995108,
        42022.15314379535,
    ]
    @test riv_storage[idx] ≈ [
        62018.708818695835,
        57857.565369606586,
        63380.96358853941,
        51729.33665390077,
        54753.02294163998,
        46767.56026545779,
        45154.0230037507,
        43141.3937287654,
        46602.85708316674,
        41961.68169876138,
    ]
    @test ssf_storage_demand[idx] ≈ [
        293459.39181716816,
        288916.24142510165,
        287348.88496321393,
        282966.7568473801,
        279360.4763743409,
        274843.2933935162,
        270349.52203453914,
        265848.9812016826,
        264341.7860038618,
        260747.35406717093,
    ]
    @test ssf_storage[idx] ≈ [
        293447.3728218427,
        288838.18266323843,
        287221.5459536154,
        282807.3893207005,
        279198.0924175473,
        274667.0012648777,
        270203.22127512813,
        265724.33233628777,
        264130.6334396194,
        260497.0820298229,
    ]
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
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
    @test sum(surfacewater_alloc) ≈ 1662.1429336566434
    @test sum(act_groundwater_abst) ≈ 373.02432487790287
    @test paddy.variables.h[[25, 42, 45]] ≈
          [43.181590971577364, 51.20409409088053, 34.473513683211834]
    @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
          [4.0442025332978435, 0.7127330271840676, 4.491063429489396]
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
          [29.24205323934821, 32.68575584329306, 39.97006803315815]
    @test reservoir.variables.storage_av ≈
          [1.8948850499097648e8, 4.279957853085982e7, 7.159979181269434e7]
    @test reservoir.variables.outflow_av ≈
          [4.839262369375678, 9.7190485747002, 58.12081274927687]
    @test reservoir.variables.demandrelease[[2, 3]] ≈ [1.182999968528626, 7.902500152587073]
end

Wflow.run_timestep!(model)

@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1776.2367788369727
    @test sum(act_groundwater_abst) ≈ 403.9133089562002
    @test paddy.variables.h[[25, 42, 45]] ≈
          [39.22749740732533, 48.04243578765524, 28.970231709713346]
    @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
          [4.722084191267045, 0.7260732965706886, 5.452294842123767]
    @test reservoir.variables.waterlevel_av ≈
          [29.2485234495118, 32.68570780540624, 39.9700676103882]
    @test reservoir.variables.storage_av ≈
          [1.8953043195283672e8, 4.279951562880186e7, 7.159979105537169e7]
    @test reservoir.variables.outflow_av ≈
          [4.84140356355074, 9.32858384747847, 54.86508305737946]
    @test reservoir.variables.demandrelease[[2, 3]] ≈ [1.182999968528626, 7.902500152587073]
end

Wflow.close_files(model; delete_output = false)

# test cyclic reservoir external inflow
tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "piave: reservoir without external negative inflow" begin
    (; reservoir) = model.routing.river_flow.boundary_conditions
    @test reservoir.boundary_conditions.external_inflow[1] == 0.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[1] == 0.0
    @test reservoir.boundary_conditions.inflow[1] ≈ 5.738423327629313
    @test reservoir.variables.storage_av[1] ≈ 1.8953043195283672e8
    @test reservoir.variables.outflow_av[1] ≈ 4.84140356355074
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
config.input.cyclic.dict["reservoir_water_inflow~external__volume_flow_rate"] =
    Wflow.InputEntry(; standard_name = "reservoir_inflow")
model = Wflow.Model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "piave: reservoir with cyclic external negative inflow" begin
    (; reservoir) = model.routing.river_flow.boundary_conditions
    @test reservoir.boundary_conditions.external_inflow[1] == -3.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[1] ≈ 3.0
    @test reservoir.boundary_conditions.inflow[1] ≈ 2.7384228457941124
    @test reservoir.variables.storage_av[1] ≈ 1.891429435839712e8
    @test reservoir.variables.outflow_av[1] ≈ 4.82162803109911
end
Wflow.close_files(model; delete_output = false)
