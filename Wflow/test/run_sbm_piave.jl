
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
        227.47934199200185,
        196.49178200511778,
        266.32968808285113,
        153.38666768554123,
        155.21478051544884,
        142.97628915389365,
        114.48721956751163,
        100.11686353860816,
        166.9766485414521,
        102.53586232954153,
    ]
    @test q_[idx] ≈ [
        230.5078482226877,
        203.14141502102652,
        274.3732493238211,
        161.41889420851422,
        166.82296005197907,
        149.48540584790365,
        124.68222889813899,
        107.93685211256276,
        169.86373756264732,
        107.16962464537252,
    ]
    @test riv_storage_demand[idx] ≈ [
        61169.05717169386,
        55765.22842825826,
        61155.7904491885,
        47857.27669499157,
        50469.24129692829,
        42471.2920166788,
        40579.945403144004,
        38334.18412908786,
        43355.320911765586,
        38430.55441951444,
    ]
    @test riv_storage[idx] ≈ [
        61787.60824041945,
        56352.115952484186,
        61377.16946584834,
        48342.337925797925,
        51107.23357577086,
        43054.500577895735,
        40932.83107946221,
        38631.246623249586,
        43240.25971175695,
        38359.735371238276,
    ]
    @test ssf_storage_demand[idx] ≈ [
        148937.39757703582,
        145104.5143277988,
        144075.95245385382,
        140693.30157751308,
        138518.72018740955,
        135438.57593621884,
        132595.94797779017,
        129797.53965593802,
        129807.77506057269,
        127817.066350789,
    ]
    @test ssf_storage[idx] ≈ [
        148925.32279961562,
        145076.77320164008,
        144018.910769367,
        140632.8033830338,
        138479.92196316362,
        135413.75092316957,
        132612.4203760569,
        129849.67878422113,
        129864.60531901068,
        127864.04770631882,
    ]
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
Wflow.run_timestep!(model)

(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    model.land.allocation.variables
(; soil) = model.land
(; river_flow) = model.routing
(; reservoir) = river_flow.boundary_conditions

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
    @test reservoir.variables.waterlevel ≈
          [29.25906414623379, 32.68607771649562, 39.970184252221905]
    @test reservoir.variables.storage ≈ [1.8959873566759476e8, 4.28e7, 7.16e7]
    @test reservoir.variables.outflow_av ≈
          [4.839239564238508, 9.691599149465155, 57.17487754699672]
    @test soil.variables.exfiltsatwater[[937, 939, 979, 1020, 1158]] ≈ [
        2.829877544764015,
        5.522561610629175,
        3.757838384021166,
        5.683601788725743,
        14.167786752295916,
    ]
    @test maximum(soil.variables.exfiltsatwater) ≈ 229.46901533940505
    @test mean(river_flow.variables.q_av) ≈ 60.222209714658796
    @test maximum(river_flow.variables.q_av) ≈ 235.54856117510604
end

Wflow.run_timestep!(model)

@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1592.1358888849306
    @test sum(act_groundwater_abst) ≈ 340.0140096928211
    @test paddy.variables.h[[25, 42, 45]] ≈
          [39.22731609310245, 48.04240591301383, 28.969111802094687]
    @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
          [4.2274274357998785, 0.71188497017256, 4.934431434115837]
    @test reservoir.variables.waterlevel ≈
          [29.250852997099564, 32.686077716495625, 39.97018425222191]
    @test reservoir.variables.storage ≈ [1.895455274212049e8, 4.28e7, 7.16e7]
    @test reservoir.variables.outflow_av ≈
          [4.841273534298489, 9.227831920778037, 52.615872475211354]
    @test soil.variables.exfiltsatwater[[937, 939, 979, 1020, 1158]] ≈ [
        3.019110189955034,
        6.058978072005602,
        4.279807292644558,
        5.994632869660209,
        14.56882001468644,
    ]
    @test maximum(soil.variables.exfiltsatwater) ≈ 211.03629198471612
    @test mean(river_flow.variables.q_av) ≈ 55.861124247540246
    @test maximum(river_flow.variables.q_av) ≈ 227.04275209379193
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
    @test reservoir.boundary_conditions.inflow[1] ≈ 5.703462938000151
    @test reservoir.variables.storage[1] ≈ 1.895455274212049e8
    @test reservoir.variables.outflow_av[1] ≈ 4.841273534298489
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
config.input.cyclic["reservoir_water__external_inflow_volume_flow_rate"] = "reservoir_inflow"
model = Wflow.Model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "piave: reservoir with cyclic external negative inflow" begin
    (; reservoir) = model.routing.river_flow.boundary_conditions
    @test reservoir.boundary_conditions.external_inflow[1] == -3.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[1] ≈ 3.0
    @test reservoir.boundary_conditions.inflow[1] ≈ 2.7034624582232096
    @test reservoir.variables.storage[1] ≈ 1.8902961043739313e8
    @test reservoir.variables.outflow_av[1] ≈ 4.821498186228839
end

# test use of observed reservoir outflow (cyclic) 
tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
config = Wflow.Config(tomlpath)
config.input.cyclic["reservoir_water__outgoing_observed_volume_flow_rate"] = "reservoir_outflow"
config.logging.loglevel = "debug"
config.logging.path_log = "log_sbm_piave_debug.txt"
config.time.endtime = DateTime(2010, 7, 3)
model = Wflow.Model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)
Wflow.close_files(model; delete_output = false)

@testset "piave: reservoir with observed outflow (cyclic)" begin
    (; reservoir) = model.routing.river_flow.boundary_conditions
    @test reservoir.boundary_conditions.external_inflow[1] == 0.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[1] ≈ 0.0
    @test reservoir.boundary_conditions.inflow[1] ≈ 5.6829979459011115
    @test reservoir.variables.storage[1] ≈ 1.9010714959114635e8
    @test reservoir.variables.outflow_av[1] ≈ 3.0
    @test reservoir.variables.outflow[1] ≈ 3.0
end

# test debug message using observed outflow for two timesteps
@testset "piave: log debug message using observed reservoir outflow" begin
    tomlpath_debug = joinpath(@__DIR__, "sbm_piave_config-debug.toml")
    open(tomlpath_debug, "w") do io
        TOML.print(io, Wflow.to_dict(config))
    end
    Wflow.run(tomlpath_debug; silent = true)
    rm(tomlpath_debug)
    path_log = Wflow.output_path(config, "log_sbm_piave_debug.txt")
    lines = readlines(path_log)
    msg = "┌ Debug: Observed outflow is used for reservoir location ids [169986]"
    @test count(contains(line, msg) for line in lines) == 2
end

@testset "water balance piave water demand" begin
    tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; land_water_balance, routing) = model.mass_balance
    (; overland_water_balance, river_water_balance, subsurface_water_balance) = routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 1e-9, land_water_balance.error)
        @test all(re -> abs(re) < 1e-9, land_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, overland_water_balance.error)
        @test all(re -> abs(re) < 1e-9, overland_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, land_water_balance.error)
        @test all(re -> abs(re) < 1e-9, land_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, routing.overland_water_balance.error)
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error)
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end
