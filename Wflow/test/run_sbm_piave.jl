@testitem "Piave with and without water demand (sbm model)" begin
    include("testing_utils.jl")
    tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    q_demand, riv_storage_demand, ssf_storage_demand = run_piave(model, 30)
    Wflow.close_files(model; delete_output = false)

    tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    q_, riv_storage, ssf_storage = run_piave(model, 30)
    Wflow.close_files(model; delete_output = false)

    idx = 1:3:28
    @test q_demand[idx] ≈ [
        227.38262076140597,
        197.77826932854217,
        267.02912980712193,
        155.28759918931823,
        155.82292184228206,
        143.33475478513975,
        115.38260642480303,
        100.99725841800637,
        166.27201893591504,
        103.1045884698738,
    ]
    @test q_[idx] ≈ [
        230.50830652411622,
        204.71391716296807,
        274.93940253146957,
        162.9967927126573,
        167.73617405644782,
        150.31152775549913,
        125.55699317343567,
        108.6868994407236,
        169.1977308700292,
        108.08518587596184,
    ]
    @test riv_storage_demand[idx] ≈ [
        61262.15586688514,
        56031.22324367828,
        61373.56490076488,
        48073.00601232394,
        50574.8858099468,
        42730.438661308515,
        40760.706597408636,
        38426.54099045719,
        43466.61231827235,
        38558.90224961908,
    ]
    @test riv_storage[idx] ≈ [
        61902.31987942963,
        56612.07374434918,
        61588.27339754675,
        48604.22242542106,
        51199.60577284767,
        43244.11730927491,
        41097.23262856786,
        38768.294030445264,
        43354.39861956778,
        38512.97533866192,
    ]
    @test ssf_storage_demand[idx] ≈ [
        148920.45317987585,
        145035.32705296585,
        143981.91015223257,
        140540.61889193894,
        138341.04974891624,
        135225.98687730395,
        132350.44182452557,
        129525.55864517286,
        129529.58758586958,
        127504.02426822705,
    ]
    @test ssf_storage[idx] ≈ [
        148907.47777948846,
        145003.35842325471,
        143918.60784373066,
        140474.26677411146,
        138293.5185650071,
        135189.54430373877,
        132355.50247198046,
        129566.85803370144,
        129569.2396631266,
        127534.96295815968,
    ]
end

@testitem "Piave water demand and allocation (sbm model)" begin
    using Statistics: mean

    tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)

    (; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
    (; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
        model.land.allocation.variables
    (; soil) = model.land
    (; river_flow) = model.routing
    (; reservoir) = river_flow.boundary_conditions

    @testset "First timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
        @test sum(surfacewater_alloc) ≈ 1706.10764866567
        @test sum(act_groundwater_abst) ≈ 388.08837400239827
        @test paddy.variables.h[[25, 42, 45]] ≈
              [43.181590971577364, 51.20409409088053, 34.473513683211834]
        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              [4.235245721225455, 0.7341536180838294, 4.691319396612657]
        @test industry.demand.demand_gross[[1, end]] ≈
              [0.2105557769536972, 0.0485190823674202]
        @test industry.demand.demand_net[[1, end]] ≈
              [0.05265098437666893, 0.012132546864449978]
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
              [29.259144530899885, 32.68607771649562, 39.970184252221905]
        @test reservoir.variables.storage ≈ [1.8959925656023118e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_av ≈
              [4.839249448770597, 9.70890592273304, 57.64196655435958]
        @test soil.variables.exfiltsatwater[[937, 939, 979, 1020, 1158]] ≈ [
            2.9401361588857626,
            5.54461271947943,
            3.76716736049869,
            5.766529601726355,
            14.050058868842683,
        ]
        @test maximum(soil.variables.exfiltsatwater) ≈ 234.1526270577704
        @test mean(river_flow.variables.q_av) ≈ 60.35958015226597
        @test maximum(river_flow.variables.q_av) ≈ 235.44682290212157
    end

    Wflow.run_timestep!(model)

    @testset "Second timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
        @test sum(surfacewater_alloc) ≈ 1591.415961657279
        @test sum(act_groundwater_abst) ≈ 337.68291766185564
        @test paddy.variables.h[[25, 42, 45]] ≈
              [39.22741049173483, 48.04241913069636, 28.96957318918968]
        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              [4.442695372669294, 0.7341536180838293, 5.022845160752851]
        @test reservoir.variables.waterlevel ≈
              [29.25110298462435, 32.686077716495625, 39.970184252221905]
        @test reservoir.variables.storage ≈ [1.8954714734036595e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_av ≈
              [4.841325917420907, 9.278008188175482, 53.34402066559018]
        @test soil.variables.exfiltsatwater[[937, 939, 979, 1020, 1158]] ≈ [
            3.1186475152972895,
            6.073982628186992,
            4.288301138492518,
            6.054701391567255,
            14.57870735626147,
        ]
        @test maximum(soil.variables.exfiltsatwater) ≈ 215.88061474014935
        @test mean(river_flow.variables.q_av) ≈ 56.22490708883267
        @test maximum(river_flow.variables.q_av) ≈ 227.21143082238987
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "Piave water demand and allocation, switch off livestock, paddy and nonpaddy" begin
    tomlpath = joinpath(@__DIR__, "sbm_piave_nonirri-demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    (; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
    (; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
        model.land.allocation.variables
    @test typeof(paddy) == Wflow.NoIrrigationPaddy
    @test typeof(nonpaddy) == Wflow.NoIrrigationNonPaddy
    @test typeof(livestock) == Wflow.NoNonIrrigationDemand
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 824.8974139426691
    @test sum(act_groundwater_abst) ≈ 115.97004946871232
    @test industry.demand.demand_gross[[1, end]] ≈ [0.2105557769536972, 0.0485190823674202]
    @test industry.demand.demand_net[[1, end]] ≈ [0.05265098437666893, 0.012132546864449978]
    @test industry.variables.returnflow[[1, end]] ≈
          [0.15790479257702827, 0.03638653550297022]
    @test domestic.demand.demand_gross[[1, end]] ≈ [0.6012673377990723, 0.0]
    @test domestic.demand.demand_net[[1, end]] ≈ [0.3802947998046875, 0.0]
    @test domestic.variables.returnflow[[1, end]] ≈ [0.2209725379943848, 0.0]
end

@testitem "Piave activate river boundary (river subsurface exchange)" begin
    using Statistics: mean

    tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.river_subsurface_exchange_head_based__flag = true
    config.input.static["river_bottom__elevation"] = "zb_river"
    config.input.static["river_water__infiltration_conductance"] = "riverbed_cond"
    config.input.static["river_water__exfiltration_conductance"] = "riverbed_cond"
    model = Wflow.Model(config)
    (; subsurface_flow) = model.routing
    (; river, recharge) = subsurface_flow.boundary_conditions

    Wflow.run_timestep!(model)

    @testset "First timestep" begin
        @test subsurface_flow.variables.head[1] ≈ 1.5649019759969596
        @test mean(subsurface_flow.variables.head) ≈ 1106.4948788348809
        @test subsurface_flow.variables.zi[1] ≈ 0.05409810125066005
        @test subsurface_flow.parameters.top[1] - subsurface_flow.variables.zi[1] ==
              subsurface_flow.variables.head[1]
        @test river.variables.flux_av[1] ≈ 37872.287583718644
        @test mean(river.variables.flux_av) ≈ -39618.370151473195
        @test recharge.variables.rate[1] ≈ -0.0002922905062717429
        @test mean(recharge.variables.rate) ≈ 0.0009271689030318317
    end

    Wflow.run_timestep!(model)

    @testset "Second timestep" begin
        @test subsurface_flow.variables.head[1] ≈ 1.563147470676254
        @test mean(subsurface_flow.variables.head) ≈ 1106.4867249666993
        @test subsurface_flow.variables.zi[1] ≈ 0.05585260657136559
        @test subsurface_flow.parameters.top[1] - subsurface_flow.variables.zi[1] ==
              subsurface_flow.variables.head[1]
        @test river.variables.flux_av[1] ≈ 45231.725584511034
        @test river.variables.flux[1] == river.variables.flux_av[1]
        @test mean(river.variables.flux_av) ≈ 129.0006905680265
        @test recharge.variables.rate[1] ≈ -0.00021663702639498745
        @test mean(recharge.variables.rate) ≈ 0.0010751049777759111
    end
end

@testitem "Piave: reservoir without external negative inflow (sbm model)" begin
    # test cyclic reservoir external inflow
    tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    (; reservoir) = model.routing.river_flow.boundary_conditions
    @test reservoir.boundary_conditions.external_inflow[1] == 0.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[1] == 0.0
    @test reservoir.boundary_conditions.inflow[1] ≈ 5.716044819323686
    @test reservoir.variables.storage[1] ≈ 1.8954714734036595e8
    @test reservoir.variables.outflow_av[1] ≈ 4.841325917420907
end

@testitem "Piave: reservoir with cyclic external negative inflow (sbm model)" begin
    tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.input.cyclic["reservoir_water__external_inflow_volume_flow_rate"] = "reservoir_inflow"
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    (; reservoir) = model.routing.river_flow.boundary_conditions
    @test reservoir.boundary_conditions.external_inflow[1] == -3.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[1] ≈ 3.0
    @test reservoir.boundary_conditions.inflow[1] ≈ 2.7160443378899912
    @test reservoir.variables.storage[1] ≈ 1.8903123036901748e8
    @test reservoir.variables.outflow_av[1] ≈ 4.8215504874860775
end

@testitem "Piave: reservoir with observed (cyclic) outflow (sbm model)" begin
    using Dates: DateTime
    # test use of observed reservoir outflow (cyclic)
    tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.input.cyclic["reservoir_water__outgoing_observed_volume_flow_rate"] = "reservoir_outflow"
    config.time.endtime = DateTime(2010, 7, 3)
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)
    Wflow.close_files(model; delete_output = false)

    (; reservoir) = model.routing.river_flow.boundary_conditions
    @test reservoir.boundary_conditions.external_inflow[1] == 0.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[1] ≈ 0.0
    @test reservoir.boundary_conditions.inflow[1] ≈ 5.6955612967900695
    @test reservoir.variables.storage[1] ≈ 1.901087558591953e8
    @test reservoir.variables.outflow_av[1] ≈ 3.0
    @test reservoir.variables.outflow[1] ≈ 3.0
end

# test debug message using observed outflow for two timesteps
@testitem "Piave: log debug message using observed reservoir outflow (sbm model)" begin
    using TOML
    using Dates: DateTime
    tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.input.cyclic["reservoir_water__outgoing_observed_volume_flow_rate"] = "reservoir_outflow"
    config.logging.loglevel = "debug"
    config.logging.path_log = "log_sbm_piave_debug.txt"
    config.time.endtime = DateTime(2010, 7, 3)
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

@testitem "Water balance Piave water demand (sbm model)" begin
    tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
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
