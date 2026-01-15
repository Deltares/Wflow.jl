@testitem "Piave with and without water demand (sbm model)" begin
    include("testing_utils.jl")
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

    idx = 1:3:28
    @test q_demand[idx] ≈ [
        227.3831999880687,
        202.19047094531086,
        277.8028528359398,
        168.31471715982528,
        176.94887099172664,
        162.37592041137452,
        140.0014053226863,
        124.77405278981702,
        194.49291677184874,
        123.23205104035759,
    ]
    @test q_[idx] ≈ [
        230.50891797513165,
        209.1351422536607,
        285.61935918748566,
        176.10753806451422,
        188.8533705755549,
        169.5089639155915,
        149.9865573749758,
        132.52624991627732,
        197.50798672292973,
        128.2211240687692,
    ]
    @test riv_storage_demand[idx] ≈ [
        61374.41368668161,
        57429.603754662334,
        63223.691456235436,
        51509.07964519735,
        54683.53384722047,
        47023.79102606856,
        45791.943657164804,
        43947.98795862534,
        47443.503274952665,
        43036.83729231518,
    ]
    @test riv_storage[idx] ≈ [
        62014.50063779726,
        58003.940275575354,
        63436.82628710128,
        52017.726843065684,
        55278.122773417505,
        47505.141419412,
        46071.90434263895,
        44216.61902610516,
        47330.7289988868,
        42960.82591181803,
    ]
    @test ssf_storage_demand[idx] ≈ [
        148902.90100638504,
        144753.45813295056,
        143305.13076818842,
        139242.78588431314,
        136244.57668144818,
        132305.65681759277,
        128513.77499731268,
        124727.24989179976,
        123626.11419601868,
        120700.05688609753,
    ]
    @test ssf_storage[idx] ≈ [
        148889.9251169022,
        144721.4870940757,
        143241.89373183213,
        139176.20841239404,
        136196.83227523384,
        132269.77925875073,
        128519.03915828765,
        124768.62986942183,
        123666.92560733219,
        120732.64759485851,
    ]
end

@testitem "Piave water demand and allocation (sbm model)" begin
    using Statistics: mean

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
              [29.259245822497967, 32.68607771649564, 39.970184252221905]
        @test reservoir.variables.storage ≈ [1.8959991292978692e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_av ≈
              [4.839261881259249, 9.717572653187853, 58.11708723477156]
        @test soil.variables.exfiltsatwater[[937, 939, 979, 1020, 1158]] ≈ [
            3.085860377722936,
            5.577976128972852,
            3.7860523390663428,
            5.881271624508976,
            13.887309622273573,
        ]
        @test maximum(soil.variables.exfiltsatwater) ≈ 238.20636973621728
        @test mean(river_flow.variables.q_av) ≈ 60.50413027502687
        @test maximum(river_flow.variables.q_av) ≈ 235.47692910316073
    end

    Wflow.run_timestep!(model)

    @testset "Second timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
        @test sum(surfacewater_alloc) ≈ 1591.4015870269836
        @test sum(act_groundwater_abst) ≈ 337.6767352329224
        @test paddy.variables.h[[25, 42, 45]] ≈
              [39.227506410553374, 48.0424361651646, 28.97022876577807]
        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              [4.442695372669294, 0.7341536180838293, 5.022845160752851]
        @test reservoir.variables.waterlevel ≈
              [29.25155309090927, 32.686077716495646, 39.970184252221905]
        @test reservoir.variables.storage ≈ [1.8955006402909216e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_av ≈
              [4.841407860671789, 9.33782883773028, 54.92001813166034]
        @test soil.variables.exfiltsatwater[[937, 939, 979, 1020, 1158]] ≈ [
            3.4947694909457154,
            6.150730906013771,
            4.334543034906983,
            6.324337963139996,
            14.292030396430834,
        ]
        @test maximum(soil.variables.exfiltsatwater) ≈ 228.5950871076442
        @test mean(river_flow.variables.q_av) ≈ 56.84871498837686
        @test maximum(river_flow.variables.q_av) ≈ 227.58250547568701
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "Piave water demand and allocation, switch off livestock, paddy and nonpaddy" begin
    tomlpath = joinpath(@__DIR__, "sbm_piave_nonirri-demand_config.toml")
    config = Wflow.Config(tomlpath)
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

@testitem "Piave: reservoir without external negative inflow (sbm model)" begin
    # test cyclic reservoir external inflow
    tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    (; reservoir) = model.routing.river_flow.boundary_conditions
    @test reservoir.boundary_conditions.external_inflow[1] == 0.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[1] == 0.0
    @test reservoir.boundary_conditions.inflow[1] ≈ 5.742102618978243
    @test reservoir.variables.storage[1] ≈ 1.895500640290987e8
    @test reservoir.variables.outflow_av[1] ≈ 4.841407860671789
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
    @test reservoir.boundary_conditions.inflow[1] ≈ 2.742102137458466
    @test reservoir.variables.storage[1] ≈ 1.8903414707267475e8
    @test reservoir.variables.outflow_av[1] ≈ 4.8216323210981376
end

@testitem "Piave: reservoir with observed (cyclic) outflow (sbm model)" begin
    using Dates: DateTime
    # test use of observed reservoir outflow (cyclic)
    tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
    config = Wflow.Config(tomlpath)
    config.input.cyclic["reservoir_water__outgoing_observed_volume_flow_rate"] = "reservoir_outflow"
    config.time.endtime = DateTime(2010, 7, 3)
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)
    Wflow.close_files(model; delete_output = false)

    (; reservoir) = model.routing.river_flow.boundary_conditions
    @test reservoir.boundary_conditions.external_inflow[1] == 0.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[1] ≈ 0.0
    @test reservoir.boundary_conditions.inflow[1] ≈ 5.721551868669384
    @test reservoir.variables.storage[1] ≈ 1.901116596885563e8
    @test reservoir.variables.outflow_av[1] ≈ 3.0
    @test reservoir.variables.outflow[1] ≈ 3.0
end

# test debug message using observed outflow for two timesteps
@testitem "Piave: log debug message using observed reservoir outflow (sbm model)" begin
    using TOML
    using Dates: DateTime
    tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
    config = Wflow.Config(tomlpath)
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
