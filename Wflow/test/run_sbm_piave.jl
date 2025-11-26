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
        227.48035999493885,
        201.66655380222255,
        277.3585133618376,
        167.06113313065387,
        172.85132307882748,
        158.4868851162215,
        133.42316714652813,
        118.79239436626797,
        189.3066342899249,
        115.74440489724931,
    ]
    @test q_[idx] ≈ [
        230.50891783820745,
        208.70476972858077,
        284.38432448798,
        175.24143075182934,
        185.6574635583777,
        165.81546406628178,
        144.3529395342142,
        126.62715698373853,
        191.90325467397034,
        122.22507259114501,
    ]
    @test riv_storage_demand[idx] ≈ [
        61395.09279571792,
        57160.65692276638,
        63092.379798845825,
        51024.87776219157,
        53982.12212292085,
        46201.9194344024,
        44610.78797109941,
        42552.63008891031,
        46560.29697212866,
        41849.31928155387,
    ]
    @test riv_storage[idx] ≈ [
        62013.48243896582,
        57826.70054420379,
        63310.699787656515,
        51620.810255230586,
        54634.346978697016,
        46652.37923131112,
        45015.6453650631,
        42995.04542365768,
        46455.15009283594,
        41789.82947520561,
    ]
    @test ssf_storage_demand[idx] ≈ [
        293460.22026866337,
        288923.6353138966,
        287378.90184287925,
        283024.21156711585,
        279450.68877513084,
        274962.15361320204,
        270498.72127685224,
        266025.07960873423,
        264560.66050100315,
        260999.95670112918,
    ]
    @test ssf_storage[idx] ≈ [
        293448.2011923609,
        288845.5752515868,
        287251.56581749313,
        282864.8960769829,
        279288.3918930026,
        274785.59093307354,
        270351.9248607254,
        265899.86716567166,
        264349.145295683,
        260749.3454407325,
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
        @test sum(surfacewater_alloc) ≈ 1662.1429336566434
        @test sum(act_groundwater_abst) ≈ 373.02432487790287
        @test paddy.variables.h[[25, 42, 45]] ≈
              [43.181590971577364, 51.20409409088053, 34.473513683211834]
        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              [4.0442025332978435, 0.7127330271840676, 4.491063429489396]
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
              [29.259245814466148, 32.686077716495625, 39.970184252221905]
        @test reservoir.variables.storage ≈ [1.8959991287774065e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_av ≈
              [4.839261880167545, 9.717572497025145, 58.105643206370274]
        @test soil.variables.exfiltsatwater[[937, 939, 979, 1020, 1158]] ≈ [
            3.0858603777229376,
            5.577976128972852,
            3.7860523390663428,
            5.881271624508977,
            13.887309622273573,
        ]
        @test maximum(soil.variables.exfiltsatwater) ≈ 238.20636973621725
        @test mean(river_flow.variables.q_av) ≈ 60.513939494438844
        @test maximum(river_flow.variables.q_av) ≈ 235.61007781125195
    end

    Wflow.run_timestep!(model)

    @testset "Second timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
        @test sum(surfacewater_alloc) ≈ 1776.237543072083
        @test sum(act_groundwater_abst) ≈ 403.9136569139513
        @test paddy.variables.h[[25, 42, 45]] ≈
              [39.227496928297654, 48.04243573106785, 28.97022876577813]
        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              [4.722084029669344, 0.7260733256347136, 5.452294821220649]
        @test reservoir.variables.waterlevel ≈
              [29.251492025015416, 32.68607771649562, 39.97018425222191]
        @test reservoir.variables.storage ≈ [1.895496683220998e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_av ≈
              [4.84140054568368, 9.325285002501143, 54.827349895897754]
        @test soil.variables.exfiltsatwater[[937, 939, 979, 1020, 1158]] ≈ [
            3.462567528286398,
            6.141947091961624,
            4.332401346606494,
            6.305582623275848,
            14.316493120446802,
        ]
        @test maximum(soil.variables.exfiltsatwater) ≈ 227.75222991987133
        @test mean(river_flow.variables.q_av) ≈ 56.77305056192323
        @test maximum(river_flow.variables.q_av) ≈ 227.42442386186528
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
    @test reservoir.boundary_conditions.inflow[1] ≈ 5.7375179314371225
    @test reservoir.variables.storage[1] ≈ 1.895496683220998e8
    @test reservoir.variables.outflow_av[1] ≈ 4.84140054568368
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
    @test reservoir.boundary_conditions.inflow[1] ≈ 2.737517449273689
    @test reservoir.variables.storage[1] ≈ 1.8903375136457917e8
    @test reservoir.variables.outflow_av[1] ≈ 4.821625016952349
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
    @test reservoir.boundary_conditions.inflow[1] ≈ 5.71697857643707
    @test reservoir.variables.storage[1] ≈ 1.9011126436313578e8
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
        @test all(e -> abs(e) < 1.e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(e -> abs(e) < 1.1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1.e-9, subsurface_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, land_water_balance.error)
        @test all(re -> abs(re) < 1e-9, land_water_balance.relative_error)
        @test all(e -> abs(e) < 1.e-9, routing.overland_water_balance.error)
        @test all(re -> abs(re) < 1.e-9, routing.overland_water_balance.relative_error)
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(e -> abs(e) < 1.2e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end
