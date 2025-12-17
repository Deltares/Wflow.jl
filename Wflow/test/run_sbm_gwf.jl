
@testitem "Run model sbm_gwf (kinematic wave routing)" begin
    using Dates: DateTime
    include("testing_utils.jl")
    tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
    config = Wflow.Config(tomlpath)

    model = Wflow.Model(config)
    (; domain) = model

    @testset "initial states and depending variables" begin
        # test if states and depending variables are consistent between soil and groundwater
        # flow models
        (; aquifer) = model.routing.subsurface_flow
        (; zi, ustorecapacity) = model.land.soil.variables
        (; land_indices) = model.domain.river.network
        @test all(
            0.001 * zi .==
            aquifer.parameters.top .- min.(aquifer.variables.head, aquifer.parameters.top),
        )
        @test all(ustorecapacity[land_indices] .== 0.0)
        @test all(
            aquifer.variables.head[land_indices] .== aquifer.parameters.top[land_indices],
        )
    end

    Wflow.run_timestep!(model)

    # test if the first timestep was written to the CSV file
    flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
    @testset "CSV output" begin
        row = csv_first_row(model.writer.csv_path)

        @test row.time == DateTime("2000-06-01T00:00:00")
        @test row.Q_av ≈ 0.01619703129434486
        @test row.head ≈ 1.6466006412729106
    end

    @testset "first timestep" begin
        sbm = model.land

        @test model.clock.iteration == 1
        @test sbm.soil.parameters.theta_s[1] ≈ 0.44999998807907104
        @test sbm.soil.variables.runoff[1] == 0.0
        @test sbm.soil.variables.soilevap[1] == 0.0
        @test sbm.soil.variables.transpiration[1] ≈ 0.30587632831650247
    end

    # run the second timestep
    Wflow.run_timestep!(model)

    @testset "second timestep" begin
        sbm = model.land
        @test sbm.soil.parameters.theta_s[1] ≈ 0.44999998807907104
        @test sbm.soil.variables.runoff[1] == 0.0
        @test sbm.soil.variables.soilevap[1] == 0.0
        @test sbm.soil.variables.transpiration[4] ≈ 0.9545461724219301
    end

    @testset "overland flow (kinematic wave)" begin
        q = model.routing.overland_flow.variables.q_av
        @test sum(q) ≈ 2.2319312569903814e-7
    end

    @testset "river domain (kinematic wave)" begin
        q = model.routing.river_flow.variables.q_av
        river = model.routing.river_flow
        @test sum(q) ≈ 0.03025208375758796
        @test q[6] ≈ 0.006665169891136884
        @test river.variables.storage[6] ≈ 4.037102927161714
        @test river.boundary_conditions.inwater[6] ≈ 0.0001323377399268941
        @test q[13] ≈ 0.0005348254420989118
        @test q[domain.river.network.order[end]] ≈ 0.006956110669004074
    end

    @testset "groundwater" begin
        gw = model.routing.subsurface_flow
        @test gw.boundaries.river.variables.stage[1] ≈ 1.212479774379469
        @test gw.aquifer.variables.head[17:21] ≈ [
            1.2499617962076572,
            1.3146958617176479,
            1.7999999523162842,
            1.6201330681295285,
            1.3726090475101977,
        ]
        @test gw.boundaries.river.variables.flux[1] ≈ -44.36320056565769
        @test gw.boundaries.drain.variables.flux[1] ≈ 0.0
        @test gw.boundaries.recharge.variables.rate[19] ≈ -0.0014241196552847502
    end

    @testset "no drains" begin
        config.model.drain__flag = false
        delete!(
            config.output.netcdf_grid.variables,
            "land_drain_water__to_subsurface_volume_flow_rate",
        )
        model = Wflow.Model(config)
        @test typeof.(Wflow.get_boundaries(model.routing.subsurface_flow.boundaries)) ==
              (Wflow.Recharge, Wflow.GwfRiver, Nothing, Nothing)
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "Run sbm_gwf (local inertial routing)" begin
    # test complete run including logging entry TOML file (not set)
    tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")

    # test local-inertial option for river flow routing
    tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
    config = Wflow.Config(tomlpath)
    config.model.river_routing = "local_inertial"

    config.input.static["river_bank_water__elevation"] = "bankfull_elevation"
    config.input.static["river_bank_water__depth"] = "bankfull_depth"

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "river domain (local inertial)" begin
        q = model.routing.river_flow.variables.q_av
        river = model.routing.river_flow
        @test sum(q) ≈ 0.021521084547539975
        @test q[6] ≈ 0.004598563634057708
        @test river.variables.storage[6] ≈ 6.169062629400154
        @test river.boundary_conditions.inwater[6] ≈ -9.510722621694914e-5
        @test q[13] ≈ 0.000389969767596223
        @test q[5] ≈ 0.004699224180900602
    end
    Wflow.close_files(model; delete_output = false)

    # test local-inertial option for river and overland flow routing
    tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
    config = Wflow.Config(tomlpath)
    config.model.river_routing = "local_inertial"
    config.model.land_routing = "local_inertial"

    config.input.static["river_bank_water__elevation"] = "bankfull_elevation"
    config.input.static["river_bank_water__depth"] = "bankfull_depth"
    config.input.static["land_surface_water_flow__ground_elevation"] = "wflow_dem"

    pop!(Dict(config.state.variables), "land_surface_water__instantaneous_volume_flow_rate")
    config.state.variables.land_surface_water__depth = "h_av_land"
    config.state.variables.land_surface_water__x_component_of_instantaneous_volume_flow_rate = "qx_land"
    config.state.variables.land_surface_water__y_component_of_instantaneous_volume_flow_rate = "qy_land"

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "river and land domain (local inertial)" begin
        q = model.routing.river_flow.variables.q_av
        @test sum(q) ≈ 0.021521072635048613
        @test q[6] ≈ 0.004598566560558798
        @test q[13] ≈ 0.00038996995865291496
        @test q[5] ≈ 0.004699203348707801
        h = model.routing.river_flow.variables.h
        @test h[6] ≈ 0.06169046479599844
        @test h[5] ≈ 0.059452526028619744
        @test h[13] ≈ 0.06377260202490082
        qx = model.routing.overland_flow.variables.qx
        qy = model.routing.overland_flow.variables.qy
        @test all(qx .== 0.0)
        @test all(qy .== 0.0)
    end
    Wflow.close_files(model; delete_output = false)

    # test with warm start
    tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
    config = Wflow.Config(tomlpath)
    config.model.cold_start__flag = false

    model = Wflow.Model(config)
    (; domain) = model

    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep warm start" begin
        sbm = model.land
        @test sbm.soil.variables.runoff[1] == 0.0
        @test sbm.soil.variables.soilevap[1] ≈ 0.28488618656022874
        @test sbm.soil.variables.transpiration[1] ≈ 1.0122634204681036
    end

    @testset "river domain warm start (kinematic wave)" begin
        q = model.routing.river_flow.variables.q_av
        river = model.routing.river_flow
        @test sum(q) ≈ 0.011927558704212266
        @test q[6] ≈ 0.002439237560893128
        @test river.variables.storage[6] ≈ 2.2297897880708595
        @test river.boundary_conditions.inwater[6] ≈ -1.4361419680049275e-5
        @test q[13] ≈ 6.915367636708339e-5
        @test q[domain.river.network.order[end]] ≈ 0.0024769176125029905
    end

    @testset "groundwater warm start" begin
        gw = model.routing.subsurface_flow
        @test gw.boundaries.river.variables.stage[1] ≈ 1.2030201719029363
        @test gw.aquifer.variables.head[17:21] ≈ [
            1.226363413451762,
            1.2844302137711752,
            1.7999999523162842,
            1.5860343731919113,
            1.2068931431179248,
        ]
        @test gw.boundaries.river.variables.flux[1] ≈ -6.472282449944586
        @test gw.boundaries.drain.variables.flux[1] ≈ 0.0
        @test gw.boundaries.recharge.variables.rate[19] ≈ -0.0014241196552847502
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "overland flow warm start (kinematic wave)" begin
    tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
    config = Wflow.Config(tomlpath)
    config.model.cold_start__flag = false

    model = Wflow.Model(config)

    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    q = model.routing.overland_flow.variables.q_av
    @test sum(q) ≈ 1.4233852635648338e-5
end

@testitem "run wflow sbm_gwf" begin
    tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
    config = Wflow.Config(tomlpath)
    config.time.endtime = "2000-06-04"
    Wflow.run(config)
end

@testitem "water balance sbm with groundwater" begin
    tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
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
        @test all(e -> abs(e) < 1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
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
        @test all(e -> abs(e) < 1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end
