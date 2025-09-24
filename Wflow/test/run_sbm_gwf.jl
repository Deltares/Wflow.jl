
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
    @test all(aquifer.variables.head[land_indices] .== aquifer.parameters.top[land_indices])
end

Wflow.run_timestep!(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2000-06-01T00:00:00")
    @test row.Q_av ≈ 0.01619703129434486
    @test row.head ≈ 1.6471323360175287
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
    @test sum(q) ≈ 0.03460420763022196
    @test q[6] ≈ 0.007809586997135836
    @test river.variables.storage[6] ≈ 4.452845039409768
    @test river.boundary_conditions.inwater[6] ≈ 0.0003627879093024724
    @test q[13] ≈ 0.000592145641193449
    @test q[domain.river.network.order[end]] ≈ 0.008317914181880935
end

@testset "groundwater" begin
    gw = model.routing.subsurface_flow
    @test gw.boundaries.river.variables.stage[1] ≈ 1.212479774379469
    @test gw.aquifer.variables.head[17:21] ≈ [
        1.284126255728794,
        1.345258485244677,
        1.7999999523162842,
        1.6224360370353534,
        1.398423446844073,
    ]
    @test gw.boundaries.river.variables.flux[1] ≈ -49.624776263594335
    @test gw.boundaries.drain.variables.flux[1] ≈ 0.0
    @test gw.boundaries.recharge.variables.rate[19] ≈ -0.0014241196552847502
end

@testset "no drains" begin
    config.model.drain__flag = false
    delete!(
        Dict(config.output.netcdf_grid.variables),
        "land_drain_water~to-subsurface__volume_flow_rate",
    )
    model = Wflow.Model(config)
    @test collect(keys(model.routing.subsurface_flow.boundaries)) == [:recharge, :river]
end

Wflow.close_files(model; delete_output = false)

# test complete run including logging entry TOML file (not set)
Wflow.run(tomlpath; silent = true)

# test local-inertial option for river flow routing
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)
config.model.river_routing = "local-inertial"

config.input.static.river_bank_water__elevation = "bankfull_elevation"
config.input.static.river_bank_water__depth = "bankfull_depth"

model = Wflow.Model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "river domain (local inertial)" begin
    q = model.routing.river_flow.variables.q_av
    river = model.routing.river_flow
    @test sum(q) ≈ 0.025242745327901266
    @test q[6] ≈ 0.005597737309703839
    @test river.variables.storage[6] ≈ 7.185727159119112
    @test river.boundary_conditions.inwater[6] ≈ 0.00013534294315862918
    @test q[13] ≈ 0.0004363336769343206
    @test q[5] ≈ 0.005899536075917488
end
Wflow.close_files(model; delete_output = false)

# test local-inertial option for river and overland flow routing
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)
config.model.river_routing = "local-inertial"
config.model.land_routing = "local-inertial"

config.input.static.river_bank_water__elevation = "bankfull_elevation"
config.input.static.river_bank_water__depth = "bankfull_depth"
config.input.static.land_surface_water_flow__ground_elevation = "wflow_dem"

pop!(Dict(config.state.variables), "land_surface_water__instantaneous_volume_flow_rate")
config.state.variables.land_surface_water__depth = "h_av_land"
config.state.variables.land_surface_water__x_component_of_instantaneous_volume_flow_rate = "qx_land"
config.state.variables.land_surface_water__y_component_of_instantaneous_volume_flow_rate = "qy_land"

model = Wflow.Model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "river and land domain (local inertial)" begin
    q = model.routing.river_flow.variables.q_av
    @test sum(q) ≈ 0.025242731947236016
    @test q[6] ≈ 0.005597739889073783
    @test q[13] ≈ 0.00043633384561159253
    @test q[5] ≈ 0.005899514829371013
    h = model.routing.river_flow.variables.h_av
    @test h[6] ≈ 0.07757871828439228
    @test h[5] ≈ 0.07502750093278375
    @test h[13] ≈ 0.07960288161114154
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

@testset "overland flow warm start (kinematic wave)" begin
    q = model.routing.overland_flow.variables.q_av
    @test sum(q) ≈ 1.4233852635648338e-5
end

@testset "river domain warm start (kinematic wave)" begin
    q = model.routing.river_flow.variables.q_av
    river = model.routing.river_flow
    @test sum(q) ≈ 0.011918302776165372
    @test q[6] ≈ 0.002435524985929997
    @test river.variables.storage[6] ≈ 2.2278805130264883
    @test river.boundary_conditions.inwater[6] ≈ -1.2985462545754242e-5
    @test q[13] ≈ 7.335056957297285e-5
    @test q[domain.river.network.order[end]] ≈ 0.0024727437195071356
end

@testset "groundwater warm start" begin
    gw = model.routing.subsurface_flow
    @test gw.boundaries.river.variables.stage[1] ≈ 1.2030201719029363
    @test gw.aquifer.variables.head[17:21] ≈ [
        1.2277387243899684,
        1.2868951594984024,
        1.7999999523162842,
        1.5901818206391656,
        1.2094134130011027,
    ]
    @test gw.boundaries.river.variables.flux[1] ≈ -6.693665350868727
    @test gw.boundaries.drain.variables.flux[1] ≈ 0.0
    @test gw.boundaries.recharge.variables.rate[19] ≈ -0.0014241196552847502
end

Wflow.close_files(model; delete_output = false)
