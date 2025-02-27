
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_gwf_model(config)
(; network) = model

Wflow.run_timestep!(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2000-06-01T00:00:00")
    @test row.Q_av ≈ 0.01619703129434486f0
    @test row.head ≈ 1.6471323360175287f0
end

@testset "first timestep" begin
    sbm = model.land

    @test model.clock.iteration == 1
    @test sbm.soil.parameters.theta_s[1] ≈ 0.44999998807907104f0
    @test sbm.soil.variables.runoff[1] == 0.0
    @test sbm.soil.variables.soilevap[1] == 0.0
    @test sbm.soil.variables.transpiration[1] ≈ 0.30587632831650247f0
end

# run the second timestep
Wflow.run_timestep!(model)

@testset "second timestep" begin
    sbm = model.land
    @test sbm.soil.parameters.theta_s[1] ≈ 0.44999998807907104f0
    @test sbm.soil.variables.runoff[1] == 0.0
    @test sbm.soil.variables.soilevap[1] == 0.0
    @test sbm.soil.variables.transpiration[4] ≈ 0.7000003898938235f0
end

@testset "overland flow (kinematic wave)" begin
    q = model.routing.overland_flow.variables.q_av
    @test sum(q) ≈ 2.2321111203610908f-7
end

@testset "river domain (kinematic wave)" begin
    q = model.routing.river_flow.variables.q_av
    river = model.routing.river_flow
    @test sum(q) ≈ 0.035425926757567935f0
    @test q[6] ≈ 0.00802617565138912f0
    @test river.variables.storage[6] ≈ 4.528690358701646f0
    @test river.boundary_conditions.inwater[6] ≈ 0.0004037674722635451f0
    @test q[13] ≈ 0.0006017024138583771f0
    @test q[network.river.order[end]] ≈ 0.008553261399338265f0
end

@testset "groundwater" begin
    gw = model.routing.subsurface_flow
    @test gw.boundaries.river.variables.stage[1] ≈ 1.2123636929067039f0
    @test gw.aquifer.variables.head[17:21] ≈ [
        1.2866380350225155f0,
        1.3477853512604643f0,
        1.7999999523162842f0,
        1.6225103807809076f0,
        1.4053590307668113f0,
    ]
    @test gw.boundaries.river.variables.flux[1] ≈ -51.32817280138138f0
    @test gw.boundaries.drain.variables.flux[1] ≈ 0.0
    @test gw.boundaries.recharge.variables.rate[19] ≈ -0.0014241196552847502f0
end

@testset "no drains" begin
    config.model.drains = false
    delete!(
        Dict(config.output.netcdf_grid.variables),
        "land_drain_water~to-subsurface__volume_flow_rate",
    )
    model = Wflow.initialize_sbm_gwf_model(config)
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

model = Wflow.initialize_sbm_gwf_model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "river domain (local inertial)" begin
    q = model.routing.river_flow.variables.q_av
    river = model.routing.river_flow
    @test sum(q) ≈ 0.025966484848150714f0
    @test q[6] ≈ 0.0057918662111618585f0
    @test river.variables.storage[6] ≈ 7.347727838567257f0
    @test river.boundary_conditions.inwater[6] ≈ 0.00017632250611970184f0
    @test q[13] ≈ 0.00044406241604129745f0
    @test q[5] ≈ 0.006109927807397471f0
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

model = Wflow.initialize_sbm_gwf_model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "river and land domain (local inertial)" begin
    q = model.routing.river_flow.variables.q_av
    @test sum(q) ≈ 0.02596647126178727f0
    @test q[6] ≈ 0.0057918687419759255f0
    @test q[13] ≈ 0.00044406258154064935f0
    @test q[5] ≈ 0.00610990650349414f0
    h = model.routing.river_flow.variables.h_av
    @test h[6] ≈ 0.07894230285870471f0
    @test h[5] ≈ 0.07635048570353754f0
    @test h[13] ≈ 0.08095204525673293f0
    qx = model.routing.overland_flow.variables.qx
    qy = model.routing.overland_flow.variables.qy
    @test all(qx .== 0.0f0)
    @test all(qy .== 0.0f0)
end
Wflow.close_files(model; delete_output = false)

# test with warm start
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)
config.model.reinit = false

model = Wflow.initialize_sbm_gwf_model(config)
(; network) = model

Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "second timestep warm start" begin
    sbm = model.land
    @test sbm.soil.variables.runoff[1] == 0.0
    @test sbm.soil.variables.soilevap[1] ≈ 0.2889306511074693f0
    @test sbm.soil.variables.transpiration[1] ≈ 0.8370726722706481f0
end

@testset "overland flow warm start (kinematic wave)" begin
    q = model.routing.overland_flow.variables.q_av
    @test sum(q) ≈ 1.4233852635648338f-5
end

@testset "river domain warm start (kinematic wave)" begin
    q = model.routing.river_flow.variables.q_av
    river = model.routing.river_flow
    @test sum(q) ≈ 0.01191742350356312f0
    @test q[6] ≈ 0.0024353072305122064f0
    @test river.variables.storage[6] ≈ 2.2277585577366357f0
    @test river.boundary_conditions.inwater[6] ≈ -1.298187928273214f-5
    @test q[13] ≈ 7.332742814063803f-5
    @test q[network.river.order[end]] ≈ 0.002472526149620472f0
end

@testset "groundwater warm start" begin
    gw = model.routing.subsurface_flow
    @test gw.boundaries.river.variables.stage[1] ≈ 1.2031171676781156f0
    @test gw.aquifer.variables.head[17:21] ≈ [
        1.2277456867225283f0,
        1.286902494792006f0,
        1.7999999523162842f0,
        1.5901747932190804f0,
        1.2094238817776854f0,
    ]
    @test gw.boundaries.river.variables.flux[1] ≈ -6.692884222603261f0
    @test gw.boundaries.drain.variables.flux[1] ≈ 0.0
    @test gw.boundaries.recharge.variables.rate[19] ≈ -0.0014241196552847502f0
end

Wflow.close_files(model; delete_output = false)
