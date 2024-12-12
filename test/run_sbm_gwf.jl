
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
    sbm = model.vertical

    @test model.clock.iteration == 1
    @test sbm.soil.parameters.theta_s[1] ≈ 0.44999998807907104f0
    @test sbm.soil.variables.runoff[1] == 0.0
    @test sbm.soil.variables.soilevap[1] == 0.0
    @test sbm.soil.variables.transpiration[1] ≈ 0.30587632831650247f0
end

# run the second timestep
Wflow.run_timestep!(model)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.soil.parameters.theta_s[1] ≈ 0.44999998807907104f0
    @test sbm.soil.variables.runoff[1] == 0.0
    @test sbm.soil.variables.soilevap[1] == 0.0
    @test sbm.soil.variables.transpiration[4] ≈ 0.7000003898938235f0
end

@testset "overland flow (kinematic wave)" begin
    q = model.lateral.land.variables.q_av
    @test sum(q) ≈ 2.2321111203610908f-7
end

@testset "river domain (kinematic wave)" begin
    q = model.lateral.river.variables.q_av
    river = model.lateral.river
    @test sum(q) ≈ 0.035468154534622556f0
    @test q[6] ≈ 0.00803825101724232f0
    @test river.variables.volume[6] ≈ 4.532124903256408f0
    @test river.boundary_conditions.inwater[6] ≈ 0.00040826084140456616f0
    @test q[13] ≈ 0.0006017024138583771f0
    @test q[network.river.order[end]] ≈ 0.00856866488665273f0
end

@testset "groundwater" begin
    gw = model.lateral.subsurface
    @test gw.river.variables.stage[1] ≈ 1.2123636929067039f0
    @test gw.flow.aquifer.variables.head[17:21] ≈ [
        1.2866380350225155f0,
        1.3477853512604643f0,
        1.7999999523162842f0,
        1.6225103807809076f0,
        1.4053590307668113f0,
    ]
    @test gw.river.variables.flux[1] ≈ -51.34674583702381f0
    @test gw.drain.variables.flux[1] ≈ 0.0
    @test gw.recharge.variables.rate[19] ≈ -0.0014241196552847502f0
end

@testset "no drains" begin
    config.model.drains = false
    delete!(Dict(config.output.lateral.subsurface), "drain")
    model = Wflow.initialize_sbm_gwf_model(config)
    @test collect(keys(model.lateral.subsurface)) == [:flow, :recharge, :river]
end

Wflow.close_files(model; delete_output = false)

# test local-inertial option for river flow routing
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)
config.model.river_routing = "local-inertial"

config.input.river_bank_water__elevation = "bankfull_elevation"
config.input.river_bank_water__depth = "bankfull_depth"

model = Wflow.initialize_sbm_gwf_model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "river domain (local inertial)" begin
    q = model.lateral.river.variables.q_av
    river = model.lateral.river
    @test sum(q) ≈ 0.02727911500112358f0
    @test q[6] ≈ 0.006111263175002127f0
    @test river.variables.volume[6] ≈ 7.6120096530771075f0
    @test river.boundary_conditions.inwater[6] ≈ 0.0002210785332342944f0
    @test q[13] ≈ 0.0004638698607639214f0
    @test q[5] ≈ 0.0064668491697542786f0
end
Wflow.close_files(model; delete_output = false)

# test local-inertial option for river and overland flow routing
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)
config.model.river_routing = "local-inertial"
config.model.land_routing = "local-inertial"

config.input.river_bank_water__elevation = "bankfull_elevation"
config.input.river_bank_water__depth = "bankfull_depth"
config.input.land_surface_water_flow__ground_elevation = "wflow_dem"

pop!(Dict(config.state.lateral.land.variables), "q")
config.state.lateral.land.variables.h_av = "h_av_land"
config.state.lateral.land.variables.qx = "qx_land"
config.state.lateral.land.variables.qy = "qy_land"

model = Wflow.initialize_sbm_gwf_model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "river and land domain (local inertial)" begin
    q = model.lateral.river.variables.q_av
    @test sum(q) ≈ 0.027286431923384962f0
    @test q[6] ≈ 0.00611309161099138f0
    @test q[13] ≈ 0.0004639786629631376f0
    @test q[5] ≈ 0.006468859889145798f0
    h = model.lateral.river.variables.h_av
    @test h[6] ≈ 0.08120137914886108f0
    @test h[5] ≈ 0.07854966203902745f0
    @test h[13] ≈ 0.08323543174453409f0
    qx = model.lateral.land.variables.qx
    qy = model.lateral.land.variables.qy
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
    sbm = model.vertical
    @test sbm.soil.variables.runoff[1] == 0.0
    @test sbm.soil.variables.soilevap[1] ≈ 0.2889306511074693f0
    @test sbm.soil.variables.transpiration[1] ≈ 0.8370726722706481f0
end

@testset "overland flow warm start (kinematic wave)" begin
    q = model.lateral.land.variables.q_av
    @test sum(q) ≈ 1.4224503548471601f-5
end

@testset "river domain warm start (kinematic wave)" begin
    q = model.lateral.river.variables.q_av
    river = model.lateral.river
    @test sum(q) ≈ 0.01191742350356312f0
    @test q[6] ≈ 0.0024353072305122064f0
    @test river.variables.volume[6] ≈ 2.2277585577366357f0
    @test river.boundary_conditions.inwater[6] ≈ -1.3042629584651168f-5
    @test q[13] ≈ 7.332742814063803f-5
    @test q[network.river.order[end]] ≈ 0.002472526149620472f0
end

@testset "groundwater warm start" begin
    gw = model.lateral.subsurface
    @test gw.river.variables.stage[1] ≈ 1.2031171676781156f0
    @test gw.flow.aquifer.variables.head[17:21] ≈ [
        1.2277456867225283f0,
        1.286902494792006f0,
        1.7999999523162842f0,
        1.5901747932190804f0,
        1.2094238817776854f0,
    ]
    @test gw.river.variables.flux[1] ≈ -6.692884222603261f0
    @test gw.drain.variables.flux[1] ≈ 0.0
    @test gw.recharge.variables.rate[19] ≈ -0.0014241196552847502f0
end

Wflow.close_files(model; delete_output = false)
