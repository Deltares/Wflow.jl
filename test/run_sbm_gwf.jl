
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_gwf_model(config)
@unpack network = model

model = Wflow.run_timestep(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2000-06-01T00:00:00")
    @test row.Q_av ≈ 0.01620324716944374f0
    @test row.head ≈ 1.6506700475116074f0
end

@testset "first timestep" begin
    sbm = model.vertical

    @test model.clock.iteration == 1
    @test sbm.theta_s[1] ≈ 0.44999998807907104f0
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.0
    @test sbm.transpiration[1] ≈ 0.0
end

# run the second timestep
model = Wflow.run_timestep(model)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.theta_s[1] ≈ 0.44999998807907104f0
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.0
    @test sbm.transpiration[4] ≈ 0.8696975782458946f0
end

@testset "overland flow (kinematic wave)" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 2.229860508650628f-7
end

@testset "river domain (kinematic wave)" begin
    q = model.lateral.river.q_av
    river = model.lateral.river
    @test sum(q) ≈ 0.035820105590583136f0
    @test q[6] ≈ 0.008130882953221714f0
    @test river.volume[6] ≈ 4.5666019378430995f0
    @test river.inwater[6] ≈ 0.00042601629503311733f0
    @test q[13] ≈ 0.0006060110827771738f0
    @test q[network.river.order[end]] ≈ 0.008668026817817552f0
end

@testset "groundwater" begin
    gw = model.lateral.subsurface
    @test gw.river.stage[1] ≈ 1.2123636929067039f0
    @test gw.flow.aquifer.head[17:21] ≈ [
        1.2880480319282992f0,
        1.3490883291679163f0,
        1.7999999523162842f0,
        1.625090613936394f0,
        1.4071376687672454f0,
    ]
    @test gw.river.flux[1] ≈ -52.12101697937548f0
    @test gw.drain.flux[1] ≈ 0.0
    @test gw.recharge.rate[19] ≈ -0.0014241196552847502f0
end

@testset "no drains" begin
    config.model.drains = false
    delete!(Dict(config.output.lateral.subsurface), "drain")
    model = Wflow.initialize_sbm_gwf_model(config)
    @test collect(keys(model.lateral.subsurface)) == [:flow, :recharge, :river]
end

Wflow.close_files(model, delete_output = false)

# test local-inertial option for river flow routing
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)
config.model.river_routing = "local-inertial"

config.input.lateral.river.bankfull_elevation = "bankfull_elevation"
config.input.lateral.river.bankfull_depth = "bankfull_depth"

model = Wflow.initialize_sbm_gwf_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "river domain (local inertial)" begin
    q = model.lateral.river.q_av
    river = model.lateral.river
    @test sum(q) ≈ 0.02761343108159262f0
    @test q[6] ≈ 0.006200879624415703f0
    @test river.volume[6] ≈ 7.682451475056573f0
    @test river.inwater[6] ≈ 0.00023950387043266292f0
    @test q[13] ≈ 0.00046741800834202964f0
    @test q[5] ≈ 0.006564513429686601f0
end
Wflow.close_files(model, delete_output = false)

# test local-inertial option for river and overland flow routing
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)
config.model.river_routing = "local-inertial"
config.model.land_routing = "local-inertial"

config.input.lateral.river.bankfull_elevation = "bankfull_elevation"
config.input.lateral.river.bankfull_depth = "bankfull_depth"
config.input.lateral.land.elevation = "wflow_dem"

pop!(Dict(config.state.lateral.land), "q")
config.state.lateral.land.h_av = "h_av_land"
config.state.lateral.land.qx = "qx_land"
config.state.lateral.land.qy = "qy_land"

model = Wflow.initialize_sbm_gwf_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "river and land domain (local inertial)" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 0.027620758141556154f0
    @test q[6] ≈ 0.006202710435741449f0
    @test q[13] ≈ 0.0004675269778568993f0
    @test q[5] ≈ 0.006566526863956613f0
    h = model.lateral.river.h_av
    @test h[6] ≈ 0.08180205294088073f0
    @test h[5] ≈ 0.07913252828380887f0
    @test h[13] ≈ 0.08383078664048482f0
    qx = model.lateral.land.qx
    qy = model.lateral.land.qy
    @test all(qx .== 0.0f0)
    @test all(qy .== 0.0f0)
end
Wflow.close_files(model, delete_output = false)

# test with warm start
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)
config.model.reinit = false

model = Wflow.initialize_sbm_gwf_model(config)
@unpack network = model

model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "second timestep warm start" begin
    sbm = model.vertical
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] ≈ 0.2889306511074693f0
    @test sbm.transpiration[1] ≈ 0.9434505457255187f0
end

@testset "overland flow warm start (kinematic wave)" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 1.4411427003142072f-5
end

@testset "river domain warm start (kinematic wave)" begin
    q = model.lateral.river.q_av
    river = model.lateral.river
    @test sum(q) ≈ 0.011323594506729224f0
    @test q[6] ≈ 0.0022573990911054567f0
    @test river.volume[6] ≈ 2.1222151581462616f0
    @test river.inwater[6] ≈ -6.857216253435187f-5
    @test q[13] ≈ 8.677840362842245f-5
    @test q[network.river.order[end]] ≈ 0.002259888027651648f0
end

@testset "groundwater warm start" begin
    gw = model.lateral.subsurface
    @test gw.river.stage[1] ≈ 1.2031171676781156f0
    @test gw.flow.aquifer.head[17:21] ≈ [
        1.2192158271216242f0,
        1.2740235737998076f0,
        1.7999999523162842f0,
        1.5873527569888766f0,
        1.2020901533167596f0,
    ]
    @test gw.river.flux[1] ≈ -6.055040670091216f0
    @test gw.drain.flux[1] ≈ 0.0
    @test gw.recharge.rate[19] ≈ -0.0014241196552847502f0
end
Wflow.close_files(model, delete_output = false)
