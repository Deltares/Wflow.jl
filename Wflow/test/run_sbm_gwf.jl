
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.Model(config)
(; domain) = model

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
    @test sbm.soil.variables.transpiration[4] ≈ 0.7000003898938235
end

@testset "overland flow (kinematic wave)" begin
    q = model.routing.overland_flow.variables.q_av
    @test sum(q) ≈ 2.2319312569903814e-7
end

@testset "river domain (kinematic wave)" begin
    q = model.routing.river_flow.variables.q_av
    river = model.routing.river_flow
    @test sum(q) ≈ 0.035425926757567935
    @test q[6] ≈ 0.00802617565138912
    @test river.variables.storage[6] ≈ 4.528690358701646
    @test river.boundary_conditions.inwater[6] ≈ 0.0004037674722635451
    @test q[13] ≈ 0.0006016241976247014
    @test q[domain.river.network.order[end]] ≈ 0.008553261399338265
end

@testset "groundwater" begin
    gw = model.routing.subsurface_flow
    @test gw.boundaries.river.variables.stage[1] ≈ 1.212479774379469
    @test gw.aquifer.variables.head[17:21] ≈ [
        1.286914528615439,
        1.3481205868263955,
        1.7999999523162842,
        1.6225103807809076,
        1.4061108676063327,
    ]
    @test gw.boundaries.river.variables.flux[1] ≈ -51.32817280138138
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
    @test sum(q) ≈ 0.025966484848150714
    @test q[6] ≈ 0.0057918662111618585
    @test river.variables.storage[6] ≈ 7.347727838567257
    @test river.boundary_conditions.inwater[6] ≈ 0.00017632250611970184
    @test q[13] ≈ 0.00044406241604129745
    @test q[5] ≈ 0.006109927807397471
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
    @test sum(q) ≈ 0.02596647126178727
    @test q[6] ≈ 0.0057918687419759255
    @test q[13] ≈ 0.00044406258154064935
    @test q[5] ≈ 0.00610990650349414
    h = model.routing.river_flow.variables.h_av
    @test h[6] ≈ 0.07894230285870471
    @test h[5] ≈ 0.07635048570353754
    @test h[13] ≈ 0.08095204525673293
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
    @test sbm.soil.variables.soilevap[1] ≈ 0.2888846391558502
    @test sbm.soil.variables.transpiration[1] ≈ 0.8370726722706481
end

@testset "overland flow warm start (kinematic wave)" begin
    q = model.routing.overland_flow.variables.q_av
    @test sum(q) ≈ 1.4233852635648338e-5
end

@testset "river domain warm start (kinematic wave)" begin
    q = model.routing.river_flow.variables.q_av
    river = model.routing.river_flow
    @test sum(q) ≈ 0.011918378125641867
    @test q[6] ≈ 0.0024355453044100587
    @test river.variables.storage[6] ≈ 2.2278920306090555
    @test river.boundary_conditions.inwater[6] ≈ -1.298187928273214e-5
    @test q[13] ≈ 7.335203306033115e-5
    @test q[domain.river.network.order[end]] ≈ 0.002472763875440307
end

@testset "groundwater warm start" begin
    gw = model.routing.subsurface_flow
    @test gw.boundaries.river.variables.stage[1] ≈ 1.2030201719029363
    @test gw.aquifer.variables.head[17:21] ≈ [
        1.2277413823642467,
        1.2868963785900465,
        1.7999999523162842,
        1.5901747023422137,
        1.2094146088822748,
    ]
    @test gw.boundaries.river.variables.flux[1] ≈ -6.693790653429019
    @test gw.boundaries.drain.variables.flux[1] ≈ 0.0
    @test gw.boundaries.recharge.variables.rate[19] ≈ -0.0014241196552847502
end

Wflow.close_files(model; delete_output = false)
