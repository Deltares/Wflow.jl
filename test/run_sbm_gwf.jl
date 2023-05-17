
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
    @test row.head ≈ 1.8416896245327632f0
end

@testset "first timestep" begin
    sbm = model.vertical

    @test model.clock.iteration == 1
    @test sbm.θₛ[1] ≈ 0.44999998807907104f0
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.0
    @test sbm.transpiration[1] ≈ 0.6117526566330049f0
end

# run the second timestep
model = Wflow.run_timestep(model)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.θₛ[1] ≈ 0.44999998807907104f0
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.0
    @test sbm.transpiration[1] ≈ 1.0122634204681036f0
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 2.2298616f-7
end

@testset "river domain" begin
    q = model.lateral.river.q_av
    river = model.lateral.river
    @test sum(q) ≈ 0.034846861707851576f0
    @test q[6] ≈ 0.007866587223009643f0
    @test river.volume[6] ≈ 4.474556355925399f0
    @test river.inwater[6] ≈ 0.00036392604840276924f0
    @test q[13] ≈ 0.0005952292000597283f0
    @test q[network.river.order[end]] ≈ 0.008387430132453135f0
end

@testset "groundwater" begin
    gw = model.lateral.subsurface
    @test gw.river.stage[1] ≈ 1.2123636929067039f0
    @test gw.flow.aquifer.head[19] ≈ 1.7999999523162842f0
    @test gw.river.flux[1] ≈ -50.46515313302901f0
    @test gw.drain.flux[1] ≈ 0.0
    @test gw.recharge.rate[19] ≈ -0.0014241196552847502f0
end

Wflow.close_files(model)

@testset "no drains" begin
    config.model.drains = false
    delete!(Dict(config.output.lateral.subsurface), "drain")
    model = Wflow.initialize_sbm_gwf_model(config)
    @test collect(keys(model.lateral.subsurface)) == [:flow, :recharge, :river]
    Wflow.close_files(model)
end
