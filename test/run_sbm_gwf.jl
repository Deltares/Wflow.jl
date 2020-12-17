
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_gwf_model(config)
@unpack network = model

model = Wflow.update_sbm_gwf(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    csv = CSV.File(model.writer.csv_path)
    row = csv[1]

    @test row.time == DateTime("2000-06-01T00:00:00")
    @test row.Q_av ≈ 0.016260519762890745
    @test row.head ≈ 1.8416896245327632
end

@testset "first timestep" begin
    sbm = model.vertical

    @test model.clock.iteration == 2
    @test sbm.altitude[1] == 1.937000036239624
    @test sbm.θₛ[1] == 0.44999998807907104
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.0
    @test sbm.transpiration[1] == 0.6117526566330049
end

# run the second timestep
model = Wflow.update_sbm_gwf(model)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.altitude[1] == 1.937000036239624
    @test sbm.θₛ[1] == 0.44999998807907104
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.0
    @test sbm.transpiration[1] == 1.0122634204681036
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 0.0
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 0.03515256649675796
    @test q[6] ≈ 0.00795266793733648
    @test q[13] ≈ 0.0005980549864213564
    @test q[network.river.order[end]] ≈ 0.008482646477950898
end

@testset "groundwater" begin
    gw = model.lateral.subsurface
    @test gw.river.stage[1] ≈ 1.21057153442702
    @test gw.flow.aquifer.head[19] ≈ 1.7999999523162842
    @test gw.river.flux[1] ≈ -50.751898489778426
    @test gw.drain.flux[1] ≈ 0.0
    @test gw.recharge.rate[19] ≈ -0.0014241196552847502
end

Wflow.close_files(model, delete_output = false)
