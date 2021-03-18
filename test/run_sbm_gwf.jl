
tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_gwf_model(config)
@unpack network = model

model = Wflow.update_sbm_gwf(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2000-06-01T00:00:00")
    @test row.Q_av ≈ 0.016260519762890745
    @test row.head ≈ 1.8416896245327632
end

@testset "first timestep" begin
    sbm = model.vertical

    @test model.clock.iteration == 2
    @test sbm.θₛ[1] == 0.44999998807907104
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.0
    @test sbm.transpiration[1] == 0.6117526566330049
end

# run the second timestep
model = Wflow.update_sbm_gwf(model)

@testset "second timestep" begin
    sbm = model.vertical
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
    @test sum(q) ≈ 0.03515257228971982
    @test q[6] ≈ 0.007952670041402076
    @test q[13] ≈ 0.0005980549969548235
    @test q[network.river.order[end]] ≈ 0.008482648782941154
end

@testset "groundwater" begin
    gw = model.lateral.subsurface
    @test gw.river.stage[1] ≈ 1.21057153442702
    @test gw.flow.aquifer.head[19] ≈ 1.7999999523162842
    @test gw.river.flux[1] ≈ -50.751898489778426
    @test gw.drain.flux[1] ≈ 0.0
    @test gw.recharge.rate[19] ≈ -0.0014241196552847502
end

Wflow.close_files(model)

benchmark = @benchmark Wflow.run(tomlpath)
trialmin = BenchmarkTools.minimum(benchmark)
println("SBM GWF Model update (run)")
print_benchmark(trialmin)

@testset "no drains" begin
    config.model.drains = false
    delete!(Dict(config.output.lateral.subsurface), "drain")
    model = Wflow.initialize_sbm_gwf_model(config)
    @test collect(keys(model.lateral.subsurface)) == [:flow, :recharge, :river]
    Wflow.close_files(model)
end
