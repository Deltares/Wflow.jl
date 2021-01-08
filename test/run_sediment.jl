using Wflow

tomlpath = joinpath(@__DIR__, "sediment_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)
@unpack network = model

model = Wflow.update(model)

@testset "first timestep sediment model (vertical)" begin
    eros = model.vertical

    @test eros.erosov[1] == 0.9
    @test model.clock.iteration == 2
    @test mean(eros.leaf_area_index) ≈ 1.7120018886212223
    @test eros.dmsand[1] == 200.0
    @test eros.dmlagg[1] == 500.0
    @test mean(eros.interception) ≈ 0.4767846753916875
    @test mean(eros.soilloss) ≈ 0.00846987746234039
end

# run the second timestep
model = Wflow.update(model)

@testset "second timestep sediment model (vertical)" begin
    eros = model.vertical

    @test mean(eros.soilloss) ≈ 0.07441259472436522
    @test mean(eros.erosclay) ≈ 0.002198312957830507
    @test mean(eros.erossand) ≈ 0.02515898751801652
    @test mean(eros.eroslagg) ≈ 0.021613617766119256
    @test mean(eros.TCsed) == 0.0
    @test mean(eros.TCsilt) ≈ 1.0988158364353527e6
    @test mean(eros.TCsand) ≈ 1.0987090622888755e6
    @test mean(eros.TCclay) ≈ 1.0992655197016734e6
    @test mean(eros.TCsilt) ≈ 1.0988158364353527e6
end

@testset "second timestep sediment model (lateral)" begin
    lat = model.lateral

    @test mean(lat.land.inlandsed) ≈ 0.0727443455767198
    @test mean(lat.land.inlandclay) ≈ 0.002190037611542877
    @test mean(lat.land.inlandsand) ≈ 0.024512058443074893
    @test mean(lat.land.olclay) ≈ 0.006212441628146587

    @test mean(lat.river.SSconc) ≈ 15.38792327918787
    @test mean(lat.river.inlandclay) ≈ 0.019393435838709512
    @test lat.river.h_riv[network.river.order[end]] ≈ 0.006103649735450745
    @test lat.river.outclay[1] ≈ 1.1654682799750028e-6
end

Wflow.close_files(model)

benchmark = @benchmark Wflow.run_simulation(tomlpath)
trialmin = BenchmarkTools.minimum(benchmark)
println("Sediment Model update (run_simulation)")
print_benchmark(trialmin)
