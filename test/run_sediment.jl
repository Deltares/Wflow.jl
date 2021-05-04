using Wflow

tomlpath = joinpath(@__DIR__, "sediment_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)
@unpack network = model

model = Wflow.update(model)

@testset "first timestep sediment model (vertical)" begin
    eros = model.vertical

    @test eros.erosov[1] ≈ 0.9f0
    @test model.clock.iteration == 2
    @test mean(eros.leaf_area_index) ≈ 1.7120018886212223f0
    @test eros.dmsand[1] == 200.0f0
    @test eros.dmlagg[1] == 500.0f0
    @test mean(eros.interception) ≈ 0.4767846753916875f0
    @test mean(eros.soilloss) ≈ 0.00846987746234039f0
end

# run the second timestep
model = Wflow.update(model)

@testset "second timestep sediment model (vertical)" begin
    eros = model.vertical

    @test mean(eros.soilloss) ≈ 0.07441259472436522f0
    @test mean(eros.erosclay) ≈ 0.002198312957830507f0
    @test mean(eros.erossand) ≈ 0.02515898751801652f0
    @test mean(eros.eroslagg) ≈ 0.021613617766119256f0
    @test mean(eros.TCsed) == 0.0
    @test mean(eros.TCsilt) ≈ 1.0988158364353527f6
    @test mean(eros.TCsand) ≈ 1.0987090622888755f6
    @test mean(eros.TCclay) ≈ 1.0992655197016734f6
    @test mean(eros.TCsilt) ≈ 1.0988158364353527f6
end

@testset "second timestep sediment model (lateral)" begin
    lat = model.lateral

    @test mean(lat.land.inlandsed) ≈ 0.0727443455767198f0
    @test mean(lat.land.inlandclay) ≈ 0.002190037611542877f0
    @test mean(lat.land.inlandsand) ≈ 0.024512058443074893f0
    @test mean(lat.land.olclay) ≈ 0.006212441628146587f0

    @test mean(lat.river.SSconc) ≈ 0.8454728220852514f0
    @test mean(lat.river.inlandclay) ≈ 0.019393435838709512f0
    @test lat.river.h_riv[network.river.order[end]] ≈ 0.006103649735450745f0
    @test lat.river.outclay[5649] ≈ 3.3447705344961993f-6
end

Wflow.close_files(model)

benchmark = @benchmark Wflow.run(tomlpath)
trialmin = BenchmarkTools.minimum(benchmark)
println("Sediment Model update (run)")
print_benchmark(trialmin)
