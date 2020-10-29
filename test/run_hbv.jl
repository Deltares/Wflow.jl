include("testing_utils.jl")

tomlpath = joinpath(@__DIR__, "hbv_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_hbv_model(config)

model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    csv = CSV.File(model.writer.csv_path)
    row = csv[1]

    @test row.time == DateTime("2000-01-01T00:00:00")
    @test row.Q ≈ 1186.2458521037543
    @test row.temp_bycoord ≈ 2.965437173843384
    @test row.temp_byindex ≈ 1.1716821193695068
    @test row.Q_1 ≈ 1127.6471572748078
    @test row.perc_33 ≈ 2.308000087738037
    @test row.perc_34 ≈ 1.8980000019073486
    @test row.perc_35 ≈ 2.7100000381469727
    @test row.perc_36 ≈ 3.818000078201294
    @test row.perc_37 ≈ 2.1440000534057617    
end

@testset "first timestep" begin
    hbv = model.vertical
    @test hbv.tt[1] ≈ 0.0
    @test model.clock.iteration == 2
    @test hbv.altitude[1] == 484.2149963378906 
    @test hbv.soilmoisture[1] == 134.35299682617188
    @test hbv.runoff[1] == 9.188024381138971
    @test hbv.soilevap[1] == 0.0
    @test hbv.snow[1] == 0.0
end

# run the second timestep
model = Wflow.update(model)

@testset "second timestep" begin
    hbv = model.vertical
    @test hbv.altitude[1] == 484.2149963378906
    @test hbv.soilmoisture[1] == 134.35299682617188
    @test hbv.runoff[1] == 4.578703319110791 
    @test hbv.soilevap[1] == 0.0
    @test hbv.snow[1] == 0.0
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 3503.798504329061
    @test q[500] ≈ 0.2601121395680727
    @test q[1566] ≈ 0.03304875788712168
    @test q[network.land.order[end]] ≈ 0.09511081839582983
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 95238.16617175122
    @test q[500] ≈ 6.6547884120641845
    @test q[88] ≈ 10.138009383376387
    @test q[network.river.order[end]] ≈ 3.2690469368626798 
end

benchmark = @benchmark Wflow.update(model)
trialmin = BenchmarkTools.minimum(benchmark)

println("HBV Model update")
print_benchmark(trialmin)

Wflow.close_files(model, delete_output = true)