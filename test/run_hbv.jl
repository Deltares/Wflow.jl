tomlpath = joinpath(@__DIR__, "hbv_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_hbv_model(config)

model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@test String(read(model.writer.csv_path)) ==
    """time,Q,temp_bycoord,temp_byindex,Q_1,perc_33,perc_34,perc_35,perc_36,perc_37
    2000-01-01T00:00:00,1186.2458521037538,2.965437173843384,1.1716821193695068,1127.6471572748078,2.308000087738037,1.8980000019073486,2.7100000381469727,3.818000078201294,2.1440000534057617
    """

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

"Prints a benchmark results just like btime"
function print_benchmark(trialmin)
    trialtime = BenchmarkTools.time(trialmin)
    trialallocs = BenchmarkTools.allocs(trialmin)
    println(
        "  ",
        BenchmarkTools.prettytime(trialtime),
        " (",
        trialallocs,
        " allocation",
        trialallocs == 1 ? "" : "s",
        ": ",
        BenchmarkTools.prettymemory(BenchmarkTools.memory(trialmin)),
        ")",
    )
end

trialmin = BenchmarkTools.minimum(benchmark)

println("HBV Model update")
print_benchmark(trialmin)

Wflow.close_files(model, delete_output = true)