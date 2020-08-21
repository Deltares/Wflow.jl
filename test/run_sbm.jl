tomlpath = joinpath(@__DIR__, "config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
@unpack network = model

nl = length(network.land.order)
index_river = filter(i -> !isequal(model.lateral.river.rivercells[i], 0), 1:nl)
frac_toriver = Wflow.fraction_runoff_toriver(
    network.land.graph,
    index_river,
    model.lateral.subsurface.βₗ,
    nl,
)

model = Wflow.update(model, frac_toriver, index_river)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@test String(read(model.writer.csv_path)) ==
      "time,Q,volume,precipitation\n2000-01-01T00:00:00,8.068415620162373,4.364435435788508e7,2.016469537990403\n"

@testset "first timestep" begin
    sbm = model.vertical

    @test sbm.tt[1] ≈ 1.2999999523162842

    @test model.clock.iteration == 2

    @test sbm.altitude[1] == 643.5469970703125
    @test sbm.θₛ[1] == 0.48343977332115173
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.0
    @test sbm.snow[1] == 0.6029994894902302
end

# run the second timestep
model = Wflow.update(model, frac_toriver, index_river)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.altitude[1] == 643.5469970703125
    @test sbm.θₛ[1] == 0.48343977332115173
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.00586565362774996
    @test sbm.snow[1] == 0.009696763863612956
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 7.005489495052358e16
    @test ssf[network.land.order[1]] ≈ 3.0449782003445332e13
    @test ssf[network.land.order[nl-100]] ≈ 7.87333555063647e11
    @test ssf[network.land.order[end]] ≈ 3.289417561401221e11
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 6.1333024054146446
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 6.804317844228025e-6
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 671.4362633660646
    @test q[4061] ≈ 0.011731014256037769
    @test q[5617] ≈ 1.0080691890749913
    @test q[network.river.order[end]] ≈ 0.0043612246315573805
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    @test res.outflow[2] ≈ 0.2174998580279266
    @test res.inflow[2] ≈ 0.004452979665488473
    @test res.volume[2] ≈ 2.7745799188635208e7
    @test res.precipitation[2] ≈ 3.0
    @test res.evaporation[2] ≈ 4.0
end

benchmark = @benchmark Wflow.update(model, frac_toriver, index_river)

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

println("Model update")
print_benchmark(trialmin)
# @profview Wflow.update(model, frac_toriver, index_river)

Wflow.close_files(model, delete_output = true)
