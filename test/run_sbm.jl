tomlpath = joinpath(@__DIR__, "config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
@unpack network = model

model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@test String(read(model.writer.csv_path)) ==
      """time,Q,volume,temp_bycoord,temp_byindex,Q_6336050,Q_6336510,Q_6836100,Q_6336500,Q_6836190,Q_6336800,Q_6336900,Q_6336930,Q_6336910,Q_6336920,Q_6136100,Q_6136500,Q_6136520,Q_6136150,Q_6136151,Q_6136160,Q_6136200,Q_6136201,Q_6136202,recharge_1
      2000-01-01T00:00:00,8.068415620162373,4.36467000689324e7,2.3279826641082764,2.3279826641082764,0.011860664734651536,0.0013409825338072177,0.007049333453420427,0.003998799925923699,0.002435666579854784,0.002912431217242948,0.157017137266319,3.6233720323561616,0.0022332698067222467,0.0029604734164055914,0.04276199866399382,0.000784601568706128,0.0013247936870766685,0.001184569211941342,1.2994428393408683,0.0023092625180203027,1.3992592952863065,5.488085457567639,0.011027748385065077,-46470.570328920454
      """
      
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
model = Wflow.update(model)

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
    @test sum(ssf) ≈ 7.005493498683713e16
    @test ssf[network.land.order[1]] ≈ 3.0449782003445332e13
    @test ssf[network.land.order[end-100]] ≈ 7.87333555063647e11
    @test ssf[network.land.order[end]] ≈ 3.289417561401221e11
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 6.1213772525470445
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 6.804317844228025e-6
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 671.3851501260151
    @test q[4061] ≈ 0.011731014256037769
    @test q[5617] ≈ 1.0080691890749913
    @test q[network.river.order[end]] ≈ 0.0043612246315573805
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    @test res.outflow[2] ≈ 0.2174998592483153
    @test res.inflow[2] ≈ 0.0003260830320098207
    @test res.volume[2] ≈ 2.776156241597619e7
    @test res.precipitation[2] ≈ 0.1765228509902954
    @test res.evaporation[2] ≈ 0.5372688174247742
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

println("Model update")
print_benchmark(trialmin)
# @profview Wflow.update(model)

Wflow.close_files(model, delete_output = true)
