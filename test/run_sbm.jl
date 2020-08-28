tomlpath = joinpath(@__DIR__, "config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
@unpack network = model

model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@test String(read(model.writer.csv_path)) ==
      """time,Q,volume,temp_bycoord,temp_byindex,actevap_6336050,actevap_6336510,actevap_6836100,actevap_6336500,actevap_6836190,actevap_6336800,actevap_6336900,actevap_6336930,actevap_6336910,actevap_6336920,actevap_6136100,actevap_6136500,actevap_6136520,actevap_6136150,actevap_6136151,actevap_6136160,actevap_6136200,actevap_6136201,actevap_6136202,recharge_1
      2000-01-01T00:00:00,8.068415620162373,4.36467000689324e7,2.3279826641082764,2.3279826641082764,0.1463467459314577,0.13163753168695083,0.09536865090657316,0.35945141852456364,0.09788319701615555,0.10258831234024719,0.12840415736172062,0.10768710651463834,0.14325774157808685,0.11276011822358965,0.0,0.04348357096162693,0.036198475225879365,0.04324215979814671,0.0803324874462755,0.01683602806488171,0.0862277163174339,0.00901282012220647,0.11086802959442675,-46470.570328920454
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
    @test res.outflow[2] ≈ 0.21749982586488784
    @test res.inflow[2] ≈ 0.0003260830320098207
    @test res.volume[2] ≈ 2.7529217030450657e7
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
