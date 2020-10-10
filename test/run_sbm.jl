tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
@unpack network = model

model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@test String(read(model.writer.csv_path)) ==
    """time,Q,volume,temp_bycoord,temp_byindex,Q_6336050,Q_6336510,Q_6836100,Q_6336500,Q_6836190,Q_6336800,Q_6336900,Q_6336930,Q_6336910,Q_6336920,Q_6136100,Q_6136500,Q_6136520,Q_6136150,Q_6136151,Q_6136160,Q_6136200,Q_6136201,Q_6136202,recharge_1
    2000-01-01T00:00:00,6.024193921311834,4.36446822871985e7,2.3279826641082764,2.3279826641082764,0.028863190668213575,0.007553259389301421,0.004587086012100417,0.007067416841214618,0.0004304131296893515,0.010064083450515805,0.0024645807204221695,0.09743429180041539,0.0019795547240291016,0.0018154706487719095,0.0012096620479210309,0.00114094172533456,0.0022239652735204685,0.002271365545812738,0.002012939917045721,3.069066166182202,0.6138784778491593,3.7635826401122445,1.0097005458225194,-0.057306898958161656
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
    @test sbm.soilevap[1] == 0.005865653627803281
    @test sbm.snow[1] == 0.009696763863612956
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 7.01804822277135e16
    @test ssf[network.land.order[1]] ≈ 3.0449782003445332e13
    @test ssf[network.land.order[end-100]] ≈ 7.880162191955529e11 
    @test ssf[network.land.order[end]] ≈ 3.289417561401221e11
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 344.3887374333076
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 2.4104063304010105e-6
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 2589.14416681887
    @test q[4061] ≈ 0.0018793927768946056
    @test q[5617] ≈ 9.065120107547822
    @test q[network.river.order[end]] ≈ 0.005616749347427413
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    @test res.outflow[2] ≈ 0.2174998592483153
    @test res.inflow[2] ≈ 0.00023469861436165696
    @test res.volume[2] ≈ 2.776155347023747e7
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

println("SBM Model update")
print_benchmark(trialmin)
# @profview Wflow.update(model)

Wflow.close_files(model, delete_output = true)
