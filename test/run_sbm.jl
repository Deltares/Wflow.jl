tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
@unpack network = model

model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@test String(read(model.writer.csv_path)) ==
    """time,Q,volume,temp_bycoord,temp_byindex,Q_6336050,Q_6336510,Q_6836100,Q_6336500,Q_6836190,Q_6336800,Q_6336900,Q_6336930,Q_6336910,Q_6336920,Q_6136100,Q_6136500,Q_6136520,Q_6136150,Q_6136151,Q_6136160,Q_6136200,Q_6136201,Q_6136202,recharge_1
    2000-01-01T00:00:00,6.007648422080027,4.364467694582078e7,2.3279826641082764,2.3279826641082764,0.025687571754422933,0.006349994321123713,0.004052191700768449,0.0025560252457928825,0.0003866259942806916,0.007977854358793413,0.0014150456117394605,0.09113824961483859,0.0015864145544560726,0.0009009598468463504,0.0009794225220044622,0.0005269666056249656,0.0018105753159215901,0.0015987117041130552,0.001445791350630152,3.0599031722043355,0.5927558144862284,3.7318971216186565,1.0036613531867922,-0.05729963601731694
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
    @test sum(q) ≈ 344.38787670362007
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 1.9195652667407746e-6
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 2581.124080145291
    @test q[4061] ≈ 0.0016541532633614464
    @test q[5617] ≈ 9.065257601055963
    @test q[network.river.order[end]] ≈ 0.002578632105794227
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
