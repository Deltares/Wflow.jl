tomlpath = joinpath(@__DIR__, "hbv_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_hbv_model(config)

model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@test String(read(model.writer.csv_path)) ==
    """time,Q,temp_bycoord,temp_byindex,Q_6336050,Q_6336510,Q_6836100,Q_6336500,Q_6836190,Q_6336800,Q_6336900,Q_6336930,Q_6336910,Q_6336920,Q_6136100,Q_6136500,Q_6136520,Q_6136150,Q_6136151,Q_6136160,Q_6136200,Q_6136201,Q_6136202,perc_1
    2000-01-01T00:00:00,896.1751154119561,2.3058688640594482,2.6774086952209473,115.93172210612857,60.615448841211,143.80408282830942,798.5125227656837,2.0690719219390803,394.62570300032155,566.5236099579497,163.63983453793338,169.0817549943326,272.4644545041471,93.15887192285689,193.0293181442796,75.95560312280786,148.01817057572595,92.67276627775239,101.7418069076657,230.92470049632183,160.2590220254171,46.85273198585401,0.39999999999999913
    """

@testset "first timestep" begin
    hbv = model.vertical
    @test hbv.tt[1] ≈ 1.2999999523162842
    @test model.clock.iteration == 2
    @test hbv.altitude[1] == 643.5469970703125
    @test hbv.soilmoisture[1] == 260.0
    @test hbv.runoff[1] == 52.05567251345292
    @test hbv.soilevap[1] == 0.0
    @test hbv.snow[1] == 0.07686296902988943
end

# run the second timestep
model = Wflow.update(model)

@testset "second timestep" begin
    hbv = model.vertical
    @test hbv.altitude[1] == 643.5469970703125
    @test hbv.soilmoisture[1] == 260.0
    @test hbv.runoff[1] == 0.33563371915719026
    @test hbv.soilevap[1] == 0.0
    @test hbv.snow[1] == 0.0
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 790.3699216973902
    @test q[500] ≈ 0.0034896147099560204
    @test q[1566] ≈ 0.0032776071559327265
    @test q[network.land.order[end]] ≈ 0.005884784829253675
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 30779.05923341059
    @test q[500] ≈ 5.959907702746569
    @test q[88] ≈ 1.3929443368514935
    @test q[network.river.order[end]] ≈ 3.6500006913919645
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