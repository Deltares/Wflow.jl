
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
@unpack network = model

model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    csv = CSV.File(model.writer.csv_path)
    row = csv[1]

    @test row.time == DateTime("2000-01-01T00:00:00")
    @test row.Q ≈ 7.815587720999557
    @test row.volume ≈ 4.364358112965045e7
    @test row.temp_bycoord ≈ 2.3279826641082764
    @test row.temp_byindex ≈ 2.3279826641082764
    @test row.Q_6336050 ≈ 0.023884907014593577
    @test row.Q_6336510 ≈ 0.012411069390755603
    @test row.Q_6836100 ≈ 0.0048516436912648285
    @test row.Q_6336500 ≈ 0.01183784263707513
    @test row.Q_6836190 ≈ 0.0005303931706115078
    @test row.Q_6336800 ≈ 0.013467698526109831
    @test row.Q_6336900 ≈ 0.0034112570102534556
    @test row.Q_6336930 ≈ 0.09773955491544312
    @test row.Q_6336910 ≈ 0.002147610203630289
    @test row.Q_6336920 ≈ 0.0026493983494987644
    @test row.Q_6136100 ≈ 0.0008708139436150139
    @test row.Q_6136500 ≈ 0.000729148906480041
    @test row.Q_6136520 ≈ 0.002155395279410574
    @test row.Q_6136150 ≈ 0.0022298330217711665
    @test row.Q_6136151 ≈ 0.0031045029537777745
    @test row.Q_6136160 ≈ 3.3424757084395074
    @test row.Q_6136200 ≈ 1.3585171944202
    @test row.Q_6136201 ≈ 5.942459422118659
    @test row.Q_6136202 ≈ 1.6809973032577248
    @test row.recharge_1 ≈ -0.02741617005700012
end

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
    @test sbm.soilevap[1] == 0.0021875758041212754
    @test sbm.snow[1] == 0.003959971215481306
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 6.366945307450936e16
    @test ssf[network.land.order[1]] ≈ 3.0449782003445332e13
    @test ssf[network.land.order[end-100]] ≈ 7.855917592164835e11 
    @test ssf[network.land.order[end]] ≈ 2.1612469198596365e11
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 319.54430079719464
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 1.2912357005050196e-5
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 2810.2801278933775
    @test q[4061] ≈ 0.0016281281135651462
    @test q[5617] ≈  7.333524491261813
    @test q[network.river.order[end]] ≈  0.006523133873784731
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    @test res.outflow[2] ≈ 0.2174998592483153
    @test res.inflow[2] ≈ 0.00022354279860698295
    @test res.volume[2] ≈  2.776155238441744e7
    @test res.precipitation[2] ≈ 0.1765228509902954
    @test res.evaporation[2] ≈ 0.5372688174247742
end

benchmark = @benchmark Wflow.update(model)
trialmin = BenchmarkTools.minimum(benchmark)

println("SBM Model update")
print_benchmark(trialmin)
# @profview Wflow.update(model)

Wflow.close_files(model, delete_output = false)
