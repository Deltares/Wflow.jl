
tomlpath = joinpath(@__DIR__, "hbv_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_hbv_model(config)
@unpack network = model

model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2000-01-01T00:00:00")
    @test row.Q ≈ 517.9761666607884f0
    @test row.temp_bycoord ≈ 2.965437173843384f0
    @test row.temp_byindex ≈ 1.1716821193695068f0
    @test row.Q_1 ≈ 501.6252231440609f0
    @test row.perc_33 ≈ 2.308000087738037f0
    @test row.perc_34 ≈ 1.8980000019073486f0
    @test row.perc_35 ≈ 2.7100000381469727f0
    @test row.perc_36 ≈ 3.818000078201294f0
    @test row.perc_37 ≈ 2.1440000534057617f0
end

@testset "first timestep" begin
    hbv = model.vertical
    @test hbv.tt[4377] ≈ 0.0
    @test model.clock.iteration == 2
    @test hbv.soilmoisture[4377] ≈ 134.35299682617188f0
    @test hbv.runoff[4377] ≈ 7.406898120121746f0
    @test hbv.soilevap[4377] == 0.0
    @test hbv.snow[4377] == 0.0
end

# run the second timestep
model = Wflow.update(model)

@testset "second timestep" begin
    hbv = model.vertical
    @test hbv.soilmoisture[4377] ≈ 134.35299682617188f0
    @test hbv.runoff[4377] ≈ 4.3533463f0
    @test hbv.soilevap[4377] == 0.0
    @test hbv.snow[4377] == 0.0
end

@testset "overland domain" begin
    q = model.lateral.land.q_av
    land = model.lateral.land
    @test sum(q) ≈ 3278.9423039643552f0
    @test q[10354] ≈ 0.2348646912841589f0
    @test land.volume[10354] ≈ 2057.7314802432425f0
    @test land.inwater[10354] ≈ 0.027351853491789667f0
    @test q[network.land.order[end]] ≈ 0.296794455575309f0
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 51943.24072462076f0
    @test q[651] ≈ 5.693405508134953f0
    @test q[1056] ≈ 8.839306661545436f0
    @test q[network.river.order[end]] ≈ 363.5647479049738f0
end

Wflow.close_files(model, delete_output = false)
