
tomlpath = joinpath(@__DIR__, "flextopo_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_flextopo_model(config)
(; network) = model

model = Wflow.run_timestep(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2010-01-01T00:00:00")
    @test row.Q_16 ≈ 0.0000370817920883859f0
    @test row.percentageH_16 ≈ 0.340213f0
    @test row.QfP_1011 ≈ 0.0f0
    @test row.kfW_10 ≈ 0.136334478855133f0
    @test row.kfH_10 ≈ 0.0454448275268077f0
    @test row.Qs_503 ≈ 0.1575031550601124f0
    @test row.Si_16 ≈ 0.0f0
    @test row.Ss_16 ≈ 29.89794354323524f0
    @test row.Ea_16 ≈ 0.0110627892408438f0
end

@testset "first timestep" begin
    flextopo = model.vertical
    @test flextopo.tt[3500] ≈ 1.3f0
    @test model.clock.iteration == 1
    @test flextopo.rootzonestorage[3500] ≈
          [147.11238663084805f0, 79.08369375691255f0, 79.23637697443984f0]
    @test flextopo.runoff[3500] ≈ 0.19008129369467497f0
    @test flextopo.rootevap[3500] ≈
          [0.009225917980074883f0, 0.009225917980074883f0, 0.009225917980074883f0]
    @test flextopo.snow[3500] == 0.0
end

# run the second timestep
model = Wflow.run_timestep(model)

@testset "second timestep" begin
    flextopo = model.vertical
    @test flextopo.rootzonestorage[3500] ≈
          [146.73222385533154f0, 78.55190222246516f0, 78.85739340140256f0]
    @test flextopo.runoff[3500] ≈ 0.18887692333975817f0
    @test flextopo.rootevap[3500] ≈
          [0.38016277551651f0, 0.38016277551651f0, 0.38016277551651f0]
    @test flextopo.snow[3500] == 0.0f0
end

@testset "overland domain" begin
    q = model.lateral.land.q_av
    land = model.lateral.land
    @test sum(q) ≈ 202.1162497378859f0
    @test q[10354] ≈ 0.001624108000162103f0
    @test land.volume[10354] ≈ 66.41398729880544f0
    @test land.inwater[10354] ≈ 0.001658182220814705f0
    @test q[network.land.order[end]] ≈ 0.0033337667005815565f0
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 108.45732176546284f0
    @test q[651] ≈ 0.014755293483684526f0
    @test q[1056] ≈ 0.0030933658056867823f0
    @test q[network.river.order[end]] ≈ 0.0005541137960384742f0
end

#change the config to use other functions for the storages for several classes:
config = Wflow.Config(tomlpath)
config["model"]["select_snow"] = ["common_snow_no_storage"]
config["model"]["select_interception"] =
    ["interception_no_storage", "interception_overflow", "interception_overflow"]
config["model"]["select_hortonponding"] =
    ["hortonponding_no_storage", "hortonponding", "hortonponding_no_storage"]
config["model"]["select_hortonrunoff"] =
    ["hortonrunoff_no_storage", "hortonrunoff", "hortonrunoff_no_storage"]
config["model"]["select_rootzone"] =
    ["rootzone_storage", "rootzone_storage", "rootzone_no_storage"]
config["model"]["select_fast"] = ["fast_no_storage", "fast_storage", "fast_storage"]
config["model"]["select_slow"] = ["common_slow_storage"]

model = Wflow.initialize_flextopo_model(config)
(; network) = model

model = Wflow.run_timestep(model)

@testset "first timestep" begin
    flextopo = model.vertical
    @test flextopo.rootzonestorage[3500] ≈
          [147.11238663084805f0, 79.08369375691255f0, 0.0f0]
    @test flextopo.runoff[3500] ≈ 0.19008129369467497f0
    @test flextopo.rootevap[3500] ≈ [0.009225917980074883f0, 0.009225917980074883f0, 0.0f0]
    @test flextopo.snow[3500] == 0.0f0
end

@testset "overland domain" begin
    q = model.lateral.land.q_av
    land = model.lateral.land
    @test sum(q) ≈ 69.51434665542678f0
    @test q[10354] ≈ 0.0009376627801121855f0
    @test land.volume[10354] ≈ 63.46997463442048f0
    @test land.inwater[10354] ≈ 0.0016722689680105701f0
    @test q[network.land.order[end]] ≈ 0.0006175366185364747f0
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 22.236556345676156f0
    @test q[651] ≈ 0.00035906990499781237f0
    @test q[1056] ≈ 0.0002647264689318855f0
    @test q[network.river.order[end]] ≈ 0.000015647838277869928f0
end

Wflow.close_files(model, delete_output = false)
