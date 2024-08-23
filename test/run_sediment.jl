using Wflow

tomlpath = joinpath(@__DIR__, "sediment_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)
@unpack network = model

model = Wflow.run_timestep(model)

@testset "first timestep sediment model (vertical)" begin
    eros = model.vertical

    @test eros.erosov[1] ≈ 0.9f0
    @test model.clock.iteration == 1
    @test mean(eros.leaf_area_index) ≈ 1.7120018886212223f0
    @test eros.dmsand[1] == 200.0f0
    @test eros.dmlagg[1] == 500.0f0
    @test mean(eros.interception) ≈ 0.4767846753916875f0
    @test mean(eros.soilloss) ≈ 0.008596682196335555f0
end

# run the second timestep
model = Wflow.run_timestep(model)

@testset "second timestep sediment model (vertical)" begin
    eros = model.vertical

    @test mean(eros.soilloss) ≈ 0.07601393657280235f0
    @test mean(eros.erosclay) ≈ 0.0022388720384961766f0
    @test mean(eros.erossand) ≈ 0.02572046244407882f0
    @test mean(eros.eroslagg) ≈ 0.022126541806118796f0
    @test mean(eros.TCsed) == 0.0
    @test mean(eros.TCsilt) ≈ 1.0988158364353527f6
    @test mean(eros.TCsand) ≈ 1.0987090622888755f6
    @test mean(eros.TCclay) ≈ 1.0992655197016734f6
    @test mean(eros.TCsilt) ≈ 1.0988158364353527f6
end

@testset "second timestep sediment model (lateral)" begin
    lat = model.lateral

    @test mean(lat.land.inlandsed) ≈ 0.07463801685030906f0
    @test mean(lat.land.inlandclay) ≈ 0.0022367786781657497f0
    @test mean(lat.land.inlandsand) ≈ 0.02519222037812127f0
    @test mean(lat.land.olclay) ≈ 0.006443036462118322f0

    @test mean(lat.river.SSconc) ≈ 0.8259993252994058f0
    @test mean(lat.river.inlandclay) ≈ 0.01980468760667709f0
    @test lat.river.h_riv[network.river.order[end]] ≈ 0.006103649735450745f0
    @test lat.river.outclay[5649] ≈ 2.359031898208781f-9
end

@testset "Exchange and grid location sediment" begin
    @test Wflow.exchange(model.vertical, :n) == 0
    @test Wflow.exchange(model.vertical, :erosk) == 1
    @test Wflow.exchange(model.vertical, :leaf_area_index) == 1
    @test Wflow.grid_location(model.vertical, :n) == "none"
    @test Wflow.grid_location(model.vertical, :erosk) == "node"
    @test Wflow.grid_location(model.vertical, :leaf_area_index) == "node"
    land = model.lateral.land
    @test Wflow.exchange(land, :n) == 0
    @test Wflow.exchange(land, :soilloss) == 1
    @test Wflow.exchange(land, :inlandsed) == 1
    @test Wflow.grid_location(land, :n) == "none"
    @test Wflow.grid_location(land, :soilloss) == "node"
    @test Wflow.grid_location(land, :inlandsed) == "node"
end

@testset "Model type equals sediment model" begin
    @test isa(model.type, Wflow.SedimentModel)
end

@testset "Set initial conditions from state file" begin
    # Given

    config_ = Wflow.Config(tomlpath)
    setproperty!(config_, :reinit, false)

    # When

    model_ = Wflow.initialize_sediment_model(config_)
    model__ = Wflow.set_states(model_)

    # Then

    @test model_ === model__
end

@testset "Overland flow sediment type equals Float" begin
    # Given

    number_of_cells = 100
    river = rand(Bool, number_of_cells)

    # When

    overland_flow_sediment = Wflow.init_overland_flow_sediment(river, number_of_cells)

    # Then

    @test eltype(overland_flow_sediment.rivcell) == Float
end

@testset "Overland flow sediment size" begin
    # Given

    number_of_cells = 100
    river = rand(Bool, number_of_cells)

    # When

    overland_flow_sediment = Wflow.init_overland_flow_sediment(river, number_of_cells)

    # Then

    @test length(overland_flow_sediment.rivcell) == number_of_cells
end

@testset "Reinit false" begin
    # Given

    config_ = Wflow.Config(tomlpath)

    # When

    model_ = Wflow.initialize_sediment_model(config_)
    model_.config.model.reinit = false

    # Then

    @test_throws ErrorException Wflow.set_states(model_)
end

Wflow.close_files(model)
