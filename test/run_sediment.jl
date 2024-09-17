using Wflow

tomlpath = joinpath(@__DIR__, "sediment_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)
(; network) = model

Wflow.run_timestep!(model)

@testset "first timestep sediment model (vertical)" begin
    eros = model.vertical

    @test eros.hydrometeo_forcing.q_land[1] ≈ 0.0f0
    @test eros.overland_flow_erosion.parameters.usle_k[1] ≈ 0.026510488241910934
    @test eros.overland_flow_erosion.variables.amount[1] ≈ 0.9f0
    @test model.clock.iteration == 1
    #@test mean(eros.leaf_area_index) ≈ 1.7120018886212223f0
    #@test eros.dmsand[1] == 200.0f0
    #@test eros.dmlagg[1] == 500.0f0
    #@test mean(eros.interception) ≈ 0.4767846753916875f0
    @test mean(eros.soil_erosion.variables.amount) ≈ 0.008596682196335555f0
end

# # run the second timestep
# model = Wflow.run_timestep(model)

# @testset "second timestep sediment model (vertical)" begin
#     eros = model.vertical

#     @test mean(eros.soil_erosion.variables.amount) ≈ 0.07601393657280235f0
#     @test mean(eros.soil_erosion.variables.clay) ≈ 0.0022388720384961766f0
#     @test mean(eros.soil_erosion.variables.sand) ≈ 0.02572046244407882f0
#     @test mean(eros.soil_erosion.variables.lagg) ≈ 0.022126541806118796f0
#     #@test mean(eros.TCsed) == 0.0
#     #@test mean(eros.TCsilt) ≈ 1.0988158364353527f6
#     #@test mean(eros.TCsand) ≈ 1.0987090622888755f6
#     #@test mean(eros.TCclay) ≈ 1.0992655197016734f6
#     #@test mean(eros.TCsilt) ≈ 1.0988158364353527f6
# end

# @testset "second timestep sediment model (lateral)" begin
#     lat = model.lateral

#     @test mean(lat.land.inlandsed) ≈ 0.07463801685030906f0
#     @test mean(lat.land.inlandclay) ≈ 0.0022367786781657497f0
#     @test mean(lat.land.inlandsand) ≈ 0.02519222037812127f0
#     @test mean(lat.land.olclay) ≈ 0.006443036462118322f0

#     @test mean(lat.river.SSconc) ≈ 0.8259993252994058f0
#     @test mean(lat.river.inlandclay) ≈ 0.01980468760667709f0
#     @test lat.river.h_riv[network.river.order[end]] ≈ 0.006103649735450745f0
#     @test lat.river.outclay[5649] ≈ 2.359031898208781f-9
# end

@testset "Exchange and grid location sediment" begin
    @test Wflow.exchange(model.vertical.n) == false
    @test Wflow.exchange(model.vertical.erosk) == true
    @test Wflow.exchange(model.vertical.leaf_area_index) == true
    @test Wflow.grid_loc(model.vertical, :n) == "none"
    @test Wflow.grid_loc(model.vertical, :erosk) == "node"
    @test Wflow.grid_loc(model.vertical, :leaf_area_index) == "node"
    land = model.lateral.land
    @test Wflow.exchange(land.n) == false
    @test Wflow.exchange(land.soilloss) == true
    @test Wflow.exchange(land.inlandsed) == true
    @test Wflow.grid_loc(land, :n) == "none"
    @test Wflow.grid_loc(land, :soilloss) == "node"
    @test Wflow.grid_loc(land, :inlandsed) == "node"
end

Wflow.close_files(model)
