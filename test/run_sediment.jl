using Wflow

tomlpath = joinpath(@__DIR__, "sediment_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)
@unpack network = model

model = Wflow.run_timestep(model)

@testset "first timestep sediment model (vertical)" begin
    eros = model.vertical

    @test eros.hydrometeo_forcing.precipitation[1] ≈ 4.086122035980225f0
    @test eros.hydrometeo_forcing.q_land[1] ≈ 0.0f0
    @test eros.overland_flow_erosion.parameters.usle_k[1] ≈ 0.026510488241910934f0
    @test eros.overland_flow_erosion.parameters.usle_c[1] ≈ 0.014194443821907043f0
    @test eros.overland_flow_erosion.parameters.answers_k[1] ≈ 0.9f0
    @test eros.overland_flow_erosion.variables.amount[1] ≈ 0.0f0
    @test eros.rainfall_erosion.variables.amount[1] ≈ 0.00027245577922893746f0
    @test model.clock.iteration == 1
    #@test mean(eros.leaf_area_index) ≈ 1.7120018886212223f0
    #@test mean(eros.interception) ≈ 0.4767846753916875f0
    @test mean(eros.overland_flow_erosion.variables.amount) ≈ 0.00861079076689589f0
    @test mean(eros.rainfall_erosion.variables.amount) ≈ 0.00016326203201620437f0
    @test mean(eros.soil_erosion.variables.amount) ≈ 0.008774052798912092f0
end

# run the second timestep
model = Wflow.run_timestep(model)

@testset "second timestep sediment model (vertical)" begin
    eros = model.vertical

    @test mean(eros.soil_erosion.variables.amount) ≈ 0.07765800489746684f0
    @test mean(eros.soil_erosion.variables.clay) ≈ 0.002287480354866626f0
    @test mean(eros.soil_erosion.variables.silt) ≈ 0.0036164773521184896f0
    @test mean(eros.soil_erosion.variables.sand) ≈ 0.026301393837924607f0
    @test mean(eros.soil_erosion.variables.lagg) ≈ 0.022577957752547836f0
    @test mean(eros.soil_erosion.variables.sagg) ≈ 0.022874695590802723f0
end

@testset "second timestep sediment model (lateral)" begin
    land = model.lateral.land
    river = model.lateral.river

    @test land.transport_capacity.parameters.dm_sand[1] == 200.0f0
    @test land.transport_capacity.parameters.dm_lagg[1] == 500.0f0

    @test mean(land.transport_capacity.boundary_conditions.q) ≈ 0.006879398771052133f0
    @test mean(land.transport_capacity.variables.silt) ≈ 1.0988158364353527f6
    @test mean(land.transport_capacity.variables.sand) ≈ 1.0987090622888755f6
    @test mean(land.transport_capacity.variables.clay) ≈ 1.0992655197016734f6

    @test mean(land.to_river.variables.amount) ≈ 0.07624135182616738f0
    @test sum(land.to_river.variables.clay) ≈ 114.42704329506047f0
    @test sum(land.to_river.variables.sand) ≈ 1289.4785484597958f0
    @test mean(land.sediment_flux.variables.clay) ≈ 0.006578791733506439f0

    @test mean(river.hydrometeo_forcing.q_river) ≈ 0.6975180562953642f0
    @test river.hydrometeo_forcing.waterlevel_river[network.river.order[end]] ≈
          0.006103649735450745f0
    @test mean(river.geometry.width) ≈ 22.628250814095523f0

    @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642f0
    @test mean(river.transport_capacity.variables.amount) ≈ 0.4458019733090582f0
    @test mean(river.potential_erosion.variables.bed) ≈ 307.18492138827116f0

    @test sum(river.sediment_flux.boundary_conditions.erosion_land_clay) ≈
          114.42704329506047f0
    @test sum(river.sediment_flux.boundary_conditions.erosion_land_sand) ≈
          1289.4785484597958f0
    @test mean(river.sediment_flux.boundary_conditions.transport_capacity) ≈
          0.4458019733090582f0
    @test mean(river.sediment_flux.variables.amount) ≈ 0.4333483865969662f0
    @test mean(river.sediment_flux.variables.erosion) ≈ 0.019077695621351014f0
    @test mean(river.sediment_flux.variables.deposition) ≈ 0.6941274181387916f0
    @test river.sediment_flux.variables.clay[5649] ≈ 2.840979764480952f-9

    @test mean(river.concentrations.variables.suspended) ≈ 0.8260083257660087f0
end

# @testset "Exchange and grid location sediment" begin
#     @test Wflow.exchange(model.vertical, :n) == 0
#     @test Wflow.exchange(model.vertical, :erosk) == 1
#     @test Wflow.exchange(model.vertical, :leaf_area_index) == 1
#     @test Wflow.grid_location(model.vertical, :n) == "none"
#     @test Wflow.grid_location(model.vertical, :erosk) == "node"
#     @test Wflow.grid_location(model.vertical, :leaf_area_index) == "node"
#     land = model.lateral.land
#     @test Wflow.exchange(land, :n) == 0
#     @test Wflow.exchange(land, :soilloss) == 1
#     @test Wflow.exchange(land, :inlandsed) == 1
#     @test Wflow.grid_location(land, :n) == "none"
#     @test Wflow.grid_location(land, :soilloss) == "node"
#     @test Wflow.grid_location(land, :inlandsed) == "node"
# end

Wflow.close_files(model)
