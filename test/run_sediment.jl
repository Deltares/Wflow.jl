using Wflow

tomlpath = joinpath(@__DIR__, "sediment_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)
(; network) = model

Wflow.run_timestep!(model)

@testset "first timestep sediment model (land part)" begin
    eros = model.land

    @test eros.atmospheric_forcing.precipitation[1] ≈ 4.086122035980225
    @test eros.hydrological_forcing.q_land[1] ≈ 0.0
    @test eros.overland_flow_erosion.parameters.usle_k[1] ≈ 0.026510488241910934
    @test eros.overland_flow_erosion.parameters.usle_c[1] ≈ 0.014194443821907043
    @test eros.overland_flow_erosion.parameters.answers_k[1] ≈ 0.8999999761581421
    @test eros.overland_flow_erosion.variables.amount[1] ≈ 0.0
    @test eros.rainfall_erosion.variables.amount[1] ≈ 0.00027245577922893746
    @test model.clock.iteration == 1
    @test mean(eros.overland_flow_erosion.variables.amount) ≈ 0.00861079076689589
    @test mean(eros.rainfall_erosion.variables.amount) ≈ 0.00016326203201620437
    @test mean(eros.soil_erosion.variables.amount) ≈ 0.008774052798912092
end

# run the second timestep
Wflow.run_timestep!(model)

@testset "second timestep sediment model (land part)" begin
    eros = model.land

    @test mean(eros.soil_erosion.variables.amount) ≈ 0.07765800489746684
    @test mean(eros.soil_erosion.variables.clay) ≈ 0.002287480354866626
    @test mean(eros.soil_erosion.variables.silt) ≈ 0.0036164773521184896
    @test mean(eros.soil_erosion.variables.sand) ≈ 0.026301393837924607
    @test mean(eros.soil_erosion.variables.lagg) ≈ 0.022577957752547836
    @test mean(eros.soil_erosion.variables.sagg) ≈ 0.022874695590802723
end

@testset "second timestep sediment model (routing)" begin
    land = model.routing.overland_flow
    river = model.routing.river_flow

    @test land.transport_capacity.parameters.dm_sand[1] == 200.0
    @test land.transport_capacity.parameters.dm_lagg[1] == 500.0

    @test mean(land.transport_capacity.boundary_conditions.q) ≈ 0.006879398771052133
    @test mean(land.transport_capacity.variables.silt) ≈ 1.0988158364353527f6
    @test mean(land.transport_capacity.variables.sand) ≈ 1.0987090622888755f6
    @test mean(land.transport_capacity.variables.clay) ≈ 1.0992655197016734f6

    @test mean(land.to_river.variables.amount) ≈ 0.0762386279230294
    @test sum(land.to_river.variables.clay) ≈ 114.42704329506047
    @test sum(land.to_river.variables.sand) ≈ 1289.4173249850346
    @test mean(land.sediment_flux.variables.clay) ≈ 0.006578791733506439

    @test mean(river.hydrological_forcing.q_river) ≈ 0.6975180562953642
    @test river.hydrological_forcing.waterlevel_river[network.river.order[end]] ≈
          0.006103649735450745
    @test mean(river.geometry.width) ≈ 22.628250814095523

    @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
    @test mean(river.transport_capacity.variables.amount) ≈ 0.4458019733090582
    @test mean(river.potential_erosion.variables.bed) ≈ 307.18492138827116

    @test sum(river.sediment_flux.boundary_conditions.erosion_land_clay) ≈
          114.42704329506047
    @test sum(river.sediment_flux.boundary_conditions.erosion_land_sand) ≈ 1289.417324985034
    @test mean(river.sediment_flux.boundary_conditions.transport_capacity) ≈
          0.4458019733090582
    @test mean(river.sediment_flux.variables.amount) ≈ 0.4333483865969662
    @test mean(river.sediment_flux.variables.erosion) ≈ 0.01907770837584654
    @test mean(river.sediment_flux.variables.deposition) ≈ 0.6941033131514452
    @test river.sediment_flux.variables.clay[5649] ≈ 2.840979764480952e-9

    @test mean(river.concentrations.variables.suspended) ≈ 0.8260075500146108
end

Wflow.close_files(model)

### Test the sediment model with a different configuration file ###
tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)
(; network) = model

Wflow.run_timestep!(model)

@testset "first timestep sediment model eurosem (land part)" begin
    eros = model.land

    @test eros.atmospheric_forcing.precipitation[1] ≈ 4.086122035980225
    @test eros.hydrological_forcing.interception[1] ≈ 0.6329902410507202
    @test eros.hydrological_forcing.q_land[1] ≈ 0.0
    @test eros.rainfall_erosion.parameters.soil_detachability[1] ≈ 2.0
    @test eros.rainfall_erosion.parameters.eurosem_exponent[1] ≈ 2.0
    @test eros.overland_flow_erosion.parameters.usle_c[1] ≈ 0.014194443821907043
    @test eros.overland_flow_erosion.variables.amount[1] ≈ 0.0
    @test eros.rainfall_erosion.variables.amount[1] ≈ 0.01232301374083337
    @test mean(eros.overland_flow_erosion.variables.amount) ≈ 0.00861079076689589
    @test mean(eros.rainfall_erosion.variables.amount) ≈ 0.0014726364432116048
    @test mean(eros.soil_erosion.variables.amount) ≈ 0.010083427210107495
end

# run the second timestep
Wflow.run_timestep!(model)

@testset "second timestep sediment model engelund (routing)" begin
    land = model.routing.overland_flow
    river = model.routing.river_flow

    @test river.transport_capacity.parameters.d50[1] == 0.05000000074505806
    @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
    @test mean(river.transport_capacity.variables.amount) ≈ 0.14184859055736687

    @test mean(river.concentrations.variables.suspended) ≈ 0.24788001458305775
end

Wflow.close_files(model)

### Test land only model configuration and transport capacity ###

tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
config = Wflow.Config(tomlpath)
# Update config to run only the land model
config.model.run_river_model = false
# Use govers equation for land transport capacity
config.model.land_transport = "govers"

model = Wflow.initialize_sediment_model(config)
(; network) = model

# run the first and second timestep
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "second timestep sediment model govers" begin
    eros = model.land
    land = model.routing.overland_flow

    @test mean(eros.soil_erosion.variables.amount) ≈ 0.0776983847440198
    @test mean(land.transport_capacity.parameters.c_govers) ≈ 0.16393911236592437
    @test mean(land.transport_capacity.variables.amount) ≈ 1.0988158364353527f6
    @test mean(land.to_river.variables.amount) ≈ 0.07708434959918917
end

Wflow.close_files(model)

tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
config = Wflow.Config(tomlpath)
# Update config to run only the land model
config.model.run_river_model = false
# Use yalin equation for land transport capacity
config.model.land_transport = "yalin"

model = Wflow.initialize_sediment_model(config)
(; network) = model

# run the first and second timestep
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "second timestep sediment model yalin" begin
    eros = model.land
    land = model.routing.overland_flow

    @test mean(eros.soil_erosion.variables.amount) ≈ 0.0776983847440198
    @test mean(land.transport_capacity.parameters.d50) ≈ 0.001534350291334408
    @test mean(land.transport_capacity.variables.amount) ≈ 1.0988158364353527f6
    @test mean(land.to_river.variables.amount) ≈ 0.07759383356462951
end

Wflow.close_files(model)

### Test all river transport capacity ###

tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
config = Wflow.Config(tomlpath)
# Use yang equation for river transport capacity
config.model.river_transport = "yang"

model = Wflow.initialize_sediment_model(config)
(; network) = model

# run the first and second timestep
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "second timestep sediment model yang (routing)" begin
    river = model.routing.river_flow

    @test river.transport_capacity.parameters.d50[1] == 0.05000000074505806
    @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
    @test mean(river.transport_capacity.variables.amount) ≈ 39.959093179632234

    @test mean(river.concentrations.variables.suspended) ≈ 0.004036947448899768
end

Wflow.close_files(model)

tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
config = Wflow.Config(tomlpath)
# Use kodatie equation for river transport capacity
config.model.river_transport = "kodatie"

model = Wflow.initialize_sediment_model(config)
(; network) = model

# run the first and second timestep
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "second timestep sediment model kodatie (routing)" begin
    river = model.routing.river_flow

    @test river.transport_capacity.parameters.a_kodatie[1] == 2829.6
    @test river.transport_capacity.parameters.b_kodatie[1] == 3.646
    @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
    @test mean(river.transport_capacity.variables.amount) ≈ 30.332588671299625

    @test mean(river.concentrations.variables.suspended) ≈ 54.75548200245741
end

Wflow.close_files(model)

tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
config = Wflow.Config(tomlpath)
# Use molinas equation for river transport capacity
config.model.river_transport = "molinas"

model = Wflow.initialize_sediment_model(config)
(; network) = model

# run the first and second timestep
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "second timestep sediment model molinas (routing)" begin
    river = model.routing.river_flow

    @test river.transport_capacity.parameters.d50[1] == 0.05000000074505806
    @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
    @test mean(river.transport_capacity.variables.amount) ≈ 350.6483600591209

    @test mean(river.concentrations.variables.suspended) ≈ 884.4255449829093
end

Wflow.close_files(model)