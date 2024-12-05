using Wflow

tomlpath = joinpath(@__DIR__, "sediment_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)
(; network) = model

Wflow.run_timestep!(model)

@testset "first timestep sediment model (vertical)" begin
    eros = model.vertical

    @test eros.atmospheric_forcing.precipitation[1] ≈ 4.086122035980225f0
    @test eros.hydrological_forcing.q_land[1] ≈ 0.0f0
    @test eros.overland_flow_erosion.parameters.usle_k[1] ≈ 0.026510488241910934f0
    @test eros.overland_flow_erosion.parameters.usle_c[1] ≈ 0.014194443821907043f0
    @test eros.overland_flow_erosion.parameters.answers_k[1] ≈ 0.9f0
    @test eros.overland_flow_erosion.variables.amount[1] ≈ 0.0f0
    @test eros.rainfall_erosion.variables.amount[1] ≈ 0.00027245577922893746f0
    @test model.clock.iteration == 1
    @test mean(eros.overland_flow_erosion.variables.amount) ≈ 0.00861079076689589f0
    @test mean(eros.rainfall_erosion.variables.amount) ≈ 0.00016326203201620437f0
    @test mean(eros.soil_erosion.variables.amount) ≈ 0.008774052798912092f0
end

# run the second timestep
Wflow.run_timestep!(model)

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

    @test mean(river.hydrological_forcing.q_river) ≈ 0.6975180562953642f0
    @test river.hydrological_forcing.waterlevel_river[network.river.order[end]] ≈
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

Wflow.close_files(model)

### Test the sediment model with a different configuration file ###
tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sediment_model(config)
(; network) = model

Wflow.run_timestep!(model)

@testset "first timestep sediment model eurosem (vertical)" begin
    eros = model.vertical

    @test eros.atmospheric_forcing.precipitation[1] ≈ 4.086122035980225f0
    @test eros.hydrological_forcing.interception[1] ≈ 0.6329902410507202f0
    @test eros.hydrological_forcing.q_land[1] ≈ 0.0f0
    @test eros.rainfall_erosion.parameters.soil_detachability[1] ≈ 2.0f0
    @test eros.rainfall_erosion.parameters.eurosem_exponent[1] ≈ 2.0f0
    @test eros.overland_flow_erosion.parameters.usle_c[1] ≈ 0.014194443821907043f0
    @test eros.overland_flow_erosion.variables.amount[1] ≈ 0.0f0
    @test eros.rainfall_erosion.variables.amount[1] ≈ 0.01232301374083337f0
    @test mean(eros.overland_flow_erosion.variables.amount) ≈ 0.00861079076689589f0
    @test mean(eros.rainfall_erosion.variables.amount) ≈ 0.0014726364432116048f0
    @test mean(eros.soil_erosion.variables.amount) ≈ 0.010083427210107495f0
end

# run the second timestep
Wflow.run_timestep!(model)

@testset "second timestep sediment model engelund (lateral)" begin
    land = model.lateral.land
    river = model.lateral.river

    @test river.transport_capacity.parameters.d50[1] == 0.05000000074505806f0
    @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642f0
    @test mean(river.transport_capacity.variables.amount) ≈ 0.14184859055736687f0

    @test mean(river.concentrations.variables.suspended) ≈ 0.24788001458305775f0
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
    eros = model.vertical
    land = model.lateral.land

    @test mean(eros.soil_erosion.variables.amount) ≈ 0.0776983847440198f0
    @test mean(land.transport_capacity.parameters.c_govers) ≈ 0.16393911236592437f0
    @test mean(land.transport_capacity.variables.amount) ≈ 1.0988158364353527f6
    @test mean(land.to_river.variables.amount) ≈ 0.07708434959918917f0
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
    eros = model.vertical
    land = model.lateral.land

    @test mean(eros.soil_erosion.variables.amount) ≈ 0.0776983847440198f0
    @test mean(land.transport_capacity.parameters.d50) ≈ 0.001534350291334408f0
    @test mean(land.transport_capacity.variables.amount) ≈ 1.0988158364353527f6
    @test mean(land.to_river.variables.amount) ≈ 0.07759391658531742f0
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

@testset "second timestep sediment model yang (lateral)" begin
    river = model.lateral.river

    @test river.transport_capacity.parameters.d50[1] == 0.05000000074505806f0
    @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642f0
    @test mean(river.transport_capacity.variables.amount) ≈ 39.959093179632234f0

    @test mean(river.concentrations.variables.suspended) ≈ 0.004036947448899768f0
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

@testset "second timestep sediment model kodatie (lateral)" begin
    river = model.lateral.river

    @test river.transport_capacity.parameters.a_kodatie[1] == 2829.6
    @test river.transport_capacity.parameters.b_kodatie[1] == 3.646
    @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642f0
    @test mean(river.transport_capacity.variables.amount) ≈ 30.332588671299625f0

    @test mean(river.concentrations.variables.suspended) ≈ 54.75483621753514f0
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

@testset "second timestep sediment model molinas (lateral)" begin
    river = model.lateral.river

    @test river.transport_capacity.parameters.d50[1] == 0.05000000074505806f0
    @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642f0
    @test mean(river.transport_capacity.variables.amount) ≈ 350.6483600591209f0

    @test mean(river.concentrations.variables.suspended) ≈ 884.4249630198293f0
end

Wflow.close_files(model)