@testitem "Run sediment" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_config.toml")
    config = Wflow.Config(tomlpath)

    model = Wflow.Model(config)
    (; domain) = model

    @testset "initial default states" begin
        river = model.routing.river_flow
        @test mean(river.sediment_flux.variables.store_clay) == 0.0
        @test mean(river.sediment_flux.variables.leftover_clay) == 0.0
        @test mean(river.sediment_flux.variables.clay_rate) == 0.0

        @test mean(river.sediment_flux.variables.store_silt) == 0.0
        @test mean(river.sediment_flux.variables.leftover_silt) == 0.0
        @test mean(river.sediment_flux.variables.silt_rate) == 0.0

        @test mean(river.sediment_flux.variables.store_sand) == 0.0
        @test mean(river.sediment_flux.variables.leftover_sand) == 0.0
        @test mean(river.sediment_flux.variables.sand_rate) == 0.0
    end

    Wflow.run_timestep!(model)

    @testset "first timestep sediment model (land part)" begin
        eros = model.land

        @test eros.atmospheric_forcing.precipitation[1] ≈ 4.086122035980225
        @test eros.hydrological_forcing.q_land[1] ≈ 0.0
        @test eros.overland_flow_erosion.parameters.usle_k[1] ≈ 0.026510488241910934
        @test eros.overland_flow_erosion.parameters.usle_c[1] ≈ 0.014194443821907043
        @test eros.overland_flow_erosion.parameters.answers_overland_flow_factor[1] ≈
              0.8999999761581421
        @test eros.overland_flow_erosion.variables.soil_erosion_rate[1] ≈ 0.0
        @test eros.rainfall_erosion.variables.soil_erosion_rate[1] ≈ 0.00027245577922893746
        @test model.clock.iteration == 1
        @test mean(eros.overland_flow_erosion.variables.soil_erosion_rate) ≈
              0.00861079076689589
        @test mean(eros.rainfall_erosion.variables.soil_erosion_rate) ≈
              0.00016326203201620437
        @test mean(eros.soil_erosion.variables.soil_erosion_rate) ≈ 0.008774052798912092
    end

    # run the second timestep
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model (land part)" begin
        eros = model.land

        @test mean(eros.soil_erosion.variables.soil_erosion_rate) ≈ 0.07765800489746684
        @test mean(eros.soil_erosion.variables.clay_erosion_rate) ≈ 0.002287480354866626
        @test mean(eros.soil_erosion.variables.silt_erosion_rate) ≈ 0.0036164773521184896
        @test mean(eros.soil_erosion.variables.sand_erosion_rate) ≈ 0.026301393837924607
        @test mean(eros.soil_erosion.variables.lagg_erosion_rate) ≈ 0.022577957752547836
        @test mean(eros.soil_erosion.variables.sagg_erosion_rate) ≈ 0.022874695590802723
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

        @test mean(land.to_river.variables.sediment_rate) ≈ 0.07623845927222056
        @test sum(land.to_river.variables.clay_rate) ≈ 114.42704329506047
        @test sum(land.to_river.variables.sand_rate) ≈ 1289.411874381068
        @test mean(land.sediment_flux.variables.clay) ≈ 0.006578791733506439

        @test mean(river.hydrological_forcing.q_river) ≈ 0.6975180562953642
        @test river.hydrological_forcing.waterlevel_river[domain.river.network.order[end]] ≈
              0.006103649735450745
        @test mean(domain.river.parameters.flow_width) ≈ 22.628250814095523

        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              0.4458019733090582
        @test mean(river.potential_erosion.variables.bed) ≈ 307.3559271135137

        @test sum(river.sediment_flux.boundary_conditions.erosion_land_clay) ≈
              114.42704329506047
        @test sum(river.sediment_flux.boundary_conditions.erosion_land_sand) ≈
              1289.411874381068
        @test mean(river.sediment_flux.boundary_conditions.transport_capacity) ≈
              0.4458019733090582
        @test mean(river.sediment_flux.variables.sediment_rate) ≈ 0.43330810531561187
        @test mean(river.sediment_flux.variables.erosion) ≈ 0.01894491899394169
        @test mean(river.sediment_flux.variables.deposition) ≈ 0.6939690338938546
        @test river.sediment_flux.variables.clay_rate[5649] ≈ 2.840979764480952e-9

        @test mean(river.concentrations.variables.suspended) ≈ 0.8261052189488919
    end

    Wflow.close_files(model)
end

@testitem "Run sediment with warm state" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_config.toml")
    config = Wflow.Config(tomlpath)
    config.model.cold_start__flag = false
    model = Wflow.Model(config)
    @testset "initial warm states" begin
        river = model.routing.river_flow
        @test mean(river.sediment_flux.variables.store_clay) ≈ 0.06881654169660432
        @test mean(river.sediment_flux.variables.leftover_clay) ≈ 1.8985203566870876e-7
        @test mean(river.sediment_flux.variables.clay_rate) ≈ 0.761820269217149

        @test mean(river.sediment_flux.variables.store_silt) ≈ 0.13222698520947598
        @test mean(river.sediment_flux.variables.leftover_silt) ≈ 3.0309355418150914e-9
        @test mean(river.sediment_flux.variables.silt_rate) ≈ 0.4054607471933968

        @test mean(river.sediment_flux.variables.store_sand) ≈ 0.6762573229681987
        @test mean(river.sediment_flux.variables.leftover_sand) ≈ 7.890080963256109e-10
        @test mean(river.sediment_flux.variables.sand_rate) ≈ 0.005085932599331321
    end

    Wflow.close_files(model)
end

@testitem "Run sediment other configuration" begin
    using Statistics: mean
    ### Test the sediment model with a different configuration file ###
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)

    model = Wflow.Model(config)
    (; domain) = model

    Wflow.run_timestep!(model)

    @testset "first timestep sediment model eurosem (land part)" begin
        eros = model.land

        @test eros.atmospheric_forcing.precipitation[1] ≈ 4.086122035980225
        @test eros.hydrological_forcing.interception[1] ≈ 0.6329902410507202
        @test eros.hydrological_forcing.q_land[1] ≈ 0.0
        @test eros.rainfall_erosion.parameters.soil_detachability[1] ≈ 2.0
        @test eros.rainfall_erosion.parameters.eurosem_exponent[1] ≈ 2.0
        @test eros.overland_flow_erosion.parameters.usle_c[1] ≈ 0.014194443821907043
        @test eros.overland_flow_erosion.variables.soil_erosion_rate[1] ≈ 0.0
        @test eros.rainfall_erosion.variables.soil_erosion_rate[1] ≈ 0.01232301374083337
        @test mean(eros.overland_flow_erosion.variables.soil_erosion_rate) ≈
              0.00861079076689589
        @test mean(eros.rainfall_erosion.variables.soil_erosion_rate) ≈
              0.0014726364432116048
        @test mean(eros.soil_erosion.variables.soil_erosion_rate) ≈ 0.010083427210107495
    end

    # run the second timestep
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model engelund (routing)" begin
        land = model.routing.overland_flow
        river = model.routing.river_flow

        @test river.transport_capacity.parameters.d50[1] == 0.05000000074505806
        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              0.1418728167951134

        @test mean(river.concentrations.variables.suspended) ≈ 0.24791810261189964
    end

    Wflow.close_files(model)
end

@testitem "Run sediment land only" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    # Update config to run only the land model
    config.model.run_river_model__flag = false
    # Use govers equation for land transport capacity
    config.model.land_transport = "govers"

    model = Wflow.Model(config)

    # run the first and second timestep
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model govers" begin
        eros = model.land
        land = model.routing.overland_flow

        @test mean(eros.soil_erosion.variables.soil_erosion_rate) ≈ 0.0776983847440198
        @test mean(land.transport_capacity.parameters.c_govers) ≈ 0.16393911236592437
        @test mean(land.transport_capacity.variables.sediment_transport_capacity) ≈
              1.0990864706347766e6
        @test mean(land.to_river.variables.sediment_rate) ≈ 0.07708434959918917
    end

    Wflow.close_files(model)
end

@testitem "Run sediment yalin transport capacity" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    # Update config to run only the land model
    config.model.run_river_model__flag = false
    # Use yalin equation for land transport capacity
    config.model.land_transport = "yalin"

    model = Wflow.Model(config)

    # run the first and second timestep
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model yalin" begin
        eros = model.land
        land = model.routing.overland_flow

        @test mean(eros.soil_erosion.variables.soil_erosion_rate) ≈ 0.0776983847440198
        @test mean(land.transport_capacity.parameters.d50) ≈ 0.001534350291334408
        @test mean(land.transport_capacity.variables.sediment_transport_capacity) ≈
              1.0988158364353527f6
        @test mean(land.to_river.variables.sediment_rate) ≈ 0.07759383356462951
    end

    Wflow.close_files(model)
end

### Test all river transport capacity ###
@testitem "Run sediment yang transport capacity" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    # Use yang equation for river transport capacity
    config.model.river_transport = "yang"

    model = Wflow.Model(config)

    # run the first and second timestep
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model yang (routing)" begin
        river = model.routing.river_flow

        @test river.transport_capacity.parameters.d50[1] == 0.05000000074505806
        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              39.955369987101946
        @test mean(river.concentrations.variables.suspended) ≈ 0.004036949009419181
    end

    Wflow.close_files(model)
end

@testitem "Run sediment kodatie transport capacity" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    # Use kodatie equation for river transport capacity
    config.model.river_transport = "kodatie"

    model = Wflow.Model(config)

    # run the first and second timestep
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model kodatie (routing)" begin
        river = model.routing.river_flow

        @test river.transport_capacity.parameters.a_kodatie[1] == 2829.6
        @test river.transport_capacity.parameters.b_kodatie[1] == 3.646
        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              30.332588671299625

        @test mean(river.concentrations.variables.suspended) ≈ 54.7559835316139
    end

    Wflow.close_files(model)
end

@testitem "Run sediment molinas transport capacity" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    # Use molinas equation for river transport capacity
    config.model.river_transport = "molinas"

    model = Wflow.Model(config)

    # run the first and second timestep
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model molinas (routing)" begin
        river = model.routing.river_flow

        @test river.transport_capacity.parameters.d50[1] == 0.05000000074505806
        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              350.7538564169241

        @test mean(river.concentrations.variables.suspended) ≈ 884.749798262824
    end

    Wflow.close_files(model)
end

@testitem "run wflow sediment" begin
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    Wflow.run(tomlpath; silent = true)
end
