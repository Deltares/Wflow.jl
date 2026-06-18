@testitem "Run sediment" begin
    using Statistics: mean
    dt = 86400.0
    tomlpath = joinpath(@__DIR__, "sediment_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    model = Wflow.Model(config)
    (; domain) = model

    @testset "initial default states" begin
        (; river_flow) = model.routing
        @test mean(river_flow.sediment_flux.variables.store_clay) == 0.0
        @test mean(river_flow.sediment_flux.variables.leftover_clay) == 0.0
        @test mean(river_flow.sediment_flux.variables.clay_rate) == 0.0

        @test mean(river_flow.sediment_flux.variables.store_silt) == 0.0
        @test mean(river_flow.sediment_flux.variables.leftover_silt) == 0.0
        @test mean(river_flow.sediment_flux.variables.silt_rate) == 0.0

        @test mean(river_flow.sediment_flux.variables.store_sand) == 0.0
        @test mean(river_flow.sediment_flux.variables.leftover_sand) == 0.0
        @test mean(river_flow.sediment_flux.variables.sand_rate) == 0.0
    end

    Wflow.run_timestep!(model)

    @testset "first timestep sediment model (land part)" begin
        (; land) = model

        @test land.atmospheric_forcing.precipitation[1] ≈ 4.7293079120141486e-8
        @test land.hydrological_forcing.q_land[1] ≈ 0.0
        @test land.overland_flow_erosion.parameters.usle_k[1] ≈ 0.026510488241910934
        @test land.overland_flow_erosion.parameters.usle_c[1] ≈ 0.014194443821907043
        @test land.overland_flow_erosion.parameters.answers_overland_flow_factor[1] ≈
              0.8999999761581421
        @test land.overland_flow_erosion.variables.soil_erosion_rate[1] ≈ 0.0
        @test land.rainfall_erosion.variables.soil_erosion_rate[1] ≈ 3.1534233707052943e-6
        @test model.clock.iteration == 1
        @test mean(land.overland_flow_erosion.variables.soil_erosion_rate) ≈
              9.966193017240613e-5
        @test mean(land.rainfall_erosion.variables.soil_erosion_rate) ≈
              1.8896068520394022e-6
        @test mean(land.soil_erosion.variables.soil_erosion_rate) ≈ 0.0001015515370244455
    end

    # run the second timestep
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model (land part)" begin
        (; variables) = model.land.soil_erosion

        @test mean(variables.soil_erosion_rate) ≈ 0.0008988195011280884
        @test mean(variables.clay_erosion_rate) ≈ 2.6475467070215576e-5
        @test mean(variables.silt_erosion_rate) ≈ 4.185737676063066e-5
        @test mean(variables.sand_erosion_rate) ≈ 0.00030441428053153477
        @test mean(variables.large_aggregates_erosion_rate) ≈ 0.00026131895546930363
        @test mean(variables.small_aggregates_erosion_rate) ≈ 0.0002647534211898463
    end

    @testset "second timestep sediment model (routing)" begin
        (; overland_flow, river_flow) = model.routing

        @test overland_flow.transport_capacity.parameters.median_diameter_sand[1] ==
              0.00019999999999999998
        @test overland_flow.transport_capacity.parameters.median_diameter_large_aggregates[1] ==
              0.0005

        @test mean(overland_flow.transport_capacity.boundary_conditions.q) ≈
              0.006879398771052133
        @test mean(overland_flow.transport_capacity.variables.silt) ≈
              Float32(12717.775884668434)
        @test mean(overland_flow.transport_capacity.variables.sand) ≈
              Float32(12716.54007278791)
        @test mean(overland_flow.transport_capacity.variables.clay) ≈
              Float32(12722.980552102701)

        @test mean(overland_flow.to_river.variables.sediment_rate) ≈ 0.0008823895749099599
        @test sum(overland_flow.to_river.variables.clay_rate) ≈ 1.324387075174311
        @test sum(overland_flow.to_river.variables.sand_rate) ≈ 14.923748546077173
        @test mean(overland_flow.sediment_flux.variables.clay) ≈ 7.614342284150971e-5

        @test mean(river_flow.hydrological_forcing.q_river) ≈ 0.6975180562953642
        @test river_flow.hydrological_forcing.waterlevel_river[domain.river.network.order[end]] ≈
              0.006103649735450745
        @test mean(domain.river.parameters.flow_width) ≈ 22.628250814095523

        @test mean(river_flow.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river_flow.transport_capacity.variables.sediment_transport_capacity) ≈
              0.0051597450614474325
        @test mean(river_flow.potential_erosion.variables.bed) ≈ 3.5573602675175193

        @test sum(river_flow.sediment_flux.boundary_conditions.erosion_land_clay) ≈
              1.324387075174311
        @test sum(river_flow.sediment_flux.boundary_conditions.erosion_land_sand) ≈
              14.9237485460772
        @test mean(river_flow.sediment_flux.boundary_conditions.transport_capacity) ≈
              0.0051597450614474325
        @test mean(river_flow.sediment_flux.variables.sediment_rate) ≈ 0.005015140116935809
        @test mean(river_flow.sediment_flux.variables.erosion) ≈ 0.00021926989781640333
        @test mean(river_flow.sediment_flux.variables.deposition) ≈ 0.008032049003401103
        @test river_flow.sediment_flux.variables.clay_rate[5649] ≈ 3.2881710237048055e-11

        @test mean(river_flow.concentrations.variables.suspended) ≈ 0.000826105218948891
    end

    Wflow.close_files(model)
end

@testitem "Run sediment with warm state" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.cold_start__flag = false
    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)
    @testset "initial warm states" begin
        river = model.routing.river_flow
        @test mean(river.sediment_flux.variables.store_clay) ≈ 68.81654169660432
        @test mean(river.sediment_flux.variables.leftover_clay) ≈ 0.00018985203566870877
        @test mean(river.sediment_flux.variables.clay_rate) ≈ 0.008817364227050335

        @test mean(river.sediment_flux.variables.store_silt) ≈ 132.22698520947598
        @test mean(river.sediment_flux.variables.leftover_silt) ≈ 3.0309355418150915e-6
        @test mean(river.sediment_flux.variables.silt_rate) ≈ 0.004692832722145796

        @test mean(river.sediment_flux.variables.store_sand) ≈ 676.2573229681987
        @test mean(river.sediment_flux.variables.leftover_sand) ≈ 7.890080963256109e-7
        @test mean(river.sediment_flux.variables.sand_rate) ≈ 5.8864960640408804e-5
    end

    Wflow.close_files(model)
end

@testitem "Run sediment other configuration" begin
    using Statistics: mean
    include("testing_utils.jl")
    ### Test the sediment model with a different configuration file ###
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    model = Wflow.Model(config)
    (; domain, clock) = model
    dt = Wflow.tosecond(clock.dt)

    Wflow.run_timestep!(model)

    @testset "First timestep: forcing" begin
        (; atmospheric_forcing, hydrological_forcing) = model.land

        @test atmospheric_forcing.precipitation[1] ≈ 4.7293079120141486e-8
        @test hydrological_forcing.interception[1] ≈ 7.326275938087039e-9
        @test hydrological_forcing.q_land[1] ≈ 0.0

        @test test_means(atmospheric_forcing, Dict(:precipitation => 2.332765022942163e-8))
        @test test_means(
            hydrological_forcing,
            Dict(
                :interception => 5.518207945452663e-9,
                :waterlevel_land => 2.669102444877212e-5,
                :q_land => 0.0008427802189709429,
                :waterlevel_river => 0.0,
                :q_river => 0.0,
            ),
        )
    end

    @testset "First timestep: rainfall erosion" begin
        (; rainfall_erosion) = model.land

        @test rainfall_erosion.variables.soil_erosion_rate[1] ≈ 0.1426274738522381
        @test mean(rainfall_erosion.variables.soil_erosion_rate) ≈ 0.01659596342237728

        @test test_means(
            rainfall_erosion.boundary_conditions,
            Dict(
                :waterlevel => 2.669102444877212e-5,
                :interception => 5.518207945452663e-9,
                :precipitation => 2.332765022942163e-8,
            ),
        )
        @test test_means(
            rainfall_erosion.parameters,
            Dict(
                :canopy_gap_fraction => 0.1,
                :canopy_height => 10.246395046934293,
                :soil_detachability => 0.0019100919738700582,
                :soilcover_fraction => 0.013012199593097453,
                :eurosem_exponent => 2000.0,
            ),
        )
        @test test_means(
            rainfall_erosion.variables,
            Dict(:soil_erosion_rate => 0.01659596342237728),
        )
    end

    @testset "First timestep: soil erosion" begin
        (; soil_erosion) = model.land

        @test test_means(
            soil_erosion.boundary_conditions,
            Dict(
                :overland_flow_erosion => 9.966193017240612e-5,
                :rainfall_erosion => 0.01659596342237728,
            ),
        )
        @test test_means(
            soil_erosion.parameters,
            Dict(
                :clay_fraction => 0.043157121493024087,
                :sand_fraction => 0.19049260624776837,
                :small_aggregates_fraction => 0.42007166542989277,
                :silt_fraction => 0.05922145470619702,
                :large_aggregates_fraction => 0.28705715197829496,
            ),
        )
        @test test_means(
            soil_erosion.variables,
            Dict(
                :large_aggregates_erosion_rate => 0.004917999180253473,
                :sand_erosion_rate => 0.00427611093783281,
                :soil_erosion_rate => 0.01669562535254969,
                :silt_erosion_rate => 0.0009091033581761446,
                :clay_erosion_rate => 0.0006038299771647817,
                :small_aggregates_erosion_rate => 0.005988581897554994,
            ),
        )
    end

    @testset "First timestep: overland flow erosion" begin
        (; overland_flow_erosion) = model.land

        @test overland_flow_erosion.parameters.usle_c[1] ≈ 0.014194443821907043
        @test overland_flow_erosion.variables.soil_erosion_rate[1] ≈ 0.0

        @test test_means(
            overland_flow_erosion.boundary_conditions,
            Dict(:q => 0.0008427802189709429),
        )
        @test test_means(
            overland_flow_erosion.parameters,
            Dict(
                :answers_overland_flow_factor => 0.8999999761581421,
                :usle_k => 0.010302870616816031,
                :usle_c => 0.1264203163327229,
            ),
        )
        @test test_means(
            overland_flow_erosion.variables,
            Dict(:soil_erosion_rate => 9.966193017240612e-5),
        )
    end

    @testset "First timestep: overland flow" begin
        (; overland_flow) = model.routing

        # Forcing
        @test test_means(
            overland_flow.hydrological_forcing,
            Dict(
                :interception => 5.518207945452663e-9,
                :waterlevel_land => 2.669102444877212e-5,
                :q_land => 0.0008427802189709429,
                :waterlevel_river => 0.0,
                :q_river => 0.0,
            ),
        )

        # Transport capacity
        @test test_means(
            overland_flow.transport_capacity.boundary_conditions,
            Dict(:waterlevel => 2.669102444877212e-5, :q => 0.0008427802189709429),
        )
        @test test_means(
            overland_flow.transport_capacity.parameters,
            Dict(
                :density => 2650.0,
                :median_diameter_clay => 2.0e-6,
                :median_diameter_small_aggregates => 2.9999999999999997e-5,
                :median_diameter_large_aggregates => 0.0005,
                :median_diameter_silt => 9.999999999999999e-6,
                :median_diameter_sand => 0.00019999999999999998,
            ),
        )
        @test test_means(
            overland_flow.transport_capacity.variables,
            Dict(
                :silt => 12713.808955655715,
                :clay => 12714.316406971457,
                :small_aggregates => 12713.724398208968,
                :large_aggregates => 12713.68470719566,
                :sand => 12713.688492720412,
                :sediment_transport_capacity => 63569.222960752224,
            ),
        )

        # Sediment flux
        @test test_means(
            overland_flow.sediment_flux.boundary_conditions,
            Dict(
                :erosion_clay => 0.0006038299771647817,
                :transport_capacity_clay => 12714.316406971457,
                :erosion_silt => 0.0009091033581761446,
                :erosion_large_aggregates => 0.004917999180253473,
                :transport_capacity_silt => 12713.808955655715,
                :erosion_small_aggregates => 0.005988581897554994,
                :transport_capacity_sand => 12713.688492720412,
                :transport_capacity_small_aggregates => 12713.724398208968,
                :transport_capacity_large_aggregates => 12713.68470719566,
                :erosion_sand => 0.00427611093783281,
            ),
        )
        @test test_means(
            overland_flow.sediment_flux.variables,
            Dict(
                :silt => 2.4005496245140833e-5,
                :sediment_rate => 0.00040047351477098884,
                :deposition_silt => 0.0009091033581761447,
                :deposition_large_aggregates => 0.004917999180253473,
                :deposition_sand => 0.00427611093783281,
                :small_aggregates => 0.00012229895178867518,
                :large_aggregates => 0.00010422878536180842,
                :clay => 1.8798050245173526e-5,
                :deposition => 0.016695625350982204,
                :deposition_clay => 0.0006038299771647818,
                :deposition_small_aggregates => 0.0059885818975549945,
                :sand => 0.00013114223113019097,
            ),
        )

        # To river
        @test test_means(
            overland_flow.to_river.boundary_conditions,
            Dict(
                :deposition_sand => 0.00427611093783281,
                :deposition_small_aggregates => 0.0059885818975549945,
                :deposition_silt => 0.0009091033581761447,
                :deposition_large_aggregates => 0.004917999180253473,
                :deposition_clay => 0.0006038299771647818,
            ),
        )
        @test test_means(
            overland_flow.to_river.variables,
            Dict(
                :small_aggregates_rate => 0.000649151785940205,
                :large_aggregates_rate => 0.0005168185844448506,
                :sand_rate => 0.00042180961215772613,
                :sediment_rate => 0.0017542224635858756,
                :clay_rate => 6.688821117271656e-5,
                :silt_rate => 9.955426987037696e-5,
            ),
        )
    end

    @testset "`First timestep: river erosion" begin
        (; river_flow) = model.routing

        # Forcing
        @test test_means(
            river_flow.hydrological_forcing,
            Dict(
                :interception => 0.0,
                :waterlevel_land => 0.0,
                :q_land => 0.0,
                :waterlevel_river => 0.05777413367369656,
                :q_river => 0.1528727680637441,
            ),
        )

        # Transport capacity
        @test test_means(
            river_flow.transport_capacity.boundary_conditions,
            Dict(:waterlevel => 0.05777413367369656, :q => 0.1528727680637441),
        )
        @test test_means(
            river_flow.transport_capacity.parameters,
            Dict(:density => 2650.0, :d50 => 0.00011542705629386542),
        )
        @test test_means(
            river_flow.transport_capacity.variables,
            Dict(:sediment_transport_capacity => 0.00028865756033258057),
        )

        # Potential erosion
        @test test_means(
            river_flow.potential_erosion.boundary_conditions,
            Dict(:waterlevel => 0.05777413367369656),
        )
        @test test_means(
            river_flow.potential_erosion.parameters,
            Dict(:d50 => 0.00011542705629386542),
        )
        @test test_means(
            river_flow.potential_erosion.variables,
            Dict(:bed => 1.613938333854238, :bank => 0.021936887334262873),
        )

        # Sediment flux
        @test test_means(
            river_flow.sediment_flux.boundary_conditions,
            Dict(
                :q => 0.1528727680637441,
                :waterlevel => 0.05777413367369656,
                :erosion_land_small_aggregates => 0.005747662232011684,
                :erosion_land_large_aggregates => 0.004575969323280933,
                :erosion_land_sand => 0.00373474929809679,
                :potential_erosion_river_bank => 0.021936887334262873,
                :transport_capacity => 0.00028865756033258057,
                :erosion_land_silt => 0.0008814645963589342,
                :erosion_land_clay => 0.0005922356734602861,
                :potential_erosion_river_bed => 1.613938333854238,
            ),
        )
        @test test_means(
            river_flow.sediment_flux.parameters,
            Dict(
                :reservoir_area => 1865.8546640141467,
                :sand_fraction => 0.2251547318080376,
                :gravel_fraction => 0.05000000074505806,
                :clay_fraction => 0.17458002217884722,
                :silt_fraction => 0.5502652340921862,
                :reservoir_trapping_efficiency => 1.0,
                :median_diameter_gravel => 0.002,
                :median_diameter_sand => 0.00019999999999999998,
                :median_diameter_small_aggregates => 2.9999999999999997e-5,
                :median_diameter_large_aggregates => 0.0005,
                :median_diameter_clay => 2.0e-6,
                :median_diameter_silt => 9.999999999999999e-6,
            ),
        )
        @test test_means(
            river_flow.sediment_flux.variables,
            Dict(
                :deposition => 0.015539610938657309,
                :erosion => 7.529860784132161e-6,
                :leftover_gravel => 0.0,
                :leftover_clay => 3.916043968626082e-6,
                :leftover_large_aggregates => 1.9500782501449444e-12,
                :leftover_small_aggregates => 3.741274625469958e-12,
                :leftover_sand => 8.753169499245311e-13,
                :leftover_silt => 4.643328043695461e-13,
                :small_aggregates_rate => 1.2796702580071101e-5,
                :sand_rate => 1.0713110703154987e-6,
                :sediment_rate => 0.0002812198523148606,
                :silt_rate => 7.567864551321369e-5,
                :large_aggregates_rate => 1.377556978587296e-7,
                :clay_rate => 0.0001915354374534016,
                :gravel_rate => 0.0,
                :store_clay => 51.26674527072368,
                :store_gravel => 0.03252899943571069,
                :store_large_aggregates => 395.36374953147066,
                :store_small_aggregates => 496.5980168458057,
                :store_sand => 322.77996911494876,
                :store_silt => 76.5813753376071,
            ),
        )

        # Concentrations
        @test test_means(
            river_flow.concentrations.boundary_conditions,
            Dict(
                :waterlevel => 0.05777413367369656,
                :q => 0.1528727680637441,
                :silt => 7.567864551321369e-5,
                :clay => 0.0001915354374534016,
                :small_aggregates => 1.2796702580071101e-5,
                :sand => 1.0713110703154987e-6,
                :large_aggregates => 1.377556978587296e-7,
                :gravel => 0.0,
            ),
        )
        @test test_means(
            river_flow.concentrations.parameters,
            Dict(
                :median_diameter_clay => 2.0e-6,
                :median_diameter_gravel => 0.002,
                :median_diameter_small_aggregates => 2.9999999999999997e-5,
                :median_diameter_large_aggregates => 0.0005,
                :median_diameter_silt => 9.999999999999999e-6,
                :median_diameter_sand => 0.00019999999999999998,
            ),
        )
        @test test_means(
            river_flow.concentrations.variables,
            Dict(
                :suspended => 0.0001563141343612164,
                :bed => 1.1188707120088609e-13,
                :total => 0.00015631413447310348,
            ),
        )
    end

    # run the second timestep
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model engelund (routing)" begin
        (; river_flow) = model.routing

        @test river_flow.transport_capacity.parameters.d50[1] == 5.0000000745058064e-5
        @test mean(river_flow.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river_flow.transport_capacity.variables.sediment_transport_capacity) ≈
              0.001642046490684183

        @test mean(river_flow.concentrations.variables.suspended) ≈ 0.000276598816325157
    end

    Wflow.close_files(model)
end

@testitem "Run sediment land only" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    # Update config to run only the land model
    config.model.run_river_model__flag = false
    # Use govers equation for land transport capacity
    config.model.land_transport = "govers"

    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)

    # run the first and second timestep
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model govers" begin
        (; land, routing) = model
        (; overland_flow) = routing

        @test mean(land.soil_erosion.variables.soil_erosion_rate) ≈ 0.0013652659345724893
        @test mean(overland_flow.transport_capacity.parameters.c_govers) ≈
              0.16393911236592437
        @test mean(overland_flow.transport_capacity.variables.sediment_transport_capacity) ≈
              12720.908224939538
        @test mean(overland_flow.to_river.variables.sediment_rate) ≈ 0.0009306672286819593
    end

    Wflow.close_files(model)
end

@testitem "Run sediment yalin transport capacity" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    # Update config to run only the land model
    config.model.run_river_model__flag = false
    # Use yalin equation for land transport capacity
    config.model.land_transport = "yalin"

    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)

    # run the first and second timestep
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model yalin" begin
        # eros = model.land
        # land = model.routing.overland_flow
        (; land, routing) = model
        (; overland_flow) = routing

        @test mean(land.soil_erosion.variables.soil_erosion_rate) ≈ 0.0013652659345724893
        @test mean(overland_flow.transport_capacity.parameters.d50) ≈ 1.534350291334408e-6
        @test mean(overland_flow.transport_capacity.variables.sediment_transport_capacity) ≈
              Float32(12717.775884668434)
        @test mean(overland_flow.to_river.variables.sediment_rate) ≈ 0.0009377521191159138
    end

    Wflow.close_files(model)
end

### Test all river transport capacity ###
@testitem "Run sediment yang transport capacity" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    # Use yang equation for river transport capacity
    config.model.river_transport = "yang"

    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)

    # run the first and second timestep
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model yang (routing)" begin
        river = model.routing.river_flow

        @test river.transport_capacity.parameters.d50[1] == 5.0000000745058064e-5
        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              0.46244641188775404
        @test mean(river.concentrations.variables.suspended) ≈ 0.0037597697024228617
    end

    Wflow.close_files(model)
end

@testitem "Run sediment kodatie transport capacity" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    # Use kodatie equation for river transport capacity
    config.model.river_transport = "kodatie"

    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)

    # run the first and second timestep
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model kodatie (routing)" begin
        river = model.routing.river_flow

        @test river.transport_capacity.parameters.a_kodatie[1] == 2829.6
        @test river.transport_capacity.parameters.b_kodatie[1] == 3.646
        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              0.35107162814004217

        @test mean(river.concentrations.variables.suspended) ≈ 0.054748501397664226
    end

    Wflow.close_files(model)
end

@testitem "Run sediment molinas transport capacity" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    # Use molinas equation for river transport capacity
    config.model.river_transport = "molinas"

    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)

    # run the first and second timestep
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model molinas (routing)" begin
        river = model.routing.river_flow

        @test river.transport_capacity.parameters.d50[1] == 5.0000000745058064e-5
        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              4.059651115936618

        @test mean(river.concentrations.variables.suspended) ≈ 0.871274558152215
    end
    Wflow.close_files(model)
end

@testitem "run wflow sediment" begin
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    Wflow.run(tomlpath; silent = true)
end
