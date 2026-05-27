@testitem "Run sediment" begin
    using Wflow: to_SI, Unit, TON_PER_DT, MM_PER_DT
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
        @test mean(variables.lagg_erosion_rate) ≈ 0.00026131895546930363
        @test mean(variables.sagg_erosion_rate) ≈ 0.0002647534211898463
    end

    @testset "second timestep sediment model (routing)" begin
        (; overland_flow, river_flow) = model.routing

        @test overland_flow.transport_capacity.parameters.dm_sand[1] ==
              0.00019999999999999998
        @test overland_flow.transport_capacity.parameters.dm_lagg[1] == 0.0005

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

        @test rainfall_erosion.parameters.eurosem_exponent[1] ≈ 2.0
        @test rainfall_erosion.variables.soil_erosion_rate[1] ≈ 0.00014262747385223806
        @test mean(rainfall_erosion.variables.soil_erosion_rate) ≈ 1.7044403277912093e-5

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
                :canopygapfraction => 0.1,
                :canopyheight => 10.246395046934293,
                :soil_detachability => 0.0019100919738700582,
                :soilcover_fraction => 0.013012199593097453,
                :eurosem_exponent => 2.0,
            ),
        )
        @test test_means(
            rainfall_erosion.variables,
            Dict(:soil_erosion_rate => 1.7044403277912093e-5),
        )
    end

    @testset "First timestep: soil erosion" begin
        (; soil_erosion) = model.land

        @test mean(soil_erosion.variables.soil_erosion_rate) ≈ 0.00011670633345031822

        @test test_means(
            soil_erosion.boundary_conditions,
            Dict(
                :overland_flow_erosion => 9.966193017240612e-5,
                :rainfall_erosion => 1.7044403277912093e-5,
            ),
        )
        @test test_means(
            soil_erosion.parameters,
            Dict(
                :clay_fraction => 0.043157121493024087,
                :sand_fraction => 0.19049260624776837,
                :sagg_fraction => 0.42007166542989277,
                :silt_fraction => 0.05922145470619702,
                :lagg_fraction => 0.28705715197829496,
            ),
        )
        @test test_means(
            soil_erosion.variables,
            Dict(
                :lagg_erosion_rate => 3.401781185687422e-5,
                :sand_erosion_rate => 3.913126488559666e-5,
                :soil_erosion_rate => 0.00011670633345031824,
                :silt_erosion_rate => 5.509212657233483e-6,
                :clay_erosion_rate => 3.4635082509624603e-6,
                :sagg_erosion_rate => 3.458453571269738e-5,
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
                :dm_clay => 2.0e-6,
                :dm_sagg => 2.9999999999999997e-5,
                :dm_lagg => 0.0005,
                :dm_silt => 9.999999999999999e-6,
                :dm_sand => 0.00019999999999999998,
            ),
        )
        @test test_means(
            overland_flow.transport_capacity.variables,
            Dict(
                :silt => 12713.808955655715,
                :clay => 12714.316406971457,
                :sagg => 12713.724398208968,
                :lagg => 12713.68470719566,
                :sand => 12713.688492720412,
                :sediment_transport_capacity => 63569.222960752224,
            ),
        )

        # Sediment flux
        @test test_means(
            overland_flow.sediment_flux.boundary_conditions,
            Dict(
                :erosion_clay => 3.46350825096246e-6,
                :transport_capacity_clay => 12714.316406971457,
                :erosion_silt => 5.509212657233483e-6,
                :erosion_lagg => 3.4017811856874225e-5,
                :transport_capacity_silt => 12713.808955655715,
                :erosion_sagg => 3.458453571269738e-5,
                :transport_capacity_sand => 12713.688492720412,
                :transport_capacity_sagg => 12713.724398208968,
                :transport_capacity_lagg => 12713.68470719566,
                :erosion_sand => 3.913126488559666e-5,
            ),
        )
        @test test_means(
            overland_flow.sediment_flux.variables,
            Dict(
                :silt => 1.3929340336677336e-5,
                :sediment_rate => 0.0003063635802931295,
                :deposition_silt => 5.509212657233483e-6,
                :deposition_lagg => 3.4017811856874225e-5,
                :deposition_sand => 3.913126488559666e-5,
                :sagg => 8.701802333723823e-5,
                :lagg => 8.779632070055166e-5,
                :clay => 8.771741406134715e-6,
                :deposition => 0.0001167063333633642,
                :deposition_clay => 3.46350825096246e-6,
                :deposition_sagg => 3.458453571269738e-5,
                :sand => 0.00010884815451252768,
            ),
        )

        # To river
        @test test_means(
            overland_flow.to_river.boundary_conditions,
            Dict(
                :deposition_sand => 3.913126488559666e-5,
                :deposition_sagg => 3.458453571269738e-5,
                :deposition_silt => 5.509212657233483e-6,
                :deposition_lagg => 3.4017811856874225e-5,
                :deposition_clay => 3.46350825096246e-6,
            ),
        )
        @test test_means(
            overland_flow.to_river.variables,
            Dict(
                :sagg_rate => 2.8893430878564412e-5,
                :lagg_rate => 2.910125752615814e-5,
                :sand_rate => 3.491356425429223e-5,
                :sediment_rate => 0.00010045351037913171,
                :clay_rate => 2.895528492414304e-6,
                :silt_rate => 4.649729227702621e-6,
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
                :erosion_land_sagg => 0.000255825655895618,
                :erosion_land_lagg => 0.0002576657761865142,
                :erosion_land_sand => 0.000309128587482301,
                :potential_erosion_river_bank => 0.021936887334262873,
                :transport_capacity => 0.00028865756033258057,
                :erosion_land_silt => 4.116922058904866e-5,
                :erosion_land_clay => 2.5637331850607295e-5,
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
                :dm_gravel => 0.002,
                :dm_sand => 0.00019999999999999998,
                :dm_sagg => 2.9999999999999997e-5,
                :dm_lagg => 0.0005,
                :dm_clay => 2.0e-6,
                :dm_silt => 9.999999999999999e-6,
            ),
        )
        @test test_means(
            river_flow.sediment_flux.variables,
            Dict(
                :deposition => 0.0009204350874843538,
                :erosion => 3.1008560815713806e-5,
                :leftover_gravel => 0.0,
                :leftover_clay => 3.916043968626082e-6,
                :leftover_lagg => 1.9500782501449444e-12,
                :leftover_sagg => 3.741274625469958e-12,
                :leftover_sand => 8.753169499245311e-13,
                :leftover_silt => 4.643328043695461e-13,
                :sagg_rate => 1.3279650822056091e-5,
                :sand_rate => 3.2815978684007007e-6,
                :sediment_rate => 0.00027668604048794676,
                :silt_rate => 0.00011697422451263947,
                :lagg_rate => 1.2393223742490294e-6,
                :clay_rate => 0.00014191124491060145,
                :gravel_rate => 0.0,
                :store_clay => 2.6169325235485176,
                :store_gravel => 0.13395698621708418,
                :store_lagg => 22.262323062512886,
                :store_sagg => 22.103336669377652,
                :store_sand => 27.135750167760686,
                :store_silt => 5.27329214923135,
            ),
        )

        # Concentrations
        @test test_means(
            river_flow.concentrations.boundary_conditions,
            Dict(
                :waterlevel => 0.05777413367369656,
                :q => 0.1528727680637441,
                :silt => 0.00011697422451263947,
                :clay => 0.00014191124491060145,
                :sagg => 1.3279650822056091e-5,
                :sand => 3.2815978684007007e-6,
                :lagg => 1.2393223742490294e-6,
                :gravel => 0.0,
            ),
        )
        @test test_means(
            river_flow.concentrations.parameters,
            Dict(
                :dm_clay => 2.0e-6,
                :dm_gravel => 0.002,
                :dm_sagg => 2.9999999999999997e-5,
                :dm_lagg => 0.0005,
                :dm_silt => 9.999999999999999e-6,
                :dm_sand => 0.00019999999999999998,
            ),
        )
        @test test_means(
            river_flow.concentrations.variables,
            Dict(
                :suspended => 0.00014755963695059,
                :bed => 1.145940858844025e-7,
                :total => 0.00014767423103647439,
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

        @test mean(river_flow.concentrations.variables.suspended) ≈ 0.00024791810261189963
    end

    Wflow.close_files(model)
end

@testitem "Run sediment land only" begin
    using Statistics: mean
    using Wflow: to_SI, TON_PER_DT
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

        @test mean(land.soil_erosion.variables.soil_erosion_rate) ≈ 0.0008992868604631915
        @test mean(overland_flow.transport_capacity.parameters.c_govers) ≈
              0.16393911236592437
        @test mean(overland_flow.transport_capacity.variables.sediment_transport_capacity) ≈
              12720.908224939538
        @test mean(overland_flow.to_river.variables.sediment_rate) ≈ 0.0008921799722128373
    end

    Wflow.close_files(model)
end

@testitem "Run sediment yalin transport capacity" begin
    using Statistics: mean
    using Wflow: to_SI, TON_PER_DT, MM
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

        @test mean(land.soil_erosion.variables.soil_erosion_rate) ≈ 0.0008992868604631915
        @test mean(overland_flow.transport_capacity.parameters.d50) ≈ 1.534350291334408e-6
        @test mean(overland_flow.transport_capacity.variables.sediment_transport_capacity) ≈
              Float32(12717.775884668434)
        @test mean(overland_flow.to_river.variables.sediment_rate) ≈ 0.0008980767762197222
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
        @test mean(river.concentrations.variables.suspended) ≈ 0.000004036949009419181
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

        @test mean(river.concentrations.variables.suspended) ≈ 0.0547559835316139
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

        @test mean(river.concentrations.variables.suspended) ≈ 0.8847497982628241
    end

    Wflow.close_files(model)
end

@testitem "run wflow sediment" begin
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    Wflow.run(tomlpath; silent = true)
end
