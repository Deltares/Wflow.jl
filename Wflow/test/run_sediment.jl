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

        @test eros.atmospheric_forcing.precipitation[1] ≈
              to_SI(4.086122035980225, MM_PER_DT; dt_val = dt)
        @test eros.hydrological_forcing.q_land[1] ≈ 0.0
        @test eros.overland_flow_erosion.parameters.usle_k[1] ≈ 0.026510488241910934
        @test eros.overland_flow_erosion.parameters.usle_c[1] ≈ 0.014194443821907043
        @test eros.overland_flow_erosion.parameters.answers_overland_flow_factor[1] ≈
              0.8999999761581421
        @test eros.overland_flow_erosion.variables.soil_erosion_rate[1] ≈ 0.0
        @test eros.rainfall_erosion.variables.soil_erosion_rate[1] ≈
              to_SI(0.00027245577922893746, TON_PER_DT; dt_val = dt)
        @test model.clock.iteration == 1
        @test mean(eros.overland_flow_erosion.variables.soil_erosion_rate) ≈
              to_SI(0.00861079076689589, TON_PER_DT; dt_val = dt)
        @test mean(eros.rainfall_erosion.variables.soil_erosion_rate) ≈
              to_SI(0.00016326203201620437, TON_PER_DT; dt_val = dt)
        @test mean(eros.soil_erosion.variables.soil_erosion_rate) ≈
              to_SI(0.008774052798912092, TON_PER_DT; dt_val = dt)
    end

    # run the second timestep
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model (land part)" begin
        eros = model.land

        @test mean(eros.soil_erosion.variables.soil_erosion_rate) ≈
              to_SI(0.07765800489746684, TON_PER_DT; dt_val = dt)
        @test mean(eros.soil_erosion.variables.clay_erosion_rate) ≈
              to_SI(0.002287480354866626, TON_PER_DT; dt_val = dt)
        @test mean(eros.soil_erosion.variables.silt_erosion_rate) ≈
              to_SI(0.003616477352118489, TON_PER_DT; dt_val = dt)
        @test mean(eros.soil_erosion.variables.sand_erosion_rate) ≈
              to_SI(0.026301393837924607, TON_PER_DT; dt_val = dt)
        @test mean(eros.soil_erosion.variables.lagg_erosion_rate) ≈
              to_SI(0.022577957752547836, TON_PER_DT; dt_val = dt)
        @test mean(eros.soil_erosion.variables.sagg_erosion_rate) ≈
              to_SI(0.022874695590802723, TON_PER_DT; dt_val = dt)
    end

    @testset "second timestep sediment model (routing)" begin
        land = model.routing.overland_flow
        river = model.routing.river_flow

        @test land.transport_capacity.parameters.dm_sand[1] == to_SI(200.0, Unit(; μm = 1))
        @test land.transport_capacity.parameters.dm_lagg[1] == to_SI(500.0, Unit(; μm = 1))

        @test mean(land.transport_capacity.boundary_conditions.q) ≈ 0.006879398771052133
        @test mean(land.transport_capacity.variables.silt) ≈
              Float32(to_SI(1.0988158364353527e6, TON_PER_DT; dt_val = dt))
        @test mean(land.transport_capacity.variables.sand) ≈
              Float32(to_SI(1.0987090622888755e6, TON_PER_DT; dt_val = dt))
        @test mean(land.transport_capacity.variables.clay) ≈
              Float32(to_SI(1.0992655197016734e6, TON_PER_DT; dt_val = dt))

        @test mean(land.to_river.variables.sediment_rate) ≈
              to_SI(0.0762386279230294, TON_PER_DT; dt_val = dt)
        @test sum(land.to_river.variables.clay_rate) ≈
              to_SI(114.42704329506047, TON_PER_DT; dt_val = dt)
        @test sum(land.to_river.variables.sand_rate) ≈
              to_SI(1289.4173249850346, TON_PER_DT; dt_val = dt)
        @test mean(land.sediment_flux.variables.clay) ≈
              to_SI(0.006578791733506439, TON_PER_DT; dt_val = dt)

        @test mean(river.hydrological_forcing.q_river) ≈ 0.6975180562953642
        @test river.hydrological_forcing.waterlevel_river[domain.river.network.order[end]] ≈
              0.006103649735450745
        @test mean(domain.river.parameters.flow_width) ≈ 22.628250814095523

        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              to_SI(0.4458019733090582, TON_PER_DT; dt_val = dt)
        @test mean(river.potential_erosion.variables.bed) ≈
              to_SI(307.3559271135137, TON_PER_DT; dt_val = dt)

        @test sum(river.sediment_flux.boundary_conditions.erosion_land_clay) ≈
              to_SI(114.42704329506047, TON_PER_DT; dt_val = dt)
        @test sum(river.sediment_flux.boundary_conditions.erosion_land_sand) ≈
              to_SI(1289.417324985034, TON_PER_DT; dt_val = dt)
        @test mean(river.sediment_flux.boundary_conditions.transport_capacity) ≈
              to_SI(0.4458019733090582, TON_PER_DT; dt_val = dt)
        @test mean(river.sediment_flux.variables.sediment_rate) ≈
              to_SI(0.43330815120929567, TON_PER_DT; dt_val = dt)
        @test mean(river.sediment_flux.variables.erosion) ≈
              to_SI(0.018944871787785745, TON_PER_DT; dt_val = dt)
        @test mean(river.sediment_flux.variables.deposition) ≈
              to_SI(0.6939704797633529, TON_PER_DT; dt_val = dt)
        @test river.sediment_flux.variables.clay_rate[5649] ≈
              to_SI(2.840979764480952e-9, TON_PER_DT; dt_val = dt)

        @test mean(river.concentrations.variables.suspended) ≈
              to_SI(0.8261173586131326, Unit(; g = 1, m = -3))
    end

    Wflow.close_files(model)
end

@testitem "Run sediment with warm state" begin
    using Statistics: mean
    using Wflow: to_SI, Unit, TON_PER_DT
    TON = Unit(; t = 1)
    tomlpath = joinpath(@__DIR__, "sediment_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.cold_start__flag = false
    model = Wflow.Model(config)
    dt = Wflow.tosecond(model.clock.dt)
    @testset "initial warm states" begin
        river = model.routing.river_flow
        @test mean(river.sediment_flux.variables.store_clay) ≈
              to_SI(0.06881654169660432, TON)
        @test mean(river.sediment_flux.variables.leftover_clay) ≈
              to_SI(1.8985203566870876e-7, TON)
        @test mean(river.sediment_flux.variables.clay_rate) ≈
              to_SI(0.761820269217149, TON_PER_DT; dt_val = dt)

        @test mean(river.sediment_flux.variables.store_silt) ≈
              to_SI(0.13222698520947598, TON)
        @test mean(river.sediment_flux.variables.leftover_silt) ≈
              to_SI(3.0309355418150914e-9, TON)
        @test mean(river.sediment_flux.variables.silt_rate) ≈
              to_SI(0.4054607471933968, TON_PER_DT; dt_val = dt)

        @test mean(river.sediment_flux.variables.store_sand) ≈
              to_SI(0.6762573229681987, TON)
        @test mean(river.sediment_flux.variables.leftover_sand) ≈
              to_SI(7.890080963256109e-10, TON)
        @test mean(river.sediment_flux.variables.sand_rate) ≈
              to_SI(0.005085932599331321, TON_PER_DT; dt_val = dt)
    end

    Wflow.close_files(model)
end

@testitem "Run sediment other configuration" begin
    using Statistics: mean
    using Wflow: to_SI, Unit, MM_PER_DT, TON_PER_DT, MM
    GRAM_PER_JOULE = Unit(; g = 1, J = -1)
    GRAM_PER_M3 = Unit(; g = 1, m = -3)
    μM = Unit(; μm = 1)
    TON = Unit(; t = 1)
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

        @test atmospheric_forcing.precipitation[1] ≈
              to_SI(4.086122035980225, MM_PER_DT; dt_val = dt)
        @test hydrological_forcing.interception[1] ≈
              to_SI(0.6329902410507202, MM_PER_DT; dt_val = dt)
        @test hydrological_forcing.q_land[1] ≈ 0.0

        @test test_means(
            atmospheric_forcing,
            Dict(:precipitation => to_SI(2.015508979822029, MM_PER_DT; dt_val = dt)),
        )
        @test test_means(
            hydrological_forcing,
            Dict(
                :interception => to_SI(0.4767731664871101, MM_PER_DT; dt_val = dt),
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
        @test rainfall_erosion.variables.soil_erosion_rate[1] ≈
              to_SI(0.01232301374083337, TON_PER_DT; dt_val = dt)
        @test mean(rainfall_erosion.variables.soil_erosion_rate) ≈
              to_SI(0.0014726364432116048, TON_PER_DT; dt_val = dt)

        @test test_means(
            rainfall_erosion.boundary_conditions,
            Dict(
                :waterlevel => 2.669102444877212e-5,
                :interception => to_SI(0.4767731664871101, MM_PER_DT; dt_val = dt),
                :precipitation => to_SI(2.015508979822029, MM_PER_DT; dt_val = dt),
            ),
        )
        @test test_means(
            rainfall_erosion.parameters,
            Dict(
                :canopygapfraction => 0.1,
                :canopyheight => 10.246395046934293,
                :soil_detachability => to_SI(1.9100919738700581, GRAM_PER_JOULE),
                :soilcover_fraction => 0.013012199593097453,
                :eurosem_exponent => 2.0,
            ),
        )
        @test test_means(
            rainfall_erosion.variables,
            Dict(
                :soil_erosion_rate => to_SI(0.0014726364432116048, TON_PER_DT; dt_val = dt),
            ),
        )
    end

    @testset "First timestep: soil erosion" begin
        (; soil_erosion) = model.land

        @test mean(soil_erosion.variables.soil_erosion_rate) ≈
              to_SI(0.010083427210107495, TON_PER_DT; dt_val = dt)

        @test test_means(
            soil_erosion.boundary_conditions,
            Dict(
                :overland_flow_erosion =>
                    to_SI(0.008610790766895889, TON_PER_DT; dt_val = dt),
                :rainfall_erosion => to_SI(0.0014726364432116048, TON_PER_DT; dt_val = dt),
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
                :lagg_erosion_rate => to_SI(0.0029391389444339325, TON_PER_DT; dt_val = dt),
                :sand_erosion_rate => to_SI(0.0033809412861155515, TON_PER_DT; dt_val = dt),
                :soil_erosion_rate => to_SI(0.010083427210107496, TON_PER_DT; dt_val = dt),
                :silt_erosion_rate =>
                    to_SI(0.00047599597358497295, TON_PER_DT; dt_val = dt),
                :clay_erosion_rate => to_SI(0.0002992471128831566, TON_PER_DT; dt_val = dt),
                :sagg_erosion_rate => to_SI(0.0029881038855770546, TON_PER_DT; dt_val = dt),
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
            Dict(
                :soil_erosion_rate => to_SI(0.008610790766895889, TON_PER_DT; dt_val = dt),
            ),
        )
    end

    @testset "First timestep: overland flow" begin
        (; overland_flow) = model.routing

        # Forcing
        @test test_means(
            overland_flow.hydrological_forcing,
            Dict(
                :interception => to_SI(0.4767731664871101, MM_PER_DT; dt_val = dt),
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
                :dm_clay => to_SI(2.0, μM),
                :dm_sagg => to_SI(30.0, μM),
                :dm_lagg => to_SI(500.0, μM),
                :dm_silt => to_SI(10.0, μM),
                :dm_sand => to_SI(200.0, μM),
            ),
        )
        @test test_means(
            overland_flow.transport_capacity.variables,
            Dict(
                :silt => to_SI(1.098473093768654e6, TON_PER_DT; dt_val = dt),
                :clay => to_SI(1.098516937562334e6, TON_PER_DT; dt_val = dt),
                :sagg => to_SI(1.098465788005255e6, TON_PER_DT; dt_val = dt),
                :lagg => to_SI(1.098462358701705e6, TON_PER_DT; dt_val = dt),
                :sand => to_SI(1.0984626857710436e6, TON_PER_DT; dt_val = dt),
                :sediment_transport_capacity =>
                    to_SI(5.492380863808992e6, TON_PER_DT; dt_val = dt),
            ),
        )

        # Sediment flux
        @test test_means(
            overland_flow.sediment_flux.boundary_conditions,
            Dict(
                :erosion_clay => to_SI(0.00029924711288315653, TON_PER_DT; dt_val = dt),
                :transport_capacity_clay =>
                    to_SI(1.098516937562334e6, TON_PER_DT; dt_val = dt),
                :erosion_silt => to_SI(0.00047599597358497295, TON_PER_DT; dt_val = dt),
                :erosion_lagg => to_SI(0.002939138944433933, TON_PER_DT; dt_val = dt),
                :transport_capacity_silt =>
                    to_SI(1.098473093768654e6, TON_PER_DT; dt_val = dt),
                :erosion_sagg => to_SI(0.002988103885577053, TON_PER_DT; dt_val = dt),
                :transport_capacity_sand =>
                    to_SI(1.0984626857710436e6, TON_PER_DT; dt_val = dt),
                :transport_capacity_sagg =>
                    to_SI(1.098465788005255e6, TON_PER_DT; dt_val = dt),
                :transport_capacity_lagg =>
                    to_SI(1.098462358701705e6, TON_PER_DT; dt_val = dt),
                :erosion_sand => to_SI(0.0033809412861155515, TON_PER_DT; dt_val = dt),
            ),
        )
        @test test_means(
            overland_flow.sediment_flux.variables,
            Dict(
                :silt => to_SI(0.0012034950050889219, TON_PER_DT; dt_val = dt),
                :sediment_rate => to_SI(0.02646981333732639, TON_PER_DT; dt_val = dt),
                :deposition_silt => to_SI(0.00047599597358497295, TON_PER_DT; dt_val = dt),
                :deposition_lagg => to_SI(0.002939138944433933, TON_PER_DT; dt_val = dt),
                :deposition_sand => to_SI(0.0033809412861155515, TON_PER_DT; dt_val = dt),
                :sagg => to_SI(0.007518357216337384, TON_PER_DT; dt_val = dt),
                :lagg => to_SI(0.007585602108527664, TON_PER_DT; dt_val = dt),
                :clay => to_SI(0.0007578784574900395, TON_PER_DT; dt_val = dt),
                :deposition => to_SI(0.010083427202594667, TON_PER_DT; dt_val = dt),
                :deposition_clay => to_SI(0.00029924711288315653, TON_PER_DT; dt_val = dt),
                :deposition_sagg => to_SI(0.002988103885577053, TON_PER_DT; dt_val = dt),
                :sand => to_SI(0.009404480549882391, TON_PER_DT; dt_val = dt),
            ),
        )

        # To river
        @test test_means(
            overland_flow.to_river.boundary_conditions,
            Dict(
                :deposition_sand => to_SI(0.0033809412861155515, TON_PER_DT; dt_val = dt),
                :deposition_sagg => to_SI(0.002988103885577053, TON_PER_DT; dt_val = dt),
                :deposition_silt => to_SI(0.00047599597358497295, TON_PER_DT; dt_val = dt),
                :deposition_lagg => to_SI(0.002939138944433933, TON_PER_DT; dt_val = dt),
                :deposition_clay => to_SI(0.00029924711288315653, TON_PER_DT; dt_val = dt),
            ),
        )
        @test test_means(
            overland_flow.to_river.variables,
            Dict(
                :sagg_rate => to_SI(0.002496392427907965, TON_PER_DT; dt_val = dt),
                :lagg_rate => to_SI(0.0025143486502600633, TON_PER_DT; dt_val = dt),
                :sand_rate => to_SI(0.0030165319515708493, TON_PER_DT; dt_val = dt),
                :sediment_rate => to_SI(0.00867918329675698, TON_PER_DT; dt_val = dt),
                :clay_rate => to_SI(0.0002501736617445959, TON_PER_DT; dt_val = dt),
                :silt_rate => to_SI(0.00040173660527350653, TON_PER_DT; dt_val = dt),
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
            Dict(:density => 2650.0, :d50 => to_SI(0.11542705629386542, MM)),
        )
        @test test_means(
            river_flow.transport_capacity.variables,
            Dict(
                :sediment_transport_capacity =>
                    to_SI(0.024940013212734955, TON_PER_DT; dt_val = dt),
            ),
        )

        # Potential erosion
        @test test_means(
            river_flow.potential_erosion.boundary_conditions,
            Dict(:waterlevel => 0.05777413367369656),
        )
        @test test_means(
            river_flow.potential_erosion.parameters,
            Dict(:d50 => to_SI(0.11542705629386542, MM)),
        )
        @test test_means(
            river_flow.potential_erosion.variables,
            Dict(
                :bed => to_SI(139.44427204500616, TON_PER_DT; dt_val = dt),
                :bank => to_SI(1.8953470656803122, TON_PER_DT; dt_val = dt),
            ),
        )

        # Sediment flux
        @test test_means(
            river_flow.sediment_flux.boundary_conditions,
            Dict(
                :q => 0.1528727680637441,
                :waterlevel => 0.05777413367369656,
                :erosion_land_sagg => to_SI(0.0221033366693814, TON_PER_DT; dt_val = dt),
                :erosion_land_lagg => to_SI(0.022262323062514834, TON_PER_DT; dt_val = dt),
                :erosion_land_sand => to_SI(0.026708709958470814, TON_PER_DT; dt_val = dt),
                :potential_erosion_river_bank =>
                    to_SI(1.8953470656803122, TON_PER_DT; dt_val = dt),
                :transport_capacity => to_SI(0.024940013212734955, TON_PER_DT; dt_val = dt),
                :erosion_land_silt => to_SI(0.003557020658893805, TON_PER_DT; dt_val = dt),
                :erosion_land_clay => to_SI(0.0022150654718924695, TON_PER_DT; dt_val = dt),
                :potential_erosion_river_bed =>
                    to_SI(139.44427204500616, TON_PER_DT; dt_val = dt),
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
                :dm_gravel => to_SI(2000.0, μM),
                :dm_sand => to_SI(200.0, μM),
                :dm_sagg => to_SI(30.0, μM),
                :dm_lagg => to_SI(500.0, μM),
                :dm_clay => to_SI(2.0, μM),
                :dm_silt => to_SI(10.0, μM),
            ),
        )
        @test test_means(
            river_flow.sediment_flux.variables,
            Dict(
                :deposition => to_SI(0.07952559155864818, TON_PER_DT; dt_val = dt),
                :erosion => to_SI(0.0026791396544776722, TON_PER_DT; dt_val = dt),
                :leftover_gravel => 0.0,
                :leftover_clay => to_SI(3.916043968626082e-9, TON),
                :leftover_lagg => to_SI(1.9500782501449444e-15, TON),
                :leftover_sagg => to_SI(3.741274625469958e-15, TON),
                :leftover_sand => to_SI(8.753169499245311e-16, TON),
                :leftover_silt => to_SI(4.643328043695461e-16, TON),
                :sagg_rate => to_SI(0.0011473618310256493, TON_PER_DT; dt_val = dt),
                :sand_rate => to_SI(0.00028353005582982047, TON_PER_DT; dt_val = dt),
                :sediment_rate => to_SI(0.023905673898158595, TON_PER_DT; dt_val = dt),
                :silt_rate => to_SI(0.010106572997892053, TON_PER_DT; dt_val = dt),
                :lagg_rate => to_SI(0.00010707745313511598, TON_PER_DT; dt_val = dt),
                :clay_rate => to_SI(0.012261131560275957, TON_PER_DT; dt_val = dt),
                :gravel_rate => 0.0,
                :store_clay => to_SI(0.0026169325235485177, TON),
                :store_gravel => to_SI(0.0001339569862170842, TON),
                :store_lagg => to_SI(0.022262323062512884, TON),
                :store_sagg => to_SI(0.022103336669377652, TON),
                :store_sand => to_SI(0.027135750167760687, TON),
                :store_silt => to_SI(0.00527329214923135, TON),
            ),
        )

        # Concentrations
        @test test_means(
            river_flow.concentrations.boundary_conditions,
            Dict(
                :waterlevel => 0.05777413367369656,
                :q => 0.1528727680637441,
                :silt => to_SI(0.010106572997892053, TON_PER_DT; dt_val = dt),
                :clay => to_SI(0.012261131560275957, TON_PER_DT; dt_val = dt),
                :sagg => to_SI(0.0011473618310256493, TON_PER_DT; dt_val = dt),
                :sand => to_SI(0.00028353005582982047, TON_PER_DT; dt_val = dt),
                :lagg => to_SI(0.00010707745313511598, TON_PER_DT; dt_val = dt),
                :gravel => 0.0,
            ),
        )
        @test test_means(
            river_flow.concentrations.parameters,
            Dict(
                :dm_clay => to_SI(2.0, μM),
                :dm_gravel => to_SI(2000.0, μM),
                :dm_sagg => to_SI(30.0, μM),
                :dm_lagg => to_SI(500.0, μM),
                :dm_silt => to_SI(10.0, μM),
                :dm_sand => to_SI(200.0, μM),
            ),
        )
        @test test_means(
            river_flow.concentrations.variables,
            Dict(
                :suspended => to_SI(0.14755963695058993, GRAM_PER_M3; dt_val = dt),
                :bed => to_SI(0.0001145940858844022, GRAM_PER_M3; dt_val = dt),
                :total => to_SI(0.14767423103647434, GRAM_PER_M3; dt_val = dt),
            ),
        )
    end

    # run the second timestep
    Wflow.run_timestep!(model)

    @testset "second timestep sediment model engelund (routing)" begin
        land = model.routing.overland_flow
        river = model.routing.river_flow

        @test river.transport_capacity.parameters.d50[1] == to_SI(0.05000000074505806, MM)
        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              to_SI(0.1418728167951134, TON_PER_DT; dt_val = dt)

        @test mean(river.concentrations.variables.suspended) ≈
              to_SI(0.24791810261189964, GRAM_PER_M3)
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
        eros = model.land
        land = model.routing.overland_flow

        @test mean(eros.soil_erosion.variables.soil_erosion_rate) ≈
              to_SI(0.0776983847440198, TON_PER_DT; dt_val = dt)
        @test mean(land.transport_capacity.parameters.c_govers) ≈ 0.16393911236592437
        @test mean(land.transport_capacity.variables.sediment_transport_capacity) ≈
              to_SI(1.0990864706347766e6, TON_PER_DT; dt_val = dt)
        @test mean(land.to_river.variables.sediment_rate) ≈
              to_SI(0.07708434959918917, TON_PER_DT; dt_val = dt)
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
        eros = model.land
        land = model.routing.overland_flow

        @test mean(eros.soil_erosion.variables.soil_erosion_rate) ≈
              to_SI(0.0776983847440198, TON_PER_DT; dt_val = dt)
        @test mean(land.transport_capacity.parameters.d50) ≈ to_SI(0.001534350291334408, MM)
        @test mean(land.transport_capacity.variables.sediment_transport_capacity) ≈
              Float32(to_SI(1.0988158364353527e6, TON_PER_DT; dt_val = dt))
        @test mean(land.to_river.variables.sediment_rate) ≈
              to_SI(0.07759383356462951, TON_PER_DT; dt_val = dt)
    end

    Wflow.close_files(model)
end

### Test all river transport capacity ###
@testitem "Run sediment yang transport capacity" begin
    using Statistics: mean
    using Wflow: to_SI, Unit, MM, TON_PER_DT
    GRAM_PER_M3 = Unit(; g = 1, m = -3)
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

        @test river.transport_capacity.parameters.d50[1] == to_SI(0.05000000074505806, MM)
        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              to_SI(39.955369987101946, TON_PER_DT; dt_val = dt)
        @test mean(river.concentrations.variables.suspended) ≈
              to_SI(0.004036949009419181, GRAM_PER_M3)
    end

    Wflow.close_files(model)
end

@testitem "Run sediment kodatie transport capacity" begin
    using Statistics: mean
    using Wflow: to_SI, Unit, TON_PER_DT
    GRAM_PER_M3 = Unit(; g = 1, m = -3)
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
              to_SI(30.332588671299625, TON_PER_DT; dt_val = dt)

        @test mean(river.concentrations.variables.suspended) ≈
              to_SI(54.7559835316139, GRAM_PER_M3)
    end

    Wflow.close_files(model)
end

@testitem "Run sediment molinas transport capacity" begin
    using Statistics: mean
    using Wflow: to_SI, Unit, MM, TON_PER_DT
    GRAM_PER_M3 = Unit(; g = 1, m = -3)
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

        @test river.transport_capacity.parameters.d50[1] == to_SI(0.05000000074505806, MM)
        @test mean(river.transport_capacity.boundary_conditions.q) ≈ 0.6975180562953642
        @test mean(river.transport_capacity.variables.sediment_transport_capacity) ≈
              to_SI(350.7538564169241, TON_PER_DT; dt_val = dt)

        @test mean(river.concentrations.variables.suspended) ≈
              to_SI(884.749798262824, GRAM_PER_M3)
    end

    Wflow.close_files(model)
end

@testitem "run wflow sediment" begin
    tomlpath = joinpath(@__DIR__, "sediment_eurosem_engelund_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    Wflow.run(tomlpath; silent = true)
end
