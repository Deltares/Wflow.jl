@testitem "unit: Rainfall erosion (Eurosem, Answers)" begin
    dt = 86400.0

    precip = 3.609333601262834e-8
    area = 551370.3190576562

    # Eurosem
    interception = 5.79804570103685e-10
    waterlevel = 0.0
    soil_detachability = 0.002
    eurosem_exponent = 2000.0
    canopy_height = 0.0
    canopy_gap_fraction = 0.1
    soilcover_fraction = 0.0
    @test Wflow.rainfall_erosion_eurosem(
        precip,
        interception,
        waterlevel,
        soil_detachability,
        eurosem_exponent,
        canopy_height,
        canopy_gap_fraction,
        soilcover_fraction,
        area,
        dt,
    ) ≈ 0.05697861124142303

    # Answers
    usle_k = 0.1
    usle_c = 0.01
    answers_rainfall_factor = 0.108
    @test Wflow.rainfall_erosion_answers(
        precip,
        usle_k,
        usle_c,
        answers_rainfall_factor,
        area,
    ) ≈ 4.654493615217501e-6
end

@testitem "unit: Overland flow erosion (Answers)" begin
    dt = 86400.0
    overland_flow = 300.0
    usle_k = 0.1
    usle_c = 0.01
    answers_overland_flow_factor = 0.9
    slope = 0.5
    area = 551370.3190576562

    @test Wflow.overland_flow_erosion_answers(
        overland_flow,
        usle_k,
        usle_c,
        answers_overland_flow_factor,
        slope,
        area,
    ) ≈ 89.66035496171442
end

@testitem "unit: River flow erosion (Julian Torres)" begin
    waterlevel = 0.0085
    d50 = 5.0e-6
    width = 150.0
    length = 635.0
    slope = 1.5e-3
    dt = 86400.0

    bed, bank = Wflow.river_erosion_julian_torres(waterlevel, d50, width, length, slope, dt)
    @test bed ≈ 18.819752773111976
    @test bank == 0.0
end

@testitem "unit: ReservoirModel deposition (Camp)" begin
    dt = 86400.0

    input = 2.3148148148148147e-5
    q = 3.0
    waterlevel = 0.57
    res_trapping_efficiency = 0.75
    dm = 2.9999999999999997e-5
    slope = 2e-4

    # Case: limited deposition, dm < dsuspf
    res_area = 1.48e4

    @test Wflow.reservoir_deposition_camp(
        input,
        q,
        waterlevel,
        res_area,
        res_trapping_efficiency,
        dm,
        slope,
    ) ≈ 1.1733796296296296e-5

    # Case: non-limited deposition, dm < dsuspf
    res_area = 1.48e6

    @test Wflow.reservoir_deposition_camp(
        input,
        q,
        waterlevel,
        res_area,
        res_trapping_efficiency,
        dm,
        slope,
    ) ≈ 2.3148148148148147e-5

    # Case non-limited deposition dm > dsuspf
    dm = 0.00039999999999999996

    @test Wflow.reservoir_deposition_camp(
        input,
        q,
        waterlevel,
        res_area,
        res_trapping_efficiency,
        dm,
        slope,
    ) ≈ 2.3148148148148147e-5
end

@testitem "unit: Transport capacity (Govers, Yalin)" begin
    q = 100.0
    waterlevel = 1.0
    density = 2650.0
    slope = 0.15
    width = 600.0
    reservoirs = false
    rivers = false
    dt = 86400.0

    # Govers
    c_govers = 0.000505
    n_govers = 4.27

    @test Wflow.transport_capacity_govers(
        q,
        waterlevel,
        c_govers,
        n_govers,
        density,
        slope,
        width,
        reservoirs,
        rivers,
        dt,
    ) ≈ 1.1088214231707866e8

    # Yalin
    d50 = 0.0001

    @test Wflow.transport_capacity_yalin(
        q,
        waterlevel,
        density,
        d50,
        slope,
        width,
        reservoirs,
        rivers,
        dt,
    ) ≈ 1.3682410405973883e6
end

@testitem "unit: differentiation (Yalin)" begin
    waterlevel = 30.0
    density = 3020.0
    median_diameter_clay = 2.0e-6
    median_diameter_silt = 9.999999999999999e-6
    median_diameter_sand = 0.00019999999999999998
    median_diameter_small_aggregates = 2.9999999999999997e-5
    median_diameter_large_aggregates = 0.0005
    slope = 0.15

    @test Wflow.transportability_yalin_differentiation(
        waterlevel,
        density,
        median_diameter_clay,
        median_diameter_silt,
        median_diameter_sand,
        median_diameter_small_aggregates,
        median_diameter_large_aggregates,
        slope,
    ) ≈ 2.3511712003217816e7

    q = 30.0
    dm = 2.9999999999999997e-5
    slope = 0.25
    width = 600.0
    reservoirs = false
    rivers = false
    dtot = 0.5
    dt = 86400.0

    @test Wflow.transport_capacity_yalin_differentiation(
        q,
        waterlevel,
        density,
        dm,
        slope,
        width,
        reservoirs,
        rivers,
        dtot,
        dt,
    ) ≈ 1.3194390560823308e16
end

@testitem "unit: Transport capacity (Bagnold, Engelund, Kodatie, Yang, Molinas)" begin
    q = 30.0
    waterlevel = 3.0
    density = 2650.0
    width = 600.0
    length = 635.0
    slope = 0.25
    dt = 86400.0
    d50 = 0.0001

    # Bagnold
    c_bagnold = 1.75e-5
    e_bagnold = 1.4
    @test Wflow.transport_capacity_bagnold(
        q,
        waterlevel,
        c_bagnold,
        e_bagnold,
        width,
        length,
        dt,
    ) ≈ 0.0024513377981088234

    # Engelund
    @test Wflow.transport_capacity_engelund(
        q,
        waterlevel,
        density,
        d50,
        width,
        length,
        slope,
        dt,
    ) ≈ 24.13093433870848

    # Kodatie
    a_kodatie = 2829.6
    b_kodatie = 3.646
    c_kodatie = 0.406
    d_kodatie = 0.412
    @test Wflow.transport_capacity_kodatie(
        q,
        waterlevel,
        a_kodatie,
        b_kodatie,
        c_kodatie,
        d_kodatie,
        width,
        length,
        slope,
        dt,
    ) ≈ 0.00821374503643595

    # Yang
    @test Wflow.transport_capacity_yang(
        q,
        waterlevel,
        density,
        d50,
        width,
        length,
        slope,
        dt,
    ) ≈ 44331.193721554555

    # Molinas
    @test Wflow.transport_capacity_molinas(
        q,
        waterlevel,
        density,
        d50,
        width,
        length,
        slope,
        dt,
    ) ≈ 2.832911300898929e-5
end

@testitem "unit: update SedimentConcentrationsRiverModel" begin
    dt = 86400.0

    n_river_cells = 1
    sediment_concentrations_model = Wflow.SedimentConcentrationsRiverModel(;
        n_river_cells,
        boundary_conditions = Wflow.SedimentConcentrationsRiverBC(;
            n_river_cells,
            q = [2.5],
            waterlevel = [0.2],
            clay = [0.0],
            silt = [0.0],
            sand = [1.1574074074074074e-12],
            small_aggregates = [0.0],
            large_aggregates = [2.3148148148148147e-14],
            gravel = [3.472222222222222e-16],
        ),
        parameters = Wflow.SedimentConcentrationsRiverParameters(;
            median_diameter_clay = [2.0e-6],
            median_diameter_silt = [9.999999999999999e-6],
            median_diameter_sand = [0.00019999999999999998], # dm < dsuspf
            median_diameter_small_aggregates = [2.9999999999999997e-5],
            median_diameter_large_aggregates = [0.0005], # dsuspf < dm < dbedf
            median_diameter_gravel = [0.002], # dbedf < dm
        ),
    )
    parameters = Wflow.RiverParameters(; slope = [1e-3])

    Wflow.update_river_sediment_concentration_model!(
        sediment_concentrations_model,
        parameters,
        dt,
    )
    @test sediment_concentrations_model.variables.suspended[1] ≈ 4.675925925925926e-13
    @test sediment_concentrations_model.variables.bed[1] ≈ 4.768518518518557e-15
    @test sediment_concentrations_model.variables.total[1] ≈ 4.723611111111112e-13
end

@testitem "unit: update SedimentRiverTransportModel" begin
    using Graphs: DiGraph
    dt = 86400.0

    function get_objects()
        n_river_cells = 4
        sediment_flux = Wflow.SedimentRiverTransportModel(;
            n_river_cells,
            boundary_conditions = Wflow.SedimentRiverTransportBC(;
                n_river_cells,
                waterlevel = fill(0.0105, n_river_cells),
                q = fill(0.002, n_river_cells),
                transport_capacity = [
                    5.787037037037036e-9,
                    5.787037037037036e-9,
                    5.787037037037036e-9,
                    2.662037037037037e-10,
                ],
                erosion_land_clay = fill(6.655092592592592e-11, n_river_cells),
                erosion_land_silt = fill(6.655092592592592e-11, n_river_cells),
                erosion_land_sand = fill(6.655092592592592e-11, n_river_cells),
                erosion_land_small_aggregates = fill(6.655092592592592e-11, n_river_cells),
                erosion_land_large_aggregates = fill(6.655092592592592e-11, n_river_cells),
                potential_erosion_river_bed = fill(1.1574074074074074e-10, n_river_cells),
                potential_erosion_river_bank = fill(2.3148148148148147e-10, n_river_cells),
            ),
            parameters = Wflow.SedimentRiverTransportParameters(;
                clay_fraction = fill(0.15, n_river_cells),
                silt_fraction = fill(0.25, n_river_cells),
                sand_fraction = fill(0.35, n_river_cells),
                gravel_fraction = fill(0.45, n_river_cells),
                median_diameter_clay = fill(2.0e-6, n_river_cells),
                median_diameter_silt = fill(9.999999999999999e-6, n_river_cells),
                median_diameter_sand = fill(0.00019999999999999998, n_river_cells),
                median_diameter_small_aggregates = fill(
                    2.9999999999999997e-5,
                    n_river_cells,
                ),
                median_diameter_large_aggregates = fill(0.0005, n_river_cells),
                median_diameter_gravel = fill(0.002, n_river_cells),
                reservoir_outlet = [true, false, false, false],
                reservoir_area = fill(6e5, n_river_cells),
                reservoir_trapping_efficiency = fill(0.5, n_river_cells),
            ),
            variables = Wflow.SedimentRiverTransportVariables(;
                n_river_cells,
                leftover_clay = fill(5.5e-6, n_river_cells),
                store_gravel = fill(5.3e-6, n_river_cells),
            ),
        )

        order = collect(1:n_river_cells)
        graph = DiGraph(n_river_cells)
        domain = Wflow.DomainRiver(;
            network = Wflow.NetworkRiver(; order, graph),
            parameters = Wflow.RiverParameters(;
                slope = fill(1e-3, n_river_cells),
                flow_width = fill(4.2, n_river_cells),
                flow_length = [785.0, 785.0, 785.0, 7850.0],
                reservoir_coverage = [false, true, false, false],
            ),
        )
        return sediment_flux, domain, graph, order, n_river_cells
    end

    function perform_tests(variables)
        @test variables.sediment_rate ≈
              [0.0, 3.96412037037037e-10, 8.238374498676117e-11, 1.3287692996519204e-10]
        @test variables.clay_rate ≈
              [0.0, 1.3020833333333334e-10, 8.238374462152822e-11, 6.499415052645262e-11]
        @test variables.silt_rate ≈
              [0.0, 6.655092592592592e-11, 3.6523294471267653e-19, 3.3219232491298004e-11]
        @test variables.sand_rate ≈
              [0.0, 6.655092592592592e-11, 0.0, 1.4443144561433946e-12]
        @test variables.small_aggregates_rate ≈
              [0.0, 6.655092592592592e-11, 0.0, 3.3219232491298004e-11]
        @test variables.large_aggregates_rate ≈ [0.0, 6.655092592592592e-11, 0.0, 0.0]
        @test variables.gravel_rate ≈ [0.0, 0.0, 0.0, 0.0]
        @test variables.deposition ≈
              [8.744212962962962e-10, 0.0, 7.920375513095351e-10, 1.3020833333333334e-10]
        @test variables.erosion ≈ [4.780092592592592e-10, 0.0, 4.780092592592592e-10, 0.0]
        @test variables.leftover_clay ≈ [0.0, 0.0, 0.0, 5.634505394514494e-6]
        @test variables.leftover_silt ≈ [0.0, 0.0, 0.0, 2.8798583127518524e-6]
        @test variables.leftover_sand ≈ [0.0, 0.0, 0.0, 1.2521123098921127e-7]
        @test variables.leftover_small_aggregates ≈ [0.0, 0.0, 0.0, 2.8798583127518524e-6]
        @test variables.leftover_large_aggregates ≈ zeros(4)
        @test variables.leftover_gravel ≈ zeros(4)
        @test variables.store_clay ≈ [1.575e-5, 0.0, 8.63204446469996e-6, 0.0]
        @test variables.store_silt ≈ [1.325e-5, 0.0, 1.3249999968443874e-5, 0.0]
        @test variables.store_sand ≈ [1.625e-5, 0.0, 1.625e-5, 5.5e-6]
        @test variables.store_small_aggregates ≈ [5.75e-6, 0.0, 5.75e-6, 0.0]
        @test variables.store_large_aggregates ≈ [5.75e-6, 0.0, 5.75e-6, 5.75e-6]
        @test variables.store_gravel ≈ [1.88e-5, 5.3e-6, 1.88e-5, 5.3e-6]
    end

    sediment_flux, domain, graph, order, n_river_cells = get_objects()

    input_particles =
        Wflow.compute_sediment_input.(Ref(sediment_flux), Ref(graph), dt, order)
    input_particles_expected = [
        1.3020833333333334e-10,
        6.655092592592592e-11,
        6.655092592592592e-11,
        6.655092592592592e-11,
        6.655092592592592e-11,
        0.0,
    ]
    @test all(x -> collect(x) ≈ input_particles_expected, input_particles)

    sediment_need =
        max.(
            sediment_flux.boundary_conditions.transport_capacity .-
            sum(input_particles_expected),
            0.0,
        )
    sediment_need[2] = 0.0
    @test sediment_need ≈ [5.390625e-9, 0.0, 5.390625e-9, 0.0]

    store_sediment = sum(
        getfield(sediment_flux.variables, name) for
        name in propertynames(sediment_flux.variables) if startswith(string(name), "store")
    )
    @test all(≈(5.3e-6), store_sediment)

    erosion_particles =
        Wflow.compute_direct_river_erosion.(
            Ref(sediment_flux),
            sediment_need,
            store_sediment,
            dt,
            order,
        )
    erosion_particles_expected = [
        5.208333333333333e-11,
        8.680555555555555e-11,
        1.2152777777777776e-10,
        0.0,
        0.0,
        1.5625e-10,
    ]
    @test collect.(erosion_particles) ≈
          [erosion_particles_expected, zeros(6), erosion_particles_expected, zeros(6)]

    store_erosion =
        Wflow.compute_store_erosion!.(
            Ref(sediment_flux.variables),
            sediment_need,
            dt,
            order,
        )
    store_erosion_expected = [0.0, 0.0, 0.0, 0.0, 0.0, 6.134259259259259e-11]
    @test collect.(store_erosion) ≈
          [store_erosion_expected, zeros(6), store_erosion_expected, zeros(6)]

    erosion_particles = [erosion_particles[i] .+ store_erosion[i] for i in 1:n_river_cells]
    @. sediment_flux.variables.erosion = sum(erosion_particles)
    erosion_expected = 4.780092592592592e-10
    @test sediment_flux.variables.erosion ≈ [erosion_expected, 0.0, erosion_expected, 0.0]

    # Index 1: reservoir outlet
    deposition_particles_1 = Wflow.compute_reservoir_deposition(
        sediment_flux,
        domain.parameters,
        input_particles[1][1:6],
        erosion_particles[1],
        1,
    )
    deposition_particles_1_expected = [
        1.8229166666666664e-10,
        1.5335648148148147e-10,
        1.8807870370370368e-10,
        6.655092592592592e-11,
        6.655092592592592e-11,
        2.1759259259259257e-10,
    ]
    @test collect(deposition_particles_1) ≈ deposition_particles_1_expected

    # Index 2: reservoir coverage
    deposition_particles_2 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Index 3: natural deposition in the river
    excess_sediment = max(
        input_particles[3][1] - sediment_flux.boundary_conditions.transport_capacity[3],
        0.0,
    )
    @test excess_sediment == 0.0
    deposition_particles_3 = Wflow.compute_natural_deposition(
        sediment_flux,
        domain.parameters,
        input_particles[3][1:6],
        erosion_particles[3],
        3,
    )
    deposition_particles_3_expected = [
        9.990792204513842e-11,
        1.5335648111624853e-10,
        1.8807870370370368e-10,
        6.655092592592592e-11,
        6.655092592592592e-11,
        2.1759259259259257e-10,
    ]
    @test collect(deposition_particles_3) ≈ deposition_particles_3_expected

    # Index 4: deposition in the river from transport capacity exceeding
    excess_sediment = max(
        sum(input_particles[4]) - sediment_flux.boundary_conditions.transport_capacity[4],
        0.0,
    )
    @test excess_sediment ≈ 1.3020833333333334e-10
    deposition_particles_4 = Wflow.compute_transport_capacity_deposition(
        excess_sediment,
        input_particles[4][1:6],
        erosion_particles[4],
    )
    deposition_particles_4_expected =
        [0.0, 0.0, 6.36574074074074e-11, 0.0, 6.655092592592592e-11, 0.0]
    @test collect(deposition_particles_4) ≈ deposition_particles_4_expected

    deposition_particles = [
        deposition_particles_1,
        deposition_particles_2,
        deposition_particles_3,
        deposition_particles_4,
    ]

    fwaterout =
        Wflow.water_outflow_fraction.(
            sediment_flux.boundary_conditions.waterlevel,
            sediment_flux.boundary_conditions.q,
            domain.parameters.flow_width,
            domain.parameters.flow_length,
            dt,
        )
    @test fwaterout ≈ [1.0, 1.0, 1.0, 0.49915507604315607]

    Wflow.update_variables!.(
        Ref(sediment_flux.variables),
        input_particles,
        erosion_particles,
        deposition_particles,
        fwaterout,
        dt,
        order,
    )

    perform_tests(sediment_flux.variables)

    # Perform calculations again by directly calling update!
    sediment_flux, domain, _, _, _ = get_objects()
    Wflow.update_sediment_river_transport_model!(sediment_flux, domain, dt)
    perform_tests(sediment_flux.variables)
end
