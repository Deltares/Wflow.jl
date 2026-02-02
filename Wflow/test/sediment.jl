@testitem "unit: Rainfall erosion (Eurosem, Answers)" begin
    using Wflow: Unit, to_SI, from_SI, TON_PER_DT, MM_PER_DT
    dt = 86400.0
    GRAM_PER_J = Unit(; g = 1, J = -1)

    precip = to_SI(3.11846423149108887, MM_PER_DT; dt_val = dt)
    area = 551370.3190576562

    # Eurosem
    interception = to_SI(0.05009511485695839, MM_PER_DT; dt_val = dt)
    waterlevel = 0.0
    soil_detachability = to_SI(2.0, GRAM_PER_J)
    eurosem_exponent = 2.0
    canopyheight = 0.0
    canopygapfraction = 0.1
    soilcover_fraction = 0.0
    @test Wflow.rainfall_erosion_eurosem(
        precip,
        interception,
        waterlevel,
        soil_detachability,
        eurosem_exponent,
        canopyheight,
        canopygapfraction,
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
    ) ≈ to_SI(0.0004021482483547921, TON_PER_DT; dt_val = dt)
end

@testitem "unit: Overland flow erosion (Answers)" begin
    using Wflow: TON_PER_DT, to_SI
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
    ) ≈ to_SI(7746.6546686921265, TON_PER_DT; dt_val = dt)
end

@testitem "unit: River flow erosion (Julian Torres)" begin
    using Wflow: to_SI, TON_PER_DT, MM
    waterlevel = 0.0085
    d50 = to_SI(0.005, MM)
    width = 150.0
    length = 635.0
    slope = 1.5e-3
    dt = 86400.0

    bed, bank = Wflow.river_erosion_julian_torres(waterlevel, d50, width, length, slope, dt)
    @test bed ≈ to_SI(1626.026639596875, TON_PER_DT; dt_val = dt)
    @test bank == 0.0
end

@testitem "unit: Reservoir deposition (Camp)" begin
    using Wflow: to_SI, Unit, TON_PER_DT
    M3_PER_DT = Unit(; m = 3, dt = -1)
    μM = Unit(; μm = 1)
    dt = 86400.0

    input = to_SI(0.002, TON_PER_DT; dt_val = dt)
    q = 3.0
    waterlevel = 0.57
    res_trapping_efficiency = 0.75
    dm = to_SI(30.0, μM)
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
    ) ≈ to_SI(0.0010138, TON_PER_DT; dt_val = dt)

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
    ) ≈ to_SI(0.002, TON_PER_DT; dt_val = dt)

    # Case non-limited deposition dm > dsuspf
    dm = to_SI(400.0, μM)

    @test Wflow.reservoir_deposition_camp(
        input,
        q,
        waterlevel,
        res_area,
        res_trapping_efficiency,
        dm,
        slope,
    ) ≈ to_SI(0.002, TON_PER_DT; dt_val = dt)
end

@testitem "unit: Transport capacity (Govers, Yalin)" begin
    using Wflow: to_SI, TON_PER_DT, MM

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
    ) ≈ to_SI(9.580217096195597e9, TON_PER_DT; dt_val = dt)

    # Yalin
    d50 = to_SI(0.1, MM)

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
    ) ≈ to_SI(1.1821602590761435e8, TON_PER_DT; dt_val = dt)
end

@testitem "unit: differentiation (Yalin)" begin
    using Wflow: Unit, to_SI, TON_PER_DT
    μM = Unit(; μm = 1)

    waterlevel = 30.0
    density = 3020.0
    dm_clay = to_SI(2.0, μM)
    dm_silt = to_SI(10.0, μM)
    dm_sand = to_SI(200.0, μM)
    dm_sagg = to_SI(30.0, μM)
    dm_lagg = to_SI(500.0, μM)
    slope = 0.15

    @test Wflow.transportability_yalin_differentiation(
        waterlevel,
        density,
        dm_clay,
        dm_silt,
        dm_sand,
        dm_sagg,
        dm_lagg,
        slope,
    ) ≈ 2.3511712003217816e7

    q = 30.0
    dm = to_SI(30.0, μM)
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
    ) ≈ to_SI(1.1399953444551338e18, TON_PER_DT; dt_val = dt)
end

@testitem "unit: Transport capacity (Bagnold, Engelund, Kodatie, Yang, Molinas)" begin
    using Wflow: to_SI, TON_PER_DT, MM

    q = 30.0
    waterlevel = 3.0
    density = 2650.0
    width = 600.0
    length = 635.0
    slope = 0.25
    dt = 86400.0
    d50 = to_SI(0.1, MM)

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
    ) ≈ to_SI(0.21179558575660237, TON_PER_DT; dt_val = dt)

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
    ) ≈ to_SI(2084.912726864413, TON_PER_DT; dt_val = dt)

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
    ) ≈ to_SI(0.709667571148066, TON_PER_DT; dt_val = dt)

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
    ) ≈ to_SI(3.830215137542314e6, TON_PER_DT; dt_val = dt)

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
    ) ≈ to_SI(0.002447635363976675, TON_PER_DT; dt_val = dt)
end

@testitem "unit: update SedimentConcentrationsRiverModel" begin
    using Wflow: to_SI, Unit
    GRAM_PER_M3 = Unit(; g = 1, m = -3)
    μM = Unit(; μm = 1)

    n = 1
    model = Wflow.SedimentConcentrationsRiverModel(;
        n,
        boundary_conditions = Wflow.SedimentConcentrationsRiverBC(;
            n,
            q = [2.5],
            waterlevel = [0.2],
            clay = [to_SI(0.0, GRAM_PER_M3)],
            silt = [to_SI(0.0, GRAM_PER_M3)],
            sand = [to_SI(1.0e-10, GRAM_PER_M3)],
            sagg = [to_SI(0.0, GRAM_PER_M3)],
            lagg = [to_SI(2.0e-12, GRAM_PER_M3)],
            gravel = [to_SI(3.0e-14, GRAM_PER_M3)],
        ),
        parameters = Wflow.SedimentConcentrationsRiverParameters(;
            dm_clay = [to_SI(2.0, μM)],
            dm_silt = [to_SI(10.0, μM)],
            dm_sand = [to_SI(200.0, μM)], # dm < dsuspf
            dm_sagg = [to_SI(30.0, μM)],
            dm_lagg = [to_SI(500.0, μM)], # dsuspf < dm < dbedf
            dm_gravel = [to_SI(2000.0, μM)], # dbedf < dm
        ),
    )
    parameters = Wflow.RiverParameters(; slope = [1e-3])
    dt = 86400.0

    Wflow.update!(model, parameters, dt)
    @test model.variables.suspended[1] ≈ to_SI(4.675925925925926e-10, GRAM_PER_M3)
    @test model.variables.bed[1] ≈ to_SI(4.768518518518557e-12, GRAM_PER_M3)
    @test model.variables.total[1] ≈ to_SI(4.723611111111112e-10, GRAM_PER_M3)
end

@testitem "unit: update SedimentRiverTransportModel" begin
    using Wflow: to_SI, Unit, TON_PER_DT
    using Graphs: DiGraph
    μM = Unit(; μm = 1)
    TON = Unit(; t = 1)
    dt = 86400.0

    function get_objects()
        n = 4
        model = Wflow.SedimentRiverTransportModel(;
            n,
            boundary_conditions = Wflow.SedimentRiverTransportBC(;
                n,
                waterlevel = fill(0.0105, n),
                q = fill(0.002, n),
                transport_capacity = to_SI(
                    [5e-7, 5e-7, 5e-7, 2.3e-8],
                    TON_PER_DT;
                    dt_val = dt,
                ),
                erosion_land_clay = fill(to_SI(5.75e-9, TON_PER_DT; dt_val = dt), n),
                erosion_land_silt = fill(to_SI(5.75e-9, TON_PER_DT; dt_val = dt), n),
                erosion_land_sand = fill(to_SI(5.75e-9, TON_PER_DT; dt_val = dt), n),
                erosion_land_sagg = fill(to_SI(5.75e-9, TON_PER_DT; dt_val = dt), n),
                erosion_land_lagg = fill(to_SI(5.75e-9, TON_PER_DT; dt_val = dt), n),
                potential_erosion_river_bed = fill(to_SI(1e-8, TON_PER_DT; dt_val = dt), n),
                potential_erosion_river_bank = fill(
                    to_SI(2e-8, TON_PER_DT; dt_val = dt),
                    n,
                ),
            ),
            parameters = Wflow.SedimentRiverTransportParameters(;
                clay_fraction = fill(0.15, n),
                silt_fraction = fill(0.25, n),
                sand_fraction = fill(0.35, n),
                gravel_fraction = fill(0.45, n),
                dm_clay = fill(to_SI(2.0, μM), n),
                dm_silt = fill(to_SI(10.0, μM), n),
                dm_sand = fill(to_SI(200.0, μM), n),
                dm_sagg = fill(to_SI(30.0, μM), n),
                dm_lagg = fill(to_SI(500.0, μM), n),
                dm_gravel = fill(to_SI(2000.0, μM), n),
                reservoir_outlet = [true, false, false, false],
                reservoir_area = fill(6e5, n),
                reservoir_trapping_efficiency = fill(0.5, n),
            ),
            variables = Wflow.SedimentRiverTransportVariables(;
                n,
                leftover_clay = fill(to_SI(5.5e-9, TON), n),
                store_gravel = fill(to_SI(5.3e-9, TON), n),
            ),
        )

        order = collect(1:n)
        graph = DiGraph(n)
        domain = Wflow.DomainRiver(;
            network = Wflow.NetworkRiver(; order, graph),
            parameters = Wflow.RiverParameters(;
                slope = fill(1e-3, n),
                flow_width = fill(4.2, n),
                flow_length = [785.0, 785.0, 785.0, 7850.0],
                reservoir_coverage = [false, true, false, false],
            ),
        )
        return model, domain, graph, order, n
    end

    function perform_tests(variables)
        @test variables.sediment_rate ≈ to_SI(
            [0.0, 3.425e-8, 7.117955566856165e-9, 1.1480566748992592e-8],
            TON_PER_DT;
            dt_val = dt,
        )
        @test variables.clay_rate ≈ to_SI(
            [0.0, 1.125e-8, 7.117955535300039e-9, 5.615494605485506e-9],
            TON_PER_DT;
            dt_val = dt,
        )
        @test variables.silt_rate ≈ to_SI(
            [0.0, 5.75e-9, 3.1556126423175254e-17, 2.8701416872481476e-9],
            TON_PER_DT;
            dt_val = dt,
        )
        @test variables.sand_rate ≈
              to_SI([0.0, 5.75e-9, 0.0, 1.247887690107893e-10], TON_PER_DT; dt_val = dt)
        @test variables.sagg_rate ≈
              to_SI([0.0, 5.75e-9, 0.0, 2.8701416872481476e-9], TON_PER_DT; dt_val = dt)
        @test variables.lagg_rate ≈ to_SI([0.0, 5.75e-9, 0.0, 0.0], TON_PER_DT; dt_val = dt)
        @test variables.gravel_rate ≈ to_SI([0.0, 0.0, 0.0, 0.0], TON_PER_DT; dt_val = dt)
        @test variables.deposition ≈ to_SI(
            [7.555e-8, 0.0, 6.843204443314384e-8, 1.125e-8],
            TON_PER_DT;
            dt_val = dt,
        )
        @test variables.erosion ≈
              to_SI([4.13e-8, 0.0, 4.13e-8, 0.0], TON_PER_DT; dt_val = dt)
        @test variables.leftover_clay ≈ to_SI([0.0, 0.0, 0.0, 5.634505394514494e-9], TON)
        @test variables.leftover_silt ≈ to_SI([0.0, 0.0, 0.0, 2.8798583127518526e-9], TON)
        @test variables.leftover_sand ≈ to_SI([0.0, 0.0, 0.0, 1.2521123098921128e-10], TON)
        @test variables.leftover_sagg ≈ to_SI([0.0, 0.0, 0.0, 2.8798583127518526e-9], TON)
        @test variables.leftover_lagg ≈ zeros(4)
        @test variables.leftover_gravel ≈ zeros(4)
        @test variables.store_clay ≈ to_SI([1.575e-8, 0.0, 8.63204446469996e-9, 0.0], TON)
        @test variables.store_silt ≈ to_SI([1.325e-8, 0.0, 1.3249999968443875e-8, 0.0], TON)
        @test variables.store_sand ≈ to_SI([1.625e-8, 0.0, 1.625e-8, 5.5e-9], TON)
        @test variables.store_sagg ≈ to_SI([5.75e-9, 0.0, 5.75e-9, 0.0], TON)
        @test variables.store_lagg ≈ to_SI([5.75e-9, 0.0, 5.75e-9, 5.75e-9], TON)
        @test variables.store_gravel ≈ to_SI([1.88e-8, 5.3e-9, 1.88e-8, 5.3e-9], TON)
    end

    model, domain, graph, order, n = get_objects()

    input_particles = Wflow.calculate_sediment_input.(Ref(model), Ref(graph), dt, order)
    input_particles_expected =
        to_SI([1.125e-8, 5.75e-9, 5.75e-9, 5.75e-9, 5.75e-9, 0.0], TON_PER_DT; dt_val = dt)
    @test all(x -> collect(x) ≈ input_particles_expected, input_particles)

    sediment_need =
        max.(
            model.boundary_conditions.transport_capacity .- sum(input_particles_expected),
            0.0,
        )
    sediment_need[2] = 0.0
    @test sediment_need ≈ to_SI([4.6575e-7, 0.0, 4.6575e-7, 0.0], TON_PER_DT; dt_val = dt)

    store_sediment = sum(
        getfield(model.variables, name) for
        name in propertynames(model.variables) if startswith(string(name), "store")
    )
    @test all(≈(to_SI(5.3e-9, TON)), store_sediment)

    erosion_particles =
        Wflow.calculate_direct_river_erosion.(
            Ref(model),
            sediment_need,
            store_sediment,
            dt,
            order,
        )
    erosion_particles_expected =
        to_SI([4.5e-9, 7.5e-9, 1.05e-8, 0.0, 0.0, 1.35e-8], TON_PER_DT; dt_val = dt)
    @test collect.(erosion_particles) ≈
          [erosion_particles_expected, zeros(6), erosion_particles_expected, zeros(6)]

    store_erosion =
        Wflow.calculate_store_erosion!.(Ref(model.variables), sediment_need, dt, order)
    store_erosion_expected =
        to_SI([0.0, 0.0, 0.0, 0.0, 0.0, 5.3e-9], TON_PER_DT; dt_val = dt)
    @test collect.(store_erosion) ≈
          [store_erosion_expected, zeros(6), store_erosion_expected, zeros(6)]

    erosion_particles = [erosion_particles[i] .+ store_erosion[i] for i in 1:n]
    @. model.variables.erosion = sum(erosion_particles)
    erosion_expected = to_SI(4.13e-8, TON_PER_DT; dt_val = dt)
    @test model.variables.erosion ≈ [erosion_expected, 0.0, erosion_expected, 0.0]

    # Index 1: reservoir outlet
    deposition_particles_1 = Wflow.calculate_reservoir_deposition(
        model,
        domain.parameters,
        input_particles[1][1:6],
        erosion_particles[1],
        1,
    )
    deposition_particles_1_expected = to_SI(
        [1.575e-8, 1.325e-8, 1.625e-8, 5.75e-9, 5.75e-9, 1.88e-8],
        TON_PER_DT;
        dt_val = dt,
    )
    @test collect(deposition_particles_1) ≈ deposition_particles_1_expected

    # Index 2: reservoir coverage
    deposition_particles_2 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Index 3: natural deposition in the river
    excess_sediment =
        max(input_particles[3][1] - model.boundary_conditions.transport_capacity[3], 0.0)
    @test excess_sediment == 0.0
    deposition_particles_3 = Wflow.calculate_natural_deposition(
        model,
        domain.parameters,
        input_particles[3][1:6],
        erosion_particles[3],
        3,
    )
    deposition_particles_3_expected = to_SI(
        [8.63204446469996e-9, 1.3249999968443875e-8, 1.625e-8, 5.75e-9, 5.75e-9, 1.88e-8],
        TON_PER_DT;
        dt_val = dt,
    )
    @test collect(deposition_particles_3) ≈ deposition_particles_3_expected

    # Index 4: deposition in the river from transport capacity exceeding
    excess_sediment =
        max(sum(input_particles[4]) - model.boundary_conditions.transport_capacity[4], 0.0)
    @test excess_sediment ≈ to_SI(1.125e-8, TON_PER_DT; dt_val = dt)
    deposition_particles_4 = Wflow.calculate_transport_capacity_deposition(
        excess_sediment,
        input_particles[4][1:6],
        erosion_particles[4],
    )
    deposition_particles_4_expected =
        to_SI([0.0, 0.0, 5.5e-9, 0.0, 5.75e-9, 0.0], TON_PER_DT; dt_val = dt)
    @test collect(deposition_particles_4) ≈ deposition_particles_4_expected

    deposition_particles = [
        deposition_particles_1,
        deposition_particles_2,
        deposition_particles_3,
        deposition_particles_4,
    ]

    fwaterout =
        Wflow.water_outflow_fraction.(
            model.boundary_conditions.waterlevel,
            model.boundary_conditions.q,
            domain.parameters.flow_width,
            domain.parameters.flow_length,
            dt,
        )
    @test fwaterout ≈ [1.0, 1.0, 1.0, 0.49915507604315607]

    Wflow.update_variables!.(
        Ref(model.variables),
        input_particles,
        erosion_particles,
        deposition_particles,
        fwaterout,
        dt,
        order,
    )

    perform_tests(model.variables)

    # Perform calculations again by directly calling update!
    model, domain, _, _, _ = get_objects()
    Wflow.update!(model, domain, dt)
    perform_tests(model.variables)
end
