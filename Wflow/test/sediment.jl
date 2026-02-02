@testitem "unit: Rainfall erosion (Eurosem, Answers)" begin
    precip = 3.11846423149108887
    area = 551370.3190576562
    dt = 86400.0

    # Eurosem
    interception = 0.05009511485695839
    waterlevel = 0.0
    soil_detachability = 2.0
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
    ) ≈ 0.004922952011258949

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
        dt,
    ) ≈ 0.0004021482483547921
end

@testitem "unit: Overland flow erosion (Answers)" begin
    overland_flow = 300.0
    usle_k = 0.1
    usle_c = 0.01
    answers_overland_flow_factor = 0.9
    slope = 0.5
    area = 551370.3190576562
    dt = 86400.0

    @test Wflow.overland_flow_erosion_answers(
        overland_flow,
        usle_k,
        usle_c,
        answers_overland_flow_factor,
        slope,
        area,
        dt,
    ) ≈ 7746.6546686921265
end

@testitem "unit: River flow erosion (Julian Torres)" begin
    waterlevel = 0.0085
    d50 = 0.005
    width = 150.0
    length = 635.0
    slope = 1.5e-3
    dt = 86400.0

    bed, bank = Wflow.river_erosion_julian_torres(waterlevel, d50, width, length, slope, dt)
    @test bed ≈ 1626.026639596875
    @test bank == 0.0
end

@testitem "unit: Reservoir deposition (Camp)" begin
    input = 0.002
    q = 3.0
    waterlevel = 0.57
    res_area = 1.48e6
    res_trapping_efficiency = 0.75
    dm = 30.0
    slope = 2e-4

    @test Wflow.reservoir_deposition_camp(
        input,
        q,
        waterlevel,
        res_area,
        res_trapping_efficiency,
        dm,
        slope,
    ) ≈ 0.002
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
    ) ≈ 9.580217096195597e9

    # Yalin
    d50 = 0.1

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
    ) ≈ 1.1821602590761435e8
end

@testitem "unit: differentiation (Yalin)" begin
    waterlevel = 30.0
    density = 3020.0
    dm_clay = 2.0
    dm_silt = 10.0
    dm_sand = 200.0
    dm_sagg = 30.0
    dm_lagg = 500.0
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
    dm = 30.0
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
    ) ≈ 1.1399953444551338e18
end

@testitem "unit: Transport capacity (Bagnold, Engelund, Kodatie, Yang, Molinas)" begin
    q = 30.0
    waterlevel = 3.0
    density = 2650.0
    width = 600.0
    length = 635.0
    slope = 0.25
    dt = 86400.0
    d50 = 0.1

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
    ) ≈ 0.21179558575660237

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
    ) ≈ 2084.912726864413

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
    ) ≈ 0.709667571148066

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
    ) ≈ 3.830215137542314e6

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
    ) ≈ 0.002447635363976675
end

@testitem "unit: update SedimentConcentrationsRiverModel" begin
    n = 1
    model = Wflow.SedimentConcentrationsRiverModel(;
        n,
        boundary_conditions = Wflow.SedimentConcentrationsRiverBC(;
            n,
            q = [2.5],
            waterlevel = [0.2],
            clay = [0.0],
            silt = [0.0],
            sand = [1.0e-10],
            sagg = [0.0],
            lagg = [2.0e-12],
            gravel = [3.0e-14],
        ),
        parameters = Wflow.SedimentConcentrationsRiverParameters(;
            dm_clay = [2.0],
            dm_silt = [10.0],
            dm_sand = [200.0], # dm < dsuspf
            dm_sagg = [30.0],
            dm_lagg = [500.0], # dsuspf < dm < dbedf
            dm_gravel = [2000.0], # dbedf < dm
        ),
    )
    parameters = Wflow.RiverParameters(; slope = [1e-3])
    dt = 86400.0

    Wflow.update!(model, parameters, dt)
    @test model.variables.suspended[1] ≈ 4.675925925925926e-10
    @test model.variables.bed[1] ≈ 4.768518518518557e-12
    @test model.variables.total[1] ≈ 4.723611111111112e-10
end

@testitem "unit: update SedimentRiverTransportModel" begin
    using Graphs: DiGraph

    function get_objects()
        n = 4
        model = Wflow.SedimentRiverTransportModel(;
            n,
            boundary_conditions = Wflow.SedimentRiverTransportBC(;
                n,
                waterlevel = fill(0.0105, n),
                q = fill(0.002, n),
                transport_capacity = [5e-7, 5e-7, 5e-7, 2.3e-8],
                erosion_land_clay = fill(5.75e-9, n),
                erosion_land_silt = fill(5.75e-9, n),
                erosion_land_sand = fill(5.75e-9, n),
                erosion_land_sagg = fill(5.75e-9, n),
                erosion_land_lagg = fill(5.75e-9, n),
                potential_erosion_river_bed = fill(1e-8, n),
                potential_erosion_river_bank = fill(2e-8, n),
            ),
            parameters = Wflow.SedimentRiverTransportParameters(;
                clay_fraction = fill(0.15, n),
                silt_fraction = fill(0.25, n),
                sand_fraction = fill(0.35, n),
                gravel_fraction = fill(0.45, n),
                dm_clay = fill(2.0, n),
                dm_silt = fill(10.0, n),
                dm_sand = fill(200.0, n),
                dm_sagg = fill(30.0, n),
                dm_lagg = fill(500.0, n),
                dm_gravel = fill(2000.0, n),
                reservoir_outlet = [true, false, false, false],
                reservoir_area = fill(6e5, n),
                reservoir_trapping_efficiency = fill(0.5, n),
            ),
            variables = Wflow.SedimentRiverTransportVariables(;
                n,
                leftover_clay = fill(5.5e-9, n),
                store_gravel = fill(5.3e-9, n),
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
        @test variables.sediment_rate ≈
              [0.0, 3.425e-8, 7.117955566856165e-9, 1.1480566748992592e-8]
        @test variables.clay_rate ≈
              [0.0, 1.125e-8, 7.117955535300039e-9, 5.615494605485506e-9]
        @test variables.silt_rate ≈
              [0.0, 5.75e-9, 3.1556126423175254e-17, 2.8701416872481476e-9]
        @test variables.sand_rate ≈ [0.0, 5.75e-9, 0.0, 1.247887690107893e-10]
        @test variables.sagg_rate ≈ [0.0, 5.75e-9, 0.0, 2.8701416872481476e-9]
        @test variables.lagg_rate ≈ [0.0, 5.75e-9, 0.0, 0.0]
        @test variables.gravel_rate ≈ [0.0, 0.0, 0.0, 0.0]
        @test variables.deposition ≈ [7.555e-8, 0.0, 6.843204443314384e-8, 1.125e-8]
        @test variables.erosion ≈ [4.13e-8, 0.0, 4.13e-8, 0.0]
        @test variables.leftover_clay ≈ [0.0, 0.0, 0.0, 5.634505394514494e-9]
        @test variables.leftover_silt ≈ [0.0, 0.0, 0.0, 2.8798583127518526e-9]
        @test variables.leftover_sand ≈ [0.0, 0.0, 0.0, 1.2521123098921128e-10]
        @test variables.leftover_sagg ≈ [0.0, 0.0, 0.0, 2.8798583127518526e-9]
        @test variables.leftover_lagg ≈ zeros(4)
        @test variables.leftover_gravel ≈ zeros(4)
        @test variables.store_clay ≈ [1.575e-8, 0.0, 8.63204446469996e-9, 0.0]
        @test variables.store_silt ≈ [1.325e-8, 0.0, 1.3249999968443875e-8, 0.0]
        @test variables.store_sand ≈ [1.625e-8, 0.0, 1.625e-8, 5.5e-9]
        @test variables.store_sagg ≈ [5.75e-9, 0.0, 5.75e-9, 0.0]
        @test variables.store_lagg ≈ [5.75e-9, 0.0, 5.75e-9, 5.75e-9]
        @test variables.store_gravel ≈ [1.88e-8, 5.3e-9, 1.88e-8, 5.3e-9]
    end

    dt = 86400.0
    model, domain, graph, order, n = get_objects()

    input_particles = Wflow.compute_sediment_input.(Ref(model), Ref(graph), order)
    input_particles_expected = [1.125e-8, 5.75e-9, 5.75e-9, 5.75e-9, 5.75e-9, 0.0]
    @test all(x -> collect(x) ≈ input_particles_expected, input_particles)

    sediment_need =
        max.(
            model.boundary_conditions.transport_capacity .- sum(input_particles_expected),
            0.0,
        )
    sediment_need[2] = 0.0
    @test sediment_need ≈ [4.6575e-7, 0.0, 4.6575e-7, 0.0]

    store_sediment = sum(
        getfield(model.variables, name) for
        name in propertynames(model.variables) if startswith(string(name), "store")
    )
    @test all(≈(5.3e-9), store_sediment)

    erosion_particles =
        Wflow.compute_direct_river_erosion.(
            Ref(model),
            sediment_need,
            store_sediment,
            order,
        )
    erosion_particles_expected = [4.5e-9, 7.5e-9, 1.05e-8, 0.0, 0.0, 1.35e-8]
    @test collect.(erosion_particles) ≈
          [erosion_particles_expected, zeros(6), erosion_particles_expected, zeros(6)]

    store_erosion =
        Wflow.compute_store_erosion!.(Ref(model.variables), sediment_need, order)
    store_erosion_expected = [0.0, 0.0, 0.0, 0.0, 0.0, 5.3e-9]
    @test collect.(store_erosion) ≈
          [store_erosion_expected, zeros(6), store_erosion_expected, zeros(6)]

    erosion_particles = [erosion_particles[i] .+ store_erosion[i] for i in 1:n]
    @. model.variables.erosion = sum(erosion_particles)
    erosion_expected = 4.13e-8
    @test model.variables.erosion ≈ [erosion_expected, 0.0, erosion_expected, 0.0]

    # Index 1: reservoir outlet
    deposition_particles_1 = Wflow.compute_reservoir_deposition(
        model,
        domain.parameters,
        input_particles[1][1:6],
        erosion_particles[1],
        1,
    )
    deposition_particles_1_expected =
        [1.575e-8, 1.325e-8, 1.625e-8, 5.75e-9, 5.75e-9, 1.88e-8]
    @test collect(deposition_particles_1) ≈ deposition_particles_1_expected

    # Index 2: reservoir coverage
    deposition_particles_2 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Index 3: natural deposition in the river
    excess_sediment =
        max(input_particles[3][1] - model.boundary_conditions.transport_capacity[3], 0.0)
    @test excess_sediment == 0.0
    deposition_particles_3 = Wflow.compute_natural_deposition(
        model,
        domain.parameters,
        input_particles[3][1:6],
        erosion_particles[3],
        3,
    )
    deposition_particles_3_expected =
        [8.63204446469996e-9, 1.3249999968443875e-8, 1.625e-8, 5.75e-9, 5.75e-9, 1.88e-8]
    @test collect(deposition_particles_3) ≈ deposition_particles_3_expected

    # Index 4: deposition in the river from transport capacity exceeding
    excess_sediment =
        max(sum(input_particles[4]) - model.boundary_conditions.transport_capacity[4], 0.0)
    @test excess_sediment ≈ 1.125e-8
    deposition_particles_4 = Wflow.compute_transport_capacity_deposition(
        excess_sediment,
        input_particles[4][1:6],
        erosion_particles[4],
    )
    deposition_particles_4_expected = [0.0, 0.0, 5.5e-9, 0.0, 5.75e-9, 0.0]
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
        order,
    )

    perform_tests(model.variables)

    # Perform calculations again by directly calling update!
    model, domain, _, _, _ = get_objects()
    Wflow.update!(model, domain, dt)
    perform_tests(model.variables)
end
