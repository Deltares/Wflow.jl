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
