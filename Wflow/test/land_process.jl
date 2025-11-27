@testitem "vertical processes" begin
    using Dates
    using Wflow: Unit, to_SI, MM_PER_DAY, MM, ABSOLUTE_DEGREES
    DEGREES = Unit(; degC = 1)
    MM_PER_DT = Unit(; mm = 1, dt = -1)
    MM_PER_DEGREE_PER_DT = Unit(; mm = 1, degC = -1, dt = -1)
    CM = Unit(; cm = 1)
    dt = 86400.0

    @test all(
        isapprox.(
            Wflow.rainfall_interception_gash(
                to_SI(3.0, MM),
                0.11,
                0.24,
                to_SI(18.0, MM_PER_DT; dt_val = dt),
                1.5,
                to_SI(4.0, MM_PER_DT; dt_val = dt),
                dt,
            ),
            (
                to_SI(13.568, MM_PER_DT; dt_val = dt),
                to_SI(4.0, MM_PER_DT; dt_val = dt),
                to_SI(0.432, MM_PER_DT; dt_val = dt),
                1.5,
            ),
        ),
    )
    @test all(
        isapprox.(
            Wflow.rainfall_interception_modrut(
                to_SI(8.6, MM_PER_DT; dt_val = dt),
                to_SI(3.8, MM_PER_DT; dt_val = dt),
                to_SI(1.5, MM; dt_val = dt),
                0.45,
                to_SI(2.8, MM; dt_val = dt),
                dt,
            ),
            (
                to_SI(3.87, MM_PER_DT; dt_val = dt),
                to_SI(3.8, MM_PER_DT; dt_val = dt),
                to_SI(0.387, MM_PER_DT; dt_val = dt),
                to_SI(2.043, MM),
            ),
        ),
    )
    @test Wflow.head_brooks_corey(0.25, 0.6, 0.15, 10.5, to_SI(-10.0, CM)) ≈
          to_SI(-90.6299820833844, CM)
    h3_high_low = (to_SI(-300.0, MM), to_SI(-600.0, MM))
    @test Wflow.feddes_h3(h3_high_low..., to_SI(3.5, MM_PER_DAY)) ≈ to_SI(-412.5, MM)
    @test Wflow.feddes_h3(h3_high_low..., to_SI(0.5, MM_PER_DAY)) == to_SI(-600.0, MM)
    @test Wflow.feddes_h3(h3_high_low..., to_SI(6.0, MM_PER_DAY)) == to_SI(-300.0, MM)
    hs = (
        to_SI(0.0, CM),
        to_SI(-10.0, CM),
        to_SI(-100.0, CM),
        to_SI(-300.0, CM),
        to_SI(-15000.0, CM),
    )
    @test Wflow.rwu_reduction_feddes(hs..., 0.0) == 0.0
    @test Wflow.rwu_reduction_feddes(hs..., 1.0) == 1.0
    @test Wflow.rwu_reduction_feddes(
        to_SI(-90.0, CM),
        to_SI(-10.0, CM),
        to_SI(-100.0, CM),
        to_SI(-412.5, CM),
        to_SI(-15000.0, CM),
        0.0,
    ) ≈ 0.8888888888888888
    @test Wflow.rwu_reduction_feddes(
        to_SI(-350.0, CM),
        to_SI(-10.0, CM),
        to_SI(-100.0, CM),
        to_SI(-412.5, CM),
        to_SI(-15000.0, CM),
        0.0,
    ) == 1.0
    @test Wflow.rwu_reduction_feddes(
        to_SI(-12000.0, CM),
        to_SI(-10.0, CM),
        to_SI(-100.0, CM),
        to_SI(-412.5, CM),
        to_SI(-15000.0, CM),
        0.0,
    ) ≈ 0.20565552699228792
    @test Wflow.rwu_reduction_feddes(
        to_SI(-16000.0, CM),
        to_SI(-10.0, CM),
        to_SI(-100.0, CM),
        to_SI(-412.5, CM),
        to_SI(-15000.0, CM),
        0.0,
    ) == 0.0
    @test all(
        isapprox.(
            Wflow.infiltration(
                to_SI(27.5, MM_PER_DT; dt_val = dt),
                0.2,
                to_SI(50.0, MM_PER_DT; dt_val = dt),
                to_SI(5.0, MM_PER_DT; dt_val = dt),
                to_SI(23.5, MM; dt_val = dt),
                1.0,
                dt,
            ),
            (to_SI(23.5, MM_PER_DT; dt_val = dt), to_SI(0.5, MM_PER_DT; dt_val = dt)),
        ),
    )
    @test all(
        isapprox.(
            Wflow.unsatzone_flow_layer(
                to_SI(43.5, MM),
                to_SI(256.0, MM_PER_DT; dt_val = dt),
                to_SI(135.0, MM),
                12.6,
                dt,
            ),
            (
                to_SI(43.49983744545384, MM),
                to_SI(0.00016255454615829025, MM_PER_DT; dt_val = dt),
            ),
        ),
    )
    @test all(
        isapprox.(
            Wflow.unsatzone_flow_sbm(
                to_SI(67.0, MM),
                to_SI(950.0, MM),
                to_SI(600.0, MM),
                to_SI(250.0, MM_PER_DT; dt_val = dt),
                to_SI(200.0, MM),
                0.6,
                0.1,
                dt,
            ),
            (
                to_SI(19.142857142857146, MM),
                to_SI(47.857142857142854, MM_PER_DT; dt_val = dt),
            ),
        ),
    )
    @test all(
        isapprox.(
            Wflow.precipitation_hbv(
                to_SI(30.1, MM_PER_DT; dt_val = dt),
                to_SI(0.54, ABSOLUTE_DEGREES),
                to_SI(2.0, DEGREES),
                to_SI(0.0, ABSOLUTE_DEGREES),
            ),
            (to_SI(6.923, MM_PER_DT; dt_val = dt), to_SI(23.177, MM_PER_DT; dt_val = dt)),
        ),
    )
    @test all(
        isapprox.(
            Wflow.snowpack_hbv(
                to_SI(201.5, MM),
                to_SI(15.0, MM),
                to_SI(6.923, MM_PER_DT; dt_val = dt),
                to_SI(23.177, MM_PER_DT; dt_val = dt),
                to_SI(0.54, ABSOLUTE_DEGREES),
                to_SI(0.0, ABSOLUTE_DEGREES),
                to_SI(2.5, MM_PER_DEGREE_PER_DT; dt_val = dt),
                0.10,
                dt,
            ),
            (
                to_SI(207.073, MM),
                to_SI(20.7073, MM),
                to_SI(227.7803, MM),
                to_SI(1.35, MM_PER_DT; dt_val = dt),
                to_SI(18.8197, MM_PER_DT; dt_val = dt),
            ),
        ),
    )
    @test Wflow.scurve(2.0, 0.0, 3.0, 2.5) ≈ 0.3325863502664285
    @test all(
        isapprox.(
            Wflow.glacier_hbv(
                0.35,
                to_SI(500.0, MM),
                to_SI(9.5, MM),
                to_SI(5.0, ABSOLUTE_DEGREES),
                to_SI(0.0, ABSOLUTE_DEGREES),
                to_SI(3.4, MM_PER_DEGREE_PER_DT; dt_val = dt),
                0.2,
                to_SI(8.0, MM_PER_DT; dt_val = dt),
                dt,
            ),
            (
                to_SI(8.835, MM),
                to_SI(1.9, MM_PER_DT; dt_val = dt),
                to_SI(484.9, MM),
                to_SI(17.0, MM_PER_DT; dt_val = dt),
            ),
        ),
    )
end
