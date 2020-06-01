@testset "vertical processes" begin
    @test all(isapprox.(
        Wflow.rainfall_interception_gash(3.0, 0.11, 0.24, 18.0, 1.5, maxevap = 4.0),
        (13.568000000000001, 4.0, 0.432, 1.5),
    ))
    @test all(isapprox.(
        Wflow.acttransp_unsat_sbm(
            300.0,
            55.0,
            400.0,
            3.6,
            1.85,
            10.5,
            300.0,
            0.6,
            0.1,
            10.0,
            false,
        ),
        (55.0, 1.85, 3.6),
    ))
    @test all(isapprox.(
        Wflow.infiltration(27.5, 0.2, 0.038, 8.9, 50.0, 5.0, 23.5, false, false),
        (23.5, 19.14814814814815, 4.351851851851852, 22.0, 5.5, 0.5),
    ))
    @test all(isapprox.(
        Wflow.unsatzone_flow_layer(43.5, 256.0, 135.0, 12.6),
        (43.49983744545384, 0.00016255454615829025),
    ))
    @test all(isapprox.(
        Wflow.unsatzone_flow_sbm(67.0, 950.0, 600.0, 250.0, 200.0, 0.6, 0.1),
        (19.142857142857146, 47.857142857142854),
    ))
    @test all(isapprox.(
        Wflow.snowpack_hbv(201.5, 15.0, 30.1, 0.54, 2.0, 0.0, 0.0, 2.5, 0.10),
        (207.073, 20.707300000000004, 1.35, 18.819699999999997, 6.923),
    ))
    @test Wflow.scurve(2.0, a = 0.0, b = 3.0, c = 2.5) ≈ 0.3325863502664285
end
