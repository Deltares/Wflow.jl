using Dates
@testset "vertical processes" begin
    @test all(
        isapprox.(
            Wflow.rainfall_interception_gash(3.0, 0.11, 0.24, 18.0, 1.5, 4.0),
            (13.568000000000001, 4.0, 0.432, 1.5),
        ),
    )
    @test all(
        isapprox.(
            Wflow.rainfall_interception_modrut(8.6, 3.8, 1.5, 0.45, 2.8),
            (4.343, 3.87, 0.387, 0.0, 3.8, 2.043),
        ),
    )
    head = Wflow.head_brooks_corey(0.25, 0.6, 0.15, 10.5, -10.0)
    @test head ≈ -90.6299820833844
    h3 = Wflow.feddes_h3(-300.0, -600.0, 3.5, Second(86400))
    @test h3 ≈ -412.5
    alpha = Wflow.rwu_reduction_feddes(head, -10.0, -100.0, h3, -15000.0, 0.0)
    @test alpha ≈ 0.8958886898153823
    @test all(
        isapprox.(
            Wflow.infiltration(27.5, 0.2, 0.038, 8.9, 50.0, 5.0, 23.5, false, false),
            (23.5, 19.14814814814815, 4.351851851851852, 22.0, 5.5, 0.5, 1.0),
        ),
    )
    @test all(
        isapprox.(
            Wflow.unsatzone_flow_layer(43.5, 256.0, 135.0, 12.6),
            (43.49983744545384, 0.00016255454615829025),
        ),
    )
    @test all(
        isapprox.(
            Wflow.unsatzone_flow_sbm(67.0, 950.0, 600.0, 250.0, 200.0, 0.6, 0.1),
            (19.142857142857146, 47.857142857142854),
        ),
    )
    @test all(
        isapprox.(
            Wflow.snowpack_hbv(201.5, 15.0, 30.1, 0.54, 2.0, 0.0, 0.0, 2.5, 0.10),
            (207.073, 20.707300000000004, 1.35, 18.819699999999997, 6.923),
        ),
    )
    @test Wflow.scurve(2.0, 0.0, 3.0, 2.5) ≈ 0.3325863502664285
    @test all(
        isapprox.(
            Wflow.glacier_hbv(0.35, 500.0, 9.5, 5.0, 0.0, 3.4, 0.2, Second(Day(1))),
            (8.835, 1.9, 484.9, 17.0),
        ),
    )
end
