@testitem "vertical processes" begin
    using Dates
    @test all(
        isapprox.(
            Wflow.rainfall_interception_gash(3.0, 0.11, 0.24, 18.0, 1.5, 4.0),
            (13.568000000000001, 4.0, 0.432, 1.5),
        ),
    )
    @test all(
        isapprox.(
            Wflow.rainfall_interception_modrut(8.6, 3.8, 1.5, 0.45, 2.8),
            (3.87, 3.8, 0.387, 2.043),
        ),
    )
    @test Wflow.head_brooks_corey(0.25, 0.6, 0.15, 10.5, -10.0) ≈ -90.6299820833844
    @test Wflow.feddes_h3(-300.0, -600.0, 3.5, 86400.0) ≈ -412.5
    @test Wflow.feddes_h3(-300.0, -600.0, 0.5, 86400.0) == -600.0
    @test Wflow.feddes_h3(-300.0, -600.0, 6.0, 86400.0) == -300.0
    @test Wflow.rwu_reduction_feddes(0.0, -10.0, -100.0, -300.0, -15000.0, 0.0) == 0.0
    @test Wflow.rwu_reduction_feddes(0.0, -10.0, -100.0, -300.0, -15000.0, 1.0) == 1.0
    @test Wflow.rwu_reduction_feddes(-90.0, -10.0, -100.0, -412.5, -15000.0, 0.0) ≈
          0.8888888888888888
    @test Wflow.rwu_reduction_feddes(-350.0, -10.0, -100.0, -412.5, -15000.0, 0.0) == 1.0
    @test Wflow.rwu_reduction_feddes(-12000.0, -10.0, -100.0, -412.5, -15000.0, 0.0) ≈
          0.20565552699228792
    @test Wflow.rwu_reduction_feddes(-16000.0, -10.0, -100.0, -412.5, -15000.0, 0.0) == 0.0
    @test all(isapprox.(Wflow.infiltration(27.5, 0.2, 50.0, 5.0, 23.5, 1.0), (23.5, 0.5)))
    @test all(
        isapprox.(
            Wflow.unsatzone_flow_layer(43.5, 256.0, 135.0, 12.6),
            (43.49983744545384, 0.00016255454615829025),
        ),
    )
    @test all(
        isapprox.(
            Wflow.precipitation_hbv(30.1, 0.54, 2.0, 0.0),
            (6.923, 23.177000000000003),
        ),
    )
    @test all(
        isapprox.(
            Wflow.snowpack_hbv(201.5, 15.0, 6.923, 23.177, 0.54, 0.0, 2.5, 0.10),
            (207.073, 20.707300000000004, 227.7803, 1.35, 18.819699999999997),
        ),
    )
    @test Wflow.scurve(2.0, 0.0, 3.0, 2.5) ≈ 0.3325863502664285
    @test all(
        isapprox.(
            Wflow.glacier_hbv(0.35, 500.0, 9.5, 5.0, 0.0, 3.4, 0.2, 8.0),
            (8.835, 1.9, 484.9, 17.0),
        ),
    )
end
