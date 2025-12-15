@testitem "unit: rainfall_intercepiton_gash" begin
    # Case cmax == 0
    cmax = 0
    e_r = 0.11
    canopy_gap_fraction = 0.24
    precipitation = 18.0
    canopy_storage_in = 1.5
    max_evaporation = 4.0
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            cmax,
            e_r,
            canopy_gap_fraction,
            precipitation,
            canopy_storage_in,
            max_evaporation,
        )
    @test throughfall == precipitation
    @test interception == 0.0
    @test stem_flow == 0.0
    @test canopy_storage_in == canopy_storage_out

    # Case cmax > 0, large_storms == true, interception > max_evaporation
    cmax = 3.0
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            cmax,
            e_r,
            canopy_gap_fraction,
            precipitation,
            canopy_storage_in,
            max_evaporation,
        )
    @test throughfall ≈ 13.568
    @test interception ≈ 4.0
    @test stem_flow ≈ 0.432
    @test canopy_storage_in == canopy_storage_out

    # Case cmax > 0, large_storms == false, interception > max_evaporation
    precipitation = 1.0
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            cmax,
            e_r,
            canopy_gap_fraction,
            precipitation,
            canopy_storage_in,
            max_evaporation,
        )
    @test throughfall ≈ 0.24
    @test interception ≈ 0.736
    @test stem_flow ≈ 0.024
    @test canopy_storage_in == canopy_storage_out
end

@testitem "unit: rainfall_interception_modrut (modified Rutter)" begin
    # Case: canopy_gap_fraction < inv(1.1), potential_evaporation < canopy_storage (after precipitation)
    precipitation = 8.6
    potential_evaporation = 3.8
    canopy_storage = 1.5
    canopy_gap_fraction = 0.45
    cmax = 2.8
    throughfall, canopy_evap, stemflow, canopy_storage = Wflow.rainfall_interception_modrut(
        precipitation,
        potential_evaporation,
        canopy_storage,
        canopy_gap_fraction,
        cmax,
    )
    @test throughfall ≈ 3.87
    @test canopy_evap ≈ 3.8
    @test stemflow ≈ 0.387
    @test canopy_storage ≈ 2.043

    # Case: canopy_gap_fraction > inv(1.1), potential_evaporation > canopy_storage
    precipitation = 1.0
    canopy_gap_fraction = 0.95
    throughfall, canopy_evap, stemflow, canopy_storage = Wflow.rainfall_interception_modrut(
        precipitation,
        potential_evaporation,
        canopy_storage,
        canopy_gap_fraction,
        cmax,
    )
    @test throughfall ≈ 0.95
    @test canopy_evap ≈ 2.043
    @test stemflow ≈ 0.05
    @test canopy_storage ≈ 0.0
end

@testitem "unit: Brooks-Corey soil hydraulic model" begin
    # Case par_lambda > 0
    vwc = 0.25
    theta_s = 0.6
    theta_r = 0.15
    c = 10.5
    hb = -10.0
    h = Wflow.head_brooks_corey(vwc, theta_s, theta_r, c, hb)
    @test h ≈ -90.6299820833844
    @test Wflow.vwc_brooks_corey(h, hb, theta_s, theta_r, c) ≈ vwc + theta_r

    # Case par_lambda < 0
    c = 2.0
    h = Wflow.head_brooks_corey(vwc, theta_s, theta_r, c, hb)
    @test h == hb
    @test Wflow.vwc_brooks_corey(h, hb, theta_s, theta_r, c) ≈ theta_s
end

@testitem "unit: other" begin
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
