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
    # Case canopy_gap_fraction < inv(1.1), potential_evaporation < canopy_storage (after precipitation)
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

    # Case canopy_gap_fraction > inv(1.1), potential_evaporation > canopy_storage
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

@testitem "unit: precipitation_hbv" begin
    ## Case tti > 0.0
    precipitation = 30.1
    temperature = 0.54
    tti = 2.0
    tt = 0.0
    snow_precip, liquid_precip =
        Wflow.precipitation_hbv(precipitation, temperature, tti, tt)
    @test snow_precip ≈ 6.923
    @test liquid_precip ≈ 23.177

    ## Case tti == 0
    # Case temperature > tt
    tti = 0.0
    snow_precip, liquid_precip =
        Wflow.precipitation_hbv(precipitation, temperature, tti, tt)
    @test snow_precip == 0.0
    @test liquid_precip == precipitation

    # Case temperate < tt
    temperature = -1.0
    snow_precip, liquid_precip =
        Wflow.precipitation_hbv(precipitation, temperature, tti, tt)
    @test snow_precip == precipitation
    @test liquid_precip == 0.0
end

@testitem "unit: snowpack_hbv" begin
    snow_storage = 201.5
    snow_water = 15.0
    snow_precip = 6.923
    liquid_precip = 23.177
    temperature = 0.54
    ttm = 0.0
    cfmax = 2.5
    whc = 0.10
    # Case temperature > ttm
    snow_storage_new, snow_water, snow_water_equivalent, snow_melt, runoff =
        Wflow.snowpack_hbv(
            snow_storage,
            snow_water,
            snow_precip,
            liquid_precip,
            temperature,
            ttm,
            cfmax,
            whc,
        )
    @test snow_storage_new ≈ 207.073
    @test snow_water ≈ 20.7073
    @test snow_water_equivalent ≈ 227.7803
    @test snow_melt ≈ 1.35
    @test runoff ≈ 18.8197

    # Case temperature < ttm
    temperature = -0.5
    snow_storage_new, snow_water, snow_water_equivalent, snow_melt, runoff =
        Wflow.snowpack_hbv(
            snow_storage,
            snow_water,
            snow_precip,
            liquid_precip,
            temperature,
            ttm,
            cfmax,
            whc,
        )
    @test snow_storage_new ≈ 208.4855
    @test snow_water ≈ 20.84855
    @test snow_water_equivalent ≈ 229.33405
    @test snow_melt ≈ 0.0
    @test runoff ≈ 22.97325
end

@testitem "unit: glacier_hbv" begin
    glacier_frac = 0.35
    glacier_store = 500.0
    snow_storage = 9.5
    temperature = 5.0
    ttm = 0.0
    cfmax = 3.4
    g_sifrac = 0.2
    max_snow_to_glacier = 8.0
    snow_storage, snow_to_glacier, glacier_storage, glacier_melt = Wflow.glacier_hbv(
        glacier_frac,
        glacier_store,
        snow_storage,
        temperature,
        ttm,
        cfmax,
        g_sifrac,
        max_snow_to_glacier,
    )
    @test snow_storage ≈ 8.835
    @test snow_to_glacier ≈ 1.9
    @test glacier_storage ≈ 484.9
    @test glacier_melt ≈ 17.0
end

@testitem "unit: infiltration" begin
    potential_infiltration = 27.5
    pathfrac = 0.2
    infiltcapsoil = 50.0
    infiltcappath = 5.0
    ustorecapacity = 23.5
    f_infilt_reduction = 1.0

    infiltsoilpath, infiltexcess = Wflow.infiltration(
        potential_infiltration,
        pathfrac,
        infiltcapsoil,
        infiltcappath,
        ustorecapacity,
        f_infilt_reduction,
    )
    @test infiltsoilpath == ustorecapacity
    @test infiltexcess ≈ 0.5
end

@testitem "unit: unsatzone_flow_layer" begin
    kv_z = 256.0
    l_sat = 135.0
    c = 12.6

    # Case usd > 0
    usd = 43.5
    usd_new, sum_ast = Wflow.unsatzone_flow_layer(usd, kv_z, l_sat, c)
    usd_new = 43.49983744545384
    sum_ast = 0.00016255454615829025

    # Case usd == 0
    usd = 0
    usd_new, sum_ast = Wflow.unsatzone_flow_layer(usd, kv_z, l_sat, c)
    @test usd_new == 0.0
    @test sum_ast == 0.0
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

@testitem "unit: Feddes root water uptake" begin
    h3_high = -300.0
    h3_low = -600.0
    dt = 86400.0

    # Case tpot_daily < 1.0
    tpot = 0.5
    @test Wflow.feddes_h3(h3_high, h3_low, tpot, dt) == h3_low

    # Case 1.0 < tpot_daily < 5.0
    tpot = 3.0
    @test Wflow.feddes_h3(h3_high, h3_low, tpot, dt) ≈ (h3_high + h3_low) / 2

    # Case tpot_daily > 5.0
    tpot = 7.5
    @test Wflow.feddes_h3(h3_high, h3_low, tpot, dt) ≈ h3_high

    h1 = -10.0
    h2 = -100.0
    h3 = -300.0
    h4 = -15000.0

    ## Case alpha == 0.0
    alpha = 0.0

    # Case h < h4
    h = -16000.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 0.0

    # Case h3 < h < h4
    h = -1000.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.4 / 1.47

    # Case h2 < h < h3
    h = -150.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0

    # Case h1 < h < h2
    h = -50.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 4 / 9

    # Case h > h1
    h = -5.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 0.0

    ## Case alpha ≠ 0.0
    alpha = 0.5

    # Case h < h4
    h = -16000.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 0.0

    # Case h3 < h < h4
    h = -1000.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.4 / 1.47

    # Case h2 < h < h3
    h = -150.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0

    # Case h1 < h < h2
    h = -50.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0

    # Case h > h1
    h = -5.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0
end

@testitem "unit: soil_temperature" begin
    tsoil_prev = 1.0
    w_soil = 2.0
    temperature = 1.5
    @test Wflow.soil_temperature(tsoil_prev, w_soil, temperature) ≈ 2.0
end

@testitem "unit: infiltration_reduction_factor" begin
    tsoil = 0.1
    cf_soil = 0.3

    # Case model_snow && soil_infiltration_reduction
    modelsnow = true
    soil_infiltration_reduction = true
    @test Wflow.infiltration_reduction_factor(
        tsoil,
        cf_soil;
        modelsnow,
        soil_infiltration_reduction,
    ) ≈ 0.8325096069489968

    # Case !(model_snow && soil_infiltration_reduction)
    soil_infiltration_reduction = false
    @test Wflow.infiltration_reduction_factor(
        tsoil,
        cf_soil;
        modelsnow,
        soil_infiltration_reduction,
    ) == 1.0
end

@testitem "unit: soil_evaporation_unsaturated_store" begin
    potential_soilevaporation = 0.302
    ustorelayerdepth = 1.23
    ustorelayerthickness = 100.0
    zi = 300.0
    theta_effective = 0.241

    # Case n_unsatlayers == 0
    n_unsatlayers = 0
    @test Wflow.soil_evaporation_unsatured_store(
        potential_soilevaporation,
        ustorelayerdepth,
        ustorelayerthickness,
        n_unsatlayers,
        zi,
        theta_effective,
    ) == 0.0

    # case n_unsatlayers == 1
    n_unsatlayers = 1
    @test Wflow.soil_evaporation_unsatured_store(
        potential_soilevaporation,
        ustorelayerdepth,
        ustorelayerthickness,
        n_unsatlayers,
        zi,
        theta_effective,
    ) ≈ 0.005137759336099585

    # Case n_unsatlayers > 1
    n_unsatlayers = 2
    @test Wflow.soil_evaporation_unsatured_store(
        potential_soilevaporation,
        ustorelayerdepth,
        ustorelayerthickness,
        n_unsatlayers,
        zi,
        theta_effective,
    ) ≈ 0.015413278008298757
end

@testitem "unit: soil_evaporation_saturated_store" begin
    potential_soilevaporation = 0.125
    layerthickness = 100.0
    zi = 300.0
    theta_effective = 0.32205961644649506

    # Case n_unsatlayers ∈ (0, 1)
    n_unsatlayers = 0
    @test Wflow.soil_evaporation_satured_store(
        potential_soilevaporation,
        n_unsatlayers,
        layerthickness,
        zi,
        theta_effective,
    ) ≈ -64.41192328929901

    # Case n_unsatlayers ∉ (0, 1)
    n_unsatlayers = 2
    @test Wflow.soil_evaporation_satured_store(
        potential_soilevaporation,
        n_unsatlayers,
        layerthickness,
        zi,
        theta_effective,
    ) == 0.0
end

@testitem "unit: actual_infiltration_soil_path" begin
    potential_infiltration = 1.627
    actinfilt = 1.627
    pathfrac = 0.1
    infiltcapsoil = 228.596
    infiltcappath = 5.0
    f_infiltration_reduction = 0.9

    # Case actinfilt > 0
    actinfilt = 1.627
    actinfiltsoil, actinfiltpath = Wflow.actual_infiltration_soil_path(
        potential_infiltration,
        actinfilt,
        pathfrac,
        infiltcapsoil,
        infiltcappath,
        f_infiltration_reduction,
    )
    @test actinfiltsoil ≈ 1.4643
    @test actinfiltpath ≈ 0.1627

    # Case actinfilt == 0
    actinfilt = 0
    actinfiltsoil, actinfiltpath = Wflow.actual_infiltration_soil_path(
        potential_infiltration,
        actinfilt,
        pathfrac,
        infiltcapsoil,
        infiltcappath,
        f_infiltration_reduction,
    )
    @test actinfiltsoil == 0.0
    @test actinfiltpath == 0.0
end
