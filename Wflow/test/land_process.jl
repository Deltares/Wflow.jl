@testitem "unit: rainfall_intercepiton_gash" begin
    # Case maximum_canopy_storage == 0
    maximum_canopy_storage = 0
    evaporation_to_precipitation_ratio = 0.11
    canopy_gap_fraction = 0.24
    precipitation = 18.0
    canopy_storage_in = 1.5
    max_evaporation = 4.0
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            maximum_canopy_storage,
            evaporation_to_precipitation_ratio,
            canopy_gap_fraction,
            precipitation,
            canopy_storage_in,
            max_evaporation,
        )
    @test throughfall == precipitation
    @test interception == 0.0
    @test stem_flow == 0.0
    @test canopy_storage_in == canopy_storage_out

    # Case maximum_canopy_storage > 0, large_storms == true, interception > max_evaporation
    maximum_canopy_storage = 3.0
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            maximum_canopy_storage,
            evaporation_to_precipitation_ratio,
            canopy_gap_fraction,
            precipitation,
            canopy_storage_in,
            max_evaporation,
        )
    @test throughfall ≈ 13.568
    @test interception ≈ 4.0
    @test stem_flow ≈ 0.432
    @test canopy_storage_in == canopy_storage_out

    # Case maximum_canopy_storage > 0, large_storms == false, interception > max_evaporation
    precipitation = 1.0
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            maximum_canopy_storage,
            evaporation_to_precipitation_ratio,
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
    maximum_canopy_storage = 2.8
    throughfall, canopy_evap, stemflow, canopy_storage = Wflow.rainfall_interception_modrut(
        precipitation,
        potential_evaporation,
        canopy_storage,
        canopy_gap_fraction,
        maximum_canopy_storage,
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
        maximum_canopy_storage,
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
    glacier_fraction = 0.35
    glacier_store = 500.0
    snow_storage = 9.5
    temperature = 5.0
    ttm = 0.0
    cfmax = 3.4
    snow_to_ice_fraction = 0.2
    maximum_snow_to_ice_rate = 8.0
    snow_storage, snow_to_glacier, glacier_storage, glacier_melt = Wflow.glacier_hbv(
        glacier_fraction,
        glacier_store,
        snow_storage,
        temperature,
        ttm,
        cfmax,
        snow_to_ice_fraction,
        maximum_snow_to_ice_rate,
    )
    @test snow_storage ≈ 8.835
    @test snow_to_glacier ≈ 1.9
    @test glacier_storage ≈ 484.9
    @test glacier_melt ≈ 17.0
end

@testitem "unit: infiltration" begin
    potential_infiltration = 27.5
    compacted_soil_area_fraction = 0.2
    infiltration_capacity_soil = 50.0
    infiltration_capacity_compacted_soil = 5.0
    unsaturated_store_capacity = 23.5
    f_infilt_reduction = 1.0

    infiltration, infiltration_excess = Wflow.infiltration(
        potential_infiltration,
        compacted_soil_area_fraction,
        infiltration_capacity_soil,
        infiltration_capacity_compacted_soil,
        unsaturated_store_capacity,
        f_infilt_reduction,
    )
    @test infiltration == unsaturated_store_capacity
    @test infiltration_excess ≈ 0.5
end

@testitem "unit: unsatzone_flow_layer" begin
    kv_z = 256.0
    l_sat = 135.0
    brooks_corey_exponent = 12.6

    # Case usd > 0
    usd = 43.5
    usd_new, sum_ast =
        Wflow.unsatzone_flow_layer(usd, kv_z, l_sat, brooks_corey_exponent)
    usd_new = 43.49983744545384
    sum_ast = 0.00016255454615829025

    # Case usd == 0
    usd = 0
    usd_new, sum_ast =
        Wflow.unsatzone_flow_layer(usd, kv_z, l_sat, brooks_corey_exponent)
    @test usd_new == 0.0
    @test sum_ast == 0.0
end

@testitem "unit: Brooks-Corey soil hydraulic model" begin
    # Case par_lambda > 0
    volumetric_water_content = 0.25
    theta_s = 0.6
    theta_r = 0.15
    brooks_corey_exponent = 10.5
    air_entry_pressure = -10.0
    h = Wflow.head_brooks_corey(
        volumetric_water_content,
        theta_s,
        theta_r,
        brooks_corey_exponent,
        air_entry_pressure,
    )
    @test h ≈ -90.6299820833844
    @test Wflow.vwc_brooks_corey(
        h,
        air_entry_pressure,
        theta_s,
        theta_r,
        brooks_corey_exponent,
    ) ≈
          volumetric_water_content + theta_r

    # Case par_lambda < 0
    brooks_corey_exponent = 2.0
    h = Wflow.head_brooks_corey(
        volumetric_water_content,
        theta_s,
        theta_r,
        brooks_corey_exponent,
        air_entry_pressure,
    )
    @test h == air_entry_pressure
    @test Wflow.vwc_brooks_corey(
        h,
        air_entry_pressure,
        theta_s,
        theta_r,
        brooks_corey_exponent,
    ) ≈ theta_s
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
    soil_surface_temperature = 0.1
    cf_soil = 0.3

    # Case model_snow && soil_infiltration_reduction
    modelsnow = true
    soil_infiltration_reduction = true
    @test Wflow.infiltration_reduction_factor(
        soil_surface_temperature,
        cf_soil;
        modelsnow,
        soil_infiltration_reduction,
    ) ≈ 0.8325096069489968

    # Case !(model_snow && soil_infiltration_reduction)
    soil_infiltration_reduction = false
    @test Wflow.infiltration_reduction_factor(
        soil_surface_temperature,
        cf_soil;
        modelsnow,
        soil_infiltration_reduction,
    ) == 1.0
end

@testitem "unit: soil_evaporation_unsaturated_store" begin
    potential_soilevaporation = 0.302
    unsaturated_layer_depth = 1.23
    unsaturated_layer_thickness = 100.0
    water_table_depth = 300.0
    theta_effective = 0.241

    # Case n_unsatlayers == 0
    n_unsatlayers = 0
    @test Wflow.soil_evaporation_unsatured_store(
        potential_soilevaporation,
        unsaturated_layer_depth,
        unsaturated_layer_thickness,
        n_unsatlayers,
        water_table_depth,
        theta_effective,
    ) == 0.0

    # case n_unsatlayers == 1
    n_unsatlayers = 1
    @test Wflow.soil_evaporation_unsatured_store(
        potential_soilevaporation,
        unsaturated_layer_depth,
        unsaturated_layer_thickness,
        n_unsatlayers,
        water_table_depth,
        theta_effective,
    ) ≈ 0.005137759336099585

    # Case n_unsatlayers > 1
    n_unsatlayers = 2
    @test Wflow.soil_evaporation_unsatured_store(
        potential_soilevaporation,
        unsaturated_layer_depth,
        unsaturated_layer_thickness,
        n_unsatlayers,
        water_table_depth,
        theta_effective,
    ) ≈ 0.015413278008298757
end

@testitem "unit: soil_evaporation_saturated_zone_store" begin
    potential_soilevaporation = 0.125
    layerthickness = 100.0
    water_table_depth = 300.0
    theta_effective = 0.32205961644649506

    # Case n_unsatlayers ∈ (0, 1)
    n_unsatlayers = 0
    @test Wflow.soil_evaporation_satured_store(
        potential_soilevaporation,
        n_unsatlayers,
        layerthickness,
        water_table_depth,
        theta_effective,
    ) ≈ -64.41192328929901

    # Case n_unsatlayers ∉ (0, 1)
    n_unsatlayers = 2
    @test Wflow.soil_evaporation_satured_store(
        potential_soilevaporation,
        n_unsatlayers,
        layerthickness,
        water_table_depth,
        theta_effective,
    ) == 0.0
end

@testitem "unit: actual_infiltration_soil_path" begin
    potential_infiltration = 1.627
    actual_infiltration = 1.627
    compacted_soil_area_fraction = 0.1
    infiltration_capacity_soil = 228.596
    infiltration_capacity_compacted_soil = 5.0
    f_infiltration_reduction = 0.9

    # Case actual_infiltration > 0
    actual_infiltration = 1.627
    actual_infiltration_soil, actual_infiltration_compacted_soil =
        Wflow.actual_infiltration_soil_path(
            potential_infiltration,
            actual_infiltration,
            compacted_soil_area_fraction,
            infiltration_capacity_soil,
            infiltration_capacity_compacted_soil,
            f_infiltration_reduction,
        )
    @test actual_infiltration_soil ≈ 1.4643
    @test actual_infiltration_compacted_soil ≈ 0.1627

    # Case actual_infiltration == 0
    actual_infiltration = 0
    actual_infiltration_soil, actual_infiltration_compacted_soil =
        Wflow.actual_infiltration_soil_path(
            potential_infiltration,
            actual_infiltration,
            compacted_soil_area_fraction,
            infiltration_capacity_soil,
            infiltration_capacity_compacted_soil,
            f_infiltration_reduction,
        )
    @test actual_infiltration_soil == 0.0
    @test actual_infiltration_compacted_soil == 0.0
end
