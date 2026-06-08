@testitem "unit: rainfall_intercepiton_gash" begin
    include("testing_utils.jl")
    dt = 86400.0

    # Case maximum_canopy_storage == 0
    maximum_canopy_storage = 0
    evaporation_to_precipitation_ratio = 0.11
    canopy_gap_fraction = 0.24
    precipitation = 2.0833333333333333e-7
    canopy_storage_in = 0.0015
    max_evaporation = 4.6296296296296295e-8
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            maximum_canopy_storage,
            evaporation_to_precipitation_ratio,
            canopy_gap_fraction,
            precipitation,
            canopy_storage_in,
            max_evaporation,
            dt,
        )
    @test throughfall == precipitation
    @test interception == 0.0
    @test stem_flow == 0.0
    @test canopy_storage_in == canopy_storage_out

    # Case maximum_canopy_storage > 0, large_storms == true, interception > max_evaporation
    maximum_canopy_storage = 0.003
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            maximum_canopy_storage,
            evaporation_to_precipitation_ratio,
            canopy_gap_fraction,
            precipitation,
            canopy_storage_in,
            max_evaporation,
            dt,
        )
    @test throughfall ≈ 1.5703703703703703e-7
    @test interception ≈ 4.6296296296296295e-8
    @test stem_flow ≈ 5.0e-9
    @test canopy_storage_in == canopy_storage_out

    # Case maximum_canopy_storage > 0, large_storms == false, interception > max_evaporation
    precipitation = 1.1574074074074074e-8
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            maximum_canopy_storage,
            evaporation_to_precipitation_ratio,
            canopy_gap_fraction,
            precipitation,
            canopy_storage_in,
            max_evaporation,
            dt,
        )
    @test throughfall ≈ 2.7777777777777776e-9
    @test interception ≈ 8.518518518518518e-9
    @test stem_flow ≈ 2.7777777777777777e-10
    @test canopy_storage_in == canopy_storage_out
end

@testitem "unit: rainfall_interception_modrut (modified Rutter)" begin
    dt = 86400.0

    # Case canopy_gap_fraction < inv(1.1), potential_evaporation < canopy_storage (after precipitation)
    precipitation = 9.953703703703703e-8
    potential_evaporation = 4.398148148148148e-8
    canopy_storage = 0.0015
    canopy_gap_fraction = 0.45
    maximum_canopy_storage = 0.0028
    throughfall, canopy_evap, stemflow, canopy_storage = Wflow.rainfall_interception_modrut(
        precipitation,
        potential_evaporation,
        canopy_storage,
        canopy_gap_fraction,
        maximum_canopy_storage,
        dt,
    )
    @test throughfall ≈ 4.4791666666666666e-8
    @test canopy_evap ≈ 4.398148148148148e-8
    @test stemflow ≈ 4.479166666666667e-9
    @test canopy_storage ≈ 0.002043

    # Case canopy_gap_fraction > inv(1.1), potential_evaporation > canopy_storage
    precipitation = 1.1574074074074074e-8
    canopy_gap_fraction = 0.95
    throughfall, canopy_evap, stemflow, canopy_storage = Wflow.rainfall_interception_modrut(
        precipitation,
        potential_evaporation,
        canopy_storage,
        canopy_gap_fraction,
        maximum_canopy_storage,
        dt,
    )
    @test throughfall ≈ 1.099537037037037e-8
    @test canopy_evap ≈ 2.3645833333333334e-8
    @test stemflow ≈ 5.787037037037037e-10
    @test canopy_storage ≈ 0.0
end

@testitem "unit: precipitation_hbv" begin
    include("testing_utils.jl")

    ## Case temperature_interval_snowfall > 0.0
    precipitation = 3.4837962962962964e-7
    temperature = 273.69
    temperature_interval_snowfall = 2.0
    temperature_threshold_snowfall = 273.15
    snow_precip, liquid_precip = Wflow.precipitation_hbv(
        precipitation,
        temperature,
        temperature_interval_snowfall,
        temperature_threshold_snowfall,
    )
    @test snow_precip ≈ 8.012731481481482e-8
    @test liquid_precip ≈ 2.682523148148148e-7

    ## Case temperature_interval_snowfall == 0
    # Case temperature > tt
    temperature_interval_snowfall = 0.0
    snow_precip, liquid_precip = Wflow.precipitation_hbv(
        precipitation,
        temperature,
        temperature_interval_snowfall,
        temperature_threshold_snowfall,
    )
    @test snow_precip == 0.0
    @test liquid_precip == precipitation

    # Case temperate < tt
    temperature = 272.15
    snow_precip, liquid_precip = Wflow.precipitation_hbv(
        precipitation,
        temperature,
        temperature_interval_snowfall,
        temperature_threshold_snowfall,
    )
    @test snow_precip == precipitation
    @test liquid_precip == 0.0
end

@testitem "unit: snowpack_hbv" begin
    dt = 86400.0

    snow_storage = 0.2015
    snow_water = 0.015
    snow_precip = 8.012731481481482e-8
    liquid_precip = 2.682523148148148e-7
    temperature = 273.69
    temperature_threshold_melt = 273.15
    degree_day_factor = 2.8935185185185185e-8
    water_holding_capacity = 0.10
    # Case temperature > temperature_threshold_melt
    snow_storage_new, snow_water, snow_water_equivalent, snow_melt, runoff =
        Wflow.snowpack_hbv(
            snow_storage,
            snow_water,
            snow_precip,
            liquid_precip,
            temperature,
            temperature_threshold_melt,
            degree_day_factor,
            water_holding_capacity,
            dt,
        )
    @test snow_storage_new ≈ 0.207073
    @test snow_water ≈ 0.0207073
    @test snow_water_equivalent ≈ 0.22778030000000002
    @test snow_melt ≈ 1.5625e-8
    @test runoff ≈ 2.1782060185185186e-7

    # Case temperature < temperature_threshold_melt
    temperature = 272.65
    snow_storage_new, snow_water, snow_water_equivalent, snow_melt, runoff =
        Wflow.snowpack_hbv(
            snow_storage,
            snow_water,
            snow_precip,
            liquid_precip,
            temperature,
            temperature_threshold_melt,
            degree_day_factor,
            water_holding_capacity,
            dt,
        )
    @test snow_storage_new ≈ 0.20848550000000002
    @test snow_water ≈ 0.02084855
    @test snow_water_equivalent ≈ 0.22933404999999998
    @test snow_melt ≈ 0.0
    @test runoff ≈ 2.658940972222222e-7
end

@testitem "unit: glacier_hbv" begin
    dt = 86400.0

    glacier_fraction = 0.35
    glacier_store = 0.5
    snow_storage = 0.0095
    temperature = 278.15
    temperature_threshold_melt = 273.15
    degree_day_factor = 3.935185185185185e-8
    snow_to_ice_fraction = 2.3148148148148148e-6
    maximum_snow_to_ice_rate = 9.259259259259259e-8
    snow_storage, snow_to_glacier, glacier_storage, glacier_melt = Wflow.glacier_hbv(
        glacier_fraction,
        glacier_store,
        snow_storage,
        temperature,
        temperature_threshold_melt,
        degree_day_factor,
        snow_to_ice_fraction,
        maximum_snow_to_ice_rate,
        dt,
    )
    @test snow_storage ≈ 0.008835
    @test snow_to_glacier ≈ 2.199074074074074e-8
    @test glacier_storage ≈ 0.4849
    @test glacier_melt ≈ 1.9675925925925924e-7
end

@testitem "unit: infiltration" begin
    dt = 86400.0

    potential_infiltration = 3.18287037037037e-7
    compacted_soil_area_fraction = 0.2
    infiltration_capacity_soil = 5.787037037037037e-7
    infiltration_capacity_compacted_soil = 5.787037037037037e-8
    unsaturated_store_capacity = 0.0235
    f_infilt_reduction = 1.0

    infiltration, infiltration_excess = Wflow.infiltration(
        potential_infiltration,
        compacted_soil_area_fraction,
        infiltration_capacity_soil,
        infiltration_capacity_compacted_soil,
        unsaturated_store_capacity,
        f_infilt_reduction,
        dt,
    )
    @test infiltration_excess ≈ 5.787037037037037e-9
end

@testitem "unit: unsatzone_flow_layer" begin
    dt = 86400.0

    kv_z = 2.962962962962963e-6
    l_sat = 0.135
    c = 12.6

    # Case usd > 0
    usd = 0.043500000000000004
    usd_new, sum_ast = Wflow.unsatzone_flow_layer(usd, kv_z, l_sat, c, dt)
    usd_new = 0.04349983744545384
    sum_ast = 1.6255454615829024e-7

    # Case usd == 0
    usd = 0
    usd_new, sum_ast = Wflow.unsatzone_flow_layer(usd, kv_z, l_sat, c, dt)
    @test usd_new == 0.0
    @test sum_ast == 0.0
end

@testitem "unit: Brooks-Corey soil hydraulic model" begin

    # Case par_lambda > 0
    volumetric_water_content = 0.25
    theta_s = 0.6
    theta_r = 0.15
    c = 10.5
    air_entry_pressure = -0.1
    h = Wflow.head_brooks_corey(
        volumetric_water_content,
        theta_s,
        theta_r,
        c,
        air_entry_pressure,
    )
    @test h ≈ -0.9062998208338441
    @test Wflow.vwc_brooks_corey(h, air_entry_pressure, theta_s, theta_r, c) ≈
          volumetric_water_content + theta_r

    # Case par_lambda < 0
    c = 2.0
    h = Wflow.head_brooks_corey(
        volumetric_water_content,
        theta_s,
        theta_r,
        c,
        air_entry_pressure,
    )
    @test h == air_entry_pressure
    @test Wflow.vwc_brooks_corey(h, air_entry_pressure, theta_s, theta_r, c) ≈ theta_s
end

@testitem "unit: Feddes root water uptake" begin
    h3_high = -3.0
    h3_low = -6.0

    # Case tpot_daily < 1.0
    tpot = 5.787037037037037e-9
    @test Wflow.feddes_h3(h3_high, h3_low, tpot) == h3_low

    # Case 1.0 < tpot_daily < 5.0
    tpot = 3.472222222222222e-8
    @test Wflow.feddes_h3(h3_high, h3_low, tpot) ≈ (h3_high + h3_low) / 2

    # Case tpot_daily > 5.0
    tpot = 8.680555555555556e-8
    @test Wflow.feddes_h3(h3_high, h3_low, tpot) ≈ h3_high

    h1 = -0.1
    h2 = -1.0
    h3 = -3.0
    h4 = -150.0

    ## Case alpha == 0.0
    alpha = 0.0

    # Case h < h4
    h = -160.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 0.0

    # Case h3 < h < h4
    h = -10.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.4 / 1.47

    # Case h2 < h < h3
    h = -1.5
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0

    # Case h1 < h < h2
    h = -0.5
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 4 / 9

    # Case h > h1
    h = -0.05
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 0.0

    ## Case alpha ≠ 0.0
    alpha = 0.5

    # Case h < h4
    h = -160.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 0.0

    # Case h3 < h < h4
    h = -10.0
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.4 / 1.47

    # Case h2 < h < h3
    h = -1.5
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0

    # Case h1 < h < h2
    h = -0.5
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0

    # Case h > h1
    h = -0.05
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0
end

@testitem "unit: soil_temperature" begin
    tsoil_prev = 274.15
    w_soil = 2.0
    temperature = 274.65
    @test Wflow.soil_temperature(tsoil_prev, w_soil, temperature) ≈ 275.15
end

@testitem "unit: infiltration_reduction_factor" begin
    soil_surface_temperature = 273.25
    cf_soil = 0.3
    modelsnow = true

    # Case model_snow && soil_infiltration_reduction
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
    potential_soilevaporation = 3.49537037037037e-9
    unsaturated_layer_depth = 0.00123
    unsaturated_layer_thickness = 0.1
    water_table_depth = 0.3
    theta_effective = 0.241

    # Case n_unsatlayers == 0
    n_unsatlayers = 0
    @test Wflow.soil_evaporation_unsaturated_store(
        potential_soilevaporation,
        unsaturated_layer_depth,
        unsaturated_layer_thickness,
        n_unsatlayers,
        water_table_depth,
        theta_effective,
    ) == 0.0

    # case n_unsatlayers == 1
    n_unsatlayers = 1
    @test Wflow.soil_evaporation_unsaturated_store(
        potential_soilevaporation,
        unsaturated_layer_depth,
        unsaturated_layer_thickness,
        n_unsatlayers,
        water_table_depth,
        theta_effective,
    ) ≈ 5.946480713078224e-11

    # Case n_unsatlayers > 1
    n_unsatlayers = 2
    @test Wflow.soil_evaporation_unsaturated_store(
        potential_soilevaporation,
        unsaturated_layer_depth,
        unsaturated_layer_thickness,
        n_unsatlayers,
        water_table_depth,
        theta_effective,
    ) ≈ 1.783944213923467e-10
end

@testitem "unit: soil_evaporation_saturated_store" begin
    dt = 86400.0

    potential_soilevaporation = 1.4467592592592592e-9
    layerthickness = 0.1
    water_table_depth = 0.3
    theta_effective = 0.32205961644649506

    # Case n_unsatlayers ∈ (0, 1)
    n_unsatlayers = 0
    @test Wflow.soil_evaporation_saturated_store(
        potential_soilevaporation,
        n_unsatlayers,
        layerthickness,
        water_table_depth,
        theta_effective,
        dt,
    ) ≈ -7.455083714039237e-7

    # Case n_unsatlayers ∉ (0, 1)
    n_unsatlayers = 2
    @test Wflow.soil_evaporation_saturated_store(
        potential_soilevaporation,
        n_unsatlayers,
        layerthickness,
        water_table_depth,
        theta_effective,
        dt,
    ) == 0.0
end

@testitem "unit: actual_infiltration_soil_path" begin
    potential_infiltration = 1.883101851851852e-8
    actual_infiltration = 1.883101851851852e-8
    compacted_soil_area_fraction = 0.1
    infiltration_capacity_soil = 2.645787037037037e-6
    infiltration_capacity_compacted_soil = 5.787037037037037e-8
    f_infiltration_reduction = 0.9

    # Case actual_infiltration > 0
    actual_infiltration = 1.883101851851852e-8
    actual_infiltration_soil, actual_infiltration_compacted_soil =
        Wflow.actual_infiltration_soil_path(
            potential_infiltration,
            actual_infiltration,
            compacted_soil_area_fraction,
            infiltration_capacity_soil,
            infiltration_capacity_compacted_soil,
            f_infiltration_reduction,
        )
    @test actual_infiltration_soil ≈ 1.6947916666666665e-8
    @test actual_infiltration_compacted_soil ≈ 1.883101851851852e-9

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
