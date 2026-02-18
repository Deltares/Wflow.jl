@testitem "unit: rainfall_intercepiton_gash" begin
    using Wflow: to_SI
    include("testing_utils.jl")
    dt = 86400.0

    # Case cmax == 0
    cmax = 0
    e_r = 0.11
    canopy_gap_fraction = 0.24
    precipitation = to_SI(18.0, "atmosphere_water__precipitation_volume_flux"; dt_val = dt)
    canopy_storage_in = to_SI(1.5, "vegetation_canopy_water__depth")
    max_evaporation =
        to_SI(4.0, "land_surface_water__potential_evaporation_volume_flux"; dt_val = dt)
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            cmax,
            e_r,
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

    # Case cmax > 0, large_storms == true, interception > max_evaporation
    cmax = to_SI(3.0, "vegetation_water__storage_capacity")
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            cmax,
            e_r,
            canopy_gap_fraction,
            precipitation,
            canopy_storage_in,
            max_evaporation,
            dt,
        )
    @test throughfall ≈
          to_SI(13.568, "vegetation_canopy_water__throughfall_volume_flux"; dt_val = dt)
    @test interception ≈
          to_SI(4.0, "vegetation_canopy_water__interception_volume_flux"; dt_val = dt)
    @test stem_flow ≈
          to_SI(0.432, "vegetation_canopy_water__stemflow_volume_flux"; dt_val = dt)
    @test canopy_storage_in == canopy_storage_out

    # Case cmax > 0, large_storms == false, interception > max_evaporation
    precipitation = to_SI(1.0, "atmosphere_water__precipitation_volume_flux"; dt_val = dt)
    throughfall, interception, stem_flow, canopy_storage_out =
        Wflow.rainfall_interception_gash(
            cmax,
            e_r,
            canopy_gap_fraction,
            precipitation,
            canopy_storage_in,
            max_evaporation,
            dt,
        )
    @test throughfall ≈
          to_SI(0.24, "vegetation_canopy_water__throughfall_volume_flux"; dt_val = dt)
    @test interception ≈
          to_SI(0.736, "vegetation_canopy_water__interception_volume_flux"; dt_val = dt)
    @test stem_flow ≈
          to_SI(0.024, "vegetation_canopy_water__stemflow_volume_flux"; dt_val = dt)
    @test canopy_storage_in == canopy_storage_out
end

@testitem "unit: rainfall_interception_modrut (modified Rutter)" begin
    using Wflow: to_SI
    dt = 86400.0

    # Case canopy_gap_fraction < inv(1.1), potential_evaporation < canopy_storage (after precipitation)
    precipitation = to_SI(8.6, "atmosphere_water__precipitation_volume_flux"; dt_val = dt)
    potential_evaporation =
        to_SI(3.8, "land_surface_water__potential_evaporation_volume_flux"; dt_val = dt)
    canopy_storage = to_SI(1.5, "vegetation_canopy_water__depth")
    canopy_gap_fraction = 0.45
    cmax = to_SI(2.8, "vegetation_water__storage_capacity")
    throughfall, canopy_evap, stemflow, canopy_storage = Wflow.rainfall_interception_modrut(
        precipitation,
        potential_evaporation,
        canopy_storage,
        canopy_gap_fraction,
        cmax,
        dt,
    )
    @test throughfall ≈
          to_SI(3.87, "vegetation_canopy_water__throughfall_volume_flux"; dt_val = dt)
    @test canopy_evap ≈
          to_SI(3.8, "vegetation_canopy_water__interception_volume_flux"; dt_val = dt)
    @test stemflow ≈
          to_SI(0.387, "vegetation_canopy_water__stemflow_volume_flux"; dt_val = dt)
    @test canopy_storage ≈ to_SI(2.043, "vegetation_canopy_water__depth")

    # Case canopy_gap_fraction > inv(1.1), potential_evaporation > canopy_storage
    precipitation = to_SI(1.0, "atmosphere_water__precipitation_volume_flux"; dt_val = dt)
    canopy_gap_fraction = 0.95
    throughfall, canopy_evap, stemflow, canopy_storage = Wflow.rainfall_interception_modrut(
        precipitation,
        potential_evaporation,
        canopy_storage,
        canopy_gap_fraction,
        cmax,
        dt,
    )
    @test throughfall ≈
          to_SI(0.95, "vegetation_canopy_water__throughfall_volume_flux"; dt_val = dt)
    @test canopy_evap ≈
          to_SI(2.043, "vegetation_canopy_water__interception_volume_flux"; dt_val = dt)
    @test stemflow ≈
          to_SI(0.05, "vegetation_canopy_water__stemflow_volume_flux"; dt_val = dt)
    @test canopy_storage ≈ 0.0
end

@testitem "unit: precipitation_hbv" begin
    using Wflow: Unit, to_SI
    dt = 86400.0

    ## Case tti > 0.0
    precipitation = to_SI(30.1, "atmosphere_water__precipitation_volume_flux"; dt_val = dt)
    temperature = to_SI(0.54, "atmosphere_air__temperature")
    tti = 2.0
    tt = to_SI(0.0, "atmosphere_air__snowfall_temperature_threshold")
    snow_precip, liquid_precip =
        Wflow.precipitation_hbv(precipitation, temperature, tti, tt)
    @test snow_precip ≈
          to_SI(6.923, "atmosphere_water__precipitation_volume_flux"; dt_val = dt)
    @test liquid_precip ≈
          to_SI(23.177, "atmosphere_water__precipitation_volume_flux"; dt_val = dt)

    ## Case tti == 0
    # Case temperature > tt
    tti = 0.0
    snow_precip, liquid_precip =
        Wflow.precipitation_hbv(precipitation, temperature, tti, tt)
    @test snow_precip == 0.0
    @test liquid_precip == precipitation

    # Case temperate < tt
    temperature = to_SI(-1.0, "atmosphere_air__temperature")
    snow_precip, liquid_precip =
        Wflow.precipitation_hbv(precipitation, temperature, tti, tt)
    @test snow_precip == precipitation
    @test liquid_precip == 0.0
end

@testitem "unit: snowpack_hbv" begin
    using Wflow: to_SI, Unit, MM, ABSOLUTE_DEGREES, MM_PER_DT
    DDF = Unit(; mm = 1, degC = -1, dt = -1)
    dt = 86400.0

    snow_storage = to_SI(201.5, "snowpack_dry_snow__leq_depth")
    snow_water = to_SI(15.0, "snowpack_liquid_water__depth")
    snow_precip = to_SI(6.923, "atmosphere_water__precipitation_volume_flux"; dt_val = dt)
    liquid_precip =
        to_SI(23.177, "atmosphere_water__precipitation_volume_flux"; dt_val = dt)
    temperature = to_SI(0.54, "atmosphere_air__temperature")
    ttm = to_SI(0.0, "snowpack__melting_temperature_threshold")
    cfmax = to_SI(2.5, "snowpack__degree_day_coefficient"; dt_val = dt)
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
            dt,
        )
    @test snow_storage_new ≈ to_SI(207.073, "snowpack_dry_snow__leq_depth")
    @test snow_water ≈ to_SI(20.7073, "snowpack_liquid_water__depth")
    @test snow_water_equivalent ≈ to_SI(227.7803, "snowpack__leq_depth")
    @test snow_melt ≈ to_SI(1.35, "snowpack_meltwater__volume_flux"; dt_val = dt)
    @test runoff ≈ to_SI(18.8197, "snowpack_water__runoff_volume_flux"; dt_val = dt)

    # Case temperature < ttm
    temperature = to_SI(-0.5, "atmosphere_air__temperature")
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
            dt,
        )
    @test snow_storage_new ≈ to_SI(208.4855, "snowpack_dry_snow__leq_depth")
    @test snow_water ≈ to_SI(20.84855, "snowpack_liquid_water__depth")
    @test snow_water_equivalent ≈ to_SI(229.33405, "snowpack__leq_depth")
    @test snow_melt ≈ 0.0
    @test runoff ≈ to_SI(22.97325, "snowpack_water__runoff_volume_flux"; dt_val = dt)
end

@testitem "unit: glacier_hbv" begin
    using Wflow: to_SI
    dt = 86400.0

    glacier_frac = 0.35
    glacier_store = to_SI(500.0, "glacier_ice__leq_depth")
    snow_storage = to_SI(9.5, "snowpack_dry_snow__leq_depth")
    temperature = to_SI(5.0, "atmosphere_air__temperature")
    ttm = to_SI(0.0, "glacier_ice__melting_temperature_threshold")
    cfmax = to_SI(3.4, "glacier_ice__degree_day_coefficient"; dt_val = dt)
    g_sifrac = to_SI(
        0.2,
        "glacier_firn_accumulation__snowpack_dry_snow_leq_depth_fraction";
        dt_val = dt,
    )
    max_snow_to_glacier = to_SI(8.0, "soil_surface_water__runoff_volume_flux"; dt_val = dt)
    snow_storage, snow_to_glacier, glacier_storage, glacier_melt = Wflow.glacier_hbv(
        glacier_frac,
        glacier_store,
        snow_storage,
        temperature,
        ttm,
        cfmax,
        g_sifrac,
        max_snow_to_glacier,
        dt,
    )
    @test snow_storage ≈ to_SI(8.835, "snowpack_dry_snow__leq_depth")
    @test snow_to_glacier ≈
          to_SI(1.9, "soil_surface_water__runoff_volume_flux"; dt_val = dt)
    @test glacier_storage ≈ to_SI(484.9, "glacier_ice__leq_depth")
    @test glacier_melt ≈ to_SI(17.0, "glacier_ice__melt_volume_flux"; dt_val = dt)
end

@testitem "unit: infiltration" begin
    using Wflow: to_SI
    dt = 86400.0

    potential_infiltration =
        to_SI(27.5, "soil_water__infiltration_volume_flux"; dt_val = dt)
    pathfrac = 0.2
    infiltcapsoil = to_SI(50.0, "soil_water__infiltration_volume_flux"; dt_val = dt)
    infiltcappath =
        to_SI(5.0, "compacted_soil_surface_water__infiltration_capacity"; dt_val = dt)
    ustorecapacity = to_SI(23.5, "soil_layer_water_unsaturated_zone__depth")
    f_infilt_reduction = 1.0

    infiltsoilpath, infiltexcess = Wflow.infiltration(
        potential_infiltration,
        pathfrac,
        infiltcapsoil,
        infiltcappath,
        ustorecapacity,
        f_infilt_reduction,
        dt,
    )
    @test infiltexcess ≈ to_SI(0.5, "soil_water__infiltration_volume_flux"; dt_val = dt)
end

@testitem "unit: unsatzone_flow_layer" begin
    using Wflow: to_SI
    dt = 86400.0

    kv_z = to_SI(
        256.0,
        "soil_layer_water__vertical_saturated_hydraulic_conductivity";
        dt_val = dt,
    )
    l_sat = to_SI(135.0, "soil_water_saturated_zone__depth")
    c = 12.6

    # Case usd > 0
    usd = to_SI(43.5, "soil_layer_water_unsaturated_zone__depth")
    usd_new, sum_ast = Wflow.unsatzone_flow_layer(usd, kv_z, l_sat, c, dt)
    usd_new = to_SI(43.49983744545384, "soil_layer_water_unsaturated_zone__depth")
    sum_ast = to_SI(0.00016255454615829025, "soil_layer_water_unsaturated_zone__depth")

    # Case usd == 0
    usd = 0
    usd_new, sum_ast = Wflow.unsatzone_flow_layer(usd, kv_z, l_sat, c, dt)
    @test usd_new == 0.0
    @test sum_ast == 0.0
end

@testitem "unit: Brooks-Corey soil hydraulic model" begin
    using Wflow: to_SI

    # Case par_lambda > 0
    vwc = 0.25
    theta_s = 0.6
    theta_r = 0.15
    c = 10.5
    hb = to_SI(-10.0, "soil_water__air_entry_pressure_head")
    h = Wflow.head_brooks_corey(vwc, theta_s, theta_r, c, hb)
    @test h ≈ to_SI(-90.6299820833844, "soil_water__air_entry_pressure_head")
    @test Wflow.vwc_brooks_corey(h, hb, theta_s, theta_r, c) ≈ vwc + theta_r

    # Case par_lambda < 0
    c = 2.0
    h = Wflow.head_brooks_corey(vwc, theta_s, theta_r, c, hb)
    @test h == hb
    @test Wflow.vwc_brooks_corey(h, hb, theta_s, theta_r, c) ≈ theta_s
end

@testitem "unit: Feddes root water uptake" begin
    using Wflow: to_SI
    dt = 86400.0

    h3_high = to_SI(-300.0, "soil_water__air_entry_pressure_head")
    h3_low = to_SI(-600.0, "soil_water__air_entry_pressure_head")

    # Case tpot_daily < 1.0
    tpot = to_SI(0.5, "soil_water__transpiration_volume_flux"; dt_val = dt)
    @test Wflow.feddes_h3(h3_high, h3_low, tpot) == h3_low

    # Case 1.0 < tpot_daily < 5.0
    tpot = to_SI(3.0, "soil_water__transpiration_volume_flux"; dt_val = dt)
    @test Wflow.feddes_h3(h3_high, h3_low, tpot) ≈ (h3_high + h3_low) / 2

    # Case tpot_daily > 5.0
    tpot = to_SI(7.5, "soil_water__transpiration_volume_flux"; dt_val = dt)
    @test Wflow.feddes_h3(h3_high, h3_low, tpot) ≈ h3_high

    h1 = to_SI(-10.0, "soil_water__air_entry_pressure_head")
    h2 = to_SI(-100.0, "soil_water__air_entry_pressure_head")
    h3 = to_SI(-300.0, "soil_water__air_entry_pressure_head")
    h4 = to_SI(-15000.0, "soil_water__air_entry_pressure_head")

    ## Case alpha == 0.0
    alpha = 0.0

    # Case h < h4
    h = to_SI(-16000.0, "soil_water__air_entry_pressure_head")
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 0.0

    # Case h3 < h < h4
    h = to_SI(-1000.0, "soil_water__air_entry_pressure_head")
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.4 / 1.47

    # Case h2 < h < h3
    h = to_SI(-150.0, "soil_water__air_entry_pressure_head")
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0

    # Case h1 < h < h2
    h = to_SI(-50.0, "soil_water__air_entry_pressure_head")
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 4 / 9

    # Case h > h1
    h = to_SI(-5.0, "soil_water__air_entry_pressure_head")
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 0.0

    ## Case alpha ≠ 0.0
    alpha = 0.5

    # Case h < h4
    h = to_SI(-16000.0, "soil_water__air_entry_pressure_head")
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 0.0

    # Case h3 < h < h4
    h = to_SI(-1000.0, "soil_water__air_entry_pressure_head")
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.4 / 1.47

    # Case h2 < h < h3
    h = to_SI(-150.0, "soil_water__air_entry_pressure_head")
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0

    # Case h1 < h < h2
    h = to_SI(-50.0, "soil_water__air_entry_pressure_head")
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0

    # Case h > h1
    h = to_SI(-5.0, "soil_water__air_entry_pressure_head")
    @test Wflow.rwu_reduction_feddes(h, h1, h2, h3, h4, alpha) ≈ 1.0
end

@testitem "unit: soil_temperature" begin
    using Wflow: to_SI

    tsoil_prev = to_SI(1.0, "soil_surface__temperature")
    w_soil = 2.0
    temperature = to_SI(1.5, "atmosphere_air__temperature")
    @test Wflow.soil_temperature(tsoil_prev, w_soil, temperature) ≈
          to_SI(2.0, "soil_surface__temperature")
end

@testitem "unit: infiltration_reduction_factor" begin
    using Wflow: to_SI

    tsoil = to_SI(0.1, "soil_surface__temperature")
    cf_soil = 0.3
    modelsnow = true

    # Case model_snow && soil_infiltration_reduction
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
    using Wflow: to_SI
    dt = 86400.0

    potential_soilevaporation =
        to_SI(0.302, "land_surface_water__potential_evaporation_volume_flux"; dt_val = dt)
    ustorelayerdepth = to_SI(1.23, "soil_layer_water_unsaturated_zone__depth")
    ustorelayerthickness = to_SI(100.0, "soil__thickness")
    zi = to_SI(300.0, "soil_water_saturated_zone_top__depth")
    theta_effective = 0.241

    # Case n_unsatlayers == 0
    n_unsatlayers = 0
    @test Wflow.soil_evaporation_unsaturated_store(
        potential_soilevaporation,
        ustorelayerdepth,
        ustorelayerthickness,
        n_unsatlayers,
        zi,
        theta_effective,
    ) == 0.0

    # case n_unsatlayers == 1
    n_unsatlayers = 1
    @test Wflow.soil_evaporation_unsaturated_store(
        potential_soilevaporation,
        ustorelayerdepth,
        ustorelayerthickness,
        n_unsatlayers,
        zi,
        theta_effective,
    ) ≈ to_SI(
        0.005137759336099585,
        "land_surface_water__potential_evaporation_volume_flux";
        dt_val = dt,
    )

    # Case n_unsatlayers > 1
    n_unsatlayers = 2
    @test Wflow.soil_evaporation_unsaturated_store(
        potential_soilevaporation,
        ustorelayerdepth,
        ustorelayerthickness,
        n_unsatlayers,
        zi,
        theta_effective,
    ) ≈ to_SI(
        0.015413278008298757,
        "land_surface_water__potential_evaporation_volume_flux";
        dt_val = dt,
    )
end

@testitem "unit: soil_evaporation_saturated_store" begin
    using Wflow: to_SI
    dt = 86400.0

    potential_soilevaporation =
        to_SI(0.125, "land_surface_water__potential_evaporation_volume_flux"; dt_val = dt)
    layerthickness = to_SI(100.0, "soil__thickness")
    zi = to_SI(300.0, "soil_water_saturated_zone_top__depth")
    theta_effective = 0.32205961644649506

    # Case n_unsatlayers ∈ (0, 1)
    n_unsatlayers = 0
    @test Wflow.soil_evaporation_saturated_store(
        potential_soilevaporation,
        n_unsatlayers,
        layerthickness,
        zi,
        theta_effective,
        dt,
    ) ≈ to_SI(
        -64.41192328929901,
        "land_surface_water__potential_evaporation_volume_flux";
        dt_val = dt,
    )

    # Case n_unsatlayers ∉ (0, 1)
    n_unsatlayers = 2
    @test Wflow.soil_evaporation_saturated_store(
        potential_soilevaporation,
        n_unsatlayers,
        layerthickness,
        zi,
        theta_effective,
        dt,
    ) == 0.0
end

@testitem "unit: actual_infiltration_soil_path" begin
    using Wflow: to_SI
    dt = 86400.0

    potential_infiltration =
        to_SI(1.627, "soil_water__infiltration_volume_flux"; dt_val = dt)
    actinfilt = to_SI(1.627, "soil_water__infiltration_volume_flux"; dt_val = dt)
    pathfrac = 0.1
    infiltcapsoil = to_SI(228.596, "soil_water__infiltration_volume_flux"; dt_val = dt)
    infiltcappath =
        to_SI(5.0, "compacted_soil_surface_water__infiltration_capacity"; dt_val = dt)
    f_infiltration_reduction = 0.9

    # Case actinfilt > 0
    actinfilt = to_SI(1.627, "soil_water__infiltration_volume_flux"; dt_val = dt)
    actinfiltsoil, actinfiltpath = Wflow.actual_infiltration_soil_path(
        potential_infiltration,
        actinfilt,
        pathfrac,
        infiltcapsoil,
        infiltcappath,
        f_infiltration_reduction,
    )
    @test actinfiltsoil ≈ to_SI(1.4643, "soil_water__infiltration_volume_flux"; dt_val = dt)
    @test actinfiltpath ≈
          to_SI(0.1627, "compacted_soil_surface_water__infiltration_capacity"; dt_val = dt)

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
