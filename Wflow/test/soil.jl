@testitem "unit: update_boundary_conditions!" begin
    using Wflow: to_SI
    using StaticArrays: SVector
    include("testing_utils.jl")
    n = 1
    dt = 86400.0
    model = init_sbm_soil_model(1, 1; soil_fraction = [0.3])
    atmospheric_forcing = Wflow.AtmosphericForcing(;
        n,
        potential_evaporation = [
            to_SI(
                1e-3,
                "land_surface_water__potential_evaporation_volume_flux";
                dt_val = dt,
            ),
        ],
    )

    interception = Wflow.GashInterceptionModel(;
        n,
        parameters = Wflow.GashParameters(;
            e_r = [0.1],
            vegetation_parameter_set = Wflow.VegetationParameters(;
                leaf_area_index = nothing,
                storage_wood = nothing,
                kext = nothing,
                storage_specific_leaf = nothing,
                canopygapfraction = [0.1],
                cmax = [to_SI(10.0, "vegetation_water__storage_capacity")],
                rootingdepth = [to_SI(10.0, "vegetation_root__depth")],
                kc = [0.1],
            ),
        ),
        variables = Wflow.InterceptionVariables(;
            n,
            canopy_potevap = [
                to_SI(
                    4e-3,
                    "land_surface_water__potential_evaporation_volume_flux";
                    dt_val = dt,
                ),
            ],
            interception_rate = [
                to_SI(
                    2e-3,
                    "vegetation_canopy_water__interception_volume_flux";
                    dt_val = dt,
                ),
            ],
        ),
    )
    runoff = Wflow.OpenWaterRunoff(;
        n,
        variables = Wflow.OpenWaterRunoffVariables(;
            n,
            runoff_land = [
                to_SI(6e-4, "soil_surface_water__runoff_volume_flux"; dt_val = dt),
            ],
            runoff_river = [
                to_SI(7e-3, "soil_surface_water__runoff_volume_flux"; dt_val = dt),
            ],
        ),
        boundary_conditions = Wflow.OpenWaterRunoffBC(;
            n,
            water_flux_surface = [
                to_SI(7.5e-3, "soil_surface_water__runoff_volume_flux"; dt_val = dt),
            ],
        ),
    )
    demand = Wflow.Demand(;
        domestic = Wflow.NoNonIrrigationDemand(n),
        industry = Wflow.NoNonIrrigationDemand(n),
        livestock = Wflow.NoNonIrrigationDemand(n),
        paddy = Wflow.Paddy(;
            n,
            parameters = Wflow.PaddyParameters(;
                irrigation_efficiency = [],
                maximum_irrigation_rate = [],
                irrigation_areas = [true],
                irrigation_trigger = [],
                h_min = [],
                h_opt = [],
                h_max = [],
            ),
            variables = Wflow.PaddyVariables(;
                n,
                h = [to_SI(1e-5, "vegetation_water__storage_capacity")],
                evaporation = [
                    to_SI(
                        4e-3,
                        "land_surface_water__potential_evaporation_volume_flux";
                        dt_val = dt,
                    ),
                ],
            ),
        ),
        nonpaddy = Wflow.NoIrrigationNonPaddy(n),
        variables = Wflow.DemandVariables(; n),
    )

    allocation = Wflow.AllocationLand(;
        n,
        parameters = Wflow.AllocationLandParameters(; frac_sw_used = [], areas = []),
        variables = Wflow.AllocationLandVariables(;
            n,
            irri_alloc = [
                to_SI(2e-2, "land__allocated_irrigation_water_volume_flux"; dt_val = dt),
            ],
        ),
    )

    external_models = (; interception, runoff, demand, allocation)

    Wflow.update_boundary_conditions!(model, atmospheric_forcing, external_models, dt)

    @test model.boundary_conditions.potential_transpiration[1] ≈
          to_SI(2e-3, "soil_water__transpiration_volume_flux"; dt_val = dt)
    @test model.boundary_conditions.potential_soilevaporation[1] ≈ to_SI(
        2.9e-4,
        "land_surface_water__potential_evaporation_volume_flux";
        dt_val = dt,
    )
    @test model.boundary_conditions.water_flux_surface[1] ≈
          to_SI(1.99e-2, "soil_surface_water__runoff_volume_flux"; dt_val = dt)
end

@testitem "unit: unsaturated_zone_flow!" begin
    using Wflow: to_SI
    using StaticArrays: SVector
    dt = 86400.0
    include("testing_utils.jl")
    n = 1
    N = 6
    model = init_sbm_soil_model(
        n,
        N;
        n,
        ustorelayerthickness = [
            to_SI.(
                SVector((50.0, 100.0, 50.0, 200.0, 800.0, 250.0)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector(8.3, 16.6, 7.7, 36.5, 190.7, 59.2),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        n_unsatlayers = [6],
        theta_s = [0.5],
        theta_r = [0.07],
        kv_profile = Wflow.KvExponential(
            [to_SI(2530.0, "soil_water__infiltration_volume_flux"; dt_val = dt)],
            [
                to_SI(
                    1e-3,
                    "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
                ),
            ],
        ),
        kvfrac = [SVector((0.9, 0.8, 0.9, 0.8, 0.9, 0.8))],
        c = [SVector((9.3, 9.1, 8.9, 8.7, 8.5, 8.3))],
        infiltsoilpath = [to_SI(1.5, "soil_water__infiltration_volume_flux"; dt_val = dt)],
    )
    Wflow.unsaturated_zone_flow!(model, dt)

    @test model.variables.ustorelayerdepth[1] ≈ to_SI(
        [
            8.86158736157703,
            17.080057700709837,
            7.851580246695791,
            36.02511190264408,
            187.18720614458448,
            59.05659877245674,
        ],
        "soil_layer_water_unsaturated_zone__depth",
    )
    @test model.variables.transfer[1] ≈ to_SI(
        4.437857871332029,
        "soil_water_saturated_zone_top__net_recharge_volume_flux";
        dt_val = dt,
    )
end

@testitem "unit: soil_evaporation!" begin
    using Wflow: to_SI, MM_PER_DT
    include("testing_utils.jl")
    n = 1
    N = 6
    dt = 86400.0
    model = init_sbm_soil_model(
        n,
        N;
        potential_soilevaporation = [to_SI(3.9, MM_PER_DT; dt_val = dt)],
        ustorelayerthickness = [
            to_SI.(
                SVector((50.0, 100.0, 50.0, 200.0, 800.0, 250.0)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector(8.3, 16.6, 7.7, 36.5, 190.7, 59.2),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        n_unsatlayers = [1],
        zi = [to_SI(40.0, "soil_water_saturated_zone_top__depth")],
        theta_s = [0.5],
        theta_r = [0.07],
        theta_fc = [0.25],
        act_thickl = [
            to_SI.(
                SVector(50.0, 100.0, 50.0, 200.0, 800.0, 573.3),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        drainable_waterdepth = [to_SI(100.0, "soil_water_saturated_zone__depth")],
    )

    Wflow.soil_evaporation!(model, dt)

    @test model.variables.soilevapsat[1] ≈
          to_SI(0.40360465116279065, MM_PER_DT; dt_val = dt)
    @test model.variables.soilevap[1] ≈ to_SI(2.2855813953488378, MM_PER_DT; dt_val = dt)
    @test model.variables.drainable_waterdepth[1] ≈
          to_SI(99.5963953488372, "soil_water_saturated_zone__depth")
end

@testitem "unit: transpiration!" begin
    using Wflow: to_SI
    dt = 86400.0
    include("testing_utils.jl")
    n = 1
    N = 6
    model = init_sbm_soil_model(
        n,
        N;
        h3_high = [to_SI(-400.0, "vegetation_root__feddes_critical_pressure_head_h3_high")],
        h3_low = [to_SI(-100.0, "vegetation_root__feddes_critical_pressure_head_h3_low")],
        potential_transpiration = [
            to_SI(1.2, "soil_water__transpiration_volume_flux"; dt_val = dt),
        ],
        n_unsatlayers = [6],
        zi = [to_SI(75.0, "soil_water_saturated_zone_top__depth")],
        ustorelayerthickness = [
            to_SI.(
                SVector((50.0, 100.0, 50.0, 200.0, 800.0, 250.0)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector(8.3, 16.6, 7.7, 36.5, 190.7, 59.2),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        rootingdepth = [to_SI(106.8, "vegetation_root__depth")],
        rootfraction = [SVector(0.125, 0.25, 0.125, 0.5, 0.0, 0.0)],
        act_thickl = [
            to_SI.(
                SVector(50.0, 100.0, 50.0, 200.0, 800.0, 573.3),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        sumlayers = [
            to_SI.(
                SVector(0.0, 50.0, 150.0, 200.0, 400.0, 1200.0, 1773.3),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        c = [SVector((9.3, 9.1, 8.9, 8.7, 8.5, 8.3))],
        theta_s = [0.5],
        theta_r = [0.07],
        hb = [to_SI(-10.0, "soil_water__air_entry_pressure_head")],
        h1 = [to_SI(0.0, "vegetation_root__feddes_critical_pressure_head_h1")],
        h2 = [to_SI(-100.0, "vegetation_root__feddes_critical_pressure_head_h2")],
        h3 = [to_SI(-915.0, "vegetation_root__feddes_critical_pressure_head_h1")],
        h4 = [to_SI(-15849.0, "vegetation_root__feddes_critical_pressure_head_h4")],
        alpha_h1 = [1.0],
        rootdistpar = [
            to_SI(
                -500.0,
                "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
            ),
        ],
        drainable_waterdepth = [to_SI(100.0, "soil_water_saturated_zone__depth")],
    )

    Wflow.transpiration!(model, dt)

    @test model.variables.ae_ustore[1] ≈
          to_SI(0.4479021069963573, "soil_water__transpiration_volume_flux"; dt_val = dt)
    @test model.variables.actevapsat[1] ≈
          to_SI(0.7520978930036426, "soil_water__transpiration_volume_flux"; dt_val = dt)
    @test model.variables.drainable_waterdepth[1] ≈
          to_SI(99.24790210699636, "soil_water_saturated_zone__depth")
    @test model.variables.transpiration[1] ≈
          to_SI(1.2, "soil_water__transpiration_volume_flux"; dt_val = dt)
end

@testitem "unit: capillary_flux!" begin
    using Wflow: to_SI
    include("testing_utils.jl")
    n = 1
    N = 6
    dt = 86400.0
    model = init_sbm_soil_model(
        n,
        N;
        rootingdepth = [to_SI(106.8, "vegetation_root__depth")],
        n_unsatlayers = [6],
        kv_profile = Wflow.KvExponential(
            [to_SI(2530.0, "soil_water__infiltration_volume_flux"; dt_val = dt)],
            [
                to_SI(
                    1e-3,
                    "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
                ),
            ],
        ),
        kvfrac = [SVector((0.9, 0.8, 0.9, 0.8, 0.9, 0.8))],
        zi = [to_SI(112.0, "soil_water_saturated_zone_top__depth")],
        ustorelayerthickness = [
            to_SI.(
                SVector((50.0, 100.0, 50.0, 200.0, 800.0, 250.0)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector(8.3, 16.6, 7.7, 36.5, 190.7, 59.2),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ae_ustore = [to_SI(1.25, "soil_water__transpiration_volume_flux"; dt_val = dt)],
        ustorecapacity = [to_SI(385.0, "soil_layer_water_unsaturated_zone__depth")],
        drainable_waterdepth = [to_SI(100.0, "soil_water_saturated_zone__depth")],
        cap_hmax = [
            to_SI(
                2000.0,
                "soil_water_saturated_zone_top__capillary_rise_max_water_table_depth",
            ),
        ],
        cap_n = [2.0],
        theta_s = [0.5],
        theta_r = [0.07],
    )

    Wflow.capillary_flux!(model, dt)

    @test model.variables.actcapflux[1] ≈ to_SI(
        1.11392,
        "soil_water_saturated_zone_top__capillary_volume_flux";
        dt_val = dt,
    )
end

@testitem "unit: update! SbmSoilModel" begin
    using StaticArrays: SVector
    using Wflow: to_SI, Unit, MM_PER_DT, MM
    M_PER_DT = Unit(; m = 1, dt = -1)
    dt = 86400.0
    include("testing_utils.jl")
    n = 1
    N = 6
    model = init_sbm_soil_model(
        n,
        N;
        runoff = [to_SI(5.0, MM_PER_DT; dt_val = dt)],
        zi = [to_SI(50.0, "soil_water_saturated_zone_top__depth")],
        ustorelayerthickness = [
            to_SI.(
                SVector((50.0, 100.0, 50.0, 200.0, 800.0, 250.0)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector(8.3, 16.6, 7.7, 36.5, 190.7, 59.2),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        n_unsatlayers = [6],
        vwc = [SVector((0.23, 0.24, 0.22, 0.24, 0.30, 0.42))],
        vwc_perc = [SVector((44.0, 45.0, 42.0, 47.0, 58.0, 81.0))],
        nlayers = [6],
        act_thickl = [
            to_SI.(
                SVector(50.0, 100.0, 50.0, 200.0, 800.0, 573.3),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        theta_s = [0.5],
        theta_r = [0.07],
        theta_fc = [0.25],
        rootingdepth = [to_SI(106.8, "vegetation_root__depth")],
        sumlayers = [
            to_SI.(
                SVector(0.0, 50.0, 150.0, 200.0, 400.0, 1200.0, 1773.3),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        soilthickness = [to_SI(2000.0, "soil__thickness")],
        soilwatercapacity = [to_SI(1215.5, "soil__thickness")],
        ustorecapacity = [0.0],
        satwaterdepth = [0.0],
        drainable_waterdepth = [0.0],
        total_soilwater_storage = [0.0],
        excesswater = [
            to_SI(0.25, "compacted_soil_surface_water__excess_volume_flux"; dt_val = dt),
        ],
        infiltexcess = [to_SI(0.1, MM_PER_DT; dt_val = dt)],
    )

    runoff = (;
        variables = (;
            runoff_land = [to_SI(0.03, MM_PER_DT; dt_val = dt)],
            ae_openw_l = [to_SI(0.05, MM_PER_DT; dt_val = dt)],
        )
    )
    demand = Wflow.Demand(;
        domestic = Wflow.NoNonIrrigationDemand(n),
        industry = Wflow.NoNonIrrigationDemand(n),
        livestock = Wflow.NoNonIrrigationDemand(n),
        paddy = Wflow.Paddy(;
            n,
            parameters = Wflow.PaddyParameters(;
                irrigation_efficiency = [],
                maximum_irrigation_rate = [],
                irrigation_areas = [true],
                irrigation_trigger = [],
                h_min = [],
                h_opt = [],
                h_max = [to_SI(2.0, "irrigated_paddy__max_depth")],
            ),
            variables = Wflow.PaddyVariables(;
                n,
                h = [1e-5],
                evaporation = [to_SI(4e-3, MM_PER_DT; dt_val = dt)],
            ),
        ),
        nonpaddy = Wflow.NoIrrigationNonPaddy(n),
        variables = Wflow.DemandVariables(; n),
    )

    subsurface_flow = (; exfiltwater = [to_SI(0.0025, M_PER_DT; dt_val = dt)])
    Wflow.get_exfiltwater(nt::NamedTuple) = nt.exfiltwater

    external_models = (; runoff, demand, subsurface_flow)

    Wflow.update!(model, external_models, dt)

    @test demand.paddy.variables.h[1] ≈ to_SI(2.0, "paddy_surface_water__depth")
    @test model.variables.runoff[1] ≈
          to_SI(0.88, "soil_surface_water__runoff_volume_flux"; dt_val = dt)
    @test model.variables.ustorecapacity[1] ≈ to_SI(58.0, "soil__thickness")
    @test model.variables.satwaterdepth[1] ≈
          to_SI(838.5, "soil_water_saturated_zone__depth")
    @test model.variables.drainable_waterdepth[1] ≈ to_SI(487.5, MM)
    @test model.variables.exfiltsatwater[1] ≈ to_SI(
        2.5,
        "soil_surface_water_saturated_zone__exfiltration_volume_flux";
        dt_val = dt,
    )
    @test model.variables.vwc[1] ≈
          [0.236, 0.236, 0.224, 0.2525, 0.308375, 0.4157509157509158]
    @test model.variables.vwc_perc[1] ≈ [47.2, 47.2, 44.8, 50.5, 61.675, 83.15018315018315]
    @test model.variables.vwc_root[1] ≈ 0.464689138576779
    @test model.variables.vwc_percroot[1] ≈ 92.93782771535581
    @test model.variables.total_soilwater_storage[1] ≈ to_SI(1157.5, MM)
end
