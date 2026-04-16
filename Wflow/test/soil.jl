@testitem "unit: update_bc_soil_model!" begin
    using Wflow: to_SI, MM_PER_DT
    using StaticArrays: SVector

    include("testing_utils.jl")
    n = 1
    dt = 86400.0

    soil_model = init_sbm_soil_model(1, 1; soil_fraction = [0.2397957498236932])
    atmospheric_forcing = Wflow.AtmosphericForcing(;
        n,
        potential_evaporation = [
            to_SI(
                0.5799999833106995,
                "land_surface_water__potential_evaporation_volume_flux";
                dt_val = dt,
            ),
        ],
    )

    interception = Wflow.GashInterceptionModel(;
        n,
        parameters = Wflow.GashParameters(;
            e_r = [0.25],
            vegetation_parameter_set = Wflow.VegetationParameters(;
                leaf_area_index = [1.5],
                storage_wood = [2.0e-4],
                kext = [0.67],
                storage_specific_leaf = [9.0e-5],
                canopygapfraction = [0.3487189230509198],
                cmax = [to_SI(0.5687566481997495, "vegetation_water__storage_capacity")],
                rootingdepth = [to_SI(410.0, "vegetation_root__depth")],
                kc = [1.0],
            ),
        ),
        variables = Wflow.InterceptionVariables(;
            n,
            canopy_potevap = [
                to_SI(
                    0.44091845241498073,
                    "land_surface_water__potential_evaporation_volume_flux";
                    dt_val = dt,
                ),
            ],
            interception_rate = [to_SI(0.029448986349521342, MM_PER_DT; dt_val = dt)],
        ),
    )
    runoff = Wflow.OpenWaterRunoff(;
        n,
        variables = Wflow.OpenWaterRunoffVariables(;
            n,
            runoff_land = [0.0],
            runoff_river = [0.0],
        ),
        boundary_conditions = Wflow.OpenWaterRunoffBC(;
            n,
            water_flux_surface = [
                to_SI(
                    0.24686226660437197,
                    "soil_surface_water__runoff_volume_flux";
                    dt_val = dt,
                ),
            ],
        ),
    )
    demand = Wflow.DemandModel(;
        domestic = Wflow.NoDemandModel(; n),
        industry = Wflow.NoDemandModel(; n),
        livestock = Wflow.NoDemandModel(; n),
        paddy = Wflow.NoIrrigationPaddyModel(n),
        nonpaddy = Wflow.NoIrrigationNonPaddyModel(n),
        variables = Wflow.DemandVariables(; n),
    )

    allocation = Wflow.NoAllocationLandModel(n)

    external_models = (; interception, runoff, demand, allocation)

    Wflow.update_bc_soil_model!(soil_model, atmospheric_forcing, external_models, dt)

    @test soil_model.boundary_conditions.potential_transpiration[1] ≈
          to_SI(0.41146946606545937, MM_PER_DT; dt_val = dt)
    @test soil_model.boundary_conditions.potential_soilevaporation[1] ≈
          to_SI(0.13908153089571873, MM_PER_DT; dt_val = dt)
    @test soil_model.boundary_conditions.water_flux_surface[1] ≈
          to_SI(0.24686226660437197, MM_PER_DT; dt_val = dt)
end

@testitem "unit: unsaturated_zone_flow!" begin
    using Wflow: to_SI
    using StaticArrays: SVector
    dt = 86400.0
    include("testing_utils.jl")
    n = 1
    N = 6

    soil_model = init_sbm_soil_model(
        n,
        N;
        n,
        ustorelayerthickness = [
            to_SI.(
                SVector((50.0, 5.081648613929929, NaN, NaN, NaN, NaN)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector(1.2855527211118947, 0.20814868098590805, 0.0, 0.0, 0.0, 0.0),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        n_unsatlayers = [2],
        theta_s = [0.4414711594581604],
        theta_r = [0.08942600339651108],
        kv_profile = Wflow.KvExponential(
            [to_SI(364.1460876464844, "soil_water__infiltration_volume_flux"; dt_val = dt)],
            [
                to_SI(
                    0.0033079576678574085,
                    "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
                ),
            ],
        ),
        kvfrac = [SVector((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))],
        c = [
            SVector((
                9.121646881103516,
                9.247220993041992,
                9.514554023742676,
                9.675407409667969,
                9.831438064575195,
                9.716856956481934,
            )),
        ],
        infiltsoilpath = [
            to_SI(0.02277405542765069, "soil_water__infiltration_volume_flux"; dt_val = dt),
        ],
    )

    Wflow.unsaturated_zone_flow!(soil_model, dt)

    @test soil_model.variables.ustorelayerdepth[1] ≈ to_SI(
        [1.3083267609636298, 0.20814799974448514, 0.0, 0.0, 0.0, 0.0],
        "soil_layer_water_unsaturated_zone__depth",
    )
    @test soil_model.variables.transfer[1] ≈ to_SI(
        6.968173386761924e-7,
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
    soil_model = init_sbm_soil_model(
        n,
        N;
        potential_soilevaporation = [to_SI(2.8, MM_PER_DT; dt_val = dt)],
        ustorelayerthickness = [
            to_SI.(
                SVector((50.0, 21.472680450878443, NaN, NaN, NaN, NaN)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector(1.537249298366254, 4.821326813899425e-8, 0.0, 0.0, 0.0, 0.0),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        n_unsatlayers = [2],
        zi = [to_SI(71.47268045087844, "soil_water_saturated_zone_top__depth")],
        theta_s = [0.44],
        theta_r = [0.09],
        theta_fc = [0.275],
        act_thickl = [
            to_SI.(
                SVector(50.0, 100.0, 50.0, 200.0, 800.0, 800.0),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        drainable_waterdepth = [
            to_SI(321.13323174500624, "soil_water_saturated_zone__depth"),
        ],
    )

    Wflow.soil_evaporation!(soil_model, dt)

    @test soil_model.variables.soilevapsat[1] == 0
    @test soil_model.variables.soilevap[1] ≈
          to_SI(0.2459598877386006, MM_PER_DT; dt_val = dt)
    @test soil_model.variables.drainable_waterdepth[1] ≈
          to_SI(321.13323174500624, "soil_water_saturated_zone__depth")
end

@testitem "unit: transpiration!" begin
    using Wflow: to_SI
    include("testing_utils.jl")
    dt = 86400.0
    n = 1
    N = 4

    soil_model = init_sbm_soil_model(
        n,
        N;
        h3_high = [to_SI(-400.0, "vegetation_root__feddes_critical_pressure_head_h3_high")],
        h3_low = [to_SI(-1000.0, "vegetation_root__feddes_critical_pressure_head_h3_low")],
        potential_transpiration = [
            to_SI(
                0.05153840858869719,
                "soil_water__transpiration_volume_flux";
                dt_val = dt,
            ),
        ],
        n_unsatlayers = [2],
        zi = [to_SI(106.89587841733061, "soil_water_saturated_zone_top__depth")],
        ustorelayerthickness = [
            to_SI.(
                SVector((100.0, 6.895878417330607, NaN, NaN)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector(10.9327977152876, 0.8620432154993639, 0.0, 0.0),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        rootingdepth = [to_SI(453.0, "vegetation_root__depth")],
        rootfraction = [
            SVector(0.22075055187637968, 0.6622516556291391, 0.11699779249448124, 0.0),
        ],
        act_thickl = [
            to_SI.(
                SVector(100.0, 300.0, 200.0, NaN),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        sumlayers = [
            to_SI.(
                SVector(0.0, 100.0, 400.0, 600.0, NaN),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        c = [
            SVector((
                9.53970437651816,
                10.007558316712927,
                10.603868189606647,
                10.662998826419395,
            )),
        ],
        theta_s = [0.4790319800376892],
        theta_r = [0.17089612782001495],
        hb = [to_SI(-10.0, "soil_water__air_entry_pressure_head")],
        h1 = [to_SI(0.0, "vegetation_root__feddes_critical_pressure_head_h1")],
        h2 = [to_SI(-100.0, "vegetation_root__feddes_critical_pressure_head_h2")],
        h3 = [to_SI(-1000.0, "vegetation_root__feddes_critical_pressure_head_h1")],
        h4 = [to_SI(-16000.0, "vegetation_root__feddes_critical_pressure_head_h4")],
        alpha_h1 = [1.0],
        rootdistpar = [
            to_SI(
                -500.0,
                "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
            ),
        ],
        drainable_waterdepth = [
            to_SI(72.40310797113221, "soil_water_saturated_zone__depth"),
        ],
    )

    Wflow.transpiration!(soil_model, dt)

    @test soil_model.variables.ae_ustore[1] ≈
          to_SI(0.05153840858869719, "soil_water__transpiration_volume_flux"; dt_val = dt)
    @test soil_model.variables.actevapsat[1] ≈ 0.0
    @test soil_model.variables.drainable_waterdepth[1] ≈
          to_SI(72.40310797113221, "soil_water_saturated_zone__depth")
    @test soil_model.variables.transpiration[1] ≈
          to_SI(0.05153840858869719, "soil_water__transpiration_volume_flux"; dt_val = dt)
end

@testitem "unit: capillary_flux!" begin
    using Wflow: to_SI
    include("testing_utils.jl")
    n = 1
    N = 6
    dt = 86400.0
    soil_model = init_sbm_soil_model(
        n,
        N;
        rootingdepth = [to_SI(384.1, "vegetation_root__depth")],
        n_unsatlayers = [6],
        kv_profile = Wflow.KvExponential(
            [to_SI(2018.0371, "soil_water__infiltration_volume_flux"; dt_val = dt)],
            [
                to_SI(
                    1.29274e-3,
                    "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
                ),
            ],
        ),
        kvfrac = [SVector((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))],
        zi = [to_SI(1266.33589, "soil_water_saturated_zone_top__depth")],
        ustorelayerthickness = [
            to_SI.(
                SVector((50.0, 100.0, 50.0, 200.0, 800.0, 66.33589243)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector((
                    8.874129508377953,
                    18.18721052056329,
                    8.54476050597162,
                    17.279406870498256,
                    147.78493107547982,
                    12.51769175241036,
                )),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ae_ustore = [
            to_SI(0.528773283329, "soil_water__transpiration_volume_flux"; dt_val = dt),
        ],
        ustorecapacity = [
            to_SI(313.1741519821792, "soil_layer_water_unsaturated_zone__depth"),
        ],
        drainable_waterdepth = [
            to_SI(169.02603585525694, "soil_water_saturated_zone__depth"),
        ],
        cap_hmax = [
            to_SI(
                2000.0,
                "soil_water_saturated_zone_top__capillary_rise_max_water_table_depth",
            ),
        ],
        cap_n = [2.0],
        theta_s = [0.4868114888668],
        theta_r = [0.0711537748575],
    )

    Wflow.capillary_flux!(soil_model, dt)

    @test soil_model.variables.actcapflux[1] ≈ to_SI(
        0.07115477692809025,
        "soil_water_saturated_zone_top__capillary_volume_flux";
        dt_val = dt,
    )
end

@testitem "unit: update_soil_water_storage! SbmSoilModel" begin
    using StaticArrays: SVector
    using Wflow: to_SI, Unit, MM_PER_DT, MM
    M_PER_DT = Unit(; m = 1, dt = -1)
    include("testing_utils.jl")
    dt = 86400.0
    n = 1
    N = 4
    soil_model = init_sbm_soil_model(
        n,
        N;
        runoff = [to_SI(0.0, MM_PER_DT; dt_val = dt)],
        zi = [to_SI(1244.5135404970033, "soil_water_saturated_zone_top__depth")],
        ustorelayerthickness = [
            to_SI.(
                SVector((100.0, 300.0, 800.0, 44.513540497003305)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector(
                    14.408928105784874,
                    19.461087500813772,
                    125.91094235311145,
                    7.223489814154271,
                ),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        n_unsatlayers = [4],
        vwc = [
            SVector((
                0.27533306202340196,
                0.1874167718952191,
                0.28208626066069253,
                0.48316940665245056,
            )),
        ],
        vwc_perc = [
            SVector((56.98478799206181, 38.7890394786585, 58.38247554104776, 100.0)),
        ],
        nlayers = [4],
        act_thickl = [
            to_SI.(
                SVector(100.0, 300.0, 800.0, 800.0),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        theta_s = [0.4831694066524056],
        theta_r = [0.12372369319200516],
        theta_fc = [0.28599992944511926],
        rootingdepth = [to_SI(380.0, "vegetation_root__depth")],
        sumlayers = [
            to_SI.(
                SVector(0.0, 100.0, 400.0, 1200.0, 2000.0),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        soilthickness = [to_SI(2000.0, "soil__thickness")],
        soilwatercapacity = [to_SI(718.8914269208908, "soil__thickness")],
        ustorecapacity = [
            to_SI(271.02508737207194, "soil_layer_water_unsaturated_zone__depth"),
        ],
        satwaterdepth = [to_SI(288.62821979680706, "soil_water_saturated_zone__depth")],
        drainable_waterdepth = [
            to_SI(158.32342151683937, "soil_water_saturated_zone__depth"),
        ],
        total_soilwater_storage = [448.5154852374101],
        excesswater = [
            to_SI(0.0, "compacted_soil_surface_water__excess_volume_flux"; dt_val = dt),
        ],
        infiltexcess = [0.0],
    )

    runoff = (; variables = (; runoff_land = [0.0], ae_openw_l = [0.0]))

    demand = Wflow.DemandModel(;
        domestic = Wflow.NoDemandModel(; n),
        industry = Wflow.NoDemandModel(; n),
        livestock = Wflow.NoDemandModel(; n),
        paddy = Wflow.NoIrrigationPaddyModel(n),
        nonpaddy = Wflow.NoIrrigationNonPaddyModel(n),
        variables = Wflow.DemandVariables(; n),
    )

    subsurface_flow = (; variables = (; exfiltwater = [0.0]))
    external_models = (; runoff, demand, subsurface_flow)

    Wflow.update_soil_water_storage!(soil_model, external_models, dt)

    @test soil_model.variables.runoff[1] ≈ 0.0
    @test soil_model.variables.ustorecapacity[1] ≈
          to_SI(280.33060970129986, "soil__thickness")
    @test soil_model.variables.satwaterdepth[1] ≈
          to_SI(271.55636944572655, "soil_water_saturated_zone__depth")
    @test soil_model.variables.drainable_waterdepth[1] ≈ to_SI(148.95887025738955, MM)
    @test soil_model.variables.exfiltsatwater[1] ≈ 0.0
    @test soil_model.variables.vwc[1] ≈
          [0.2678129742498539, 0.1885939848613844, 0.2811123711333945, 0.4721985172668562]
    @test soil_model.variables.vwc_perc[1] ≈
          [55.42837989378741, 39.03268341595556, 58.18091279434587, 97.72939072000436]
    @test soil_model.variables.vwc_root[1] ≈ 0.20944108733203426
    @test soil_model.variables.vwc_percroot[1] ≈ 43.34734038380604
    @test soil_model.variables.total_soilwater_storage[1] ≈ to_SI(438.56081721959094, MM)
end
