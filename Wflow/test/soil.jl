@testitem "unit: update_bc_soil_model!" begin
    using Wflow: to_SI, MM_PER_DT
    using StaticArrays: SVector

    include("testing_utils.jl")
    n = 1
    dt = 86400.0

    soil_model = init_sbm_soil_model(1, 1; soil_fraction = [0.30670667824303055])
    atmospheric_forcing = Wflow.AtmosphericForcing(;
        n,
        potential_evaporation = [
            to_SI(
                4.8000001907348635,
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
                canopygapfraction = [0.3487189230509198],
                cmax = [to_SI(0.2235326728960274, "vegetation_water__storage_capacity")],
                rootingdepth = [to_SI(390.531982421875, "vegetation_root__depth")],
                kc = [1.1469999551773071],
            ),
        ),
        variables = Wflow.InterceptionVariables(;
            n,
            canopy_potevap = [
                to_SI(
                    3.585693099611687,
                    "land_surface_water__potential_evaporation_volume_flux";
                    dt_val = dt,
                ),
            ],
            interception_rate = [to_SI(0.12944592473374397, MM_PER_DT; dt_val = dt)],
        ),
    )
    runoff = Wflow.OpenWaterRunoff(;
        n,
        variables = Wflow.OpenWaterRunoffVariables(;
            n,
            runoff_land = [0.0],
            runoff_river = [
                to_SI(
                    0.003384257254905341,
                    "soil_surface_water__runoff_volume_flux";
                    dt_val = dt,
                ),
            ],
        ),
        boundary_conditions = Wflow.OpenWaterRunoffBC(;
            n,
            water_flux_surface = [
                to_SI(7.5e-3, "soil_surface_water__runoff_volume_flux"; dt_val = dt),
            ],
        ),
    )
    demand = Wflow.DemandModel(;
        domestic = Wflow.NoDemandModel(; n),
        industry = Wflow.NoDemandModel(; n),
        livestock = Wflow.NoDemandModel(; n),
        paddy = Wflow.PaddyModel(;
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

    Wflow.update_bc_soil_model!(soil_model, atmospheric_forcing, external_models, dt)

    @test soil_model.boundary_conditions.potential_transpiration[1] ≈
          to_SI(3.456247174877943, "soil_water__transpiration_volume_flux"; dt_val = dt)
    @test soil_model.boundary_conditions.potential_soilevaporation[1] ≈ to_SI(
        1.472182114066203,
        "land_surface_water__potential_evaporation_volume_flux";
        dt_val = dt,
    )
    @test soil_model.boundary_conditions.water_flux_surface[1] ≈
          to_SI(0.02411574274509466, "soil_surface_water__runoff_volume_flux"; dt_val = dt)
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

    @test model.variables.ustorelayerdepth[1] ≈ to_SI(
        [1.3083267609636298, 0.20814799974448514, 0.0, 0.0, 0.0, 0.0],
        "soil_layer_water_unsaturated_zone__depth",
    )
    @test model.variables.transfer[1] ≈ to_SI(
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
    dt = 86400
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

    @test model.variables.ae_ustore[1] ≈
          to_SI(0.05153840858869719, "soil_water__transpiration_volume_flux"; dt_val = dt)
    @test model.variables.actevapsat[1] ≈ 0.0
    @test model.variables.drainable_waterdepth[1] ≈
          to_SI(72.40310797113221, "soil_water_saturated_zone__depth")
    @test model.variables.transpiration[1] ≈
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

    @test model.variables.actcapflux[1] ≈ to_SI(
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
    N = 6
    soil_model = init_sbm_soil_model(
        n,
        N;
        runoff = [to_SI(0.002036245, MM_PER_DT; dt_val = dt)],
        zi = [to_SI(0.000656547, "soil_water_saturated_zone_top__depth")],
        ustorelayerthickness = [
            to_SI.(
                SVector((0.000656547, NaN, NaN, NaN, NaN, NaN)),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        ustorelayerdepth = [
            to_SI.(
                SVector(1.998965544e-5, 0.0, 0.0, 0.0, 0.0, 0.0),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        n_unsatlayers = [1],
        vwc = [SVector((0.440136, 0.440140, 0.440140, 0.440140, 0.440140, 0.440140))],
        vwc_perc = [SVector((99.99912, 100.0, 100.0, 100.0, 100.0, 100.0))],
        nlayers = [6],
        act_thickl = [
            to_SI.(
                SVector(50.0, 100.0, 50.0, 200.0, 800.0, 800.0),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        theta_s = [0.44014],
        theta_r = [0.0880263],
        theta_fc = [0.26583],
        rootingdepth = [to_SI(390.39999, "vegetation_root__depth")],
        sumlayers = [
            to_SI.(
                SVector(0.0, 50.0, 150.0, 200.0, 400.0, 1200.0, 2000.0),
                Ref("soil_layer_water_unsaturated_zone__depth"),
            ),
        ],
        soilthickness = [to_SI(704.2274028, "soil__thickness")],
        soilwatercapacity = [to_SI(1215.5, "soil__thickness")],
        [to_SI(0.0001915933, "soil_layer_water_unsaturated_zone__depth")]satwaterdepth = [
            to_SI(704.227211212, "soil_water_saturated_zone__depth"),
        ],
        drainable_waterdepth = [to_SI(348.61417105, "soil_water_saturated_zone__depth")],
        total_soilwater_storage = [704.227211212],
        excesswater = [
            to_SI(
                0.047112555773,
                "compacted_soil_surface_water__excess_volume_flux";
                dt_val = dt,
            ),
        ],
        infiltexcess = [0.0],
    )

    runoff = (; variables = (; runoff_land = [0.0], ae_openw_l = [0.0]))

    demand = Wflow.DemandModel(;
        domestic = Wflow.NoDemandModel(; n),
        industry = Wflow.NoDemandModel(; n),
        livestock = Wflow.NoDemandModel(; n),
        paddy = Wflow.PaddyModel(;
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
        nonpaddy = Wflow.NoIrrigationNonPaddyModel(n),
        variables = Wflow.DemandVariables(; n),
    )

    subsurface_flow = (; variables = (; exfiltwater = [0.0]))
    external_models = (; runoff, demand, subsurface_flow)

    Wflow.update_soil_water_storage!(soil_model, external_models, dt)

    @test demand.paddy.variables.h[1] ≈ to_SI(0.047112555773, "paddy_surface_water__depth")
    @test model.variables.runoff[1] ≈ 0.0
    @test model.variables.ustorecapacity[1] ≈
          to_SI(0.00021398953802734852, "soil__thickness")
    @test model.variables.satwaterdepth[1] ≈
          to_SI(704.2271688208066, "soil_water_saturated_zone__depth")
    @test model.variables.drainable_waterdepth[1] ≈ to_SI(348.6198855572923, MM)
    @test model.variables.exfiltsatwater[1] ≈ 0.0
    @test model.variables.vwc[1] ≈
          [0.4401357762092409, 0.44014, 0.44014, 0.44014, 0.44014, 0.44014]
    @test model.variables.vwc_perc[1] ≈
          [99.99904035289701, 100.0, 100.0, 100.0, 100.0, 100.0]
    @test model.variables.vwc_root[1] ≈ 0.44013945904317786
    @test model.variables.vwc_percroot[1] ≈ 99.99987709437403
    @test model.variables.total_soilwater_storage[1] ≈ to_SI(704.227188810462, MM)
end
