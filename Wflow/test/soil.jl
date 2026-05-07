@testitem "unit: update_bc_soil_model!" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    n = 1
    soil_model = init_sbm_soil_model(1, 1; soil_fraction = [0.2397957498236932])
    atmospheric_forcing =
        Wflow.AtmosphericForcing(; n, potential_evaporation = [0.5799999833106995])

    interception = Wflow.GashInterceptionModel(;
        parameters = Wflow.GashParameters(;
            evaporation_to_precipitation_ratio = [0.25],
            vegetation_parameter_set = Wflow.VegetationParameters(;
                leaf_area_index = [1.5],
                storage_wood = [0.2],
                light_extinction_coefficient = [0.67],
                storage_specific_leaf = [0.09],
                canopy_gap_fraction = [0.2397957498236932],
                maximum_canopy_storage = [0.5687566481997495],
                rooting_depth = [410.0],
                crop_coefficient = [1.0],
            ),
        ),
        variables = Wflow.InterceptionVariables(;
            n = 1,
            canopy_potevap = [0.44091845241498073],
            interception_rate = [0.029448986349521342],
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
            water_flux_surface = [0.24686226660437197],
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

    Wflow.update_bc_soil_model!(soil_model, atmospheric_forcing, external_models)

    @test soil_model.boundary_conditions.potential_transpiration[1] ≈ 0.41146946606545937
    @test soil_model.boundary_conditions.potential_soilevaporation[1] ≈ 0.13908153089571873
    @test soil_model.boundary_conditions.water_flux_surface[1] ≈ 0.24686226660437197
end

@testitem "unit: unsaturated_zone_flow!" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    n = 1
    N = 6
    soil_model = init_sbm_soil_model(
        n,
        N;
        unsaturated_layer_thickness = [
            SVector((50.0, 5.081648613929929, NaN, NaN, NaN, NaN)),
        ],
        unsaturated_layer_depth = [
            SVector(1.2855527211118947, 0.20814868098590805, 0.0, 0.0, 0.0, 0.0),
        ],
        n_unsatlayers = [2],
        theta_s = [0.4414711594581604],
        theta_r = [0.08942600339651108],
        kv_profile = Wflow.KvExponential([364.1460876464844], [0.0033079576678574085]),
        vertical_conductivity_factor = [SVector((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))],
        brooks_corey_exponent = [
            SVector((
                9.121646881103516,
                9.247220993041992,
                9.514554023742676,
                9.675407409667969,
                9.831438064575195,
                9.716856956481934,
            )),
        ],
        infiltration_soil_and_path = [0.02277405542765069],
    )
    Wflow.unsaturated_zone_flow!(soil_model)

    @test soil_model.variables.unsaturated_layer_depth[1] ≈
          [1.3083267609636298, 0.20814799974448514, 0.0, 0.0, 0.0, 0.0]
    @test soil_model.variables.transfer[1] ≈ 6.968173386761924e-7
end

@testitem "unit: soil_evaporation!" begin
    include("testing_utils.jl")
    n = 1
    N = 6
    soil_model = init_sbm_soil_model(
        n,
        N;
        potential_soilevaporation = [2.8],
        unsaturated_layer_thickness = [
            SVector((50.0, 21.472680450878443, NaN, NaN, NaN, NaN)),
        ],
        unsaturated_layer_depth = [
            SVector(1.537249298366254, 4.821326813899425e-8, 0.0, 0.0, 0.0, 0.0),
        ],
        n_unsatlayers = [2],
        water_table_depth = [71.47268045087844],
        theta_s = [0.44],
        theta_r = [0.09],
        theta_fc = [0.275],
        actual_layer_thickness = [SVector(50.0, 100.0, 50.0, 200.0, 800.0, 800.0)],
        drainable_water_depth = [321.13323174500624],
    )

    Wflow.soil_evaporation!(soil_model)

    @test soil_model.variables.soil_evaporation_saturated[1] == 0.0
    @test soil_model.variables.soil_evaporation[1] ≈ 0.2459598877386006
    @test soil_model.variables.drainable_water_depth[1] ≈ 321.13323174500624
end

@testitem "unit: transpiration!" begin
    include("testing_utils.jl")
    n = 1
    N = 4
    soil_model = init_sbm_soil_model(
        n,
        N;
        h3_high = [-400.0],
        h3_low = [-1000.0],
        potential_transpiration = [0.05153840858869719],
        n_unsatlayers = [2],
        water_table_depth = [106.89587841733061],
        unsaturated_layer_thickness = [SVector((100.0, 6.895878417330607, NaN, NaN))],
        unsaturated_layer_depth = [SVector(10.9327977152876, 0.8620432154993639, 0.0, 0.0)],
        rooting_depth = [453.0],
        rootfraction = [
            SVector(0.22075055187637968, 0.6622516556291391, 0.11699779249448124, 0.0),
        ],
        actual_layer_thickness = [SVector(100.0, 300.0, 200.0, NaN)],
        cumulative_layer_depth = [SVector(0.0, 100.0, 400.0, 600.0, NaN)],
        brooks_corey_exponent = [
            SVector((
                9.53970437651816,
                10.007558316712927,
                10.603868189606647,
                10.662998826419395,
            )),
        ],
        theta_s = [0.4790319800376892],
        theta_r = [0.17089612782001495],
        air_entry_pressure = [-10.0],
        h1 = [0.0],
        h2 = [-100.0],
        h3 = [1000.0],
        h4 = [-16000.0],
        alpha_h1 = [1.0],
        root_distribution_parameter = [-500.0],
        drainable_water_depth = [72.40310797113221],
    )
    dt = 86400.0

    Wflow.transpiration!(soil_model, dt)

    @test soil_model.variables.actual_evaporation_unsaturated_store[1] ≈ 0.05153840858869719
    @test soil_model.variables.actual_evaporation_saturated[1] ≈ 0.0
    @test soil_model.variables.drainable_water_depth[1] ≈ 72.40310797113221
    @test soil_model.variables.transpiration[1] ≈ 0.05153840858869719
end

@testitem "unit: capillary_flux!" begin
    include("testing_utils.jl")
    n = 1
    N = 6
    soil_model = init_sbm_soil_model(
        n,
        N;
        rooting_depth = [384.1],
        n_unsatlayers = [6],
        kv_profile = Wflow.KvExponential([2018.0371], [1.29274e-3]),
        vertical_conductivity_factor = [SVector((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))],
        water_table_depth = [1266.33589],
        unsaturated_layer_thickness = [
            SVector((50.0, 100.0, 50.0, 200.0, 800.0, 66.33589243)),
        ],
        unsaturated_layer_depth = [
            SVector((
                8.874129508377953,
                18.18721052056329,
                8.54476050597162,
                17.279406870498256,
                147.78493107547982,
                12.51769175241036,
            )),
        ],
        actual_evaporation_unsaturated_store = [0.528773283329],
        unsaturated_store_capacity = [313.1741519821792],
        drainable_water_depth = [169.02603585525694],
        cap_hmax = [2000.0],
        cap_n = [2.0],
        theta_s = [0.4868114888668],
        theta_r = [0.0711537748575],
    )

    Wflow.capillary_flux!(soil_model)

    @test soil_model.variables.actual_capillary_flux[1] ≈ 0.07115477692809025
end

@testitem "unit: update_soil_water_storage! SbmSoilModel" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    n = 1
    N = 4
    soil_model = init_sbm_soil_model(
        n,
        N;
        runoff = [0.0],
        water_table_depth = [1244.5135404970033],
        unsaturated_layer_thickness = [SVector((100.0, 300.0, 800.0, 44.513540497003305))],
        unsaturated_layer_depth = [
            SVector(
                14.408928105784874,
                19.461087500813772,
                125.91094235311145,
                7.223489814154271,
            ),
        ],
        n_unsatlayers = [4],
        volumetric_water_content = [
            SVector((
                0.27533306202340196,
                0.1874167718952191,
                0.28208626066069253,
                0.48316940665245056,
            )),
        ],
        volumetric_water_content_percent = [
            SVector((56.98478799206181, 38.7890394786585, 58.38247554104776, 100.0)),
        ],
        number_of_layers = [4],
        actual_layer_thickness = [SVector(100.0, 300.0, 800.0, 800.0)],
        theta_s = [0.4831694066524056],
        theta_r = [0.12372369319200516],
        theta_fc = [0.28599992944511926],
        rooting_depth = [380.0],
        cumulative_layer_depth = [SVector(0.0, 100.0, 400.0, 1200.0, 2000.0)],
        soil_thickness = [2000.0],
        soil_water_capacity = [718.8914269208908],
        unsaturated_store_capacity = [271.02508737207194],
        saturated_water_depth = [288.62821979680706],
        drainable_water_depth = [158.32342151683937],
        total_soil_water_storage = [448.5154852374101],
        saturation_excess_water = [0.0],
        infiltration_excess = [0.0],
    )

    runoff =
        (; variables = (; runoff_land = [0.0], actual_evaporation_open_water_land = [0.0]))
    demand = Wflow.DemandModel(;
        domestic = Wflow.NoDemandModel(; n),
        industry = Wflow.NoDemandModel(; n),
        livestock = Wflow.NoDemandModel(; n),
        paddy = Wflow.NoIrrigationPaddyModel(n),
        nonpaddy = Wflow.NoIrrigationNonPaddyModel(n),
        variables = Wflow.DemandVariables(; n),
    )

    subsurface_flow = (; variables = (; exfiltration_water = [0.0]))
    external_models = (; runoff, demand, subsurface_flow)

    Wflow.update_soil_water_storage!(soil_model, external_models)

    @test soil_model.variables.runoff[1] ≈ 0.0
    @test soil_model.variables.unsaturated_store_capacity[1] ≈ 280.33060970129986
    @test soil_model.variables.saturated_water_depth[1] ≈ 271.55636944572655
    @test soil_model.variables.drainable_water_depth[1] ≈ 148.95887025738955
    @test soil_model.variables.exfiltration_saturated_water[1] ≈ 0.0
    @test soil_model.variables.volumetric_water_content[1] ≈
          [0.2678129742498539, 0.1885939848613844, 0.2811123711333945, 0.4721985172668562]
    @test soil_model.variables.volumetric_water_content_percent[1] ≈
          [55.42837989378741, 39.03268341595556, 58.18091279434587, 97.72939072000436]
    @test soil_model.variables.volumetric_water_content_root_zone[1] ≈ 0.20944108733203426
    @test soil_model.variables.volumetric_water_content_percent_root_zone[1] ≈
          43.34734038380604
    @test soil_model.variables.total_soil_water_storage[1] ≈ 438.56081721959094
end
