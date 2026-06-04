@testitem "unit: update_bc_soil_model!" begin
    using StaticArrays: SVector

    include("testing_utils.jl")
    n = 1
    dt = 86400.0

    soil_model = init_sbm_soil_model(1, 1; soil_fraction = [0.2397957498236932])
    atmospheric_forcing = Wflow.AtmosphericForcing(;
        n_cells = n,
        potential_evaporation = [6.712962769799762e-9],
    )

    interception = Wflow.GashInterceptionModel(;
        n,
        parameters = Wflow.GashParameters(;
            evaporation_to_precipitation_ratio = [0.25],
            vegetation_parameter_set = Wflow.VegetationParameters(;
                leaf_area_index = [1.5],
                storage_wood = [2.0e-4],
                light_extinction_coefficient = [0.67],
                storage_specific_leaf = [9.0e-5],
                canopy_gap_fraction = [0.3487189230509198],
                maximum_canopy_storage = [0.0005687566481997495],
                rooting_depth = [0.41000000000000003],
                crop_coefficient = [1.0],
            ),
        ),
        variables = Wflow.InterceptionVariables(;
            n,
            canopy_potevap = [5.103222828877092e-9],
            interception_rate = [3.408447494157563e-10],
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
            water_flux_surface = [2.8572021597728237e-9],
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

    @test soil_model.boundary_conditions.potential_transpiration[1] ≈ 4.762378079461335e-9
    @test soil_model.boundary_conditions.potential_soilevaporation[1] ≈
          1.6097399409226706e-9
    @test soil_model.boundary_conditions.water_flux_surface[1] ≈ 2.8572021597728237e-9
end

@testitem "unit: unsaturated_zone_flow!" begin
    using StaticArrays: SVector
    dt = 86400.0
    include("testing_utils.jl")
    n = 1
    N = 6

    soil_model = init_sbm_soil_model(
        n,
        N;
        unsaturated_layer_thickness = [
            SVector((0.05, 0.005081648613929929, NaN, NaN, NaN, NaN)),
        ],
        unsaturated_layer_depth = [
            SVector(0.0012855527211118947, 0.00020814868098590806, 0.0, 0.0, 0.0, 0.0),
        ],
        n_unsatlayers = [2],
        theta_s = [0.4414711594581604],
        theta_r = [0.08942600339651108],
        kv_profile = Wflow.KvExponential([4.21465379220468e-6], [3.3079576678574085]),
        vertical_hydraulic_conductivity_factor = [SVector((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))],
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
        infiltration = [2.635886044866978e-10],
    )

    Wflow.unsaturated_zone_flow!(soil_model, dt)

    @test soil_model.variables.unsaturated_layer_depth[1] ≈
          SVector(0.0013083267609636298, 0.00020814799974448514, 0.0, 0.0, 0.0, 0.0)
    @test soil_model.variables.transfer[1] ≈ 8.065015493937412e-15
end

@testitem "unit: soil_evaporation!" begin
    include("testing_utils.jl")
    n = 1
    N = 6
    dt = 86400.0
    soil_model = init_sbm_soil_model(
        n,
        N;
        potential_soilevaporation = [3.2407407407407403e-8],
        unsaturated_layer_thickness = [
            SVector((0.05, 0.021472680450878443, NaN, NaN, NaN, NaN)),
        ],
        unsaturated_layer_depth = [
            SVector(0.001537249298366254, 4.8213268138994254e-11, 0.0, 0.0, 0.0, 0.0),
        ],
        n_unsatlayers = [2],
        water_table_depth = [0.07147268045087844],
        theta_s = [0.44],
        theta_r = [0.09],
        theta_fc = [0.275],
        actual_layer_thickness = [SVector(0.05, 0.1, 0.05, 0.2, 0.8, 0.8)],
        drainable_water_depth = [0.32113323174500624],
    )

    Wflow.soil_evaporation!(soil_model, dt)

    @test soil_model.variables.soil_evaporation_saturated_zone[1] == 0
    @test soil_model.variables.soil_evaporation[1] ≈ 2.846757959937507e-9
    @test soil_model.variables.drainable_water_depth[1] ≈ 0.32113323174500624
end

@testitem "unit: transpiration!" begin
    include("testing_utils.jl")
    dt = 86400.0
    n = 1
    N = 4

    soil_model = init_sbm_soil_model(
        n,
        N;
        h3_high = [-4.0],
        h3_low = [-10.0],
        potential_transpiration = [5.965093586654767e-10],
        n_unsatlayers = [2],
        water_table_depth = [0.10689587841733061],
        unsaturated_layer_thickness = [SVector((0.1, 0.006895878417330607, NaN, NaN))],
        unsaturated_layer_depth = [
            SVector(0.010932797715287601, 0.000862043215499364, 0.0, 0.0),
        ],
        rooting_depth = [0.453],
        rootfraction = [
            SVector(0.22075055187637968, 0.6622516556291391, 0.11699779249448124, 0.0),
        ],
        actual_layer_thickness = [SVector(0.1, 0.3, 0.2, NaN)],
        cumulative_layer_depth = [SVector(0.0, 0.1, 0.4, 0.6, NaN)],
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
        air_entry_pressure = [-0.1],
        h1 = [0.0],
        h2 = [-1.0],
        h3 = [-10.0],
        h4 = [-160.0],
        alpha_h1 = [1.0],
        wet_root_distribution_parameter = [-500000.0],
        drainable_water_depth = [0.07240310797113221],
    )

    Wflow.transpiration!(soil_model, dt)

    @test soil_model.variables.actual_evaporation_unsaturated_store[1] ≈
          5.965093586654767e-10
    @test soil_model.variables.actual_evaporation_saturated_zone[1] ≈ 0.0
    @test soil_model.variables.drainable_water_depth[1] ≈ 0.07240310797113221
    @test soil_model.variables.transpiration[1] ≈ 5.965093586654767e-10
end

@testitem "unit: capillary_flux!" begin
    include("testing_utils.jl")
    n = 1
    N = 6
    dt = 86400.0
    soil_model = init_sbm_soil_model(
        n,
        N;
        rooting_depth = [0.38410000000000005],
        n_unsatlayers = [6],
        kv_profile = Wflow.KvExponential([2.335691087962963e-5], [1.29274]),
        vertical_hydraulic_conductivity_factor = [SVector((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))],
        water_table_depth = [1.2663358900000001],
        unsaturated_layer_thickness = [SVector((0.05, 0.1, 0.05, 0.2, 0.8, 0.06633589243))],
        unsaturated_layer_depth = [
            SVector((
                0.008874129508377954,
                0.018187210520563293,
                0.00854476050597162,
                0.017279406870498257,
                0.14778493107547983,
                0.01251769175241036,
            )),
        ],
        actual_evaporation_unsaturated_store = [6.120061149641203e-9],
        unsaturated_store_capacity = [0.3131741519821792],
        drainable_water_depth = [0.16902603585525694],
        cap_hmax = [2.0],
        cap_n = [2.0],
        theta_s = [0.4868114888668],
        theta_r = [0.0711537748575],
    )

    Wflow.capillary_flux!(soil_model, dt)

    @test soil_model.variables.actual_capillary_flux[1] ≈ 8.235506588899334e-10
end

@testitem "unit: update_soil_water_storage! SbmSoilModel" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    dt = 86400.0
    n = 1
    N = 4
    soil_model = init_sbm_soil_model(
        n,
        N;
        runoff = [0.0],
        water_table_depth = [1.2445135404970034],
        unsaturated_layer_thickness = [SVector((0.1, 0.3, 0.8, 0.044513540497003304))],
        unsaturated_layer_depth = [
            SVector(
                0.014408928105784874,
                0.01946108750081377,
                0.12591094235311145,
                0.007223489814154271,
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
        relative_volumetric_water_content = [
            SVector((56.98478799206181, 38.7890394786585, 58.38247554104776, 100.0)),
        ],
        number_of_layers = [4],
        actual_layer_thickness = [SVector(0.1, 0.3, 0.8, 0.8)],
        theta_s = [0.4831694066524056],
        theta_r = [0.12372369319200516],
        theta_fc = [0.28599992944511926],
        rooting_depth = [0.38],
        cumulative_layer_depth = [SVector(0.0, 0.1, 0.4, 1.2, 2.0)],
        soil_thickness = [2.0],
        soil_water_capacity = [0.7188914269208908],
        unsaturated_store_capacity = [0.2710250873720719],
        saturated_water_depth = [0.2886282197968071],
        drainable_water_depth = [0.15832342151683937],
        total_soil_water_storage = [448.5154852374101],
        saturation_excess_water = [0.0],
        infiltration_excess = [0.0],
    )

    runoff =
        (; variables = (; runoff_land = [0.0], actual_open_water_evaporation_land = [0.0]))

    demand = Wflow.DemandModel(;
        domestic = Wflow.NoDemandModel(; n),
        industry = Wflow.NoDemandModel(; n),
        livestock = Wflow.NoDemandModel(; n),
        paddy = Wflow.NoIrrigationPaddyModel(n),
        nonpaddy = Wflow.NoIrrigationNonPaddyModel(n),
        variables = Wflow.DemandVariables(; n),
    )

    subsurface_flow = (; variables = (; exfiltwater_average = [0.0]))
    external_models = (; runoff, demand, subsurface_flow)

    Wflow.update_soil_water_storage!(soil_model, external_models, dt)

    @test soil_model.variables.runoff[1] ≈ 0.0
    @test soil_model.variables.unsaturated_store_capacity[1] ≈ 0.28033060970129986
    @test soil_model.variables.saturated_water_depth[1] ≈ 0.27155636944572653
    @test soil_model.variables.drainable_water_depth[1] ≈ 0.14895887025738955
    @test soil_model.variables.exfiltration_saturated_water[1] ≈ 0.0
    @test soil_model.variables.volumetric_water_content[1] ≈
          [0.2678129742498539, 0.1885939848613844, 0.2811123711333945, 0.4721985172668562]
    @test soil_model.variables.relative_volumetric_water_content[1] ≈
          [55.42837989378741, 39.03268341595556, 58.18091279434587, 97.72939072000436]
    @test soil_model.variables.volumetric_water_content_root_zone[1] ≈ 0.20944108733203426
    @test soil_model.variables.relative_volumetric_water_content_root_zone[1] ≈
          43.34734038380604
    @test soil_model.variables.total_soil_water_storage[1] ≈ 0.43856081721959095
end
