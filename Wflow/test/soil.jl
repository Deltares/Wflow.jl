@testitem "unit: update_bc_soil_model!" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    n_land_cells = 1
    soil_model = init_sbm_soil_model(1, 1; soil_fraction = [0.2397957498236932])
    atmospheric_forcing = Wflow.AtmosphericForcing(;
        n_land_cells,
        potential_evaporation = [0.5799999833106995],
    )

    interception = Wflow.GashInterceptionModel(;
        parameters = Wflow.GashParameters(;
            e_r = [0.25],
            vegetation_parameter_set = Wflow.VegetationParameters(;
                leaf_area_index = [1.5],
                storage_wood = [0.2],
                kext = [0.67],
                storage_specific_leaf = [0.09],
                canopygapfraction = [0.2397957498236932],
                cmax = [0.5687566481997495],
                rootingdepth = [410.0],
                kc = [1.0],
            ),
        ),
        variables = Wflow.InterceptionVariables(;
            n_land_cells,
            canopy_potevap = [0.44091845241498073],
            interception_rate = [0.029448986349521342],
        ),
    )
    runoff = Wflow.OpenWaterRunoff(;
        n_land_cells,
        variables = Wflow.OpenWaterRunoffVariables(;
            n_land_cells,
            runoff_land = [0.0],
            runoff_river = [0.0],
        ),
        boundary_conditions = Wflow.OpenWaterRunoffBC(;
            n_land_cells,
            water_flux_surface = [0.24686226660437197],
        ),
    )
    demand = Wflow.DemandModel(;
        domestic = Wflow.NoDemandModel(; n_land_cells),
        industry = Wflow.NoDemandModel(; n_land_cells),
        livestock = Wflow.NoDemandModel(; n_land_cells),
        paddy = Wflow.NoIrrigationPaddyModel(n_land_cells),
        nonpaddy = Wflow.NoIrrigationNonPaddyModel(n_land_cells),
        variables = Wflow.DemandVariables(; n_land_cells),
    )

    allocation = Wflow.NoAllocationLandModel(n_land_cells)

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
        ustorelayerthickness = [SVector((50.0, 5.081648613929929, NaN, NaN, NaN, NaN))],
        ustorelayerdepth = [
            SVector(1.2855527211118947, 0.20814868098590805, 0.0, 0.0, 0.0, 0.0),
        ],
        n_unsatlayers = [2],
        theta_s = [0.4414711594581604],
        theta_r = [0.08942600339651108],
        kv_profile = Wflow.KvExponential([364.1460876464844], [0.0033079576678574085]),
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
        infiltsoilpath = [0.02277405542765069],
    )
    Wflow.unsaturated_zone_flow!(soil_model)

    @test soil_model.variables.ustorelayerdepth[1] ≈
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
        ustorelayerthickness = [SVector((50.0, 21.472680450878443, NaN, NaN, NaN, NaN))],
        ustorelayerdepth = [
            SVector(1.537249298366254, 4.821326813899425e-8, 0.0, 0.0, 0.0, 0.0),
        ],
        n_unsatlayers = [2],
        zi = [71.47268045087844],
        theta_s = [0.44],
        theta_r = [0.09],
        theta_fc = [0.275],
        act_thickl = [SVector(50.0, 100.0, 50.0, 200.0, 800.0, 800.0)],
        drainable_waterdepth = [321.13323174500624],
    )

    Wflow.soil_evaporation!(soil_model)

    @test soil_model.variables.soilevapsat[1] == 0.0
    @test soil_model.variables.soilevap[1] ≈ 0.2459598877386006
    @test soil_model.variables.drainable_waterdepth[1] ≈ 321.13323174500624
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
        zi = [106.89587841733061],
        ustorelayerthickness = [SVector((100.0, 6.895878417330607, NaN, NaN))],
        ustorelayerdepth = [SVector(10.9327977152876, 0.8620432154993639, 0.0, 0.0)],
        rootingdepth = [453.0],
        rootfraction = [
            SVector(0.22075055187637968, 0.6622516556291391, 0.11699779249448124, 0.0),
        ],
        act_thickl = [SVector(100.0, 300.0, 200.0, NaN)],
        sumlayers = [SVector(0.0, 100.0, 400.0, 600.0, NaN)],
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
        hb = [-10.0],
        h1 = [0.0],
        h2 = [-100.0],
        h3 = [1000.0],
        h4 = [-16000.0],
        alpha_h1 = [1.0],
        rootdistpar = [-500.0],
        drainable_waterdepth = [72.40310797113221],
    )
    dt = 86400.0

    Wflow.transpiration!(soil_model, dt)

    @test soil_model.variables.ae_ustore[1] ≈ 0.05153840858869719
    @test soil_model.variables.actevapsat[1] ≈ 0.0
    @test soil_model.variables.drainable_waterdepth[1] ≈ 72.40310797113221
    @test soil_model.variables.transpiration[1] ≈ 0.05153840858869719
end

@testitem "unit: capillary_flux!" begin
    include("testing_utils.jl")
    n = 1
    N = 6
    soil_model = init_sbm_soil_model(
        n,
        N;
        rootingdepth = [384.1],
        n_unsatlayers = [6],
        kv_profile = Wflow.KvExponential([2018.0371], [1.29274e-3]),
        kvfrac = [SVector((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))],
        zi = [1266.33589],
        ustorelayerthickness = [SVector((50.0, 100.0, 50.0, 200.0, 800.0, 66.33589243))],
        ustorelayerdepth = [
            SVector((
                8.874129508377953,
                18.18721052056329,
                8.54476050597162,
                17.279406870498256,
                147.78493107547982,
                12.51769175241036,
            )),
        ],
        ae_ustore = [0.528773283329],
        ustorecapacity = [313.1741519821792],
        drainable_waterdepth = [169.02603585525694],
        cap_hmax = [2000.0],
        cap_n = [2.0],
        theta_s = [0.4868114888668],
        theta_r = [0.0711537748575],
    )

    Wflow.capillary_flux!(soil_model)

    @test soil_model.variables.actcapflux[1] ≈ 0.07115477692809025
end

@testitem "unit: update_soil_water_storage! SbmSoilModel" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    n_land_cells = 1
    N = 4
    soil_model = init_sbm_soil_model(
        n_land_cells,
        N;
        runoff = [0.0],
        zi = [1244.5135404970033],
        ustorelayerthickness = [SVector((100.0, 300.0, 800.0, 44.513540497003305))],
        ustorelayerdepth = [
            SVector(
                14.408928105784874,
                19.461087500813772,
                125.91094235311145,
                7.223489814154271,
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
        act_thickl = [SVector(100.0, 300.0, 800.0, 800.0)],
        theta_s = [0.4831694066524056],
        theta_r = [0.12372369319200516],
        theta_fc = [0.28599992944511926],
        rootingdepth = [380.0],
        sumlayers = [SVector(0.0, 100.0, 400.0, 1200.0, 2000.0)],
        soilthickness = [2000.0],
        soilwatercapacity = [718.8914269208908],
        ustorecapacity = [271.02508737207194],
        satwaterdepth = [288.62821979680706],
        drainable_waterdepth = [158.32342151683937],
        total_soilwater_storage = [448.5154852374101],
        excesswater = [0.0],
        infiltexcess = [0.0],
    )

    runoff = (; variables = (; runoff_land = [0.0], ae_openw_l = [0.0]))
    demand = Wflow.DemandModel(;
        domestic = Wflow.NoDemandModel(; n_land_cells),
        industry = Wflow.NoDemandModel(; n_land_cells),
        livestock = Wflow.NoDemandModel(; n_land_cells),
        paddy = Wflow.NoIrrigationPaddyModel(n_land_cells),
        nonpaddy = Wflow.NoIrrigationNonPaddyModel(n_land_cells),
        variables = Wflow.DemandVariables(; n_land_cells),
    )

    subsurface_flow = (; variables = (; exfiltwater = [0.0]))
    external_models = (; runoff, demand, subsurface_flow)

    Wflow.update_soil_water_storage!(soil_model, external_models)

    @test soil_model.variables.runoff[1] ≈ 0.0
    @test soil_model.variables.ustorecapacity[1] ≈ 280.33060970129986
    @test soil_model.variables.satwaterdepth[1] ≈ 271.55636944572655
    @test soil_model.variables.drainable_waterdepth[1] ≈ 148.95887025738955
    @test soil_model.variables.exfiltsatwater[1] ≈ 0.0
    @test soil_model.variables.vwc[1] ≈
          [0.2678129742498539, 0.1885939848613844, 0.2811123711333945, 0.4721985172668562]
    @test soil_model.variables.vwc_perc[1] ≈
          [55.42837989378741, 39.03268341595556, 58.18091279434587, 97.72939072000436]
    @test soil_model.variables.vwc_root[1] ≈ 0.20944108733203426
    @test soil_model.variables.vwc_percroot[1] ≈ 43.34734038380604
    @test soil_model.variables.total_soilwater_storage[1] ≈ 438.56081721959094
end
