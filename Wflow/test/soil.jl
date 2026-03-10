@testitem "unit: update_bc_soil_model!" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    n = 1
    soil_model = init_sbm_soil_model(1, 1; soil_fraction = [0.30670667824303055])
    atmospheric_forcing =
        Wflow.AtmosphericForcing(; n, potential_evaporation = [4.8000001907348635])

    interception = Wflow.GashInterceptionModel(;
        parameters = Wflow.GashParameters(;
            e_r = [0.1],
            vegetation_parameter_set = Wflow.VegetationParameters(;
                leaf_area_index = nothing,
                storage_wood = nothing,
                kext = nothing,
                storage_specific_leaf = nothing,
                canopygapfraction = [0.3487189230509198],
                cmax = [0.2235326728960274],
                rootingdepth = [390.531982421875],
                kc = [1.1469999551773071],
            ),
        ),
        variables = Wflow.InterceptionVariables(;
            n = 1,
            canopy_potevap = [3.585693099611687],
            interception_rate = [0.12944592473374397],
        ),
    )
    runoff = Wflow.OpenWaterRunoff(;
        n,
        variables = Wflow.OpenWaterRunoffVariables(;
            n,
            runoff_land = [0.0],
            runoff_river = [0.003384257254905341],
        ),
        boundary_conditions = Wflow.OpenWaterRunoffBC(; n, water_flux_surface = [7.5e-3]),
    )
    demand = Wflow.Demand(;
        domestic = Wflow.NoDemand(; n),
        industry = Wflow.NoDemand(; n),
        livestock = Wflow.NoDemand(; n),
        paddy = Wflow.Paddy(;
            parameters = Wflow.PaddyParameters(;
                irrigation_efficiency = [],
                maximum_irrigation_rate = [],
                irrigation_areas = [true],
                irrigation_trigger = [],
                h_min = [],
                h_opt = [],
                h_max = [],
            ),
            variables = Wflow.PaddyVariables(; n, h = [1e-5], evaporation = [4e-3]),
        ),
        nonpaddy = Wflow.NoIrrigationNonPaddy(n),
        variables = Wflow.DemandVariables(; n),
    )

    allocation = Wflow.AllocationLand(;
        n,
        parameters = Wflow.AllocationLandParameters(; frac_sw_used = [], areas = []),
        variables = Wflow.AllocationLandVariables(; n, irri_alloc = [2e-2]),
    )

    external_models = (; interception, runoff, demand, allocation)

    Wflow.update_bc_soil_model!(soil_model, atmospheric_forcing, external_models)

    @test soil_model.boundary_conditions.potential_transpiration[1] ≈ 3.456247174877943
    @test soil_model.boundary_conditions.potential_soilevaporation[1] ≈ 1.472182114066203
    @test soil_model.boundary_conditions.water_flux_surface[1] ≈ 0.02411574274509466
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
    model = init_sbm_soil_model(
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

    Wflow.soil_evaporation!(model)

    @test model.variables.soilevapsat[1] == 0.0
    @test model.variables.soilevap[1] ≈ 0.2459598877386006
    @test model.variables.drainable_waterdepth[1] ≈ 321.13323174500624
end

@testitem "unit: transpiration!" begin
    include("testing_utils.jl")
    n = 1
    N = 6
    soil_model = init_sbm_soil_model(
        n,
        N;
        h3_high = [-400.0],
        h3_low = [-100.0],
        potential_transpiration = [0.0],
        n_unsatlayers = [5],
        zi = [586.4715307142369],
        ustorelayerthickness = [
            SVector((50.0, 100.0, 50.0, 200.0, 186.4715307142369, NaN)),
        ],
        ustorelayerdepth = [
            SVector(
                10.095952041768657,
                21.371723437347292,
                9.834201105799886,
                8.50563731464462,
                34.84891547900849,
                0.0,
            ),
        ],
        rootingdepth = [403.322998046875],
        rootfraction = [SVector(0.124, 0.25, 0.125, 0.495, 0.008, 0.0)],
        act_thickl = [SVector(50.0, 100.0, 50.0, 200.0, 200.0, NaN)],
        sumlayers = [SVector(0.0, 50.0, 150.0, 200.0, 400.0, 600.0, NaN)],
        c = [SVector((8.91, 8.24, 8.67, 8.62, 8.49, 8.43))],
        theta_s = [0.5405],
        theta_r = [0.06145],
        hb = [-10.0],
        h1 = [0.0],
        h2 = [-100.0],
        h3 = [-755.506],
        h4 = [-15849.0],
        alpha_h1 = [1.0],
        rootdistpar = [-500.0],
        drainable_waterdepth = [3.6183247],
    )
    dt = 86400.0

    Wflow.transpiration!(soil_model, dt)

    @test soil_model.variables.ae_ustore[1] ≈ 0.0
    @test soil_model.variables.actevapsat[1] ≈ 0.0
    @test soil_model.variables.drainable_waterdepth[1] ≈ 3.6183247
    @test soil_model.variables.transpiration[1] ≈ 0.0
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
    n = 1
    N = 6
    soil_model = init_sbm_soil_model(
        n,
        N;
        runoff = [0.002036245],
        zi = [0.000656547],
        ustorelayerthickness = [SVector((0.000656547, NaN, NaN, NaN, NaN, NaN))],
        ustorelayerdepth = [SVector(1.998965544e-5, 0.0, 0.0, 0.0, 0.0, 0.0)],
        n_unsatlayers = [1],
        vwc = [SVector((0.440136, 0.440140, 0.440140, 0.440140, 0.440140, 0.440140))],
        vwc_perc = [SVector((99.99912, 100.0, 100.0, 100.0, 100.0, 100.0))],
        nlayers = [6],
        act_thickl = [SVector(50.0, 100.0, 50.0, 200.0, 800.0, 800.0)],
        theta_s = [0.44014],
        theta_r = [0.0880263],
        theta_fc = [0.26583],
        rootingdepth = [390.39999],
        sumlayers = [SVector(0.0, 50.0, 150.0, 200.0, 400.0, 1200.0, 2000.0)],
        soilthickness = [2000.0],
        soilwatercapacity = [704.2274028],
        ustorecapacity = [0.0001915933],
        satwaterdepth = [704.227211212],
        drainable_waterdepth = [348.61417105],
        total_soilwater_storage = [704.227211212],
        excesswater = [0.047112555773],
        infiltexcess = [0.0],
    )

    runoff = (; variables = (; runoff_land = [0.0], ae_openw_l = [0.0]))
    demand = Wflow.Demand(;
        domestic = Wflow.NoDemand(; n),
        industry = Wflow.NoDemand(; n),
        livestock = Wflow.NoDemand(; n),
        paddy = Wflow.Paddy(;
            parameters = Wflow.PaddyParameters(;
                irrigation_efficiency = [],
                maximum_irrigation_rate = [],
                irrigation_areas = [true],
                irrigation_trigger = [],
                h_min = [],
                h_opt = [],
                h_max = [2.0],
            ),
            variables = Wflow.PaddyVariables(; n, h = [1e-5], evaporation = [4e-3]),
        ),
        nonpaddy = Wflow.NoIrrigationNonPaddy(n),
        variables = Wflow.DemandVariables(; n),
    )

    subsurface_flow = (; exfiltwater = [0.0])
    Wflow.get_exfiltwater(nt::NamedTuple) = nt.exfiltwater

    external_models = (; runoff, demand, subsurface_flow)

    Wflow.update_soil_water_storage!(soil_model, external_models)

    @test demand.paddy.variables.h[1] ≈ 0.047112555773
    @test soil_model.variables.runoff[1] ≈ 0.0
    @test soil_model.variables.ustorecapacity[1] ≈ 0.00021398953802734852
    @test soil_model.variables.satwaterdepth[1] ≈ 704.2271688208066
    @test soil_model.variables.drainable_waterdepth[1] ≈ 348.6198855572923
    @test soil_model.variables.exfiltsatwater[1] ≈ 0.0
    @test soil_model.variables.vwc[1] ≈
          [0.4401357762092409, 0.44014, 0.44014, 0.44014, 0.44014, 0.44014]
    @test soil_model.variables.vwc_perc[1] ≈
          [99.99904035289701, 100.0, 100.0, 100.0, 100.0, 100.0]
    @test soil_model.variables.vwc_root[1] ≈ 0.44013945904317786
    @test soil_model.variables.vwc_percroot[1] ≈ 99.99987709437403
    @test soil_model.variables.total_soilwater_storage[1] ≈ 704.227188810462
end
