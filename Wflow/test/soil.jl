@testitem "unit: update_boundary_conditions!" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    n = 1
    model = init_sbm_soil_model(1, 1; soil_fraction = [0.3])
    atmospheric_forcing = Wflow.AtmosphericForcing(; n, potential_evaporation = [1e-3])

    interception = Wflow.GashInterceptionModel(;
        parameters = Wflow.GashParameters(;
            e_r = [0.1],
            vegetation_parameter_set = Wflow.VegetationParameters(;
                leaf_area_index = nothing,
                storage_wood = nothing,
                kext = nothing,
                storage_specific_leaf = nothing,
                canopygapfraction = [0.1],
                cmax = [10.0],
                rootingdepth = [10.0],
                kc = [0.1],
            ),
        ),
        variables = Wflow.InterceptionVariables(;
            n = 1,
            canopy_potevap = [4e-3],
            interception_rate = [2e-3],
        ),
    )
    runoff = Wflow.OpenWaterRunoff(;
        n,
        variables = Wflow.OpenWaterRunoffVariables(;
            n,
            runoff_land = [6e-4],
            runoff_river = [7e-3],
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

    Wflow.update_boundary_conditions!(model, atmospheric_forcing, external_models)

    @test model.boundary_conditions.potential_transpiration[1] ≈ 2e-3
    @test model.boundary_conditions.potential_soilevaporation[1] ≈ 2.9e-4
    @test model.boundary_conditions.water_flux_surface[1] ≈ 1.99e-2
end

@testitem "unit: unsaturated_zone_flow!" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    n = 1
    N = 6
    model = init_sbm_soil_model(
        n,
        N;
        ustorelayerthickness = [SVector((50.0, 100.0, 50.0, 200.0, 800.0, 250.0))],
        ustorelayerdepth = [SVector(8.3, 16.6, 7.7, 36.5, 190.7, 59.2)],
        n_unsatlayers = [6],
        theta_s = [0.5],
        theta_r = [0.07],
        kv_profile = Wflow.KvExponential([2530.0], [1e-3]),
        kvfrac = [SVector((0.9, 0.8, 0.9, 0.8, 0.9, 0.8))],
        c = [SVector((9.3, 9.1, 8.9, 8.7, 8.5, 8.3))],
        infiltsoilpath = [1.5],
    )
    Wflow.unsaturated_zone_flow!(model)

    @test model.variables.ustorelayerdepth[1] ≈ [
        8.86158736157703,
        17.080057700709837,
        7.851580246695791,
        36.02511190264408,
        187.18720614458448,
        59.05659877245674,
    ]
    @show model.variables.transfer[1] ≈ 4.437857871332029
end

@testitem "unit: soil_evaporation!" begin
    include("testing_utils.jl")
    n = 1
    N = 6
    model = init_sbm_soil_model(
        n,
        N;
        potential_soilevaporation = [3.9],
        ustorelayerthickness = [SVector((50.0, 100.0, 50.0, 200.0, 800.0, 250.0))],
        ustorelayerdepth = [SVector(8.3, 16.6, 7.7, 36.5, 190.7, 59.2)],
        n_unsatlayers = [1],
        zi = [40.0],
        theta_s = [0.5],
        theta_r = [0.07],
        theta_fc = [0.25],
        act_thickl = [SVector(50.0, 100.0, 50.0, 200.0, 800.0, 573.3)],
        drainable_waterdepth = [100.0],
    )

    Wflow.soil_evaporation!(model)

    @show model.variables.soilevapsat[1] ≈ 0.40360465116279065
    @show model.variables.soilevap[1] ≈ 2.2855813953488378
    @show model.variables.drainable_waterdepth[1] ≈ 99.5963953488372
end

@testitem "unit: transpiration!" begin
    include("testing_utils.jl")
    n = 1
    N = 6
    model = init_sbm_soil_model(
        n,
        N;
        h3_high = [-400.0],
        h3_low = [-100.0],
        potential_transpiration = [1.2],
        n_unsatlayers = [6],
        zi = [75.0],
        ustorelayerthickness = [SVector((50.0, 100.0, 50.0, 200.0, 800.0, 250.0))],
        ustorelayerdepth = [SVector(8.3, 16.6, 7.7, 36.5, 190.7, 59.2)],
        rootingdepth = [106.8],
        rootfraction = [SVector(0.125, 0.25, 0.125, 0.5, 0.0, 0.0)],
        act_thickl = [SVector(50.0, 100.0, 50.0, 200.0, 800.0, 573.3)],
        sumlayers = [SVector(0.0, 50.0, 150.0, 200.0, 400.0, 1200.0, 1773.3)],
        c = [SVector((9.3, 9.1, 8.9, 8.7, 8.5, 8.3))],
        theta_s = [0.5],
        theta_r = [0.07],
        hb = [-10.0],
        h1 = [0.0],
        h2 = [-100.0],
        h3 = [-915.0],
        h4 = [-15849.0],
        alpha_h1 = [1.0],
        rootdistpar = [-500.0],
        drainable_waterdepth = [100.0],
    )
    dt = 86400.0

    Wflow.transpiration!(model, dt)

    @show model.variables.ae_ustore[1] ≈ 0.4479021069963573
    @show model.variables.actevapsat[1] ≈ 0.7520978930036426
    @show model.variables.drainable_waterdepth[1] ≈ 99.24790210699636
    @show model.variables.transpiration[1] ≈ 1.2
end

@testitem "unit: capillary_flux!" begin
    include("testing_utils.jl")
    n = 1
    N = 6
    model = init_sbm_soil_model(
        n,
        N;
        rootingdepth = [106.8],
        n_unsatlayers = [6],
        kv_profile = Wflow.KvExponential([2530.0], [1e-3]),
        kvfrac = [SVector((0.9, 0.8, 0.9, 0.8, 0.9, 0.8))],
        zi = [112.0],
        ustorelayerthickness = [SVector((50.0, 100.0, 50.0, 200.0, 800.0, 250.0))],
        ustorelayerdepth = [SVector(8.3, 16.6, 7.7, 36.5, 190.7, 59.2)],
        ae_ustore = [1.25],
        ustorecapacity = [385.0],
        drainable_waterdepth = [100.0],
        cap_hmax = [2000.0],
        cap_n = [2.0],
        theta_s = [0.5],
        theta_r = [0.07],
    )

    Wflow.capillary_flux!(model)

    @show model.variables.actcapflux[1] ≈ 1.11392
end

@testitem "unit: update! SbmSoilModel" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    n = 1
    N = 6
    model = init_sbm_soil_model(
        n,
        N;
        runoff = [5.0],
        zi = [50.0],
        ustorelayerthickness = [SVector((50.0, 100.0, 50.0, 200.0, 800.0, 250.0))],
        ustorelayerdepth = [SVector(8.3, 16.6, 7.7, 36.5, 190.7, 59.2)],
        n_unsatlayers = [6],
        vwc = [SVector((0.23, 0.24, 0.22, 0.24, 0.30, 0.42))],
        vwc_perc = [SVector((44.0, 45.0, 42.0, 47.0, 58.0, 81.0))],
        nlayers = [6],
        act_thickl = [SVector(50.0, 100.0, 50.0, 200.0, 800.0, 573.3)],
        theta_s = [0.5],
        theta_r = [0.07],
        theta_fc = [0.25],
        rootingdepth = [106.8],
        sumlayers = [SVector(0.0, 50.0, 150.0, 200.0, 400.0, 1200.0, 1773.3)],
        soilthickness = [2000.0],
        soilwatercapacity = [1215.5],
        ustorecapacity = [0.0],
        satwaterdepth = [0.0],
        drainable_waterdepth = [0.0],
        total_soilwater_storage = [0.0],
        excesswater = [0.25],
        infiltexcess = [0.1],
    )

    runoff = (; variables = (; runoff_land = [0.03], ae_openw_l = [0.05]))
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

    subsurface_flow = (; exfiltwater = [0.0025])
    Wflow.get_exfiltwater(nt::NamedTuple) = nt.exfiltwater

    external_models = (; runoff, demand, subsurface_flow)

    Wflow.update!(model, external_models)

    @test demand.paddy.variables.h[1] ≈ 2.0
    @test model.variables.runoff[1] ≈ 0.88
    @test model.variables.ustorecapacity[1] ≈ 58.0
    @test model.variables.satwaterdepth[1] ≈ 838.5
    @test model.variables.drainable_waterdepth[1] ≈ 487.5
    @test model.variables.exfiltsatwater[1] ≈ 2.5
    @test model.variables.vwc[1] ≈
          [0.236, 0.236, 0.224, 0.2525, 0.308375, 0.4157509157509158]
    @test model.variables.vwc_perc[1] ≈ [47.2, 47.2, 44.8, 50.5, 61.675, 83.15018315018315]
    @test model.variables.vwc_root[1] ≈ 0.464689138576779
    @test model.variables.vwc_percroot[1] ≈ 92.93782771535581
    @test model.variables.total_soilwater_storage[1] ≈ 1157.5
end
