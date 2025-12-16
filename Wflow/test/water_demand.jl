@testitem "unit: update_demand_gross (NonPaddy)" begin
    n_unsatlayers = 1
    rootingdepth = 396.8680114746094
    sumlayers = [0.0, 50.0, 150.0, 200.0, 400.0, 1200.0, 2000.0]
    ustorelayerthickness = [2.3346668393934804, NaN, NaN, NaN, NaN, NaN]
    ustorelayerdepth = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    hb = -10.0
    theta_s = 0.4417283535003662
    theta_r = 0.09082602709531784
    c = [
        9.195682525634766,
        9.297739028930664,
        9.597416877746582,
        9.782596588134766,
        9.974676132202148,
        9.876368522644043,
    ]
    h3 = -934.9109542889025
    maximum_irrigation_rate = 25.0
    f_infiltration_reduction = 1.0
    pathfrac = 0.0
    infiltcapsoil = 334.45526123046875
    irrigation_efficiency = 1.0
    demand_gross = 0.8604076280853505

    @test Wflow.update_demand_gross(
        n_unsatlayers,
        rootingdepth,
        sumlayers,
        ustorelayerthickness,
        ustorelayerdepth,
        hb,
        theta_s,
        theta_r,
        c,
        h3,
        maximum_irrigation_rate,
        f_infiltration_reduction,
        pathfrac,
        infiltcapsoil,
        irrigation_efficiency,
        demand_gross,
    ) ≈ 0.38958813885549093
end

@testitem "unit: update_demand_gross (Paddy)" begin
    irrigation_efficiency = 1.0
    maximum_irrigation_rate = 25.0
    h_opt = 50.0
    h_min = 20.0
    h = 15.0
    demand_gross = 3.0

    @test Wflow.update_demand_gross(
        irrigation_efficiency,
        maximum_irrigation_rate,
        h_opt,
        h_min,
        h,
        demand_gross,
    ) ≈ 25.0
end

@testitem "unit: surface_water_allocation_local" begin
    surfacewater_demand = 0.02
    area = 600_000.0
    external_inflow = 5.0
    storage = 375.0
    dt = 86400.0

    abstraction_vol,
    avail_surfacewater,
    surfacewater_demand,
    act_surfacewater_abst,
    surfacewater_alloc = Wflow.surface_water_allocation_local(
        surfacewater_demand,
        area,
        external_inflow,
        storage,
        dt,
    )
    @test abstraction_vol ≈ 12.0
    @test avail_surfacewater ≈ 288.0
    @test surfacewater_demand ≈ 0.0
    @test act_surfacewater_abst == surfacewater_alloc ≈ 0.02
end

@testitem "unit: surface_water_allocation_area!" begin end

@testitem "unit: groundwater_allocation_local" begin end

@testitem "unit: groundwater_allocation_area!" begin end
