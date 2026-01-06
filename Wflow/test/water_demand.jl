@testitem "unit: update_demand_gross! (NonPaddy)" begin
    include("testing_utils.jl")
    using StaticArrays: SVector
    n = 1
    N = 6

    model = Wflow.NonPaddy(;
        parameters = Wflow.NonPaddyParameters(;
            irrigation_efficiency = [1.0],
            maximum_irrigation_rate = [25.0],
            irrigation_areas = [true],
            irrigation_trigger = [true],
        ),
        variables = Wflow.NonPaddyVariables(; n, demand_gross = [0.8604076280853505]),
    )

    soil = simple_soil(
        n,
        N;
        # Variables
        ustorelayerthickness = [SVector(2.3346668393934804, NaN, NaN, NaN, NaN, NaN)],
        ustorelayerdepth = [SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)],
        n_unsatlayers = [3],
        h3 = [-934.9109542889025],
        f_infiltration_reduction = [0.8],
        # Parameters
        maxlayers = 7,
        sumlayers = [SVector(0.0, 50.0, 150.0, 200.0, 400.0, 1200.0, 2000.0)],
        c = [
            SVector(
                9.195682525634766,
                9.297739028930664,
                9.597416877746582,
                9.782596588134766,
                9.974676132202148,
                9.876368522644043,
            ),
        ],
        nlayers = [7],
        theta_s = [0.4417283535003662],
        theta_r = [0.09082602709531784],
        hb = [-10.0],
        infiltcapsoil = [334.45526123046875],
        pathfrac = [0.0],
        vegetation_parameter_set = Wflow.VegetationParameters(;
            rootingdepth = [396.8680114746094],
        ),
    )

    i = 1
    k = 1

    depletion, readily_available_water = Wflow.water_demand_root_zone(soil, i, k)
    @test depletion ≈ 0.38958813885549093
    @test readily_available_water ≈ 0.2002509988361609
    irri_dem_gross = depletion
    demand_gross = Wflow.apply_infiltration_capacity(model, soil, irri_dem_gross, i)
    @test demand_gross ≈ 0.38958813885549093
end

@testitem "unit: update_demand_gross! (Paddy)" begin
    variables = Wflow.PaddyVariables(; n = 1, h = [15.0])
    parameters = Wflow.PaddyParameters(;
        irrigation_efficiency = [1.0],
        maximum_irrigation_rate = [25.0],
        irrigation_areas = [true],
        irrigation_trigger = [true],
        h_min = [20.0],
        h_opt = [50.0],
        h_max = [NaN],
    )
    model = Wflow.Paddy(; parameters, variables)
    @test Wflow.calc_demand_gross(model, 1) ≈ 35.0

    Wflow.update_demand_gross!(model)
    @test only(variables.demand_gross) ≈ 25.0
end

@testitem "unit: surface_water_allocation_local!" begin
    include("testing_utils.jl")
    n = 1

    model = Wflow.AllocationLand(;
        n,
        parameters = Wflow.AllocationLandParameters(;
            frac_sw_used = [NaN],
            areas = [600_000.0],
        ),
    )

    demand_variables = Wflow.DemandVariables(; n, surfacewater_demand = [0.02])

    river = DummyRiver(;
        allocation = (;
            variables = (;
                act_surfacewater_abst_vol = zeros(n),
                act_surfacewater_abst = zeros(n),
                available_surfacewater = zeros(n),
            )
        ),
        boundary_conditions = (; external_inflow = [5.0]),
        variables = (; storage = [375.0]),
    )

    domain = Wflow.DomainLand(;
        network = Wflow.NetworkLand(; river_inds_excl_reservoir = [1]),
        parameters = Wflow.LandParameters(; area = [600_000.0]),
    )

    dt = 86400.0

    Wflow.surface_water_allocation_local!(model, demand_variables, river, domain, dt)

    @test river.allocation.variables.act_surfacewater_abst_vol |> only ≈ 12.0
    @test river.allocation.variables.available_surfacewater |> only ≈ 288.0
    @test demand_variables.surfacewater_demand |> only ≈ 0.0
    @test river.allocation.variables.act_surfacewater_abst |> only ≈ 0.02
    @test model.variables.surfacewater_alloc |> only ≈ 0.02
end

@testitem "unit: surface_water_allocation_area!" begin
    include("testing_utils.jl")

    n = 3
    model = Wflow.AllocationLand(;
        n,
        parameters = Wflow.AllocationLandParameters(; frac_sw_used = [], areas = []),
    )

    demand_variables = Wflow.DemandVariables(; n, surfacewater_demand = [0.65, 0.77, 0.331])

    river = DummyRiver(;
        allocation = (;
            variables = (;
                act_surfacewater_abst_vol = zeros(n),
                act_surfacewater_abst = zeros(n),
                available_surfacewater = zeros(n),
            )
        ),
        boundary_conditions = (;
            reservoir = (;
                boundary_conditions = (; external_inflow = [2.0, 3.0, 4.0]),
                variables = (; storage = [1.5531612276024342e8, 4.28e7, 7.16e7]),
            )
        ),
    )

    domain = Wflow.Domain(;
        land = Wflow.DomainLand(;
            network = Wflow.NetworkLand(; allocation_area_indices = [[1, 2, 3]]),
            parameters = Wflow.LandParameters(; area = [603121.7, 603121.7, 603121.7]),
        ),
        river = Wflow.DomainRiver(;
            network = Wflow.NetworkRiver(;
                allocation_area_indices = [[1, 2, 3]],
                reservoir_indices = [1, 2, 3],
            ),
            parameters = Wflow.RiverParameters(;
                cell_area = [602945.46, 602857.3, 602857.3],
            ),
        ),
    )

    dt = 86400.0

    @test Wflow.available_surface_water!(
        river.allocation.variables.available_surfacewater,
        river.boundary_conditions.reservoir,
        domain.river.network.allocation_area_indices[1],
        domain.river.network.reservoir_indices,
        dt,
    ) ≈ 2.6432180030503854e8

    Wflow.surface_water_allocation_area!(model, demand_variables, river, domain, dt)

    @test river.allocation.variables.available_surfacewater ≈
          [1.5220980030503854e8, 4.1944e7, 7.0168e7]
    @test river.allocation.variables.act_surfacewater_abst_vol ≈
          [608.1360277590558, 167.5822285897938, 280.34784035115035]
    @test river.allocation.variables.act_surfacewater_abst ≈
          [1.008608685367754, 0.27797992757124085, 0.46503184145095416]
    @test model.variables.surfacewater_alloc ≈ [0.65, 0.77, 0.331]
end

@testitem "unit: groundwater_allocation_local" begin
    groundwater_demand = 100.0
    area = 591286.4
    groundwater_storage = 315947.6

    abstraction_vol,
    available_groundwater,
    groundwater_demand,
    act_groundwater_abst,
    groundwater_alloc =
        Wflow.groundwater_allocation_local(groundwater_demand, area, groundwater_storage)

    @test abstraction_vol ≈ 59128.64
    @test available_groundwater ≈ 177832.06
    @test groundwater_demand == 0.0
    @test act_groundwater_abst ≈ 100.0
    @test groundwater_alloc ≈ 100.0
end

@testitem "unit: groundwater_allocation_area!" begin
    act_groundwater_abst_vol = zeros(3)
    act_groundwater_abst = zeros(3)
    groundwater_alloc = zeros(3)
    groundwater_demand = [23.23450, 12.261, 674.32]
    area = [603121.0, 603121.0, 603121.0]
    available_groundwater = [303505.6, 308331.8, 306516.2]
    inds_land = [1, 2, 3]

    Wflow.groundwater_allocation_area!(
        act_groundwater_abst_vol,
        act_groundwater_abst,
        groundwater_alloc,
        groundwater_demand,
        area,
        available_groundwater,
        inds_land,
    )

    @test act_groundwater_abst_vol ≈
          [141483.796500842, 143733.60374878853, 142887.23392586954]
    @test act_groundwater_abst ≈
          [234.58608886250354, 238.31636396144145, 236.91304717605513]
    @test groundwater_alloc ≈ [23.2345, 12.261, 674.32]
end
