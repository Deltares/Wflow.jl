@testitem "unit: update_demand_gross! (NonPaddy)" begin
    include("testing_utils.jl")
    using Wflow: to_SI, Unit, MM
    using StaticArrays: SVector
    n = 1
    N = 3

    MM_PER_DT = Unit(; mm = 1, dt = -1)
    CM = Unit(; cm = 1)
    dt = 86400.0

    model = Wflow.NonPaddy(;
        n,
        parameters = Wflow.NonPaddyParameters(;
            irrigation_efficiency = [1.0],
            maximum_irrigation_rate = [to_SI(25.0, MM_PER_DT; dt_val = dt)],
            irrigation_areas = [true],
            irrigation_trigger = [true],
        ),
        variables = Wflow.NonPaddyVariables(;
            n,
            demand_gross = [to_SI(0.8604076280853505, MM_PER_DT; dt_val = dt)],
        ),
    )

    soil = init_sbm_soil_model(
        n,
        N;
        # Variables
        ustorelayerthickness = [to_SI.(SVector(50.0, 100.0, 50.0), Ref(MM))],
        ustorelayerdepth = [SVector(0.0, 0.0, 0.0)],
        n_unsatlayers = [3],
        h3 = [to_SI(-934.9109542889025, CM)],
        f_infiltration_reduction = [0.8],
        # Parameters
        maxlayers = 3,
        sumlayers = [to_SI.(SVector(0.0, 50.0, 150.0, 200.0), Ref(MM))],
        c = [SVector(9.195682525634766, 9.297739028930664, 9.597416877746582)],
        nlayers = [3],
        theta_s = [0.4417283535003662],
        theta_r = [0.09082602709531784],
        hb = [to_SI(-10.0, CM)],
        infiltcapsoil = [to_SI(334.45526123046875, MM_PER_DT; dt_val = dt)],
        pathfrac = [0.0],
        vegetation_parameter_set = Wflow.VegetationParameters(;
            rootingdepth = [to_SI(150.0, MM)],
            leaf_area_index = nothing,
            canopygapfraction = [0.1],
            cmax = [0.2],
        ),
    )

    i = 1
    k = 1

    depletion, readily_available_water = Wflow.water_demand_root_zone(soil, i, k)
    @test depletion ≈ to_SI(8.343548901322073, MM)
    @test readily_available_water ≈ to_SI(4.288641862240691, MM)
    irri_dem_gross = depletion / dt
    demand_gross = Wflow.compute_demand_gross(model, soil, irri_dem_gross, i)
    @test demand_gross ≈ to_SI(8.343548901322073, MM_PER_DT; dt_val = dt)
end

@testitem "unit: update_demand_gross! (Paddy)" begin
    using Wflow: to_SI, MM, Unit
    MM_PER_DT = Unit(; mm = 1, dt = -1)
    dt = 86400.0
    n = 1

    variables = Wflow.PaddyVariables(; n, h = [to_SI(15.0, MM)])
    parameters = Wflow.PaddyParameters(;
        irrigation_efficiency = [1.0],
        maximum_irrigation_rate = [to_SI(25.0, MM_PER_DT; dt_val = dt)],
        irrigation_areas = [true],
        irrigation_trigger = [true],
        h_min = [to_SI(20.0, MM)],
        h_opt = [to_SI(50.0, MM)],
        h_max = [to_SI(80.0, MM)],
    )
    model = Wflow.Paddy(; n, parameters, variables)
    @test Wflow.compute_irrigation_depth(model, 1) ≈ to_SI(35.0, MM)

    Wflow.update_demand_gross!(model, dt)
    @test only(variables.demand_gross) ≈ to_SI(25.0, MM_PER_DT; dt_val = dt)
end

@testitem "unit: surface_water_allocation_local!" begin
    using Wflow: to_SI, Unit
    MM_PER_DT = Unit(; mm = 1, dt = -1)
    M3_PER_DT = Unit(; m = 3, dt = -1)
    dt = 86400.0
    include("testing_utils.jl")
    n = 1

    model = Wflow.AllocationLand(;
        n,
        parameters = Wflow.AllocationLandParameters(;
            frac_sw_used = [1.0],
            areas = [600_000],
        ),
    )

    demand_variables = Wflow.DemandVariables(;
        n,
        surfacewater_demand = [to_SI(0.02, MM_PER_DT; dt_val = dt)],
    )

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

    @test river.allocation.variables.act_surfacewater_abst_vol |> only ≈
          to_SI(12.0, M3_PER_DT; dt_val = dt)
    @test river.allocation.variables.available_surfacewater |> only ≈ 288.0
    @test demand_variables.surfacewater_demand |> only ≈ 0.0
    @test river.allocation.variables.act_surfacewater_abst |> only ≈
          to_SI(0.02, MM_PER_DT; dt_val = dt)
    @test model.variables.surfacewater_alloc |> only ≈ to_SI(0.02, MM_PER_DT; dt_val = dt)
end

@testitem "unit: surface_water_allocation_area!" begin
    using Wflow: to_SI, Unit
    MM_PER_DT = Unit(; mm = 1, dt = -1)
    M3_PER_DT = Unit(; m = 3, dt = -1)
    dt = 86400.0
    include("testing_utils.jl")

    n = 3
    model = Wflow.AllocationLand(;
        n,
        parameters = Wflow.AllocationLandParameters(;
            frac_sw_used = [1.0],
            areas = [600_000.0],
        ),
    )

    demand_variables = Wflow.DemandVariables(;
        n,
        surfacewater_demand = to_SI.([0.65, 0.77, 0.331], Ref(MM_PER_DT); dt_val = dt),
    )

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
          to_SI.(
        [608.1360277590558, 167.5822285897938, 280.34784035115035],
        Ref(M3_PER_DT);
        dt_val = dt,
    )
    @test river.allocation.variables.act_surfacewater_abst ≈
          to_SI.(
        [1.008608685367754, 0.27797992757124085, 0.46503184145095416],
        Ref(MM_PER_DT);
        dt_val = dt,
    )
    @test model.variables.surfacewater_alloc ≈
          to_SI.([0.65, 0.77, 0.331], Ref(MM_PER_DT); dt_val = dt)
end

@testitem "unit: groundwater_allocation_local!" begin
    using Wflow: to_SI, Unit
    MM_PER_DT = Unit(; mm = 1, dt = -1)
    M3_PER_DT = Unit(; m = 3, dt = -1)
    dt = 86400.0

    n = 1

    model = Wflow.AllocationLand(;
        n,
        parameters = Wflow.AllocationLandParameters(; frac_sw_used = [], areas = []),
    )

    demand_variables = Wflow.DemandVariables(;
        n,
        total_gross_demand = [to_SI(100.0, MM_PER_DT; dt_val = dt)],
    )

    groundwater_storage = [315947.6]

    parameters = Wflow.LandParameters(; area = [591286.4], reservoir_coverage = [false])

    Wflow.groundwater_allocation_local!(
        model,
        demand_variables,
        groundwater_storage,
        parameters,
        dt,
    )

    @test model.variables.act_groundwater_abst_vol |> only ≈
          to_SI(59128.64, M3_PER_DT; dt_val = dt)
    @test model.variables.available_groundwater |> only ≈ 177832.06
    @test demand_variables.groundwater_demand |> only == 0.0
    @test model.variables.act_groundwater_abst |> only ≈
          to_SI(100.0, MM_PER_DT; dt_val = dt)
    @test model.variables.groundwater_alloc |> only ≈ to_SI(100.0, MM_PER_DT; dt_val = dt)
end

@testitem "unit: groundwater_allocation_area!" begin
    using Wflow: to_SI, Unit
    MM_PER_DT = Unit(; mm = 1, dt = -1)
    M3_PER_DT = Unit(; m = 3, dt = -1)
    dt = 86400.0

    n = 3

    model = Wflow.AllocationLand(;
        n,
        parameters = Wflow.AllocationLandParameters(; frac_sw_used = [], areas = []),
        variables = Wflow.AllocationLandVariables(;
            n,
            available_groundwater = [303505.6, 308331.8, 306516.2],
        ),
    )

    demand_variables = Wflow.DemandVariables(;
        n,
        groundwater_demand = to_SI.(
            [23.23450, 12.261, 674.32],
            Ref(MM_PER_DT);
            dt_val = dt,
        ),
    )

    domain = Wflow.Domain(;
        land = Wflow.DomainLand(;
            parameters = Wflow.LandParameters(; area = [603121.0, 603121.0, 603121.0]),
            network = Wflow.NetworkLand(; allocation_area_indices = [[1, 2, 3]]),
        ),
        river = Wflow.DomainRiver(;
            network = Wflow.NetworkRiver(; allocation_area_indices = [[1]]),
        ),
    )

    Wflow.groundwater_allocation_area!(model, demand_variables, domain, dt)

    @test model.variables.act_groundwater_abst_vol ≈
          to_SI.(
        [141483.796500842, 143733.60374878853, 142887.23392586954],
        Ref(M3_PER_DT);
        dt_val = dt,
    )
    @test model.variables.act_groundwater_abst ≈
          to_SI.(
        [234.58608886250354, 238.31636396144145, 236.91304717605513],
        Ref(MM_PER_DT);
        dt_val = dt,
    )
    @test model.variables.groundwater_alloc ≈
          to_SI.([23.2345, 12.261, 674.32], Ref(MM_PER_DT); dt_val = dt)
end
