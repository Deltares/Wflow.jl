@testitem "unit: update_demand_gross! (NonPaddyModel)" begin
    include("testing_utils.jl")
    using StaticArrays: SVector
    n = 1
    N = 3

    dt = 86400.0

    model = Wflow.NonPaddyModel(;
        n,
        parameters = Wflow.NonPaddyParameters(;
            irrigation_efficiency = [1.0],
            maximum_irrigation_rate = [2.8935185185185185e-7],
            irrigation_areas = [true],
            irrigation_trigger = [true],
        ),
        variables = Wflow.NonPaddyVariables(; n, demand_gross = [9.958421621358223e-9]),
    )

    soil = init_sbm_soil_model(
        n,
        N;
        # Variables
        unsaturated_layer_thickness = [SVector(0.05, 0.1, 0.05)],
        unsaturated_layer_depth = [SVector(0.0, 0.0, 0.0)],
        n_unsatlayers = [3],
        h3 = [-9.349109542889025],
        f_infiltration_reduction = [0.8],
        # Parameters
        maximum_number_of_layers = 3,
        cumulative_layer_depth = [SVector(0.0, 0.05, 0.15, 0.2)],
        brooks_corey_exponent = [
            SVector(9.195682525634766, 9.297739028930664, 9.597416877746582),
        ],
        number_of_layers = [3],
        theta_s = [0.4417283535003662],
        theta_fc = [0.26063963369395],
        theta_r = [0.09082602709531784],
        air_entry_pressure = [-0.1],
        infiltration_capacity_soil = [3.87100996794524e-6],
        compacted_soil_area_fraction = [0.0],
        vegetation_parameter_set = Wflow.VegetationParameters(;
            n,
            rooting_depth = [0.15],
            leaf_area_index = nothing,
            canopy_gap_fraction = [0.1],
            maximum_canopy_storage = [0.2],
        ),
    )

    i = 1
    k = 1

    depletion, readily_available_water = Wflow.water_demand_root_zone(soil, i, k)
    @test depletion ≈ 0.008490680329931607
    @test readily_available_water ≈ 0.004435773290850226
    irri_dem_gross = depletion / dt
    demand_gross = Wflow.compute_demand_gross(model, soil, irri_dem_gross, i)
    @test demand_gross ≈ 9.827176307791211e-8
end

@testitem "unit: update_demand_gross! (PaddyModel)" begin
    dt = 86400.0
    n = 1

    variables = Wflow.PaddyVariables(; n, h = [0.015])
    parameters = Wflow.PaddyParameters(;
        irrigation_efficiency = [1.0],
        maximum_irrigation_rate = [2.8935185185185185e-7],
        irrigation_areas = [true],
        irrigation_trigger = [true],
        h_min = [0.02],
        h_opt = [0.05],
        h_max = [0.08],
    )
    paddy_model = Wflow.PaddyModel(; n, parameters, variables)
    @test Wflow.compute_irrigation_depth(paddy_model, 1) ≈ 0.035

    Wflow.update_demand_gross!(paddy_model, dt)
    @test only(variables.demand_gross) ≈ 2.8935185185185185e-7
end

@testitem "unit: surface_water_allocation_local!" begin
    dt = 86400.0
    include("testing_utils.jl")
    n = 1

    allocation_model = Wflow.AllocationLandModel(;
        n,
        parameters = Wflow.AllocationLandParameters(;
            fraction_surfacewater_used = [1.0],
            areas = [600_000],
        ),
    )

    demand_variables =
        Wflow.DemandVariables(; n, surfacewater_demand = [2.3148148148148147e-10])

    river = DummyRiver(;
        allocation = (;
            variables = (;
                actual_surfacewater_abstraction_volume = zeros(n),
                actual_surfacewater_abstraction = zeros(n),
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

    Wflow.surface_water_allocation_local!(
        allocation_model,
        demand_variables,
        river,
        domain,
        dt,
    )

    @test river.allocation.variables.actual_surfacewater_abstraction_volume |> only ≈
          0.0001388888888888889
    @test river.allocation.variables.available_surfacewater |> only ≈ 288.0
    @test demand_variables.surfacewater_demand |> only ≈ 0.0
    @test river.allocation.variables.actual_surfacewater_abstraction |> only ≈
          2.3148148148148147e-10
    @test allocation_model.variables.surfacewater_allocation |> only ≈
          2.3148148148148147e-10
end

@testitem "unit: surface_water_allocation_area!" begin
    dt = 86400.0
    include("testing_utils.jl")

    n = 3
    allocation_model = Wflow.AllocationLandModel(;
        n,
        parameters = Wflow.AllocationLandParameters(;
            fraction_surfacewater_used = [1.0],
            areas = [600_000.0],
        ),
    )

    demand_variables = Wflow.DemandVariables(;
        n,
        surfacewater_demand = [
            7.523148148148148e-9,
            8.912037037037038e-9,
            3.831018518518519e-9,
        ],
    )

    river_flow_model = DummyRiver(;
        allocation = (;
            variables = (;
                actual_surfacewater_abstraction_volume = zeros(n),
                actual_surfacewater_abstraction = zeros(n),
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
        river_flow_model.allocation.variables.available_surfacewater,
        river_flow_model.boundary_conditions.reservoir,
        domain.river.network.allocation_area_indices[1],
        domain.river.network.reservoir_indices,
        dt,
    ) ≈ 2.6432180030503854e8

    Wflow.surface_water_allocation_area!(
        allocation_model,
        demand_variables,
        river_flow_model,
        domain,
        dt,
    )

    @test river_flow_model.allocation.variables.available_surfacewater ≈
          [1.5220980030503854e8, 4.1944e7, 7.0168e7]
    @test river_flow_model.allocation.variables.actual_surfacewater_abstraction_volume ≈
          [0.007038611432396478, 0.0019396091271966873, 0.0032447666707309066]
    @test river_flow_model.allocation.variables.actual_surfacewater_abstraction ≈
          [1.1673711636200856e-8, 3.2173602728152876e-9, 5.382312979756414e-9]
    @test allocation_model.variables.surfacewater_allocation ≈
          [7.523148148148148e-9, 8.912037037037038e-9, 3.831018518518519e-9]
end

@testitem "unit: groundwater_allocation_local!" begin
    dt = 86400.0

    n = 1

    allocation_model = Wflow.AllocationLandModel(;
        n,
        parameters = Wflow.AllocationLandParameters(;
            fraction_surfacewater_used = [],
            areas = [],
        ),
    )

    demand_variables =
        Wflow.DemandVariables(; n, total_gross_demand = [1.1574074074074074e-6])

    groundwater_storage = [315947.6]

    parameters = Wflow.LandParameters(; area = [591286.4], reservoir_coverage = [false])

    Wflow.groundwater_allocation_local!(
        allocation_model,
        demand_variables,
        groundwater_storage,
        parameters,
        dt,
    )

    @test allocation_model.variables.actual_groundwater_abstraction_volume |> only ≈
          0.6843592592592592
    @test allocation_model.variables.available_groundwater |> only ≈ 177832.06
    @test demand_variables.groundwater_demand |> only == 0.0
    @test allocation_model.variables.actual_groundwater_abstraction |> only ≈
          1.1574074074074074e-6
    @test allocation_model.variables.groundwater_allocation |> only ≈ 1.1574074074074074e-6
end

@testitem "unit: groundwater_allocation_area!" begin
    dt = 86400.0

    n = 3

    allocation_model = Wflow.AllocationLandModel(;
        n,
        parameters = Wflow.AllocationLandParameters(;
            fraction_surfacewater_used = [],
            areas = [],
        ),
        variables = Wflow.AllocationLandVariables(;
            n,
            available_groundwater = [303505.6, 308331.8, 306516.2],
        ),
    )

    demand_variables = Wflow.DemandVariables(;
        n,
        groundwater_demand = [
            2.689178240740741e-7,
            1.4190972222222222e-7,
            7.80462962962963e-6,
        ],
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

    Wflow.groundwater_allocation_area!(allocation_model, demand_variables, domain, dt)

    @test allocation_model.variables.actual_groundwater_abstraction_volume ≈
          [1.6375439409819674, 1.6635833767220893, 1.653787429697564]
    @test allocation_model.variables.actual_groundwater_abstraction ≈
          [2.715116769241939e-6, 2.7582912495537205e-6, 2.742049157130268e-6]
    @test allocation_model.variables.groundwater_allocation ≈
          [2.689178240740741e-7, 1.4190972222222222e-7, 7.80462962962963e-6]
end
