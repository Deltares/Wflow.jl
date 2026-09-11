@testitem "vertical processes" begin
    using Dates
    @test all(
        isapprox.(
            Wflow.rainfall_interception_gash(3.0, 0.11, 0.24, 18.0, 1.5, 4.0),
            (13.568000000000001, 4.0, 0.432, 1.5),
        ),
    )
    @test all(
        isapprox.(
            Wflow.rainfall_interception_modrut(8.6, 3.8, 1.5, 0.45, 2.8),
            (3.87, 3.8, 0.387, 2.043),
        ),
    )
    @test Wflow.head_brooks_corey(0.25, 0.6, 0.15, 10.5, -10.0) ≈ -90.6299820833844
    @test Wflow.feddes_h3(-300.0, -600.0, 3.5, 86400.0) ≈ -412.5
    @test Wflow.feddes_h3(-300.0, -600.0, 0.5, 86400.0) == -600.0
    @test Wflow.feddes_h3(-300.0, -600.0, 6.0, 86400.0) == -300.0
    @test Wflow.rwu_reduction_feddes(0.0, -10.0, -100.0, -300.0, -15000.0, 0.0) == 0.0
    @test Wflow.rwu_reduction_feddes(0.0, -10.0, -100.0, -300.0, -15000.0, 1.0) == 1.0
    @test Wflow.rwu_reduction_feddes(-90.0, -10.0, -100.0, -412.5, -15000.0, 0.0) ≈
          0.8888888888888888
    @test Wflow.rwu_reduction_feddes(-350.0, -10.0, -100.0, -412.5, -15000.0, 0.0) == 1.0
    @test Wflow.rwu_reduction_feddes(-12000.0, -10.0, -100.0, -412.5, -15000.0, 0.0) ≈
          0.20565552699228792
    @test Wflow.rwu_reduction_feddes(-16000.0, -10.0, -100.0, -412.5, -15000.0, 0.0) == 0.0
    @test all(isapprox.(Wflow.infiltration(27.5, 0.2, 50.0, 5.0, 23.5, 1.0), (23.5, 0.5)))
    @test all(
        isapprox.(
            Wflow.unsatzone_flow_layer(43.5, 256.0, 135.0, 12.6),
            (43.49983744545384, 0.00016255454615829025),
        ),
    )
    @test all(
        isapprox.(
            Wflow.precipitation_hbv(30.1, 0.54, 2.0, 0.0),
            (6.923, 23.177000000000003),
        ),
    )
    @test all(
        isapprox.(
            Wflow.snowpack_hbv(201.5, 15.0, 6.923, 23.177, 0.54, 0.0, 2.5, 0.10),
            (207.073, 20.707300000000004, 227.7803, 1.35, 18.819699999999997),
        ),
    )
    @test Wflow.scurve(2.0, 0.0, 3.0, 2.5) ≈ 0.3325863502664285
    @test all(
        isapprox.(
            Wflow.glacier_hbv(0.35, 500.0, 9.5, 5.0, 0.0, 3.4, 0.2, 8.0),
            (8.835, 1.9, 484.9, 17.0),
        ),
    )
end

@testitem "unit: update_infiltration_fluxes" begin
    # Test with infiltration from surface water and precipitation
    potential_infiltration = 10.0
    potential_infiltration_surfacewater = 2.0
    water_flux_surface_input = 10.0
    actual_infiltration_input = 6.0 # this includes the infiltration from surface water
    infiltexcess_input = 1.0

    infilt_surfacewater,
    actual_infiltration,
    infiltexcess,
    excesswater,
    water_flux_surface = Wflow.update_infiltration_fluxes(
        potential_infiltration,
        potential_infiltration_surfacewater,
        water_flux_surface_input,
        actual_infiltration_input,
        infiltexcess_input,
    )

    @test infilt_surfacewater == 1.2 # this is the infiltration from surface water
    @test actual_infiltration ≈ 4.8 # this excludes the infiltration from surface water
    @test infilt_surfacewater + actual_infiltration ≈ actual_infiltration_input
    @test infiltexcess == 0.8
    @test water_flux_surface == 8.0
    @test excesswater ≈ 2.4

    # Test with infiltration from only precipitation, no infiltration from surface water
    potential_infiltration = 10.0
    potential_infiltration_surfacewater = 0.0
    water_flux_surface_input = 10.0
    actual_infiltration_input = 6.0
    infiltexcess_input = 1.0

    infilt_surfacewater,
    actual_infiltration,
    infiltexcess,
    excesswater,
    water_flux_surface = Wflow.update_infiltration_fluxes(
        potential_infiltration,
        potential_infiltration_surfacewater,
        water_flux_surface_input,
        actual_infiltration_input,
        infiltexcess_input,
    )

    @test infilt_surfacewater == 0.0
    @test actual_infiltration == actual_infiltration_input
    @test infilt_surfacewater + actual_infiltration == actual_infiltration_input
    @test infiltexcess == infiltexcess_input
    @test water_flux_surface == water_flux_surface_input
    @test excesswater == 3.0

    # Test with infiltration from only surface water, no infiltration from precipitation
    potential_infiltration = 10.0
    potential_infiltration_surfacewater = 10.0
    water_flux_surface_input = 10.0
    actual_infiltration_input = 6.0
    infiltexcess_input = 1.0

    infilt_surfacewater,
    actual_infiltration,
    infiltexcess,
    excesswater,
    water_flux_surface = Wflow.update_infiltration_fluxes(
        potential_infiltration,
        potential_infiltration_surfacewater,
        water_flux_surface_input,
        actual_infiltration_input,
        infiltexcess_input,
    )

    @test infilt_surfacewater == 6.0
    @test actual_infiltration == 0.0
    @test infilt_surfacewater + actual_infiltration == 6.0
    @test infiltexcess == 0.0
    @test water_flux_surface == 0.0
    @test excesswater == 0.0
end

@testitem "unit: update_overland_flow_and_depth!" begin
    using Wflow:
        KinWaveOverlandFlow,
        ManningFlowParameters,
        OverLandFlowVariables,
        FlowVariables,
        LandFlowBC,
        TimeStepping

    n = 1
    infiltration_amount = 5.0 # mm
    original_depth = 0.02 # m
    river_fraction = 0.2
    expected_water_depth =
        original_depth - ((infiltration_amount * 1e-3) / (1 - river_fraction))

    flow_vars = FlowVariables(n)
    flow_vars.q[1] = 0.0
    variables = OverLandFlowVariables(; flow = flow_vars, to_river = zeros(Float64, n))
    variables.h[1] = original_depth

    mannings_n = [0.072]
    slope = [0.01]
    parameters = ManningFlowParameters(mannings_n, slope)
    parameters.alpha[1] = 2.0

    boundary_conditions = LandFlowBC(; inwater = zeros(Float64, n))
    timestepping =
        TimeStepping(; adaptive = false, dt_fixed = 900.0, stable_timesteps = zeros(n))

    overland_flow_model =
        KinWaveOverlandFlow(; timestepping, boundary_conditions, parameters, variables)

    land_parameters = (; river_fraction = [river_fraction], surface_flow_width = [10.0])

    # Test with positive infiltration
    Wflow.update_overland_flow_and_depth!(
        overland_flow_model,
        infiltration_amount,
        land_parameters,
        1,
    )
    @test overland_flow_model.variables.h[1] ≈ expected_water_depth
    @test overland_flow_model.variables.h[1] ≈ 0.01375
    @test overland_flow_model.variables.flow.q[1] ≈ 0.011537751232883156

    # Test with zero infiltration (no update should occur)
    overland_flow_model.variables.h[1] = original_depth
    overland_flow_model.variables.flow.q[1] = 0.1
    Wflow.update_overland_flow_and_depth!(overland_flow_model, 0.0, land_parameters, 1)
    @test overland_flow_model.variables.flow.q[1] == 0.1
    @test overland_flow_model.variables.h[1] == original_depth
end