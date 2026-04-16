@testitem "Piave water demand and allocation (sbm_gwf model)" begin
    using Statistics: mean
    using Wflow: to_SI, to_SI!, Unit, MM, MM_PER_DT
    tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    dt = Wflow.tosecond(model.clock.dt)

    (; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
    (; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
        model.land.allocation.variables
    (; soil) = model.land
    (; river_flow) = model.routing
    (; reservoir) = river_flow.boundary_conditions

    @testset "piave water demand and allocation first timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc

        @test sum(surfacewater_alloc) ≈ to_SI(1789.6488381807214, MM_PER_DT; dt_val = dt)
        @test sum(act_groundwater_abst) ≈ to_SI(407.090604383911, MM_PER_DT; dt_val = dt)
        @test paddy.variables.h[[25, 42, 45]] ≈
              to_SI!([42.968584878704704, 0.0, 33.2318007065323], MM)

        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈
              to_SI!([0.0, 25.0, 0.0], MM_PER_DT; dt_val = dt)
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              to_SI!([0.0, 0.0, 4.686514627361238], MM_PER_DT; dt_val = dt)
        @test industry.demand.demand_gross[[1, end]] ≈
              to_SI!([0.2105557769536972, 0.0485190823674202], MM_PER_DT; dt_val = dt)
        @test industry.demand.demand_net[[1, end]] ≈
              to_SI!([0.05265098437666893, 0.012132546864449978], MM_PER_DT; dt_val = dt)
        @test industry.variables.returnflow[[1, end]] ≈
              to_SI!([0.15790479257702827, 0.03638653550297022], MM_PER_DT; dt_val = dt)
        @test livestock.demand.demand_gross[[1, end]] ≈
              to_SI!([9.896758274408057e-5, 6.352497439365834e-5], MM_PER_DT; dt_val = dt)
        @test livestock.demand.demand_net[[1, end]] ≈
              to_SI!([9.896758274408057e-5, 6.352497439365834e-5], MM_PER_DT; dt_val = dt)
        @test livestock.variables.returnflow[[1, end]] ≈ [0.0, 0.0]
        @test domestic.demand.demand_gross[[1, end]] ≈
              to_SI!([0.6012673377990723, 0.0], MM_PER_DT; dt_val = dt)
        @test domestic.demand.demand_net[[1, end]] ≈
              to_SI!([0.3802947998046875, 0.0], MM_PER_DT; dt_val = dt)
        @test domestic.variables.returnflow[[1, end]] ≈
              to_SI!([0.2209725379943848, 0.0], MM_PER_DT; dt_val = dt)
        @test reservoir.variables.waterlevel ≈
              [23.970590790933628, 32.68607771649563, 39.97018425222192]
        @test reservoir.variables.storage ≈ [1.5532942832524955e8, 4.28e7, 7.16e7]
        @test Wflow.get_average(reservoir.variables.outflow_av) ≈
              [3.2489121397532985, 8.556416216129914, 28.17032166161074]
        @test soil.variables.exfiltsatwater[27:31] ≈ to_SI!(
            [
                25.18951221336381,
                0.5051906077504189,
                10.146835651996671,
                6.953613376684251,
                19.4382101688517,
            ],
            MM_PER_DT;
            dt_val = dt,
        )
        @test maximum(soil.variables.exfiltsatwater) ≈
              to_SI(221.53945489105732, MM_PER_DT; dt_val = dt)
        @test soil.variables.exfiltsatwater[17] == 0.0
        @test mean(Wflow.get_average(river_flow.variables.q_av)) ≈ 30.071991490895094
        @test maximum(Wflow.get_average(river_flow.variables.q_av)) ≈ 117.48258852034441
        @test soil.variables.total_storage[7503] ≈ to_SI(472.9217078886107, MM)
        @test soil.variables.total_storage[17] ≈ to_SI(817.4107029296706, MM) # river cell
    end

    Wflow.run_timestep!(model)

    @testset "piave water demand and allocation second timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
        @test sum(surfacewater_alloc) ≈ to_SI(1646.0718643945509, MM_PER_DT; dt_val = dt)
        @test sum(act_groundwater_abst) ≈ to_SI(350.0797526292415, MM_PER_DT; dt_val = dt)
        @test paddy.variables.h[[25, 42, 45]] ≈
              to_SI([38.996467765915135, 0.0, 27.60963170377481], MM)
        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈
              to_SI([0.0, 25.0, 0.0], MM_PER_DT; dt_val = dt)
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              to_SI([0.0, 4.264347104462701, 5.022735931644931], MM_PER_DT; dt_val = dt)
        @test reservoir.variables.waterlevel ≈
              [23.964702651630343, 32.686077716495625, 39.97018425222192]
        @test reservoir.variables.storage ≈ [1.5529127318256468e8, 4.28e7, 7.16e7]
        @test Wflow.get_average(reservoir.variables.outflow_av) ≈
              [3.2489840968665207, 9.467087432483991, 38.61905891069846]
        @test soil.variables.exfiltsatwater[27:33] ≈ to_SI(
            [
                38.453465491093255,
                1.8429603097441243,
                16.540511062605493,
                11.381863905570254,
                28.73860377214116,
                18.206070056605835,
                19.592872450180813,
            ],
            MM_PER_DT;
            dt_val = dt,
        )
        @test maximum(soil.variables.exfiltsatwater) ≈
              to_SI(334.9786549671658, MM_PER_DT; dt_val = dt)
        @test soil.variables.exfiltsatwater[17] == 0.0
        @test mean(Wflow.get_average(river_flow.variables.q_av)) ≈ 36.80765857144827
        @test maximum(Wflow.get_average(river_flow.variables.q_av)) ≈ 141.84332548772514
        @test soil.variables.total_storage[7503] ≈ to_SI(463.4074094243159, MM)
        @test soil.variables.total_storage[17] ≈ to_SI(839.5331546043192, MM) # river cell
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "water balance piave water demand (sbm_gwf model)" begin
    tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; land_water_balance, routing) = model.mass_balance
    (; overland_water_balance, river_water_balance, subsurface_water_balance) = routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 0.30, land_water_balance.error)
        @test all(re -> abs(re) < 0.06, land_water_balance.relative_error)
        inds = findall(
            x -> !iszero(x),
            Wflow.saturated_thickness(model.routing.subsurface_flow),
        )
        @test length(inds) == 7502
        @test all(e -> abs(e) < 1e-9, land_water_balance.error[inds])
        @test all(re -> abs(re) < 1e-9, land_water_balance.relative_error[inds])
        @test all(e -> abs(e) < 1e-9, overland_water_balance.error)
        @test all(re -> abs(re) < 1e-9, overland_water_balance.relative_error)
        @test all(e -> abs(e) < 1.e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(e -> abs(e) < 155.0, subsurface_water_balance.error)
        @test all(re -> abs(re) <= 2.0, subsurface_water_balance.relative_error)
        @test all(e -> abs(e) < 5.6e-8, subsurface_water_balance.error[inds])
        @test all(re -> abs(re) <= 4.3e-5, subsurface_water_balance.relative_error[inds])
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 0.07, land_water_balance.error)
        @test all(re -> abs(re) < 0.04, land_water_balance.relative_error)
        inds = findall(
            x -> !iszero(x),
            Wflow.saturated_thickness(model.routing.subsurface_flow),
        )
        @test length(inds) == 7501
        @test all(e -> abs(e) < 1e-9, land_water_balance.error[inds])
        @test all(re -> abs(re) < 1e-9, land_water_balance.relative_error[inds])
        @test all(e -> abs(e) < 1.e-9, routing.overland_water_balance.error)
        @test all(re -> abs(re) < 1.e-9, routing.overland_water_balance.relative_error)
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(e -> abs(e) < 37.0, subsurface_water_balance.error)
        @test all(re -> abs(re) < 0.6, subsurface_water_balance.relative_error)
        @test all(e -> abs(e) < 5.5e-8, subsurface_water_balance.error[inds])
        @test all(re -> abs(re) <= 0.009, subsurface_water_balance.relative_error[inds])
    end
    Wflow.close_files(model; delete_output = false)
end
