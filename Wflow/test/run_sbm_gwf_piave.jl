@testitem "Piave water demand and allocation (sbm_gwf model)" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)

    (; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
    (; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
        model.land.allocation.variables
    (; soil) = model.land
    (; river_flow) = model.routing
    (; reservoir) = river_flow.boundary_conditions

    @testset "piave water demand and allocation first timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc

        @test sum(surfacewater_alloc) ≈ 2.0713528219684277e-5
        @test sum(act_groundwater_abst) ≈ 4.71169680999897e-6
        @test paddy.variables.h[[25, 42, 45]] ≈
              [0.0429685848787047, 0.0, 0.0332318007065323]

        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 2.8935185185185185e-7, 0.0]
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              [0.0, 0.0, 5.424206744631063e-8]
        @test industry.demand.demand_gross[[1, end]] ≈
              [2.4369881591863102e-9, 5.615634533266227e-10]
        @test industry.demand.demand_net[[1, end]] ≈
              [6.09386393248483e-10, 1.4042299611631919e-10]
        @test industry.variables.returnflow[[1, end]] ≈
              [1.8276017659378273e-9, 4.2114045721030344e-10]
        @test livestock.demand.demand_gross[[1, end]] ≈
              [1.1454581336120437e-12, 7.352427591858604e-13]
        @test livestock.demand.demand_net[[1, end]] ≈
              [1.1454581336120437e-12, 7.352427591858604e-13]
        @test livestock.variables.returnflow[[1, end]] ≈ [0.0, 0.0]
        @test domestic.demand.demand_gross[[1, end]] ≈ [6.9591127060077805e-9, 0.0]
        @test domestic.demand.demand_net[[1, end]] ≈ [4.401560182924624e-9, 0.0]
        @test domestic.variables.returnflow[[1, end]] ≈ [2.557552523083157e-9, 0.0]
        @test reservoir.variables.waterlevel ≈
              [23.970590790933628, 32.68607771649563, 39.97018425222192]
        @test reservoir.variables.storage ≈ [1.5532942832524955e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_average ≈
              [3.2489121397532985, 8.556416216129914, 28.17032166161074]
        @test soil.variables.exfiltsatwater[27:31] ≈ [
            2.915452802472663e-7,
            5.847113515629848e-9,
            1.1744022745366517e-7,
            8.048163630421586e-8,
            2.249792843617095e-7,
        ]
        @test maximum(soil.variables.exfiltsatwater) ≈ 2.5641140612390893e-6
        @test soil.variables.exfiltsatwater[17] == 0.0
        @test mean(river_flow.variables.q_average) ≈ 30.071991490895094
        @test maximum(river_flow.variables.q_average) ≈ 117.48258852034441
        @test soil.variables.total_storage[7503] ≈ 0.4729217078886107
        @test soil.variables.total_storage[17] ≈ 0.8174107029296707 # river cell
    end

    Wflow.run_timestep!(model)

    @testset "piave water demand and allocation second timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
        @test sum(surfacewater_alloc) ≈ 1.9051757689751745e-5
        @test sum(act_groundwater_abst) ≈ 4.0518489887643695e-6
        @test paddy.variables.h[[25, 42, 45]] ≈
              [0.03899646776591514, 0.0, 0.02760963170377481]
        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 2.8935185185185185e-7, 0.0]
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              [0.0, 4.935586926461459e-8, 5.8133517727371886e-8]
        @test reservoir.variables.waterlevel ≈
              [23.964702651630343, 32.686077716495625, 39.97018425222192]
        @test reservoir.variables.storage ≈ [1.5529127318256468e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_average ≈
              [3.2489840968665207, 9.467087432483991, 38.61905891069846]
        @test soil.variables.exfiltsatwater[27:33] ≈ [
            4.450632579987645e-7,
            2.1330559140556994e-8,
            1.9144110026163764e-7,
            1.3173453594410016e-7,
            3.326227288442264e-7,
            2.107184034329379e-7,
            2.2676935706227792e-7,
        ]
        @test maximum(soil.variables.exfiltsatwater) ≈ 3.877067765823678e-6
        @test soil.variables.exfiltsatwater[17] == 0.0
        @test mean(river_flow.variables.q_average) ≈ 36.80765857144827
        @test maximum(river_flow.variables.q_average) ≈ 141.84332548772514
        @test soil.variables.total_storage[7503] ≈ 0.4634074094243159
        @test soil.variables.total_storage[17] ≈ 0.8395331546043192 # river cell
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
        @test all(e -> abs(e) < 3.5e-9, land_water_balance.error)
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
        @test all(e -> abs(e) < 1.8e-3, subsurface_water_balance.error)
        @test all(re -> abs(re) <= 2.0, subsurface_water_balance.relative_error)
        @test all(e -> abs(e) < 6.5e-13, subsurface_water_balance.error[inds])
        @test all(re -> abs(re) <= 4.3e-5, subsurface_water_balance.relative_error[inds])
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, land_water_balance.error)
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
        @test all(e -> abs(e) < 4.3e-4, subsurface_water_balance.error)
        @test all(re -> abs(re) < 0.6, subsurface_water_balance.relative_error)
        @test all(e -> abs(e) < 6.4e-13, subsurface_water_balance.error[inds])
        @test all(re -> abs(re) <= 0.009, subsurface_water_balance.relative_error[inds])
    end
    Wflow.close_files(model; delete_output = false)
end
