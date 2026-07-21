@testitem "Piave water demand and allocation (sbm_gwf model)" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
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
        @test sum(surfacewater_alloc) ≈ 1747.401299613164
        @test sum(act_groundwater_abst) ≈ 392.3593761085565
        @test paddy.variables.h[[25, 42, 45]] ≈ [42.968584878704704, 0.0, 33.2318007065323]
        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 25.0, 0.0]
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈ [0.0, 0.0, 4.487607239115634]
        @test industry.demand.demand_gross[[1, end]] ≈
              [0.2105557769536972, 0.0485190823674202]
        @test industry.demand.demand_net[[1, end]] ≈
              [0.05265098437666893, 0.012132546864449978]
        @test industry.variables.returnflow[[1, end]] ≈
              [0.15790479257702827, 0.03638653550297022]
        @test livestock.demand.demand_gross[[1, end]] ≈
              [9.896758274408057e-5, 6.352497439365834e-5]
        @test livestock.demand.demand_net[[1, end]] ≈
              [9.896758274408057e-5, 6.352497439365834e-5]
        @test livestock.variables.returnflow[[1, end]] ≈ [0.0, 0.0]
        @test domestic.demand.demand_gross[[1, end]] ≈ [0.6012673377990723, 0.0]
        @test domestic.demand.demand_net[[1, end]] ≈ [0.3802947998046875, 0.0]
        @test domestic.variables.returnflow[[1, end]] ≈ [0.2209725379943848, 0.0]
        @test reservoir.variables.waterlevel ≈
              [23.965881805419507, 32.68607771649564, 39.97018425222193]
        @test reservoir.variables.storage ≈ [1.5529891409911808e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_av ≈
              [3.248409441617693, 6.806487195712414, 23.263459419915545]
        @test soil.variables.exfiltsatwater[27:31] ≈ [
            25.189315601597876,
            0.4952061689425148,
            10.146835651996657,
            6.953613376684237,
            19.43821387789424,
        ]
        @test maximum(soil.variables.exfiltsatwater) ≈ 219.59865048404689
        @test soil.variables.exfiltsatwater[17] == 0.0
        @test mean(river_flow.variables.q_av) ≈ 27.05663605784115
        @test maximum(river_flow.variables.q_av) ≈ 108.13408379740821
        @test soil.variables.total_storage[7503] ≈ 472.24913071204054
        @test soil.variables.total_storage[17] ≈ 817.8249788425467 # river cell
    end

    Wflow.run_timestep!(model)

    @testset "piave water demand and allocation second timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
        @test sum(surfacewater_alloc) ≈ 1886.5601313569218
        @test sum(act_groundwater_abst) ≈ 427.16036339859636
        @test paddy.variables.h[[25, 42, 45]] ≈ [38.99642235422809, 0.0, 27.60922173097248]
        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 25.0, 0.0]
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              [0.0, 3.9965040974684207, 5.448194215483002]
        @test reservoir.variables.waterlevel ≈
              [23.95601448606711, 32.68607771649562, 39.970184252221905]
        @test reservoir.variables.storage ≈ [1.552349738697147e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_av ≈
              [3.2471577503452305, 7.762152557160623, 30.551177722294874]
        @test soil.variables.exfiltsatwater[27:33] ≈ [
            38.82639108429825,
            1.8662456728851977,
            16.705432364636337,
            11.497758097330564,
            28.938941979610913,
            18.347740312089048,
            19.767393038510107,
        ]
        @test maximum(soil.variables.exfiltsatwater) ≈ 293.6652322043109
        @test soil.variables.exfiltsatwater[17] == 0.0
        @test mean(river_flow.variables.q_av) ≈ 30.3877491618561
        @test maximum(river_flow.variables.q_av) ≈ 116.84021882103436
        @test soil.variables.total_storage[7503] ≈ 462.4333376594114
        @test soil.variables.total_storage[17] ≈ 823.1079236867769 # river cell
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "water balance piave water demand (sbm_gwf model)" begin
    tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; land_water_balance, routing) = model.mass_balance
    (; overland_water_balance, river_water_balance, subsurface_water_balance) = routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 0.4, land_water_balance.error)
        @test all(re -> abs(re) < 0.0035, land_water_balance.relative_error)
        inds_swd = findall(x -> x >= 0.0, model.land.soil.variables.satwaterdepth)
        @test all(e -> abs(e) < 6.5e-06, land_water_balance.error[inds_swd])
        @test all(e -> abs(e) < 1e-7, land_water_balance.relative_error[inds_swd])
        @test all(e -> abs(e) < 1e-9, overland_water_balance.error)
        @test all(re -> abs(re) < 1e-9, overland_water_balance.relative_error)
        @test all(e -> abs(e) < 0.085, river_water_balance.error)
        @test all(re -> abs(re) < 0.17, river_water_balance.relative_error)
        inds_q = findall(x -> x > 1e-30, model.routing.river_flow.variables.q)
        @test all(e -> abs(e) < 1e-9, river_water_balance.error[inds_q])
        @test all(re -> abs(re) < 1e9, river_water_balance.relative_error[inds_q])
        @test all(e -> abs(e) < 233, subsurface_water_balance.error)
        @test all(re -> abs(re) < 0.0035, subsurface_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-7, subsurface_water_balance.error[inds_swd])
        @test all(re -> abs(re) < 1e-6, subsurface_water_balance.relative_error[inds_swd])
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 6e-6, land_water_balance.error)
        @test all(re -> abs(re) < 1.1e-7, land_water_balance.relative_error)
        @test all(e -> abs(e) < 1.e-9, routing.overland_water_balance.error)
        @test all(re -> abs(re) < 1.e-9, routing.overland_water_balance.relative_error)
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-7, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1.5e-5, subsurface_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end
