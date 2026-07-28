@testitem "Piave water demand and allocation (sbm_gwf model)" begin
    using Statistics: mean
    tomlpath = joinpath(@__DIR__, "sbm_gwf_piave_demand_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)

    (; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
    (;
        total_alloc,
        irrigation_allocation,
        non_irrigation_allocation,
        surfacewater_allocation,
        actual_groundwater_abstraction,
    ) = model.land.allocation.variables
    (; soil) = model.land
    (; river_flow) = model.routing
    (; reservoir) = river_flow.boundary_conditions

    @testset "piave water demand and allocation first timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irrigation_allocation) + sum(non_irrigation_allocation) ≈ sum_total_alloc

        @test sum(surfacewater_allocation) ≈ 2.0713528219684277e-5
        @test sum(actual_groundwater_abstraction) ≈ 4.71169680999897e-6
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
              [23.969018337050507, 32.68607771649563, 39.97018425222191]
        @test reservoir.variables.storage ≈ [1.5531923882408744e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_average ≈
              [3.248754128864373, 7.716341989778147, 23.646737744565055]
        @test soil.variables.exfiltration_saturated_water[27:31] ≈ [
            2.9154273159320574e-7,
            5.731552881278919e-9,
            1.1744022745366522e-7,
            8.048163630421585e-8,
            2.2497894376592591e-7,
        ]
        @test maximum(soil.variables.exfiltration_saturated_water) ≈ 2.5421719998228685e-6
        @test soil.variables.exfiltration_saturated_water[17] == 0.0
        @test mean(river_flow.variables.q_average) ≈ 28.022713543431458
        @test maximum(river_flow.variables.q_average) ≈ 111.3683220747618
        @test soil.variables.total_storage[7503] ≈ 0.4726724673643876
        @test soil.variables.total_storage[17] ≈ 0.8173997403014627 # river cell
    end

    Wflow.run_timestep!(model)

    @testset "piave water demand and allocation second timestep" begin
        sum_total_alloc = sum(total_alloc)
        @test sum(irrigation_allocation) + sum(non_irrigation_allocation) ≈ sum_total_alloc
        @test sum(surfacewater_allocation) ≈ 1.9052384471605346e-5
        @test sum(actual_groundwater_abstraction) ≈ 4.052118012811413e-6
        @test paddy.variables.h[[25, 42, 45]] ≈
              [0.03899640238002065, 0.0, 0.02760916845586318]
        @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
        @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 2.8935185185185185e-7, 0.0]
        @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
        @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
              [0.0, 4.935586926461459e-8, 5.8133517727371886e-8]
        @test reservoir.variables.waterlevel ≈
              [23.95962920244703, 32.68607771649562, 39.9701842522219]
        @test reservoir.variables.storage ≈ [1.5525839723185676e8, 4.28e7, 7.16e7]
        @test reservoir.variables.outflow_average ≈
              [3.2481416708540873, 7.983571794690894, 30.843454389874523]
        @test soil.variables.exfiltration_saturated_water[27:33] ≈ [
            4.450183798563528e-7,
            2.1210622200963962e-8,
            1.9141828192072716e-7,
            1.317189512235159e-7,
            3.326289992470102e-7,
            2.1072154311388192e-7,
            2.267434539195128e-7,
        ]
        @test maximum(soil.variables.exfiltration_saturated_water) ≈ 3.4165805335899716e-6
        @test soil.variables.exfiltration_saturated_water[17] == 0.0
        @test mean(river_flow.variables.q_average) ≈ 31.665942201394245
        @test maximum(river_flow.variables.q_average) ≈ 122.91999441565204
        @test soil.variables.total_storage[7503] ≈ 0.4637333311448208
        @test soil.variables.total_storage[17] ≈ 0.8290373705017975 # river cell
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
        @test all(e -> abs(e) < 1e-9, land_water_balance.error)
        @test all(re -> abs(re) < 1e-09, land_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, overland_water_balance.error)
        @test all(re -> abs(re) < 1e-9, overland_water_balance.relative_error)
        @test all(e -> abs(e) < 0.065, river_water_balance.error)
        @test all(re -> abs(re) < 0.13, river_water_balance.relative_error)
        inds_q = findall(x -> x > 1e-30, model.routing.river_flow.variables.q)
        @test all(e -> abs(e) < 1e-9, river_water_balance.error[inds_q])
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error[inds_q])
        @test all(e -> abs(e) < 1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, land_water_balance.error)
        @test all(re -> abs(re) < 1e-9, land_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, routing.overland_water_balance.error)
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end
