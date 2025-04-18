
function run_piave(model, steps)
    q = zeros(steps)
    ssf_storage = zeros(steps)
    riv_storage = zeros(steps)
    for i in 1:steps
        Wflow.run_timestep!(model)
        ssf_storage[i] = mean(model.routing.subsurface_flow.variables.storage)
        riv_storage[i] = mean(model.routing.river_flow.variables.storage)
        q[i] = model.routing.river_flow.variables.q_av[1]
    end
    return q, riv_storage, ssf_storage
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
q_demand, riv_storage_demand, ssf_storage_demand = run_piave(model, 30)
Wflow.close_files(model; delete_output = false)

tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
q_, riv_storage, ssf_storage = run_piave(model, 30)
Wflow.close_files(model; delete_output = false)

@testset "piave with and without water demand" begin
    idx = 1:3:28
    @test q_demand[idx] ≈ [
        232.21872226317748,
        206.37375194726908,
        283.4020828361598,
        172.57960383677306,
        178.07343462854928,
        163.21566221793816,
        137.21564464378784,
        122.74421137559834,
        192.62288427631623,
        122.23708908443999,
    ]
    @test q_[idx] ≈ [
        233.55677258454014,
        209.54119025449393,
        288.56407748672325,
        176.49212495508897,
        187.1757131785339,
        167.03798004492944,
        144.88980522987447,
        126.89156532854736,
        195.70222987183843,
        124.11042682514802,
    ]
    @test riv_storage_demand[idx] ≈ [
        61867.695252375735,
        57811.24186409073,
        63597.61445768724,
        51589.564013367,
        54555.583506296425,
        46665.95493213195,
        44936.20032381044,
        43067.63775611285,
        47106.75536353861,
        42497.31760385939,
    ]
    @test riv_storage[idx] ≈ [
        62227.976323430026,
        57973.316933646274,
        63732.81459276585,
        51838.45431961542,
        54936.95712553251,
        46795.33035409312,
        45107.809026071554,
        43087.08618047318,
        46933.017637935765,
        42151.73831729207,
    ]
    @test ssf_storage_demand[idx] ≈ [
        293719.7652093973,
        289154.59810634685,
        287708.2752211687,
        283268.8890246969,
        279727.7213808278,
        275096.4755249175,
        270613.46791007475,
        266244.0655551083,
        264989.9257173528,
        261402.26576472278,
    ]
    @test ssf_storage[idx] ≈ [
        293702.19922404754,
        289084.77973608236,
        287624.6254231974,
        283175.87651864055,
        279643.2523878947,
        274994.91193095397,
        270519.85147105134,
        266139.88635706156,
        264856.5571256621,
        261244.76356086627,
    ]
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
Wflow.run_timestep!(model)
sbm = model.land
(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    sbm.allocation.variables

@testset "piave water demand and allocation first timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1059.4215161438574
    @test sum(act_groundwater_abst) ≈ 187.24142226501345
    @test paddy.variables.h[[25, 42, 45]] ≈
          [43.42843518721493, 51.453583244841425, 35.52247154921834]
    @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
          [1.9688564838354994, 0.3644130691439055, 0.1118568579535997]
    @test industry.demand.demand_gross[[1, end]] ≈ [0.2105557769536972, 0.0485190823674202]
    @test industry.demand.demand_net[[1, end]] ≈ [0.05265098437666893, 0.012132546864449978]
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
end

Wflow.run_timestep!(model)
sbm = model.land
(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    sbm.allocation.variables
@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1117.8312694975748
    @test sum(act_groundwater_abst) ≈ 195.80941021945517
    @test paddy.variables.h[[25, 42, 45]] ≈
          [39.51896384567826, 48.34070462363231, 30.228023083909022]
    @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈ [0.0, 0.0, 4.566942543149944]
end

Wflow.close_files(model; delete_output = false)
