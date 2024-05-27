
function run_piave(model, steps)
    q = zeros(steps)
    ssf_vol = zeros(steps)
    riv_vol = zeros(steps)
    for i = 1:steps
        model = Wflow.run_timestep(model)
        ssf_vol[i] = mean(model.lateral.subsurface.volume)
        riv_vol[i] = mean(model.lateral.river.volume)
        q[i] = model.lateral.river.q_av[1]
    end
    return q, riv_vol, ssf_vol
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.initialize_sbm_model(config)
q_demand, riv_vol_demand, ssf_vol_demand = run_piave(model, 30)
Wflow.close_files(model, delete_output = false)

tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.initialize_sbm_model(config)
q_, riv_vol, ssf_vol = run_piave(model, 30)
Wflow.close_files(model, delete_output = false)

@testset "piave with and without water demand" begin
    idx = 1:3:28
    @test q_demand[idx] ≈ [
        218.12841486736417f0,
        191.4648019414086f0,
        268.7523049972588f0,
        159.1985853669781f0,
        163.52384700836126f0,
        147.18362470295426f0,
        119.50883477637801f0,
        102.32314623731465f0,
        173.3385130195958f0,
        102.78920312709577f0,
    ]
    @test q_[idx] ≈ [
        219.9040394356358f0,
        197.6443639739224f0,
        277.01619562761306f0,
        165.63176310805122f0,
        172.79660240551553f0,
        153.7615009625052f0,
        128.3181374980589f0,
        111.81874211064547f0,
        178.4115119906761f0,
        110.00512168415872f0,
    ]
    @test riv_vol_demand[idx] ≈ [
        60025.17307286695f0,
        55185.129338795574f0,
        61864.38747070659f0,
        48409.79063343459f0,
        51438.36870005233f0,
        43102.31357010651f0,
        41077.4912655103f0,
        38668.25040366511f0,
        44049.26734471644f0,
        38669.64811667626f0,
    ]
    @test riv_vol[idx] ≈ [
        60616.76522336617f0,
        55907.563701271094f0,
        62445.64244715067f0,
        49250.63648358728f0,
        52259.461501724596f0,
        43996.06924202663f0,
        41955.5689678035f0,
        39709.505122092494f0,
        44686.720481505414f0,
        39393.633346173156f0,
    ]
    @test ssf_vol_demand[idx] ≈ [
        244164.36452182673f0,
        239251.33252052963f0,
        237840.70731844133f0,
        233367.45954980003f0,
        229975.36237559692f0,
        225408.34640607517f0,
        221132.71178529615f0,
        217000.35371304862f0,
        216282.7491124505f0,
        212897.1776435227f0,
    ]
    @test ssf_vol[idx] ≈ [
        244138.43308252218f0,
        239156.33098526628f0,
        237706.4375652105f0,
        233210.95033208298f0,
        229790.6284027817f0,
        225218.1340839258f0,
        220921.3928333043f0,
        216784.53043204424f0,
        216074.33950813126f0,
        212679.70047345175f0,
    ]
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
sbm = model.vertical

@testset "piave water demand and allocation first timestep" begin
    sum_total_alloc = sum(sbm.waterallocation.total_alloc)
    @test sum(sbm.waterallocation.irri_alloc) + sum(sbm.waterallocation.nonirri_alloc) ≈
          sum_total_alloc
    @test sum(sbm.waterallocation.surfacewater_alloc) ≈ sum_total_alloc
    @test sum(sbm.waterallocation.act_groundwater_abst) ≈ 0.0
    @test sbm.paddy.h[[45, 76, 296]] ≈
          [34.759292485507515f0, 42.517504464353635f0, 35.83766686539591f0]
    @test sbm.paddy.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test sbm.paddy.demand_gross[[45, 76, 296]] ≈ [0.0, 0.0, 0.0]
    @test sbm.nonpaddy.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test sbm.nonpaddy.demand_gross[[10, 33, 1293]] ≈
          [3.031574420740574f0, 1.8618934872392217f0, 0.42233065562148375f0]
    @test sbm.industry.demand_gross[[1, end]] ≈ [0.2105557769536972f0, 0.0485190823674202f0]
    @test sbm.industry.demand_net[[1, end]] ≈
          [0.05265098437666893f0, 0.012132546864449978f0]
    @test sbm.industry.returnflow[[1, end]] ≈ [0.15790479257702827f0, 0.03638653550297022f0]
    @test sbm.livestock.demand_gross[[1, end]] ≈
          [9.896758274408057f-5, 6.352497439365834f-5]
    @test sbm.livestock.demand_net[[1, end]] ≈ [9.896758274408057f-5, 6.352497439365834f-5]
    @test sbm.livestock.returnflow[[1, end]] ≈ [0.0f0, 0.0f0]
    @test sbm.domestic.demand_gross[[1, end]] ≈ [0.5389957427978516f0, 0.0f0]
    @test sbm.domestic.demand_net[[1, end]] ≈ [0.33949509263038635f0, 0.0f0]
    @test sbm.domestic.returnflow[[1, end]] ≈ [0.1995004952035704f0, 0.0f0]
end

model = Wflow.run_timestep(model)
sbm = model.vertical

@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(sbm.waterallocation.total_alloc)
    @test sum(sbm.waterallocation.irri_alloc) + sum(sbm.waterallocation.nonirri_alloc) ≈
          sum_total_alloc
    @test sum(sbm.waterallocation.surfacewater_alloc) ≈ sum_total_alloc
    @test sum(sbm.waterallocation.act_groundwater_abst) ≈ 0.0
    @test sbm.paddy.h[[45, 76, 296]] ≈
          [29.320512733904195f0, 38.0260019155862f0, 30.538289526380538f0]
    @test sbm.paddy.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test sbm.paddy.demand_gross[[45, 76, 296]] ≈ [0.0f0, 0.0f0, 0.0f0]
    @test sbm.nonpaddy.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test sbm.nonpaddy.demand_gross[[10, 33, 1293]] ≈
          [3.3516483822247363f0, 0.0f0, 1.3680586501618137f0]
end

Wflow.close_files(model, delete_output = false)
