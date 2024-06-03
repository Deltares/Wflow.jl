
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
        218.49829767205395f0,
        192.26513603179683f0,
        268.8372738623575f0,
        160.22423821589794f0,
        162.1093423655203f0,
        148.22988516272903f0,
        118.46881724216179f0,
        106.06772222743426f0,
        171.3264853732719f0,
        104.6559884963633f0,
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
        60119.52320796903f0,
        55280.67738410332f0,
        61843.175499167526f0,
        48502.070698987896f0,
        51315.2679065027f0,
        43198.47343529057f0,
        40986.80096794426f0,
        38822.45987754742f0,
        43899.56858086549f0,
        38558.88933342881f0,
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
        244149.81771426898f0,
        239205.5037877743f0,
        237749.2625556856f0,
        233257.1175894683f0,
        229822.63769575363f0,
        225255.17832343685f0,
        220944.8345415463f0,
        216807.0855559181f0,
        216112.96364731618f0,
        212731.17253347262f0,
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
    @test sum(sbm.waterallocation.surfacewater_alloc) ≈ 1030.1528204311428f0
    @test sum(sbm.waterallocation.act_groundwater_abst) ≈ 184.47031930837645f0
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
    @test sum(sbm.waterallocation.surfacewater_alloc) ≈ 920.8714263051696f0
    @test sum(sbm.waterallocation.act_groundwater_abst) ≈ 156.97872686011868f0
    @test sbm.paddy.h[[45, 76, 296]] ≈
          [29.320510520992286f0, 37.98384192241392f0, 30.496788535054463f0]
    @test sbm.paddy.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test sbm.paddy.demand_gross[[45, 76, 296]] ≈ [0.0f0, 0.0f0, 0.0f0]
    @test sbm.nonpaddy.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test sbm.nonpaddy.demand_gross[[10, 33, 1293]] ≈
          [3.7240080443089862f0, 0.0f0, 1.3680586501618137f0]
end

Wflow.close_files(model, delete_output = false)
