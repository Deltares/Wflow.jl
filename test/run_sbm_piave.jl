
function run_piave(model, steps)
    q = zeros(steps)
    ssf_vol = zeros(steps)
    riv_vol = zeros(steps)
    for i in 1:steps
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
Wflow.close_files(model; delete_output = false)

tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.initialize_sbm_model(config)
q_, riv_vol, ssf_vol = run_piave(model, 30)
Wflow.close_files(model; delete_output = false)

@testset "piave with and without water demand" begin
    idx = 1:3:28
    @test q_demand[idx] ≈ [
        218.52770747627918f0,
        193.02890386240665f0,
        271.68778616960685f0,
        161.80734173386602f0,
        164.81279958487485f0,
        150.4149716788231f0,
        121.02659706677429f0,
        108.15426854851842f0,
        174.63022487569546f0,
        108.20883498755789f0,
    ]
    @test q_[idx] ≈ [
        219.9042687883228f0,
        197.73238933696658f0,
        276.99163278840734f0,
        165.74244080863346f0,
        172.99856395134805f0,
        153.91963517651158f0,
        128.52653903788593f0,
        112.02923578000491f0,
        178.19599851038708f0,
        109.99414262238396f0,
    ]
    @test riv_vol_demand[idx] ≈ [
        60125.790043761706f0,
        55476.580177102805f0,
        62195.90597660078f0,
        48963.00368970478f0,
        51871.46897847274f0,
        43794.737044826295f0,
        41710.964839545784f0,
        39623.67826731185f0,
        44719.74807290461f0,
        39494.768128540185f0,
    ]
    @test riv_vol[idx] ≈ [
        60620.42943610538f0,
        55927.61042067993f0,
        62448.21099253601f0,
        49283.04780863544f0,
        52293.53301649927f0,
        44029.56815119591f0,
        41996.885591977836f0,
        39751.39389784692f0,
        44619.74874948371f0,
        39285.722489743566f0,
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
(; paddy, nonpaddy, industry, livestock, domestic) = model.vertical.demand

@testset "piave water demand and allocation first timestep" begin
    sum_total_alloc = sum(sbm.allocation.total_alloc)
    @test sum(sbm.allocation.irri_alloc) + sum(sbm.allocation.nonirri_alloc) ≈
          sum_total_alloc
    @test sum(sbm.allocation.surfacewater_alloc) ≈ 1030.1528204311428f0
    @test sum(sbm.allocation.act_groundwater_abst) ≈ 184.47031930837645f0
    @test paddy.h[[45, 76, 296]] ≈
          [34.759292485507515f0, 42.517504464353635f0, 35.83766686539591f0]
    @test paddy.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test paddy.demand_gross[[45, 76, 296]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test nonpaddy.demand_gross[[10, 33, 1293]] ≈
          [3.031574420740574f0, 1.8618934872392217f0, 0.42233065562148375f0]
    @test industry.demand_gross[[1, end]] ≈ [0.2105557769536972f0, 0.0485190823674202f0]
    @test industry.demand_net[[1, end]] ≈ [0.05265098437666893f0, 0.012132546864449978f0]
    @test industry.returnflow[[1, end]] ≈ [0.15790479257702827f0, 0.03638653550297022f0]
    @test livestock.demand_gross[[1, end]] ≈ [9.896758274408057f-5, 6.352497439365834f-5]
    @test livestock.demand_net[[1, end]] ≈ [9.896758274408057f-5, 6.352497439365834f-5]
    @test livestock.returnflow[[1, end]] ≈ [0.0f0, 0.0f0]
    @test domestic.demand_gross[[1, end]] ≈ [0.5389957427978516f0, 0.0f0]
    @test domestic.demand_net[[1, end]] ≈ [0.33949509263038635f0, 0.0f0]
    @test domestic.returnflow[[1, end]] ≈ [0.1995004952035704f0, 0.0f0]
end

model = Wflow.run_timestep(model)
sbm = model.vertical
(; paddy, nonpaddy, industry, livestock, domestic) = model.vertical.demand
@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(sbm.allocation.total_alloc)
    @test sum(sbm.allocation.irri_alloc) + sum(sbm.allocation.nonirri_alloc) ≈
          sum_total_alloc
    @test sum(sbm.allocation.surfacewater_alloc) ≈ 940.1941924010235f0
    @test sum(sbm.allocation.act_groundwater_abst) ≈ 163.59123515939052f0
    @test paddy.h[[45, 76, 296]] ≈
          [29.50674105890283f0, 38.11966817469463f0, 30.679920457042897f0]
    @test paddy.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test paddy.demand_gross[[45, 76, 296]] ≈ [0.0f0, 0.0f0, 0.0f0]
    @test nonpaddy.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test nonpaddy.demand_gross[[10, 33, 1293]] ≈
          [3.8551909476218054f0, 0.0f0, 1.3531153828385536f0]
end

Wflow.close_files(model; delete_output = false)
