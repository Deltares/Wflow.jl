
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
        211.99857265562886f0,
        187.36611456603907f0,
        254.44400351847014f0,
        157.84674222435459f0,
        155.96495231968538f0,
        140.0909893567672f0,
        114.39284413616171f0,
        97.59179667674222f0,
        162.73486135299845f0,
        98.02322384877952f0,
    ]
    @test q_[idx] ≈ [
        214.11044869595924f0,
        194.62341296307127f0,
        262.8784653320827f0,
        165.0866387234042f0,
        168.35443098223672f0,
        149.6184845898368f0,
        125.9357923707458f0,
        108.93173996691478f0,
        168.72826201444647f0,
        104.30259811567151f0,
    ]
    @test riv_vol_demand[idx] ≈ [
        59533.69655965322f0,
        54620.40891412857f0,
        61245.609172260665f0,
        47971.64604359289f0,
        50443.46458977556f0,
        42334.82163836999f0,
        40319.38974793671f0,
        38049.68626252409f0,
        42254.82673655943f0,
        37021.20109521477f0,
    ]
    @test riv_vol[idx] ≈ [
        60170.899828551424f0,
        55441.96744406564f0,
        61876.88671129906f0,
        48934.40419781716f0,
        51529.963663951305f0,
        43536.64514794115f0,
        41396.32544232497f0,
        39318.78688860317f0,
        43039.982293072826f0,
        37921.92492566406f0,
    ]
    @test ssf_vol_demand[idx] ≈ [
        242620.5973130343f0,
        237217.57048043972f0,
        235555.69645613394f0,
        230595.59391329496f0,
        226919.45706105232f0,
        221847.62544826753f0,
        217163.81153879076f0,
        212685.36199288018f0,
        211539.88950085186f0,
        208016.27823978703f0,
    ]
    @test ssf_vol[idx] ≈ [
        242540.8197370262f0,
        237053.14368167287f0,
        235345.29805380543f0,
        230352.50258261946f0,
        226654.8064348857f0,
        221565.5289960294f0,
        216859.27766727592f0,
        212365.27882439757f0,
        211055.71212857193f0,
        207513.09867755277f0,
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
          [35.69494952679905f0, 43.37769296614162f0, 36.03789232562498f0]
    @test sbm.paddy.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test sbm.paddy.demand_gross[[45, 76, 296]] ≈ [0.0, 0.0, 0.0]
    @test sbm.nonpaddy.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test sbm.nonpaddy.demand_gross[[10, 33, 1293]] ≈
          [2.713358341229611f0, 3.1355235274375515f0, 0.5028707708193706f0]
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
          [30.145759167389734f0, 38.81110855637239f0, 30.60798017152777f0]
    @test sbm.paddy.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test sbm.paddy.demand_gross[[45, 76, 296]] ≈ [0.0f0, 0.0f0, 0.0f0]
    @test sbm.nonpaddy.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test sbm.nonpaddy.demand_gross[[10, 33, 1293]] ≈
          [3.1395876909926344f0, 3.2744316785721637f0, 1.464243282784869f0]
end

Wflow.close_files(model, delete_output = false)
