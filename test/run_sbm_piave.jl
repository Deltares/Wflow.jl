
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
model = Wflow.initialize_sbm_model(config)
q_demand, riv_storage_demand, ssf_storage_demand = run_piave(model, 30)
Wflow.close_files(model; delete_output = false)

tomlpath = joinpath(@__DIR__, "sbm_piave_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.initialize_sbm_model(config)
q_, riv_storage, ssf_storage = run_piave(model, 30)
Wflow.close_files(model; delete_output = false)

@testset "piave with and without water demand" begin
    idx = 1:3:28
    @test q_demand[idx] ≈ [
        218.52013827747507,
        193.0181407476844,
        272.4117008552742,
        161.88920870957202,
        164.81617085174602,
        150.26626520831422,
        121.20488584889729,
        108.10330603823991,
        175.13616509450137,
        108.27468535654428,
    ]
    @test q_[idx] ≈ [
        219.87655637644409,
        197.70869218367324,
        277.7116430199065,
        165.80580629513835,
        173.04068321795665,
        153.84003971725414,
        128.7204398570874,
        112.02619559817663,
        178.8188772116224,
        110.06688910621169,
    ]
    @test riv_storage_demand[idx] ≈ [
        60108.57886064378,
        55483.10543528659,
        62098.337266892515,
        48990.24401897208,
        51879.77039055423,
        43785.59183643158,
        41700.32403346594,
        39630.829260251776,
        44667.073368375095,
        39489.922591430215,
    ]
    @test riv_storage[idx] ≈ [
        60612.54149229024,
        55934.72053279234,
        62385.194284078105,
        49306.029081867746,
        52307.74862683506,
        44028.570459115086,
        41998.65360510959,
        39769.223216103164,
        44589.20697670092,
        39295.941105011,
    ]
    @test ssf_storage_demand[idx] ≈ [
        244164.31393797265,
        239239.45648121895,
        237777.24196149886,
        233306.76472066468,
        229870.4168327788,
        225328.61107592838,
        221034.78635553847,
        216907.78094949326,
        216184.82316118773,
        212805.01551129774,
    ]
    @test ssf_storage[idx] ≈ [
        244152.92920401436,
        239193.8663179015,
        237724.70328343246,
        233250.89649031824,
        229830.5621480702,
        225286.71427856764,
        221007.60584048514,
        216885.51653659056,
        216145.4457927561,
        212756.40451810253,
    ]
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.initialize_sbm_model(config)
Wflow.run_timestep!(model)
sbm = model.land
(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    sbm.allocation.variables

@testset "piave water demand and allocation first timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1030.1528204311428
    @test sum(act_groundwater_abst) ≈ 184.47031930837645
    @test paddy.variables.h[[45, 76, 296]] ≈
          [34.759292485507515, 42.517504464353635, 35.83766686539591]
    @test paddy.parameters.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[45, 76, 296]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[10, 33, 1293]] ≈
          [3.031574420740574, 1.8618934872392217, 0.42233065562148375]
    @test industry.demand.demand_gross[[1, end]] ≈ [0.2105557769536972, 0.0485190823674202]
    @test industry.demand.demand_net[[1, end]] ≈ [0.05265098437666893, 0.012132546864449978]
    @test industry.variables.returnflow[[1, end]] ≈
          [0.15790479257702827, 0.03638653550297022]
    @test livestock.demand.demand_gross[[1, end]] ≈
          [9.896758274408057e-5, 6.352497439365834e-5]
    @test livestock.demand.demand_net[[1, end]] ≈
          [9.896758274408057e-5, 6.352497439365834e-5]
    @test livestock.variables.returnflow[[1, end]] ≈ [0.0, 0.0]
    @test domestic.demand.demand_gross[[1, end]] ≈ [0.5389957427978516, 0.0]
    @test domestic.demand.demand_net[[1, end]] ≈ [0.33949509263038635, 0.0]
    @test domestic.variables.returnflow[[1, end]] ≈ [0.1995006501674652, 0.0]
end

Wflow.run_timestep!(model)
sbm = model.land
(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    sbm.allocation.variables
@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 940.1941924010235
    @test sum(act_groundwater_abst) ≈ 163.59123515939052
    @test paddy.variables.h[[45, 76, 296]] ≈
          [29.50674105890283, 38.11966817469463, 30.679920457042897]
    @test paddy.parameters.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[45, 76, 296]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[10, 33, 1293]] ≈
          [3.8551909476218054, 0.0, 1.3531153828385536]
end

Wflow.close_files(model; delete_output = false)
