
function run_piave(model, steps)
    q = zeros(steps)
    ssf_vol = zeros(steps)
    riv_vol = zeros(steps)
    for i in 1:steps
        Wflow.run_timestep!(model)
        ssf_vol[i] = mean(model.lateral.subsurface.variables.volume)
        riv_vol[i] = mean(model.lateral.river.variables.volume)
        q[i] = model.lateral.river.variables.q_av[1]
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
        218.52013823809472f0,
        193.0134951603773f0,
        272.4111837647947f0,
        161.88264628787172f0,
        164.8199089671644f0,
        150.2681168314876f0,
        121.20070337007452f0,
        108.10106381132412f0,
        175.13799714754256f0,
        108.26190463186364f0,
    ]
    @test q_[idx] ≈ [
        219.87655632704903f0,
        197.7038754807009f0,
        277.7110869134211f0,
        165.79913971520423f0,
        173.04466296857905f0,
        153.84187794694486f0,
        128.71609293239374f0,
        112.02394903669563f0,
        178.8207608992179f0,
        110.0540286256144f0,
    ]
    @test riv_vol_demand[idx] ≈ [
        60108.85346812383f0,
        55483.137044195726f0,
        62098.68355844664f0,
        48989.58100040463f0,
        51879.53087453071f0,
        43784.59321641879f0,
        41700.35166194379f0,
        39630.55231098758f0,
        44664.51932077705f0,
        39490.08957716613f0,
    ]
    @test riv_vol[idx] ≈ [
        60612.8152223446f0,
        55934.74021061718f0,
        62385.538865370385f0,
        49305.35863834837f0,
        52307.51338053202f0,
        44027.55259408914f0,
        41998.69770760942f0,
        39768.94289967264f0,
        44586.63915460806f0,
        39296.07358095679f0,
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
Wflow.run_timestep!(model)
sbm = model.vertical
(; paddy, nonpaddy, industry, livestock, domestic) = model.vertical.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    sbm.allocation.variables

@testset "piave water demand and allocation first timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1030.1528204311428f0
    @test sum(act_groundwater_abst) ≈ 184.47031930837645f0
    @test paddy.variables.h[[45, 76, 296]] ≈
          [34.759292485507515f0, 42.517504464353635f0, 35.83766686539591f0]
    @test paddy.parameters.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[45, 76, 296]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[10, 33, 1293]] ≈
          [3.031574420740574f0, 1.8618934872392217f0, 0.42233065562148375f0]
    @test industry.demand.demand_gross[[1, end]] ≈
          [0.2105557769536972f0, 0.0485190823674202f0]
    @test industry.demand.demand_net[[1, end]] ≈
          [0.05265098437666893f0, 0.012132546864449978f0]
    @test industry.variables.returnflow[[1, end]] ≈
          [0.15790479257702827f0, 0.03638653550297022f0]
    @test livestock.demand.demand_gross[[1, end]] ≈
          [9.896758274408057f-5, 6.352497439365834f-5]
    @test livestock.demand.demand_net[[1, end]] ≈
          [9.896758274408057f-5, 6.352497439365834f-5]
    @test livestock.variables.returnflow[[1, end]] ≈ [0.0f0, 0.0f0]
    @test domestic.demand.demand_gross[[1, end]] ≈ [0.5389957427978516f0, 0.0f0]
    @test domestic.demand.demand_net[[1, end]] ≈ [0.33949509263038635f0, 0.0f0]
    @test domestic.variables.returnflow[[1, end]] ≈ [0.1995004952035704f0, 0.0f0]
end

Wflow.run_timestep!(model)
sbm = model.vertical
(; paddy, nonpaddy, industry, livestock, domestic) = model.vertical.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    sbm.allocation.variables
@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 940.1941924010235f0
    @test sum(act_groundwater_abst) ≈ 163.59123515939052f0
    @test paddy.variables.h[[45, 76, 296]] ≈
          [29.50674105890283f0, 38.11966817469463f0, 30.679920457042897f0]
    @test paddy.parameters.irrigation_trigger[[45, 76, 296]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[45, 76, 296]] ≈ [0.0f0, 0.0f0, 0.0f0]
    @test nonpaddy.parameters.irrigation_trigger[[10, 33, 1293]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[10, 33, 1293]] ≈
          [3.8551909476218054f0, 0.0f0, 1.3531153828385536f0]
end

Wflow.close_files(model; delete_output = false)
