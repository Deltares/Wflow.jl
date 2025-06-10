
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
        232.2190840397006,
        202.33701987699484,
        279.149069926829,
        168.30384212407193,
        173.54980693804254,
        158.94278207941792,
        134.37527199930102,
        119.92194171416998,
        190.2080942823956,
        117.21918571350398,
    ]
    @test q_[idx] ≈ [
        233.55717697324866,
        209.84672617198166,
        287.0763139710337,
        176.5245820288135,
        187.03655355782638,
        166.8574912628909,
        145.39994910270275,
        127.58694898697834,
        193.46922956618573,
        123.39989025858517,
    ]
    @test riv_storage_demand[idx] ≈ [
        61889.76178114282,
        57407.04754054791,
        63381.35514200626,
        51249.03165593826,
        54132.20108490535,
        46380.934070560324,
        44749.90001555048,
        42776.01142061794,
        46740.83476372728,
        42129.08350382318,
    ]
    @test riv_storage[idx] ≈ [
        62246.98305925363,
        58008.199331478856,
        63635.712621391496,
        51846.0263625739,
        54857.21536572722,
        46831.03806407208,
        45209.83074250328,
        43188.938080798514,
        46682.53538313581,
        42023.94848706722,
    ]
    @test ssf_storage_demand[idx] ≈ [
        293808.760133334,
        289295.84389213583,
        287742.7480485652,
        283334.41860999924,
        279698.4490921784,
        275163.37629329006,
        270653.8736164591,
        266144.71073891403,
        264639.2237807471,
        261036.52700983256,
    ]
    @test ssf_storage[idx] ≈ [
        293790.6559832811,
        289214.03703327687,
        287610.9738061458,
        283168.07534106117,
        279529.783919687,
        274974.6702018979,
        270494.2769330302,
        266002.2359187149,
        264413.2218890027,
        260768.50688967117,
    ]
end

tomlpath = joinpath(@__DIR__, "sbm_piave_demand_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
Wflow.run_timestep!(model)

(; paddy, nonpaddy, industry, livestock, domestic) = model.land.demand
(; total_alloc, irri_alloc, nonirri_alloc, surfacewater_alloc, act_groundwater_abst) =
    model.land.allocation.variables
(; lake, reservoir) = model.routing.river_flow.boundary_conditions

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
    @test lake.variables.waterlevel_av ≈ [29.317055978102708]
    @test lake.variables.storage_av ≈ [1.899745227381052e8]
    @test lake.variables.outflow_av ≈ [4.864118606306291]
    @test reservoir.variables.storage_av ≈ [4.279958230848676e7, 7.159979217490469e7]
    @test reservoir.variables.outflow_av ≈ [9.756582435870733, 58.08691586159454]
    @test reservoir.variables.demandrelease ≈ [1.182999968528626, 7.902500152587074]
end

Wflow.run_timestep!(model)

@testset "piave water demand and allocation second timestep" begin
    sum_total_alloc = sum(total_alloc)
    @test sum(irri_alloc) + sum(nonirri_alloc) ≈ sum_total_alloc
    @test sum(surfacewater_alloc) ≈ 1831.7297114844905
    @test sum(act_groundwater_abst) ≈ 421.8983408546614
    @test paddy.variables.h[[25, 42, 45]] ≈
          [39.473869456052256, 48.291983056584876, 30.019478609449436]
    @test paddy.parameters.irrigation_trigger[[25, 42, 45]] == [1, 1, 1]
    @test paddy.variables.demand_gross[[25, 42, 45]] ≈ [0.0, 0.0, 0.0]
    @test nonpaddy.parameters.irrigation_trigger[[32, 38, 41]] == [1, 1, 1]
    @test nonpaddy.variables.demand_gross[[32, 38, 41]] ≈
          [3.6972708565911336, 0.600567004968002, 7.152276973950631]
    @test lake.variables.waterlevel_av ≈ [29.324133605862848]
    @test lake.variables.storage_av ≈ [1.9002038576599127e8]
    @test lake.variables.outflow_av ≈ [4.866466879468348]
    @test reservoir.variables.storage_av ≈ [4.2799519150418326e7, 7.159979154258619e7]
    @test reservoir.variables.outflow_av ≈ [9.35321545413373, 54.86473366042814]
    @test reservoir.variables.demandrelease ≈ [1.182999968528626, 7.9025001525870735]
end

Wflow.close_files(model; delete_output = false)
