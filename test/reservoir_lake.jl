res = Wflow.SimpleReservoir{Float64}(
    dt = 86400.0,
    demand = [52.523],
    maxrelease = [420.184],
    maxvolume = [25_000_000.0],
    volume = [1.925e7],
    totaloutflow = [0.0],
    inflow = [0.0],
    area = [1885665.353626924],
    targetfullfrac = [0.8],
    targetminfrac = [0.2425554726620697],
    precipitation = [4.2],
    evaporation = [1.5],
    actevap = [0.0],
    outflow = [NaN],
    percfull = [NaN],
    demandrelease = [NaN],
)
@testset "Update reservoir simple" begin
    Wflow.update(res, 1, 100.0, 86400.0)
    @test res.outflow[1] ≈ 91.3783714867453
    @test res.totaloutflow[1] ≈ 7.895091296454794e6
    @test res.volume[1] ≈ 2.0e7
    @test res.percfull[1] ≈ 0.80
    @test res.demandrelease[1] ≈ 52.5229994727611
    @test res.precipitation[1] ≈ 4.2
    @test res.evaporation[1] ≈ 1.5
    @test res.actevap[1] ≈ 1.5
end

@testset "Exchange and grid location reservoir" begin
    @test Wflow.exchange(res, :dt) == 0
    @test Wflow.exchange(res, :volume) == 1
    @test Wflow.exchange(res, :outflow) == 1
    @test Wflow.grid_location(res, :dt) == "none"
    @test Wflow.grid_location(res, :volume) == "node"
    @test Wflow.grid_location(res, :outflow) == "node"
end

lake = Wflow.Lake{Float64}(
    dt = 86400.0,
    lowerlake_ind = [0],
    area = [180510409.0],
    maxstorage = Wflow.maximum_storage([1], [3], [180510409.0], [missing], [missing]),
    threshold = [0.0],
    storfunc = [1],
    outflowfunc = [3],
    totaloutflow = [0.0],
    inflow = [0.0],
    b = [0.22],
    e = [2.0],
    sh = [missing],
    hq = [missing],
    storage = Wflow.initialize_storage([1], [180510409.0], [18.5], [missing]),
    waterlevel = [18.5],
    precipitation = [20.0],
    evaporation = [3.2],
    actevap = [0.0],
    outflow = [NaN],
)
@testset "Update lake" begin
    Wflow.update(lake, 1, 2500.0, 181, 86400.0)
    @test lake.outflow[1] ≈ 85.14292808113598
    @test lake.totaloutflow[1] ≈ 7.356348986210149e6
    @test lake.storage[1] ≈ 3.55111879238499e9
    @test lake.waterlevel[1] ≈ 19.672653848925634
    @test lake.precipitation[1] ≈ 20.0
    @test lake.evaporation[1] ≈ 3.2
    @test lake.actevap[1] ≈ 3.2
end

@testset "Exchange and grid location lake" begin
    @test Wflow.exchange(lake, :dt) == 0
    @test Wflow.exchange(lake, :storage) == 1
    @test Wflow.exchange(lake, :outflow) == 1
    @test Wflow.grid_location(lake, :dt) == "none"
    @test Wflow.grid_location(lake, :storage) == "node"
    @test Wflow.grid_location(lake, :outflow) == "node"
end

datadir = joinpath(@__DIR__, "data")
sh = [
    Wflow.read_sh_csv(joinpath(datadir, "input", "lake_sh_1.csv")),
    Wflow.read_sh_csv(joinpath(datadir, "input", "lake_sh_2.csv")),
]
@testset "linked lakes (HBV)" begin
    @test keys(sh[1]) == (:H, :S)
    @test typeof(values(sh[1])) == Tuple{Vector{Float},Vector{Float}}

    lake = Wflow.Lake{Float}(
        dt = 86400.0,
        lowerlake_ind = [2, 0],
        area = [472461536.0, 60851088.0],
        maxstorage = Wflow.maximum_storage(
            [2, 2],
            [2, 1],
            [472461536.0, 60851088.0],
            sh,
            [missing, Wflow.read_hq_csv(joinpath(datadir, "input", "lake_hq_2.csv"))],
        ),
        threshold = [393.7, 0.0],
        storfunc = [2, 2],
        inflow = [0.0, 0.0],
        totaloutflow = [0.0, 0.0],
        outflowfunc = [2, 1],
        b = [140.0, 0.0],
        e = [1.5, 1.5],
        sh = sh,
        hq = [missing, Wflow.read_hq_csv(joinpath(datadir, "input", "lake_hq_2.csv"))],
        waterlevel = [395.03027, 394.87833],
        precipitation = [10.0, 10.0],
        evaporation = [2.0, 2.0],
        actevap = [0.0, 0.0],
        outflow = [NaN, NaN],
        storage = Wflow.initialize_storage(
            [2, 2],
            [472461536.0, 60851088.0],
            [395.03027, 394.87833],
            sh,
        ),
    )

    Wflow.update(lake, 1, 500.0, 15, 86400.0)
    Wflow.update(lake, 2, 500.0, 15, 86400.0)
    @test lake.outflow ≈ [214.80170846121263, 236.83281600000214] atol = 1e-2
    @test lake.totaloutflow ≈ [1.855886761104877e7, 2.0462355302400187e7] atol = 1e3
    @test lake.storage ≈ [1.2737435094769483e9, 2.6019755340159863e8] atol = 1e4
    @test lake.waterlevel ≈ [395.0912274997361, 395.2101079057371] atol = 1e-2
    lake.actevap .= 0.0
    lake.totaloutflow .= 0.0
    lake.inflow .= 0.0
    Wflow.update(lake, 1, 500.0, 15, 86400.0)
    Wflow.update(lake, 2, 500.0, 15, 86400.0)
    @test lake.outflow ≈ [0.0, 239.66710359986183] atol = 1e-2
    @test lake.totaloutflow ≈ [-2.2446764487487033e7, 4.3154002238515094e7] atol = 1e3
    @test lake.storage ≈ [1.3431699662524352e9, 2.6073035986708355e8] atol = 1e4
    @test lake.waterlevel ≈ [395.239782021054, 395.21771942667266] atol = 1e-2
    @test lake.actevap ≈ [2.0, 2.0]
end

@testset "overflowing lake with sh and hq" begin
    lake = Wflow.Lake{Float}(
        dt = 86400.0,
        lowerlake_ind = [0],
        area = [200_000_000],
        maxstorage = Wflow.maximum_storage(
            [2],
            [1],
            [200_000_000],
            [Wflow.read_sh_csv(joinpath(datadir, "input", "lake_sh_2.csv"))],
            [Wflow.read_hq_csv(joinpath(datadir, "input", "lake_hq_2.csv"))],
        ),
        threshold = [0.0],
        storfunc = [2],
        inflow = [0.00],
        totaloutflow = [0.0],
        outflowfunc = [1],
        b = [0.0],
        e = [0.0],
        sh = [Wflow.read_sh_csv(joinpath(datadir, "input", "lake_sh_2.csv"))],
        hq = [Wflow.read_hq_csv(joinpath(datadir, "input", "lake_hq_2.csv"))],
        waterlevel = [397.75],
        precipitation = [10.0],
        evaporation = [2.0],
        actevap = [0.0],
        outflow = [NaN],
        storage = [410_760_000],
    )

    Wflow.update(lake, 1, 1500.0, 15, 86400.0)
    @test lake.outflow ≈ [1303.67476852] atol = 1e-2
    @test lake.totaloutflow ≈ [11.26375000e7] atol = 1e3
    @test lake.storage ≈ [4.293225e8] atol = 1e4
    @test lake.waterlevel ≈ [398.000000] atol = 1e-2
end
