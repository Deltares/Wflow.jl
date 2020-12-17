@testset "reservoir simple" begin
    res = Wflow.SimpleReservoir{Float64}(
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
        outflow = [NaN],
        percfull = [NaN],
        demandrelease = [NaN],
    )

    Wflow.update(res, 1, 100.0, 86400.0)
    @test res.outflow[1] ≈ 91.3783714867453
    @test res.totaloutflow[1] ≈ 7.895091296454794e6
    @test res.volume[1] ≈ 2.0e7
    @test res.percfull[1] ≈ 0.80
    @test res.demandrelease[1] ≈ 52.5229994727611
    @test res.precipitation[1] ≈ 4.2
    @test res.evaporation[1] ≈ 1.5
end

@testset "natural lake" begin
    lake = Wflow.NaturalLake{Float64}(
        loc_id = [1],
        lowerlake_ind = [0],
        area = [180510409.0],
        threshold = [0.0],
        storfunc = [1],
        outflowfunc = [3],
        totaloutflow = [0.0],
        inflow = [0.0],
        b = [0.22],
        e = [2.0],
        sh = [DataFrame()],
        hq = [DataFrame()],
        storage = Wflow.initialize_storage([1], [180510409.0], [18.5], [DataFrame()]),
        waterlevel = [18.5],
        precipitation = [20.0],
        evaporation = [3.2],
        outflow = [NaN],
    )

    Wflow.update(lake, 1, 2500.0, 181, 86400.0)
    @test lake.outflow[1] ≈ 85.31903276150577
    @test lake.totaloutflow[1] ≈ 7.371564430594098e6
    @test lake.storage[1] ≈ 3.551103576940606e9
    @test lake.waterlevel[1] ≈ 19.672569557695734
    @test lake.precipitation[1] ≈ 20.0
    @test lake.evaporation[1] ≈ 3.2
end

datadir = joinpath(@__DIR__, "data")
sh = [
    CSV.read(joinpath(datadir, "lake_sh_1.csv"), DataFrame, type = Float64),
    CSV.read(joinpath(datadir, "lake_sh_2.csv"), DataFrame, type = Float64),
]
@testset "linked lakes (HBV)" begin
    lake = Wflow.NaturalLake{Float64}(
        loc_id = [1, 2],
        lowerlake_ind = [2, 0],
        area = [472461536.0, 60851088.0],
        threshold = [393.7, 0.0],
        storfunc = [2, 2],
        inflow = [0.0, 0.0],
        totaloutflow = [0.0, 0.0],
        outflowfunc = [2, 1],
        b = [140.0, 0.0],
        e = [1.5, 1.5],
        sh = sh,
        hq = [
            DataFrame(),
            CSV.read(joinpath(datadir, "lake_hq_2.csv"), DataFrame, type = Float64),
        ],
        waterlevel = [395.03027, 394.87833],
        precipitation = [10.0, 10.0],
        evaporation = [2.0, 2.0],
        outflow = [NaN, NaN],
        storage = Wflow.initialize_storage([2, 2], [472461536.0, 60851088.0], [395.03027, 394.87833], sh),
    )

    Wflow.update(lake, 1, 500.0, 15, 86400.0)
    Wflow.update(lake, 2, 500.0, 15, 86400.0)
    @test lake.outflow ≈ [214.80170846121263, 236.83281600000214]
    @test lake.totaloutflow ≈ [1.855886761104877e7, 2.0462355302400187e7]
    @test lake.storage ≈ [1.2737435094769483e9, 2.6019755340159863e8]
    @test lake.waterlevel ≈ [395.0912274997361, 395.2101079057371]
    lake.totaloutflow .= 0.0
    lake.inflow .= 0.0
    Wflow.update(lake, 1, 500.0, 15, 86400.0)
    Wflow.update(lake, 2, 500.0, 15, 86400.0)
    @test lake.outflow ≈ [0.0, 239.66710359986183]
    @test lake.totaloutflow ≈ [-2.2446764487487033e7, 2.070723775102806e7]
    @test lake.storage ≈ [1.3431699662524352e9, 2.6073035986708355e8]
    @test lake.waterlevel ≈ [395.239782021054, 395.21771942667266]

end
