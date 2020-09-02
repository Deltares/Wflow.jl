@testset "reservoir simple" begin
    res = Wflow.SimpleReservoir{Float64}(
        demand = [52.523],
        maxrelease = [420.184],
        maxvolume = [25_000_000.0],
        volume = [1.925e7],
        area = [1885665.353626924],
        targetfullfrac = [0.8],
        targetminfrac = [0.2425554726620697],
        precipitation = [4.2],
        evaporation = [1.5],
    )

    Wflow.update(res, 1, 100.0, 86400.0)
    @test res.outflow[1] ≈ 91.3783714867453
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
        b = [0.22],
        e = [2.0],
        sh = [Wflow.DataFrame()],
        hq = [Wflow.DataFrame()],
        avg_waterlevel = [18.5],
        precipitation = [20.0],
        evaporation = [3.2],
    )

    Wflow.update(lake, 1, 2500.0, 181, 86400.0)
    @test lake.outflow[1] ≈ 85.31903276150577
    @test lake.storage[1] ≈ 3.551103576940606e9
    @test lake.waterlevel[1] ≈ 19.672569557695734
    @test lake.precipitation[1] ≈ 20.0
    @test lake.evaporation[1] ≈ 3.2
end
