@testset "reservoir simple" begin
    res = Wflow.SimpleReservoir{Float64}(
        demand = 52.523,
        maxrelease = 420.184,
        maxvolume = 25000000.0,
        volume = 1.925e7,
        area = 1885665.353626924,
        targetfullfrac = 0.8,
        targetminfrac = 0.2425554726620697,
    )

    res_update = Wflow.update(res, 100.0, 4.2, 1.5, 86400.0)
    @test res_update.outflow ≈ 91.3783714867453
    @test res_update.volume ≈ 2.0e7
    @test res_update.percfull ≈ 0.80
    @test res_update.demandrelease ≈ 52.5229994727611
    @test res_update.precipitation ≈ 4.2
    @test res_update.evaporation ≈ 1.5
end

@testset "natural lake" begin
    lake = Wflow.NaturalLake{Float64}(
        loc_id = 1,
        lowerlake_ind = 0,
        area = 180510409.0,
        threshold = 0.0,
        storfunc = 1,
        outflowfunc = 3,
        b = 0.22,
        e = 2.0,
        sh = Wflow.DataFrame(),
        hq = Wflow.DataFrame(),
        avg_waterlevel = 18.5,
    )

    lake_update = Wflow.update(lake, 2500.0, 20.0, 3.2, 181, 86400.0)
    @test lake_update.outflow ≈ 85.31903276150577
    @test lake_update.storage ≈ 3.551103576940606e9
    @test lake_update.waterlevel ≈ 19.672569557695734
    @test lake_update.precipitation ≈ 20.0
    @test lake_update.evaporation ≈ 3.2
end
