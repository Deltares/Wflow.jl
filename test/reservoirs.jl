@testset "reservoir simple" begin
    res = Wflow.SimpleReservoir{Float64}(
        demand = 52.523,
        maxrelease = 420.184,
        maxvolume = 25000000.0,
        volume = 1.925e7,
        area = 1885665.353626924,
        targetfullfrac = 0.8,
        targetminfrac = 0.2425554726620697,
        Δt = 86400.0,
    )

    res_update = Wflow.update(res, 100.0, 4.2, 1.5)
    @test res_update.outflow ≈ 65.33670482007858
    @test res_update.volume ≈ 2.0e7
    @test res_update.percfull ≈ 0.80
    @test res_update.demandrelease ≈ 52.5229994727611
    @test res_update.precipitation ≈ 4.2
    @test res_update.evaporation ≈ 1.5
end
