
res_bc =
    Wflow.ReservoirBC{Float64}(; inflow = [0.0], precipitation = [4.2], evaporation = [1.5])
res_params = Wflow.ReservoirParameters{Float64}(;
    demand = [52.523],
    maxrelease = [420.184],
    maxvolume = [25_000_000.0],
    area = [1885665.353626924],
    targetfullfrac = [0.8],
    targetminfrac = [0.2425554726620697],
)
res_vars = Wflow.ReservoirVariables{Float64}(;
    volume = [1.925e7],
    outflow_av = [0.0],
    actevap = [0.0],
    outflow = [NaN],
    percfull = [NaN],
    demandrelease = [NaN],
)

res = Wflow.SimpleReservoir{Float64}(;
    boundary_conditions = res_bc,
    parameters = res_params,
    variables = res_vars,
)
@testset "Update reservoir simple" begin
    Wflow.update!(res, 1, 100.0, 86400.0, 86400.0)
    res.variables.outflow_av ./= 86400.0
    @test res.variables.outflow[1] ≈ 91.3783714867453
    @test res.variables.outflow_av[1] == res.variables.outflow[1]
    @test res.variables.volume[1] ≈ 2.0e7
    @test res.variables.percfull[1] ≈ 0.80
    @test res.variables.demandrelease[1] ≈ 52.5229994727611
    @test res.boundary_conditions.precipitation[1] ≈ 4.2
    @test res.boundary_conditions.evaporation[1] ≈ 1.5
    @test res.variables.actevap[1] ≈ 1.5
end

lake_bc =
    Wflow.LakeBC{Float64}(; inflow = [0.0], precipitation = [20.0], evaporation = [3.2])
lake_params = Wflow.LakeParameters{Float64}(;
    lowerlake_ind = [0],
    area = [180510409.0],
    maxstorage = Wflow.maximum_storage([1], [3], [180510409.0], [missing], [missing]),
    threshold = [0.0],
    storfunc = [1],
    outflowfunc = [3],
    b = [0.22],
    e = [2.0],
    sh = [missing],
    hq = [missing],
)
lake_vars = Wflow.LakeVariables{Float64}(;
    outflow_av = [0.0],
    storage = Wflow.initialize_storage([1], [180510409.0], [18.5], [missing]),
    waterlevel = [18.5],
    actevap = [0.0],
    outflow = [NaN],
)

lake = Wflow.Lake{Float64}(;
    boundary_conditions = lake_bc,
    parameters = lake_params,
    variables = lake_vars,
)
@testset "Update lake" begin
    lake_p = lake.parameters
    lake_v = lake.variables
    lake_bc = lake.boundary_conditions
    Wflow.update!(lake, 1, 2500.0, 181, 86400.0, 86400.0)
    lake_v.outflow_av ./= 86400.0
    @test Wflow.waterlevel(lake_p.storfunc, lake_p.area, lake_v.storage, lake_p.sh)[1] ≈
          19.672653848925634
    @test lake_v.outflow[1] ≈ 85.14292808113598
    @test lake_v.outflow_av[1] ≈ lake_v.outflow[1]
    @test lake_v.storage[1] ≈ 3.55111879238499e9
    @test lake_v.waterlevel[1] ≈ 19.672653848925634
    @test lake_bc.precipitation[1] ≈ 20.0
    @test lake_bc.evaporation[1] ≈ 3.2
    @test lake_v.actevap[1] ≈ 3.2
end

datadir = joinpath(@__DIR__, "data")
sh = [
    Wflow.read_sh_csv(joinpath(datadir, "input", "lake_sh_1.csv")),
    Wflow.read_sh_csv(joinpath(datadir, "input", "lake_sh_2.csv")),
]
@testset "linked lakes (HBV)" begin
    @test keys(sh[1]) == (:H, :S)
    @test typeof(values(sh[1])) == Tuple{Vector{Float64}, Vector{Float64}}

    lake_params = Wflow.LakeParameters{Float64}(;
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
        outflowfunc = [2, 1],
        b = [140.0, 0.0],
        e = [1.5, 1.5],
        sh = sh,
        hq = [missing, Wflow.read_hq_csv(joinpath(datadir, "input", "lake_hq_2.csv"))],
    )
    lake_vars = Wflow.LakeVariables{Float64}(;
        outflow_av = [0.0, 0.0],
        waterlevel = [395.03027, 394.87833],
        actevap = [0.0, 0.0],
        outflow = [NaN, NaN],
        storage = Wflow.initialize_storage(
            [2, 2],
            [472461536.0, 60851088.0],
            [395.03027, 394.87833],
            sh,
        ),
    )
    lake_bc = Wflow.LakeBC{Float64}(;
        inflow = [0.0, 0.0],
        precipitation = [10.0, 10.0],
        evaporation = [2.0, 2.0],
    )

    lake = Wflow.Lake{Float64}(;
        boundary_conditions = lake_bc,
        parameters = lake_params,
        variables = lake_vars,
    )

    Wflow.update!(lake, 1, 500.0, 15, 86400.0, 86400.0)
    Wflow.update!(lake, 2, 500.0, 15, 86400.0, 86400.0)
    lake.variables.outflow_av ./= 86400.0
    lake_v = lake.variables
    lake_bc = lake.boundary_conditions
    @test lake_v.outflow ≈ [214.80170846121263, 236.83281600000214]
    @test lake_v.outflow_av ≈ lake_v.outflow
    @test lake_v.storage ≈ [1.2737435094769483e9, 2.6019755340159863e8]
    lake_v.actevap .= 0.0
    lake_v.outflow_av .= 0.0
    lake_bc.inflow .= 0.0
    Wflow.update!(lake, 1, 500.0, 15, 86400.0, 86400.0)
    Wflow.update!(lake, 2, 500.0, 15, 86400.0, 86400.0)
    lake.variables.outflow_av ./= 86400.0
    @test lake_v.outflow ≈ [-259.8005149014703, 239.66710359986183]
    @test lake_v.outflow_av ≈ [-259.8005149014703, 499.4676185013321]
    @test lake_v.storage ≈ [1.3431699662524352e9, 2.6073035986708355e8]
    @test lake_v.waterlevel ≈ [395.239782021054, 395.21771942667266]
    @test lake_v.actevap ≈ [2.0, 2.0]
end

@testset "overflowing lake with sh and hq" begin
    lake_bc = Wflow.LakeBC{Float64}(;
        inflow = [0.00],
        precipitation = [10.0],
        evaporation = [2.0],
    )
    lake_params = Wflow.LakeParameters{Float64}(;
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
        outflowfunc = [1],
        b = [0.0],
        e = [0.0],
        sh = [Wflow.read_sh_csv(joinpath(datadir, "input", "lake_sh_2.csv"))],
        hq = [Wflow.read_hq_csv(joinpath(datadir, "input", "lake_hq_2.csv"))],
    )
    lake_vars = Wflow.LakeVariables{Float64}(;
        outflow_av = [0.0],
        waterlevel = [397.75],
        actevap = [0.0],
        outflow = [NaN],
        storage = [410_760_000],
    )
    lake = Wflow.Lake{Float64}(;
        boundary_conditions = lake_bc,
        parameters = lake_params,
        variables = lake_vars,
    )

    Wflow.update!(lake, 1, 1500.0, 15, 86400.0, 86400.0)
    lake.variables.outflow_av ./= 86400.0
    lake_p = lake.parameters
    lake_v = lake.variables
    @test Wflow.waterlevel(lake_p.storfunc, lake_p.area, lake_v.storage, lake_p.sh) ≈
          [398.0] atol = 1e-2
    @test lake_v.outflow ≈ [1303.67476852]
    @test lake_v.outflow_av ≈ lake_v.outflow
    @test lake_v.storage ≈ [4.293225e8]
    @test lake_v.waterlevel ≈ [398.000000]
end
