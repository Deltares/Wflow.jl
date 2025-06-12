# Simple reservoir (outflowfunc = 4)  
res_bc = Wflow.ReservoirBC(; inflow = [0.0], precipitation = [4.2], evaporation = [1.5])
res_params = Wflow.ReservoirParameters(;
    demand = [52.523],
    maxrelease = [420.184],
    maxstorage = [25_000_000.0],
    area = [1885665.353626924],
    targetfullfrac = [0.8],
    targetminfrac = [0.2425554726620697],
    storfunc = [1],
    outflowfunc = [4],
)
res_vars =
    Wflow.ReservoirVariables(; storage = [1.925e7], waterlevel = [10.208598234556407])

res = Wflow.Reservoir(;
    boundary_conditions = res_bc,
    parameters = res_params,
    variables = res_vars,
)
@testset "Update reservoir simple (outflowfunc = 4)" begin
    Wflow.set_reservoir_vars!(res)
    Wflow.update!(res, 1, 100.0, 86400.0, 86400.0)
    Wflow.average_reservoir_vars!(res, 86400.0)
    @test res.variables.outflow[1] ≈ 91.3783714867453
    @test res.variables.outflow_av[1] == res.variables.outflow[1]
    @test res.variables.storage[1] ≈ 2.0e7
    @test res.variables.storage[1] == res.variables.storage_av[1]
    @test res.variables.percfull[1] ≈ 0.80
    @test res.variables.demandrelease[1] ≈ 52.5229994727611
    @test res.boundary_conditions.precipitation[1] ≈ 4.2
    @test res.boundary_conditions.evaporation[1] ≈ 1.5
    @test res.variables.actevap[1] ≈ 1.5
end

# Reservoir Modified Puls approach (outflowfunc = 3)  
res_bc = Wflow.ReservoirBC(; inflow = [0.0], precipitation = [20.0], evaporation = [3.2])
res_params = Wflow.ReservoirParameters(;
    area = [180510409.0],
    threshold = [0.0],
    storfunc = [1],
    outflowfunc = [3],
    b = [0.22],
    e = [2.0],
)
res_vars = Wflow.ReservoirVariables(;
    storage = Wflow.initialize_storage([1], [180510409.0], [18.5], res_params.sh),
    waterlevel = [18.5],
)

res = Wflow.Reservoir(;
    boundary_conditions = res_bc,
    parameters = res_params,
    variables = res_vars,
)
@testset "Update reservoir Modified Puls approach (outflowfunc = 3)" begin
    res_p = res.parameters
    res_v = res.variables
    res_bc = res.boundary_conditions
    Wflow.set_reservoir_vars!(res)
    Wflow.update!(res, 1, 2500.0, 86400.0, 86400.0)
    Wflow.average_reservoir_vars!(res, 86400.0)
    @test Wflow.waterlevel(
        res_p.storfunc[1],
        res_p.area[1],
        res_v.storage[1],
        res_p.sh[1],
    ) ≈ 19.672653848925634
    @test res_v.outflow[1] ≈ 85.14292808113598
    @test res_v.outflow_av ≈ res_v.outflow
    @test res_v.storage[1] ≈ 3.55111879238499e9
    @test res_v.storage ≈ res_v.storage_av
    @test res_v.waterlevel[1] ≈ 19.672653848925634
    @test res_v.waterlevel ≈ res_v.waterlevel_av
    @test res_bc.precipitation[1] ≈ 20.0
    @test res_bc.evaporation[1] ≈ 3.2
    @test res_v.actevap[1] ≈ 3.2
end

# Linked reservoirs with free weir (outflowfunc = 1)
datadir = joinpath(@__DIR__, "data")
sh = Vector{Union{Wflow.SH, Missing}}([
    Wflow.read_sh_csv(joinpath(datadir, "input", "reservoir_sh_1.csv")),
    Wflow.read_sh_csv(joinpath(datadir, "input", "reservoir_sh_2.csv")),
])
hq = Vector{Union{Wflow.HQ, Missing}}([
    missing,
    Wflow.read_hq_csv(joinpath(datadir, "input", "reservoir_hq_2.csv")),
])
@testset "Linked reservoirs with free weir (outflowfunc = 2)" begin
    @test keys(sh[1]) == (:H, :S)
    @test typeof(values(sh[1])) == Tuple{Vector{Float64}, Vector{Float64}}

    res_params = Wflow.ReservoirParameters(;
        lower_reservoir_ind = [2, 0],
        area = [472461536.0, 60851088.0],
        threshold = [393.7, NaN],
        storfunc = [2, 2],
        outflowfunc = [2, 1],
        b = [140.0, NaN],
        e = [1.5, NaN],
        sh,
        hq,
        col_index_hq = [15],
    )
    @reset res_params.maxstorage[2] = Wflow.maximum_storage(res_params, 2)

    res_vars = Wflow.ReservoirVariables(;
        waterlevel = [395.03027, 394.87833],
        storage = Wflow.initialize_storage(
            [2, 2],
            [472461536.0, 60851088.0],
            [395.03027, 394.87833],
            sh,
        ),
    )
    res_bc = Wflow.ReservoirBC(;
        inflow = [0.0, 0.0],
        precipitation = [10.0, 10.0],
        evaporation = [2.0, 2.0],
    )

    res = Wflow.Reservoir(;
        boundary_conditions = res_bc,
        parameters = res_params,
        variables = res_vars,
    )
    Wflow.set_reservoir_vars!(res)
    Wflow.update!(res, 1, 500.0, 86400.0, 86400.0)
    Wflow.update!(res, 2, 500.0, 86400.0, 86400.0)
    Wflow.average_reservoir_vars!(res, 86400.0)
    res_v = res.variables
    res_bc = res.boundary_conditions
    @test res_v.outflow ≈ [214.80170846121263, 236.83281600000214]
    @test res_v.outflow_av ≈ res_v.outflow
    @test res_v.storage ≈ [1.2737435094769483e9, 2.6019755340159863e8]
    @test res_v.storage ≈ res_v.storage_av
    Wflow.set_reservoir_vars!(res)
    Wflow.update!(res, 1, 500.0, 86400.0, 86400.0)
    Wflow.update!(res, 2, 500.0, 86400.0, 86400.0)
    Wflow.average_reservoir_vars!(res, 86400.0)
    @test res_v.outflow ≈ [-259.8005149014703, 239.66710359986183]
    @test res_v.outflow_av ≈ [-259.8005149014703, 499.4676185013321]
    @test res_v.storage ≈ [1.3431699662524352e9, 2.6073035986708355e8]
    @test res_v.storage ≈ res_v.storage
    @test res_v.waterlevel ≈ [395.239782021054, 395.21771942667266]
    @test res_v.waterlevel ≈ res_v.waterlevel_av
    @test res_v.actevap ≈ [2.0, 2.0]
end

# Overflowing reservoir with SH and HQ (outflowfunc = 1)
@testset "Overflowing reservoir with SH and HQ" begin
    res_bc =
        Wflow.ReservoirBC(; inflow = [0.0], precipitation = [10.0], evaporation = [2.0])
    sh = Vector{Union{Wflow.SH, Missing}}([
        Wflow.read_sh_csv(joinpath(datadir, "input", "reservoir_sh_2.csv")),
    ])
    hq = Vector{Union{Wflow.HQ, Missing}}([
        Wflow.read_hq_csv(joinpath(datadir, "input", "reservoir_hq_2.csv")),
    ])
    res_params = Wflow.ReservoirParameters(;
        area = [200_000_000],
        storfunc = [2],
        outflowfunc = [1],
        sh,
        hq,
        col_index_hq = [15],
    )
    @reset res_params.maxstorage[1] = Wflow.maximum_storage(res_params, 1)
    res_vars = Wflow.ReservoirVariables(; waterlevel = [397.75], storage = [410_760_000])
    res = Wflow.Reservoir(;
        boundary_conditions = res_bc,
        parameters = res_params,
        variables = res_vars,
    )
    Wflow.set_reservoir_vars!(res)
    Wflow.update!(res, 1, 1500.0, 86400.0, 86400.0)
    Wflow.average_reservoir_vars!(res, 86400.0)
    res_p = res.parameters
    res_v = res.variables
    @test Wflow.waterlevel(
        res_p.storfunc[1],
        res_p.area[1],
        res_v.storage[1],
        res_p.sh[1],
    ) ≈ 398.0 atol = 1e-2
    @test res_v.outflow ≈ [1303.67476852]
    @test res_v.outflow_av ≈ res_v.outflow
    @test res_v.storage ≈ [4.293225e8]
    @test res_v.storage ≈ res_v.storage_av
    @test res_v.waterlevel ≈ [398.000000]
    @test res_v.waterlevel ≈ res_v.waterlevel_av
end
