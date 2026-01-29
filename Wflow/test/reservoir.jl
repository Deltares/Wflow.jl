@testitem "Update reservoir simple" begin
    using Wflow: ReservoirProfileType, ReservoirOutflowType, Unit, to_SI, MM_PER_DT
    dt = 86400.0
    # Simple reservoir (outflowfunc = 4)
    n = 1
    res_bc = Wflow.ReservoirBC(;
        n,
        precipitation = [to_SI(4.2, MM_PER_DT; dt_val = dt)],
        evaporation = [to_SI(1.5, MM_PER_DT; dt_val = dt)],
    )
    res_params = Wflow.ReservoirParameters(;
        id = [1],
        demand = [52.523],
        maxrelease = [420.184],
        maxstorage = [25_000_000.0],
        area = [1885665.353626924],
        targetfullfrac = [0.8],
        targetminfrac = [0.2425554726620697],
        storfunc = [ReservoirProfileType.linear],
        outflowfunc = [ReservoirOutflowType.simple],
    )
    res_vars = Wflow.ReservoirVariables(;
        outflow_obs = [Wflow.MISSING_VALUE],
        storage = [1.925e7],
        waterlevel = [10.208598234556407],
    )

    res = Wflow.Reservoir(;
        boundary_conditions = res_bc,
        parameters = res_params,
        variables = res_vars,
    )
    @testset "Update reservoir simple (outflowfunc = 4)" begin
        Wflow.set_reservoir_vars!(res)
        Wflow.update!(res, 1, 100.0, dt)
        @test res.variables.outflow[1] ≈ 91.3783714867453
        @test res.variables.outflow_av.cumulative[1] == res.variables.outflow[1] * dt
        @test res.variables.storage[1] ≈ 2.0e7
        @test res.boundary_conditions.precipitation[1] ≈ to_SI(4.2, MM_PER_DT; dt_val = dt)
        @test res.boundary_conditions.evaporation[1] ≈ to_SI(1.5, MM_PER_DT; dt_val = dt)
        @test res.variables.actevap.cumulative[1] ≈ to_SI(1.5, MM_PER_DT; dt_val = dt) * dt
    end

    # reset storage and waterlevel and set observed outflow
    res.variables.outflow_obs[1] = 80.0
    res.variables.storage[1] = 1.925e7
    res.variables.waterlevel[1] = 10.208598234556407
    @testset "Update reservoir simple (outflowfunc = 4) with observed outflow" begin
        Wflow.set_reservoir_vars!(res)
        Wflow.update!(res, 1, 100.0, dt)
        Wflow.average_reservoir_vars!(res, 86400.0)
        @test res.variables.outflow[1] ≈ 80.0
        @test res.variables.outflow_av.average[1] == res.variables.outflow[1]
        @test res.variables.storage[1] ≈ 2.0983091296454795e7
    end
end

@testitem "Update reservoir Modified Puls approach (outflowfunc = 3)" begin
    using Wflow: to_SI, Unit, MM_PER_DT, ReservoirProfileType, ReservoirOutflowType
    # Reservoir Modified Puls approach (outflowfunc = 3)
    n = 1
    dt = 86400.0
    res_bc = Wflow.ReservoirBC(;
        n,
        precipitation = [to_SI(20.0, MM_PER_DT; dt_val = dt)],
        evaporation = [to_SI(3.2, MM_PER_DT; dt_val = dt)],
    )
    res_params = Wflow.ReservoirParameters(;
        id = [1],
        area = [180510409.0],
        threshold = [0.0],
        storfunc = [ReservoirProfileType.linear],
        outflowfunc = [ReservoirOutflowType.modified_puls],
        b = [0.22],
        e = [2.0],
    )
    waterlevel = [18.5]
    res_vars = Wflow.ReservoirVariables(;
        storage = Wflow.initialize_storage(
            res_params.storfunc,
            res_params.area,
            waterlevel,
            res_params.sh,
        ),
        waterlevel,
    )

    res = Wflow.Reservoir(;
        boundary_conditions = res_bc,
        parameters = res_params,
        variables = res_vars,
    )
    res_p = res.parameters
    res_v = res.variables
    res_bc = res.boundary_conditions
    Wflow.set_reservoir_vars!(res)
    Wflow.update!(res, 1, 2500.0, dt)
    Wflow.average_reservoir_vars!(res, dt)
    @test Wflow.waterlevel(
        res_p.storfunc[1],
        res_p.area[1],
        res_v.storage[1],
        res_p.sh[1],
    ) ≈ 19.672653848925634
    @test res_v.outflow[1] ≈ 85.14292808113598
    @test res_v.outflow_av.average ≈ res_v.outflow
    @test res_v.storage[1] ≈ 3.55111879238499e9
    @test res_v.waterlevel[1] ≈ 19.672653848925634
    @test res_bc.precipitation[1] ≈ to_SI(20.0, MM_PER_DT; dt_val = dt)
    @test res_bc.evaporation[1] ≈ to_SI(3.2, MM_PER_DT; dt_val = dt)
    @test res_v.actevap.cumulative[1] ≈ to_SI(3.2, MM_PER_DT; dt_val = dt) * dt
end

@testitem "Linked reservoirs with free weir (outflowfunc = 2)" begin
    using Wflow: to_SI, Unit, MM_PER_DT, ReservoirProfileType, ReservoirOutflowType
    dt = 86400.0
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

    @test keys(sh[1]) == (:H, :S)
    @test typeof(values(sh[1])) == Tuple{Vector{Float64}, Vector{Float64}}

    res_params = Wflow.ReservoirParameters(;
        id = [1, 2],
        lower_reservoir_ind = [2, 0],
        area = [472461536.0, 60851088.0],
        threshold = [393.7, NaN],
        storfunc = fill(ReservoirProfileType.interpolation, 2),
        outflowfunc = [ReservoirOutflowType.free_weir, ReservoirOutflowType.rating_curve],
        b = [140.0, NaN],
        e = [1.5, NaN],
        sh,
        hq,
        col_index_hq = [15],
    )
    res_params.maxstorage[2] = Wflow.maximum_storage(res_params, 2)

    waterlevel = [395.03027, 394.87833]
    res_vars = Wflow.ReservoirVariables(;
        waterlevel,
        storage = Wflow.initialize_storage(
            res_params.storfunc,
            [472461536.0, 60851088.0],
            waterlevel,
            sh,
        ),
    )
    n = 2
    res_bc = Wflow.ReservoirBC(;
        n = 2,
        precipitation = [
            to_SI(10.0, MM_PER_DT; dt_val = dt),
            to_SI(10.0, MM_PER_DT; dt_val = dt),
        ],
        evaporation = [
            to_SI(2.0, MM_PER_DT; dt_val = dt),
            to_SI(2.0, MM_PER_DT; dt_val = dt),
        ],
    )

    res = Wflow.Reservoir(;
        boundary_conditions = res_bc,
        parameters = res_params,
        variables = res_vars,
    )
    Wflow.set_reservoir_vars!(res)
    Wflow.update!(res, 1, 500.0, dt)
    Wflow.update!(res, 2, 500.0, dt)
    Wflow.average_reservoir_vars!(res, dt)
    res_v = res.variables
    res_bc = res.boundary_conditions
    @test res_v.outflow ≈ [214.80170846121263, 236.83281600000214]
    @test res_v.outflow_av.average ≈ res_v.outflow
    @test res_v.storage ≈ [1.2737435094769483e9, 2.6019755340159863e8]
    Wflow.set_reservoir_vars!(res)
    Wflow.update!(res, 1, 500.0, dt)
    Wflow.update!(res, 2, 500.0, dt)
    Wflow.average_reservoir_vars!(res, dt)
    @test res_v.outflow ≈ [-259.8005149014703, 239.66710359986183]
    @test res_v.outflow_av.average ≈ [-259.8005149014703, 499.4676185013321]
    @test res_v.storage ≈ [1.3431699662524352e9, 2.6073035986708355e8]
    @test res_v.waterlevel ≈ [395.239782021054, 395.21771942667266]
    @test res_v.actevap.cumulative ≈
          [to_SI(2.0, MM_PER_DT; dt_val = dt) * dt, to_SI(2.0, MM_PER_DT; dt_val = dt) * dt]
end

# Overflowing reservoir with SH and HQ (outflowfunc = 1)
@testitem "Overflowing reservoir with SH and HQ" begin
    using Wflow: to_SI, Unit, MM_PER_DT, ReservoirProfileType, ReservoirOutflowType
    datadir = joinpath(@__DIR__, "data")
    n = 1
    dt = 86400.0
    res_bc = Wflow.ReservoirBC(;
        n,
        precipitation = [to_SI(10.0, MM_PER_DT; dt_val = dt)],
        evaporation = [to_SI(2.0, MM_PER_DT; dt_val = dt)],
    )
    sh = Vector{Union{Wflow.SH, Missing}}([
        Wflow.read_sh_csv(joinpath(datadir, "input", "reservoir_sh_2.csv")),
    ])
    hq = Vector{Union{Wflow.HQ, Missing}}([
        Wflow.read_hq_csv(joinpath(datadir, "input", "reservoir_hq_2.csv")),
    ])
    res_params = Wflow.ReservoirParameters(;
        id = [1],
        area = [200_000_000],
        storfunc = [ReservoirProfileType.interpolation],
        outflowfunc = [ReservoirOutflowType.rating_curve],
        sh,
        hq,
        col_index_hq = [15],
    )
    res_params.maxstorage[1] = Wflow.maximum_storage(res_params, 1)
    res_vars = Wflow.ReservoirVariables(; waterlevel = [397.75], storage = [410_760_000])
    res = Wflow.Reservoir(;
        boundary_conditions = res_bc,
        parameters = res_params,
        variables = res_vars,
    )
    Wflow.set_reservoir_vars!(res)
    Wflow.update!(res, 1, 1500.0, dt)
    Wflow.average_reservoir_vars!(res, dt)
    res_p = res.parameters
    res_v = res.variables
    @test Wflow.waterlevel(
        res_p.storfunc[1],
        res_p.area[1],
        res_v.storage[1],
        res_p.sh[1],
    ) ≈ 398.0 atol = 1e-2
    @test res_v.outflow ≈ [1303.67476852]
    @test res_v.outflow_av.average ≈ res_v.outflow
    @test res_v.storage ≈ [4.293225e8]
    @test res_v.waterlevel ≈ [398.000000]
end
