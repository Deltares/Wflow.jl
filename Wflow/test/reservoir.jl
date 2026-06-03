@testitem "unit: update reservoir simple" begin
    using Wflow: ReservoirProfileType, ReservoirOutflowType
    include("testing_utils.jl")
    dt = 86400.0
    # Simple reservoir (outflow_curve_type = 4)
    n = 1
    res_bc = Wflow.ReservoirBC(;
        n,
        precipitation = [4.861111111111111e-8],
        evaporation = [1.736111111111111e-8],
    )
    res_params = Wflow.ReservoirParameters(;
        id = [1],
        demand = [52.523],
        maximum_release = [420.184],
        maximum_storage = [25_000_000.0],
        area = [1885665.353626924],
        target_full_fraction = [0.8],
        target_minimum_fraction = [0.2425554726620697],
        storage_curve_type = [ReservoirProfileType.linear],
        outflow_curve_type = [ReservoirOutflowType.simple],
    )
    res_vars = Wflow.ReservoirVariables(;
        outflow_obs = [Wflow.MISSING_VALUE],
        storage = [1.925e7],
        waterlevel = [10.208598234556407],
    )

    res = Wflow.ReservoirModel(;
        boundary_conditions = res_bc,
        parameters = res_params,
        variables = res_vars,
    )
    @testset "Update reservoir simple (outflow_curve_type = 4)" begin
        Wflow.set_reservoir_vars!(res)
        Wflow.update_reservoir_model!(res, 1, 100.0, dt)
        Wflow.average_reservoir_vars!(res, dt)
        @test res.variables.outflow[1] ≈ 91.3783714867453
        @test res.variables.outflow_cumulative[1] == res.variables.outflow_average[1] * dt
        @test res.variables.storage[1] ≈ 2.0e7
        @test res.boundary_conditions.precipitation[1] ≈ 4.861111111111111e-8
        @test res.boundary_conditions.evaporation[1] ≈ 1.736111111111111e-8
        @test res.variables.actevap_cumulative[1] ≈ 0.0014999999999999998
    end

    # reset storage and waterlevel and set observed outflow
    res.variables.outflow_obs[1] = 80.0
    res.variables.storage[1] = 1.925e7
    res.variables.waterlevel[1] = 10.208598234556407
    @testset "Update reservoir simple (outflow_curve_type = 4) with observed outflow" begin
        Wflow.set_reservoir_vars!(res)
        Wflow.update_reservoir_model!(res, 1, 100.0, dt)
        Wflow.average_reservoir_vars!(res, 86400.0)
        @test res.variables.outflow[1] ≈ 80.0
        @test res.variables.outflow_average[1] == res.variables.outflow[1]
        @test res.variables.storage[1] ≈ 2.0983091296454795e7
    end
end

@testitem "unit: update reservoir Modified Puls approach (outflow_curve_type = 3)" begin
    using Wflow: ReservoirProfileType, ReservoirOutflowType
    include("testing_utils.jl")
    # ReservoirModel Modified Puls approach (outflow_curve_type = 3)
    n = 1
    dt = 86400.0
    res_bc = Wflow.ReservoirBC(;
        n,
        precipitation = [2.3148148148148148e-7],
        evaporation = [3.7037037037037036e-8],
    )
    res_params = Wflow.ReservoirParameters(;
        id = [1],
        area = [180510409.0],
        threshold = [0.0],
        storage_curve_type = [ReservoirProfileType.linear],
        outflow_curve_type = [ReservoirOutflowType.modified_puls],
        rating_curve_coefficient = [0.22],
        rating_curve_exponent = [2.0],
    )
    waterlevel = [18.5]
    res_vars = Wflow.ReservoirVariables(;
        storage = Wflow.initialize_storage(
            res_params.storage_curve_type,
            res_params.area,
            waterlevel,
            res_params.storage_waterlevel_curve,
        ),
        waterlevel,
    )

    res = Wflow.ReservoirModel(;
        boundary_conditions = res_bc,
        parameters = res_params,
        variables = res_vars,
    )
    res_p = res.parameters
    res_v = res.variables
    res_bc = res.boundary_conditions
    Wflow.set_reservoir_vars!(res)
    Wflow.update_reservoir_model!(res, 1, 2500.0, dt)
    Wflow.average_reservoir_vars!(res, dt)
    @test Wflow.waterlevel(
        res_p.storage_curve_type[1],
        res_p.area[1],
        res_v.storage[1],
        res_p.storage_waterlevel_curve[1],
    ) ≈ 19.672653848925634
    @test res_v.outflow[1] ≈ 85.14292808113598
    @test res_v.outflow_average ≈ res_v.outflow
    @test res_v.storage[1] ≈ 3.55111879238499e9
    @test res_v.waterlevel[1] ≈ 19.672653848925634
    @test res_bc.precipitation[1] ≈ 2.3148148148148148e-7
    @test res_bc.evaporation[1] ≈ 3.7037037037037036e-8
    @test res_v.actevap_cumulative[1] ≈ 0.0032
end

@testitem "update_reservoir!" begin
    using Graphs: DiGraph, add_edge!
    include("testing_utils.jl")

    dt = 86400.0

    n = 1
    reservoir = Wflow.ReservoirModel(;
        boundary_conditions = Wflow.ReservoirBC(;
            n,
            external_inflow = [-1.0],
            inflow_overland = [0.02],
            inflow_subsurface = [0.04],
            precipitation = [5.787037037037037e-9],
            evaporation = [1.1574074074074074e-9],
        ),
        parameters = Wflow.ReservoirParameters(;
            id = [1],
            storage_curve_type = [Wflow.ReservoirProfileType.linear],
            outflow_curve_type = [Wflow.ReservoirOutflowType.simple],
            area = [6.0e4],
            threshold = [0.0],
            rating_curve_coefficient = [0.0],
            rating_curve_exponent = [0.0],
        ),
        variables = Wflow.ReservoirVariables(;
            waterlevel = [1.0],
            storage = [4.5e7],
            outflow = [3.0],
            outflow_average = [2.0],
            outflow_obs = [1.0],
            actevap_average = [1.1574074074074074e-10],
        ),
    )

    river_flow_vars = Wflow.RiverFlowVariables(; n = 2, q = [0.04, 0.04])

    graph = DiGraph(2)
    add_edge!(graph, 1, 2)
    network = Wflow.NetworkRiver(; graph, reservoir_indices = [1])

    v = 1
    dt = 1000.0

    Wflow.update_reservoir_model!(reservoir, river_flow_vars, network, v, dt)
    @test river_flow_vars.qin[2] ≈ 1.0
    @test reservoir.boundary_conditions.actual_external_abstraction_cumulative[1] ≈ 1e3
    @test reservoir.variables.storage[1] ≈ 4.4998100277777776e7
    @test reservoir.variables.waterlevel[1] ≈ 0.9683379629629354
    @test reservoir.variables.outflow[1] ≈ 1.0
end

@testitem "unit: update_reservoir_model!" begin
    using Graphs: DiGraph, add_edge!

    dt = 86400.0

    n = 1
    reservoir_model = Wflow.ReservoirModel(;
        boundary_conditions = Wflow.ReservoirBC(;
            n,
            external_inflow = [-1.0],
            inflow_overland = [0.0],
            inflow_subsurface = [0.00041241066945499203],
            precipitation = [8.101853611016715e-10],
            evaporation = [6.134257548385196e-9],
        ),
        parameters = Wflow.ReservoirParameters(;
            id = [1],
            storage_curve_type = [Wflow.ReservoirProfileType.linear],
            outflow_curve_type = [Wflow.ReservoirOutflowType.simple],
            area = [9.069779e4],
            maximum_release = [1.74],
            demand = [0.2175],
            target_minimum_fraction = [0.358469158],
            target_full_fraction = [0.83492106199],
            maximum_storage = [3.3e7],
        ),
        variables = Wflow.ReservoirVariables(;
            waterlevel = [3.0266425035195113],
            storage = [2.7450978618928656e7],
        ),
    )

    n_river = 2
    river_flow_vars = Wflow.RiverFlowVariables(;
        n = n_river,
        q = [0.00012002923701686638, 0.21747539140212965],
    )

    graph = DiGraph(2)
    add_edge!(graph, 1, 2)
    network = Wflow.NetworkRiver(; graph, reservoir_indices = [1])

    v = 1
    dt = 1000.0

    Wflow.update_reservoir_model!(reservoir_model, river_flow_vars, network, v, dt)
    @test river_flow_vars.qin[2] ≈ 0.21749985206208133
    @test reservoir_model.boundary_conditions.actual_external_abstraction_cumulative[1] ≈
          1000.0
    @test reservoir_model.variables.storage[1] ≈ 2.744976116863499e7
    @test reservoir_model.variables.waterlevel[1] ≈ 3.013219350720886
    @test reservoir_model.variables.outflow[1] ≈ 0.21749985206208133
end

@testitem "Linked reservoirs with free weir (outflow_curve_type = 2)" begin
    using Wflow: ReservoirProfileType, ReservoirOutflowType
    include("testing_utils.jl")
    dt = 86400.0
    # Linked reservoirs with free weir (outflow_curve_type = 1)
    datadir = joinpath(@__DIR__, "data")
    storage_waterlevel_curve = Vector{Union{Wflow.SH, Missing}}([
        Wflow.read_sh_csv(joinpath(datadir, "input", "reservoir_sh_1.csv")),
        Wflow.read_sh_csv(joinpath(datadir, "input", "reservoir_sh_2.csv")),
    ])
    waterlevel_discharge_curve = Vector{Union{Wflow.HQ, Missing}}([
        missing,
        Wflow.read_hq_csv(joinpath(datadir, "input", "reservoir_hq_2.csv")),
    ])

    @test keys(storage_waterlevel_curve[1]) == (:H, :S)
    @test typeof(values(storage_waterlevel_curve[1])) == Tuple{Vector{Float64}, Vector{Float64}}

    res_params = Wflow.ReservoirParameters(;
        id = [1, 2],
        lower_reservoir_ind = [2, 0],
        area = [472461536.0, 60851088.0],
        threshold = [393.7, NaN],
        storage_curve_type = fill(ReservoirProfileType.interpolation, 2),
        outflow_curve_type = [ReservoirOutflowType.free_weir, ReservoirOutflowType.rating_curve],
        rating_curve_coefficient = [140.0, NaN],
        rating_curve_exponent = [1.5, NaN],
        storage_waterlevel_curve,
        waterlevel_discharge_curve,
        col_index_hq = [15],
    )
    res_params.maximum_storage[2] = Wflow.maximum_storage(res_params, 2)

    waterlevel = [395.03027, 394.87833]
    res_vars = Wflow.ReservoirVariables(;
        waterlevel,
        storage = Wflow.initialize_storage(
            res_params.storage_curve_type,
            [472461536.0, 60851088.0],
            waterlevel,
            storage_waterlevel_curve,
        ),
    )
    n = 2
    res_bc = Wflow.ReservoirBC(;
        n = 2,
        precipitation = [1.1574074074074074e-7, 1.1574074074074074e-7],
        evaporation = [2.3148148148148148e-8, 2.3148148148148148e-8],
    )

    res = Wflow.ReservoirModel(;
        boundary_conditions = res_bc,
        parameters = res_params,
        variables = res_vars,
    )
    Wflow.set_reservoir_vars!(res)
    Wflow.update_reservoir_model!(res, 1, 500.0, dt)
    Wflow.update_reservoir_model!(res, 2, 500.0, dt)
    Wflow.average_reservoir_vars!(res, dt)
    res_v = res.variables
    res_bc = res.boundary_conditions
    @test res_v.outflow ≈ [214.80170846121263, 236.83281600000214]
    @test res_v.outflow_average ≈ res_v.outflow
    @test res_v.storage ≈ [1.2737435094769483e9, 2.6019755340159863e8]
    Wflow.set_reservoir_vars!(res)
    Wflow.update_reservoir_model!(res, 1, 500.0, dt)
    Wflow.update_reservoir_model!(res, 2, 500.0, dt)
    Wflow.average_reservoir_vars!(res, 86400.0)
    @test res_v.outflow ≈ [-259.8005149014703, 239.66710359986183]
    @test res_v.outflow_average ≈ [-259.8005149014703, 499.4676185013321]
    @test res_v.storage ≈ [1.3431699662524352e9, 2.6073035986708355e8]
    @test res_v.waterlevel ≈ [395.239782021054, 395.21771942667266]
    @test res_v.actevap_cumulative ≈ [0.002, 0.002]
end

# Overflowing reservoir with SH and HQ (outflow_curve_type = 1)
@testitem "Overflowing reservoir with SH and HQ" begin
    include("testing_utils.jl")
    using Wflow: ReservoirProfileType, ReservoirOutflowType
    datadir = joinpath(@__DIR__, "data")
    n = 1
    dt = 86400.0
    res_bc = Wflow.ReservoirBC(;
        n,
        precipitation = [1.1574074074074074e-7],
        evaporation = [2.3148148148148148e-8],
    )
    storage_waterlevel_curve = Vector{Union{Wflow.SH, Missing}}([
        Wflow.read_sh_csv(joinpath(datadir, "input", "reservoir_sh_2.csv")),
    ])
    waterlevel_discharge_curve = Vector{Union{Wflow.HQ, Missing}}([
        Wflow.read_hq_csv(joinpath(datadir, "input", "reservoir_hq_2.csv")),
    ])
    res_params = Wflow.ReservoirParameters(;
        id = [1],
        area = [200_000_000],
        storage_curve_type = [ReservoirProfileType.interpolation],
        outflow_curve_type = [ReservoirOutflowType.rating_curve],
        storage_waterlevel_curve,
        waterlevel_discharge_curve,
        col_index_hq = [15],
    )
    res_params.maximum_storage[1] = Wflow.maximum_storage(res_params, 1)
    res_vars = Wflow.ReservoirVariables(; waterlevel = [397.75], storage = [410_760_000])
    res = Wflow.ReservoirModel(;
        boundary_conditions = res_bc,
        parameters = res_params,
        variables = res_vars,
    )
    Wflow.set_reservoir_vars!(res)
    Wflow.update_reservoir_model!(res, 1, 1500.0, dt)
    Wflow.average_reservoir_vars!(res, dt)
    res_p = res.parameters
    res_v = res.variables
    @test Wflow.waterlevel(
        res_p.storage_curve_type[1],
        res_p.area[1],
        res_v.storage[1],
        res_p.storage_waterlevel_curve[1],
    ) ≈ 398.0 atol = 1e-2
    @test res_v.outflow ≈ [1303.67476852]
    @test res_v.outflow_average ≈ res_v.outflow
    @test res_v.storage ≈ [4.293225e8]
    @test res_v.waterlevel ≈ [398.000000]
end
