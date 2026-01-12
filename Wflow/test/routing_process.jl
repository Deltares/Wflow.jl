@testitem "kinematic wave overland flow" begin
    using NCDatasets: NCDataset
    using Graphs: topological_sort_by_dfs

    dt_sec = 86400.0
    ldd_MISSING_VALUE = 255

    # read the staticmaps into memory
    nc = NCDataset("data/input/staticmaps-rhine.nc")
    # helper function to get the axis order and directionality right
    read_right(nc, var) = reverse(permutedims(Array(nc[var])); dims = 2)
    ldd_2d = read_right(nc, "ldd")

    inds, _ = Wflow.active_indices(ldd_2d, ldd_MISSING_VALUE)
    n = length(inds)

    # take out only the active cells
    ldd = ldd_2d[inds]
    slope = read_right(nc, "slope")[inds]
    N = read_right(nc, "N")[inds]
    Qold = read_right(nc, "Qold")[inds]
    Bw = read_right(nc, "Bw")[inds]
    waterlevel = read_right(nc, "waterlevel")[inds]
    DCL = read_right(nc, "DCL")[inds]
    close(nc)

    # create the directed acyclic graph from the drainage direction array
    graph = Wflow.flowgraph(ldd, inds, Wflow.PCR_DIR)
    # a topological sort is used for visiting nodes in order from upstream to downstream
    toposort = topological_sort_by_dfs(graph)
    sink = toposort[end]
    @test ldd[sink] == 5  # the most downstream node must be a sink

    # calculate parameters of kinematic wave
    q = 0.000001
    beta = 0.6
    AlpPow = (2.0 / 3.0) * beta
    AlpTermR = (N ./ sqrt.(slope)) .^ beta
    P = Bw + (2.0 * waterlevel)
    alpha = AlpTermR .* P .^ AlpPow

    Q = zeros(n)
    Q = Wflow.kin_wave!(Q, graph, toposort, Qold, q, alpha, beta, DCL, dt_sec)

    @test sum(Q) ≈ 2.957806043289641e6
    @test Q[toposort[1]] ≈ 0.007260052312634069
    @test Q[toposort[n - 100]] ≈ 3945.762718338739
    @test Q[sink] ≈ 4131.101474418251
end

@testitem "kinematic wave subsurface flow" begin
    using StaticArrays: SVector
    kh_exp_profile = Wflow.KhExponential([205.5965576171875], [1.0141291422769427])

    ustorelayerthickness = fill(SVector{4}([100.0, 300.0, 119.83408703759733, NaN]), 1)
    ustorelayerdepth = fill(
        SVector{4}([0.1909439890049523, 16.27933934181815, 19.508197676020185, 0.0]),
        1,
    )
    n_unsatlayers = fill(3, 1)
    zi = fill(519.8340870375973, 1)
    theta_s = fill(0.48642662167549133, 1)
    theta_fc = fill(0.28219206182657536, 1)
    theta_r = fill(0.11939866840839386, 1)
    sumlayers = fill(SVector{5}([0.0, 100.0, 400.0, 1200.0, 2000.0]), 1)
    act_thickl = fill(SVector{4}([100.0, 300.0, 800.0, 800.0]), 1)

    variables = (; ustorelayerthickness, ustorelayerdepth, n_unsatlayers, zi)
    parameters = (; theta_s, theta_fc, theta_r, sumlayers, act_thickl)
    soil = (; parameters, variables)

    @test all(
        isapprox.(
            Wflow.kinematic_wave_ssf(
                0.0,
                25953.147860945584,
                0.5198340870375974,
                0.4346106913943182,
                0.4522336721420288,
                0.20423455984891598,
                2.0,
                1.0,
                1117.0150713112287,
                517.495693771673,
                79.62016166711079,
                kh_exp_profile,
                soil,
                1,
            ),
            (23675.53215503045, 0.7162637030758906, 0.0, 0.20423455984891598),
        ),
    )
    kh_exp_const_profile = Wflow.KhExponentialConstant(kh_exp_profile, [0.2])
    ustorelayerthickness = fill(SVector{4}([100.0, 300.0, 348.31246153148595, NaN]), 1)
    ustorelayerdepth =
        fill(SVector{4}([0.1909439890049523, 16.27933934181815, 58.42501219303608, 0.0]), 1)
    n_unsatlayers = fill(3, 1)
    zi = fill(758.8905603985703, 1)
    theta_s = fill(0.48642662167549133, 1)
    theta_fc = fill(0.28219206182657536, 1)
    theta_r = fill(0.11939866840839386, 1)
    sumlayers = fill(SVector{5}([0.0, 100.0, 400.0, 1200.0, 2000.0]), 1)
    act_thickl = fill(SVector{4}([100.0, 300.0, 800.0, 800.0]), 1)

    variables = (; ustorelayerthickness, ustorelayerdepth, n_unsatlayers, zi)
    parameters = (; theta_s, theta_fc, theta_r, sumlayers, act_thickl)
    soil = (; parameters, variables)

    @test all(
        isapprox.(
            Wflow.kinematic_wave_ssf(
                0.0,
                54175.65003911068,
                0.7588905603985703,
                0.6928420612599803,
                0.4522336721420288,
                0.20423455984891598,
                2.0,
                1.0,
                1117.0150713112287,
                517.495693771673,
                153.46698446681825,
                kh_exp_const_profile,
                soil,
                1,
            ),
            (48524.884193017664, 1.163361369433857, 0.0, 0.20423455984891598),
        ),
    )
end

@testitem "accucapacity" begin
    using Graphs: DiGraph, add_edge!
    # test based on a subset of the examples at
    # https://pcraster.geo.uu.nl/pcraster/4.3.0/documentation/pcraster_manual/sphinx/op_accucapacity.html#examples
    # of the node at (row 3, column 2) and upstream nodes
    g = DiGraph(6)
    add_edge!(g, 1, 4)
    add_edge!(g, 2, 5)
    add_edge!(g, 3, 5)
    add_edge!(g, 4, 6)
    add_edge!(g, 5, 6)
    network = (graph = g, order = [1, 2, 3, 4, 5, 6])

    # example 1, accucapacityflux
    material = Float64[0.5, 2, 2, 0.5, 2, 0.5]
    capacity = fill(1.5, 6)
    flux, new_material = Wflow.accucapacityflux_state(material, network, capacity)
    @test new_material != material
    @test new_material == [0.0, 0.5, 0.5, 0.0, 3.5, 1.5]
    @test flux == Float64[0.5, 1.5, 1.5, 1, 1.5, 1.5]
    flux_ = Wflow.accucapacityflux(material, network, capacity)
    @test flux == flux_

    # example 2, accucapacityflux
    material = fill(10.0, 6)
    capacity = Float64[2, 30, 30, 2, 30, 2]
    flux, new_material = Wflow.accucapacityflux_state(material, network, capacity)
    @test new_material != material
    @test new_material == [8.0, 0.0, 0.0, 10.0, 0.0, 40.0]
    @test flux == Float64[2, 10, 10, 2, 30, 2]
    flux_ = Wflow.accucapacityflux(material, network, capacity)
    @test flux == flux_

    # example 1, accucapacitystate
    material = Float64[0.5, 2, 2, 0.5, 2, 0.5]
    capacity = fill(1.5, 6)
    new_material = Wflow.accucapacitystate(material, network, capacity)
    @test new_material != material
    @test new_material == [0.0, 0.5, 0.5, 0.0, 3.5, 1.5]

    # example 2, accucapacitystate
    material = fill(10.0, 6)
    capacity = Float64[2, 30, 30, 2, 30, 2]
    new_material = Wflow.accucapacitystate(material, network, capacity)
    @test new_material != material
    @test new_material == [8.0, 0.0, 0.0, 10.0, 0.0, 40.0]
end

@testitem "local inertial long channel MacDonald (1997)" begin
    using QuadGK: quadgk
    using Graphs: DiGraph, add_edge!, ne
    using Statistics: mean
    using Wflow: GRAVITATIONAL_ACCELERATION

    L = 1000.0
    dx = 5.0
    n = Int(L / dx)

    # analytical solution MacDonald (1997) for channel with length L of 1000.0 m, Manning's
    # n of 0.03, constant inflow of 20.0 m3/s at upper boundary and channel width of 10.0 m
    # water depth profile h(x)
    h(x) = cbrt(4 / GRAVITATIONAL_ACCELERATION) * (1.0 + 0.5 * exp(-16.0 * (x / L - 0.5)^2))
    # spatial derivative of h(x)
    h_acc(x) =
        -cbrt(4 / GRAVITATIONAL_ACCELERATION) * 16.0 / L *
        (x / L - 0.5) *
        exp(-16 * (x / L - 0.5)^2)
    # solution for channel slope s(x)
    s(x) =
        (1.0 - 4.0 / (GRAVITATIONAL_ACCELERATION * h(x)^3)) * h_acc(x) +
        0.36 * (2 * h(x) + 10.0)^(4.0 / 3.0) / ((10.0 * h(x))^(10.0 / 3.0))

    h_a = h.([dx:dx:L;]) # water depth profile (analytical solution)
    # integrate slope to get elevation (bed level) z
    x = [dx:dx:L;]
    zb = first.([quadgk(s, xi, L; rtol = 1e-12) for xi in x])

    # initialize local inertial river flow model
    graph = DiGraph(n)
    for i in 1:n
        add_edge!(graph, i, i + 1)
    end

    dl = fill(dx, n)
    width = fill(10.0, n)
    n_river = fill(0.03, n)

    # for each edge the src and dst node is required
    nodes_at_edge = Wflow.adjacent_nodes_at_edge(graph)
    _ne = ne(graph)

    # determine z, width, length and manning's n at edges
    zb_max = fill(0.0, _ne)
    width_at_edge = fill(0.0, _ne)
    length_at_edge = fill(0.0, _ne)
    mannings_n_sq = fill(0.0, _ne)
    for i in 1:_ne
        zb_max[i] = max(zb[nodes_at_edge.src[i]], zb[nodes_at_edge.dst[i]])
        width_at_edge[i] = min(width[nodes_at_edge.dst[i]], width[nodes_at_edge.src[i]])
        length_at_edge[i] = 0.5 * (dl[nodes_at_edge.dst[i]] + dl[nodes_at_edge.src[i]])
        mannings_n =
            (
                n_river[nodes_at_edge.dst[i]] * dl[nodes_at_edge.dst[i]] +
                n_river[nodes_at_edge.src[i]] * dl[nodes_at_edge.src[i]]
            ) / (dl[nodes_at_edge.dst[i]] + dl[nodes_at_edge.src[i]])
        mannings_n_sq[i] = mannings_n * mannings_n
    end

    river_network = Wflow.NetworkRiver(;
        nodes_at_edge = Wflow.NodesAtEdge(; nodes_at_edge...),
        edges_at_node = Wflow.EdgesAtNode(;
            Wflow.adjacent_edges_at_node(graph, nodes_at_edge)...,
        ),
    )

    params_river = Wflow.RiverParameters(;
        flow_width = width,
        flow_length = dl,
        reservoir_outlet = zeros(n),
    )
    domain_river = Wflow.DomainRiver(; network = river_network, parameters = params_river)
    domain = Wflow.Domain(; river = domain_river)

    h_thresh = 1.0e-03
    froude_limit = true
    h_init = zeros(n - 1)
    push!(h_init, h_a[n])

    timestepping = Wflow.TimeStepping(; cfl = 0.7)
    parameters = Wflow.LocalInertialRiverFlowParameters(;
        n,
        ne = _ne,
        active_n = collect(1:(n - 1)),
        active_e = collect(1:_ne),
        h_thresh,
        zb_max,
        mannings_n_sq,
        mannings_n = n_river,
        flow_width_at_edge = width_at_edge,
        flow_length_at_edge = length_at_edge,
        bankfull_storage = fill(Wflow.MISSING_VALUE, n),
        bankfull_depth = fill(Wflow.MISSING_VALUE, n),
        zb,
        froude_limit,
    )

    variables = Wflow.LocalInertialRiverFlowVariables(;
        n,
        n_edges = _ne,
        q0 = zeros(_ne),
        q = zeros(_ne),
        q_av = zeros(_ne),
        q_channel_av = zeros(_ne),
        h = h_init,
        zs_max = zeros(_ne),
        zs_src = zeros(_ne),
        zs_dst = zeros(_ne),
        hf = zeros(_ne),
        a = zeros(_ne),
        r = zeros(_ne),
        storage = fill(0.0, n),
        error = zeros(n),
    )

    boundary_conditions = Wflow.RiverFlowBC(; n, reservoir = nothing)

    sw_river = Wflow.LocalInertialRiverFlow(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain = nothing,
        allocation = Wflow.NoAllocationRiver(n),
    )

    # run until steady state is reached
    epsilon = 1.0e-12
    (; flow_length) = domain_river.parameters
    while true
        sw_river.boundary_conditions.inwater[1] = 20.0
        h0 = mean(sw_river.variables.h)
        dt = Wflow.stable_timestep(sw_river, flow_length)
        Wflow.local_inertial_river_update!(sw_river, domain, dt, 86400.0, true)
        d = abs(h0 - mean(sw_river.variables.h))
        if d <= epsilon
            break
        end
    end

    # test for mean absolute error [cm]
    @test mean(abs.(sw_river.variables.h .- h_a)) * 100.0 ≈ 1.873574206931199
end
