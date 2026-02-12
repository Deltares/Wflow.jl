@testitem "unit: kinematic_wave" begin
    alpha = 2.586
    beta = 0.6
    dt = 600.0
    dx = 1061.375

    # Case !(q_in + q_prev + q_lat ≈ 0.0)
    q_in = 1.104e-6
    q_prev = 0.0
    q_lat = 1.142e-6
    @test Wflow.kinematic_wave(q_in, q_prev, q_lat, alpha, beta, dt, dx) ≈
          1.09308660753423e-6

    # Case q_in + q_prev + q_lat ≈ 0.0
    q_in = 0.0
    q_prev = 0.0
    q_lat = 0.0
    @test Wflow.kinematic_wave(q_in, q_prev, q_lat, alpha, beta, dt, dx) == 0.0
end

@testitem "unit: ssf_celerity" begin
    zi = 0.3
    theta_e = 0.274
    slope = 0.00586
    i = 1

    # Case kh_profile::KhExponential
    kh_profile = Wflow.KhExponential([24.152037048339846], [1.8001038115471601])
    @test Wflow.ssf_celerity(zi, slope, theta_e, kh_profile, i) ≈ 0.3010012323985728

    # Case kh_profile::KhExponentialConstant
    kh_profile = Wflow.KhExponentialConstant(kh_profile, [0.2])
    @test Wflow.ssf_celerity(zi, slope, theta_e, kh_profile, i) ≈ 0.3603676427614705
end

@testitem "unit: kw_ssf_newton_raphson" begin
    ssf = 754.993
    constant_term = 77.774
    celerity = 12.254
    dt = 1.0
    dx = 1103.816

    @test Wflow.kw_ssf_newton_raphson(ssf, constant_term, celerity, dt, dx) ≈
          942.5785713676884
end

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

@testitem "unit: kinematic_wave_ssf" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    n = 1
    N = 4

    soil = init_sbm_soil_model(
        n,
        N;
        # Variables
        ustorelayerthickness = [SVector(100.0, 300.0, 119.83408703759733, NaN)],
        ustorelayerdepth = [
            SVector(0.1909439890049523, 16.27933934181815, 19.508197676020185, 0.0),
        ],
        n_unsatlayers = [3],
        zi = [519.8340870375973],
        # Parameters
        maxlayers = 4,
        sumlayers = [SVector(0.0, 100.0, 400.0, 1200.0, 2000.0)],
        nlayers = [4],
        theta_s = [0.48642662167549133],
        theta_r = [0.11939866840839386],
        theta_fc = [0.28219206182657536],
        act_thickl = [SVector(100.0, 300.0, 800.0, 800.0)],
    )

    ssfin = 0.0
    ssf_prev = 25953.147860945584
    zi_prev = 0.5198340870375974
    r = 0.4346106913943182
    slope = 0.4522336721420288
    sy = 0.20423455984891598
    d = 2.0
    dt = 1.0
    dx = 1117.0150713112287
    dw = 517.495693771673
    ssfmax = 79.62016166711079
    kh_profile = Wflow.KhExponential([205.5965576171875], [1.0141291422769427])
    i = 1

    ssf, zi, exfilt, sy_d = Wflow.kinematic_wave_ssf(
        ssfin,
        ssf_prev,
        zi_prev,
        r,
        slope,
        sy,
        d,
        dt,
        dx,
        dw,
        ssfmax,
        kh_profile,
        soil,
        i,
    )
    @test ssf ≈ 22100.628024231868
    @test zi ≈ 0.7029236021516849
    @test exfilt ≈ 0.0
    @test sy_d ≈ 0.20423455984891598

    soil = init_sbm_soil_model(
        n,
        N;
        # Variables
        ustorelayerthickness = [SVector(100.0, 300.0, 348.31246153148595, NaN)],
        ustorelayerdepth = [
            SVector(0.1909439890049523, 16.27933934181815, 58.42501219303608, 0.0),
        ],
        n_unsatlayers = [3],
        zi = [758.8905603985703],
        # Parameters
        maxlayers = 4,
        sumlayers = [SVector(0.0, 100.0, 400.0, 1200.0, 2000.0)],
        nlayers = [4],
        theta_s = [0.48642662167549133],
        theta_r = [0.11939866840839386],
        theta_fc = [0.28219206182657536],
        act_thickl = [SVector(100.0, 300.0, 800.0, 800.0)],
    )

    ssfin = 0.0
    ssf_prev = 54175.65003911068
    zi_prev = 0.7588905603985703
    r = 0.6928420612599803
    slope = 0.4522336721420288
    sy = 0.20423455984891598
    d = 2.0
    dt = 1.0
    dx = 1117.0150713112287
    dw = 517.495693771673
    ssfmax = 153.46698446681825
    kh_profile = Wflow.KhExponentialConstant(kh_profile, [0.2])
    i = 1

    ssf, zi, exfilt, sy_d = Wflow.kinematic_wave_ssf(
        ssfin,
        ssf_prev,
        zi_prev,
        r,
        slope,
        sy,
        d,
        dt,
        dx,
        dw,
        ssfmax,
        kh_profile,
        soil,
        1,
    )

    @test ssf ≈ 44680.57723298823
    @test zi ≈ 1.130798471269119
    @test exfilt ≈ 0.0
    @test sy_d ≈ 0.20423455984891598
end

@testitem "unit: accucapacity" begin
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
    dt = 86400.0
    material = dt .* [0.5, 2.0, 2.0, 0.5, 2.0, 0.5]
    capacity = fill(1.5, 6)
    flux, new_material = Wflow.accucapacityflux_state(material, network, capacity, dt)
    @test new_material != material
    @test new_material == dt .* [0.0, 0.5, 0.5, 0.0, 3.5, 1.5]
    @test flux == Float64[0.5, 1.5, 1.5, 1, 1.5, 1.5]
    flux_ = Wflow.accucapacityflux(material, network, capacity, dt)
    @test flux == flux_

    # example 2, accucapacityflux
    material = dt .* fill(10.0, 6)
    capacity = Float64[2, 30, 30, 2, 30, 2]
    flux, new_material = Wflow.accucapacityflux_state(material, network, capacity, dt)
    @test new_material != material
    @test new_material == dt .* [8.0, 0.0, 0.0, 10.0, 0.0, 40.0]
    @test flux == Float64[2, 10, 10, 2, 30, 2]
    flux_ = Wflow.accucapacityflux(material, network, capacity, dt)
    @test flux == flux_
end

@testitem "unit: local_inertial_river_update!" begin
    using Wflow: to_SI, MM_PER_DT, MM
    dt = 86400.0
    n = 2
    model = Wflow.LocalInertialRiverFlow(;
        timestepping = Wflow.TimeStepping(),
        boundary_conditions = Wflow.RiverFlowBC(;
            n,
            external_inflow = [-1.0, -1.0],
            inwater = [100.0, 100.0],
            reservoir = Wflow.Reservoir(;
                boundary_conditions = Wflow.ReservoirBC(;
                    n,
                    external_inflow = [-1.0],
                    inflow_overland = [3000.0],
                    inflow_subsurface = [5000.0],
                    precipitation = [to_SI(2.0, MM_PER_DT; dt_val = dt)],
                    evaporation = [to_SI(1.0, MM_PER_DT; dt_val = dt)],
                    inflow = Wflow.AverageVector(; n),
                ),
                parameters = Wflow.ReservoirParameters(;
                    id = [1, 2],
                    storfunc = [Wflow.ReservoirProfileType.linear],
                    outflowfunc = [Wflow.ReservoirOutflowType.simple],
                    area = [5.0e6, 4.0e6],
                    maxrelease = [10.0, 10.0],
                    demand = [1.5, 1.5],
                    targetminfrac = [0.3, 0.3],
                    targetfullfrac = [0.1, 0.1],
                    maxstorage = [Inf, Inf],
                    threshold = [0.0, 0.0],
                ),
                variables = Wflow.ReservoirVariables(;
                    waterlevel = [1.0, 1.0],
                    storage = [2.5e5, 2.5e5],
                    outflow = [1.8, 1.8],
                ),
            ),
        ),
        parameters = Wflow.LocalInertialRiverFlowParameters(;
            n = n,
            ne = 2,
            active_n = [2],
            active_e = [1],
            froude_limit = true,
            h_thresh = 0.0,
            zb = [0.0, 0.0],
            zb_max = [0.5, 0.5],
            bankfull_storage = [1e3, 1e3],
            bankfull_depth = [0.5, 0.5],
            mannings_n_sq = [1e-3, 1e-3],
            mannings_n = [1e-2, 1e-2],
            flow_length_at_edge = [1000.0, 1000.0],
            flow_width_at_edge = [100.0, 100.0],
        ),
        variables = Wflow.LocalInertialRiverFlowVariables(;
            n_cells = n,
            n_edges = 2,
            h = [1.0, 2.0],
            q = [0.0, 1e-4],
        ),
        floodplain = Wflow.FloodPlain(;
            parameters = Wflow.FloodPlainParameters(;
                profile = Wflow.FloodPlainProfile(;
                    depth = [10.0, 10.0],
                    storage = [1e5 1e6; 1e5 1e6],
                    width = [100.0 100.0; 100.0 100.0],
                    a = [1e3 1e3; 1e3 1e3],
                    p = [400.0 400.0; 400.0 400.0],
                ),
                mannings_n = [0.04],
                mannings_n_sq = [1.2e-4],
                zb_max = [1.0],
            ),
            variables = Wflow.FloodPlainVariables(; n, n_edges = 1, h = [0.1, 0.2]),
        ),
        allocation = Wflow.AllocationRiver(; n),
    )
    domain = Wflow.Domain(;
        river = Wflow.DomainRiver(;
            network = Wflow.NetworkRiver(;
                nodes_at_edge = Wflow.NodesAtEdge(; src = [1], dst = [2]),
                edges_at_node = Wflow.EdgesAtNode(; src = [[1], [1]], dst = [[1], [1]]),
            ),
            parameters = Wflow.RiverParameters(;
                flow_width = [10.0, 10.0],
                flow_length = [100.0, 100.0],
            ),
        ),
        reservoir = Wflow.DomainReservoir(;
            network = Wflow.NetworkReservoir(; river_indices = [1]),
        ),
    )
    dt = 1000.0

    Wflow.update_river_channel_flow!(model, domain.river, dt)

    @test model.variables.zs_src[1] ≈ 1.0
    @test model.variables.zs_dst[1] ≈ 2.0
    @test model.variables.zs_max[1] ≈ 2.0
    @test model.variables.hf[1] ≈ 1.5
    @test model.variables.a[1] ≈ 150.0
    @test model.variables.r[1] ≈ 1.4563106796116505
    @test model.variables.q[1] ≈ -575.3037784510024
    @test model.variables.q_av.cumulative_material[1] ≈ -575303.7784510023

    Wflow.update_floodplain_flow!(model, domain.river, dt)

    @test model.floodplain.variables.hf[1] ≈ 1.0
    @test model.floodplain.variables.hf_index[1] == 1
    @test model.floodplain.variables.a[1] ≈ 90.0
    @test model.floodplain.variables.r[1] ≈ 0.2356020942408377
    @test model.floodplain.variables.q[1] ≈ -281.84014086002725
    @test model.floodplain.variables.q_av.cumulative_material[1] ≈ -281840.1408600272

    Wflow.update_boundary_conditions_reservoir!(model, domain, dt)

    @test model.boundary_conditions.reservoir.variables.storage[1] ≈ 7.391913765967477e6
    @test model.boundary_conditions.reservoir.variables.waterlevel[1] ≈ 2.428382753193495
    @test model.boundary_conditions.reservoir.variables.outflow[1] ≈ 0.00018509186397934759
    @test model.boundary_conditions.reservoir.boundary_conditions.inflow.cumulative_material[1] ≈
          7.141856080688971e6
    @test model.boundary_conditions.reservoir.variables.outflow_av.cumulative_material[1] ≈
          0.1850918639793476
    @test model.boundary_conditions.reservoir.variables.actevap.cumulative_material[1] ≈
          to_SI(0.011574074074074073, MM)
    @test model.variables.q[1] ≈ 0.00018509186397934759
    @test model.variables.q_av.cumulative_material[1] ≈ -575303.5933591384

    Wflow.update_water_depth_and_storage!(model, domain.river, dt)

    @test model.variables.storage[2] ≈ 99000.0
    @test model.variables.h[2] ≈ 99.0

    Wflow.update_water_depth_and_storage!(model.floodplain, model, domain.river, dt)

    @test model.variables.h ≈ [1.0, -79.7]
    @test model.variables.storage ≈ [0.0, -79700.0]
    @test model.floodplain.variables.storage ≈ [0.0, 178700.0]
end

@testitem "unit: update_directional_flow!" begin
    n = 3
    land = Wflow.LocalInertialOverlandFlow(;
        timestepping = Wflow.TimeStepping(),
        boundary_conditions = Wflow.LocalInertialOverlandFlowBC(; n),
        parameters = Wflow.LocalInertialOverlandFlowParameters(;
            n,
            ywidth = fill(900.0, n),
            xwidth = [250.0, 300.0, 450.0],
            zx_max = [750.0, 900.0, 800.0],
            theta = 1.0,
            h_thresh = 1e-3,
            zy_max = [800.0, 800.0, 800.0],
            mannings_n_sq = [0.06, 0.06, 0.06],
            z = [800.0, 800.0, 800.0],
            froude_limit = true,
        ),
        variables = Wflow.LocalInertialOverlandFlowVariables(;
            n,
            qx0 = [1e-3, 2e-3, 3e-3],
            h = [0.03, 0.02, 0.05],
        ),
    )
    domain = Wflow.Domain(;
        land = Wflow.DomainLand(;
            network = Wflow.NetworkLand(;
                edge_indices = Wflow.EdgeConnectivity(; n = 1, xu = [2], xd = [3]),
            ),
            parameters = Wflow.LandParameters(;
                x_length = fill(600.0, n),
                y_length = fill(900.0, n),
            ),
        ),
    )
    i = 1
    dt = 60.0
    is_x_direction = true

    Wflow.update_directional_flow!(land, domain, i, dt, is_x_direction)
    @test land.variables.qx_av.cumulative_material[1] ≈ 26493.90166029366
end

@testitem "unit: local_inertial_update_water_depth!" begin
    n = 2
    land = Wflow.LocalInertialOverlandFlow(;
        timestepping = Wflow.TimeStepping(),
        variables = Wflow.LocalInertialOverlandFlowVariables(;
            n,
            qx = [0.1, 0.3],
            qy = [0.25, 0.15],
            storage = [1000.0, 1250.0],
        ),
        boundary_conditions = Wflow.LocalInertialOverlandFlowBC(; n, runoff = [0.2, 0.3]),
        parameters = Wflow.LocalInertialOverlandFlowParameters(;
            n,
            xwidth = [600.0],
            ywidth = [900.0],
            theta = 1.0,
            h_thresh = 1e-3,
            zx_max = [800.0],
            zy_max = [800.0],
            mannings_n_sq = [0.06],
            z = [700.0],
            froude_limit = true,
        ),
    )
    river = Wflow.LocalInertialRiverFlow(;
        timestepping = Wflow.TimeStepping(),
        boundary_conditions = Wflow.RiverFlowBC(;
            n,
            external_inflow = [-0.2, -0.1],
            reservoir = nothing,
        ),
        parameters = Wflow.LocalInertialRiverFlowParameters(;
            n,
            ne = 2,
            active_n = [1, 1],
            active_e = [1, 1],
            froude_limit = true,
            h_thresh = 1e-3,
            zb = [500.0, 500.0],
            zb_max = [400.0, 400.0],
            bankfull_storage = [500.0, 500.0],
            bankfull_depth = [1.0, 1.0],
            mannings_n_sq = [9.0e-3, 9.0e-3],
            mannings_n = [0.03, 0.03],
            flow_length_at_edge = [800.0, 800.0],
            flow_width_at_edge = [30.0, 30.0],
        ),
        variables = Wflow.LocalInertialRiverFlowVariables(;
            n_cells = n,
            n_edges = 2,
            q = [0.03, 0.04],
            h = [1.0, 1.0],
            storage = [30.0e3, 25e3],
        ),
        floodplain = nothing,
        allocation = Wflow.NoAllocationRiver(1),
    )
    domain = Wflow.Domain(;
        land = Wflow.DomainLand(;
            network = Wflow.NetworkLand(;
                river_indices = [1],
                edge_indices = Wflow.EdgeConnectivity(; n, xd = [2, 1], yd = [2, 1]),
            ),
            parameters = Wflow.LandParameters(;
                x_length = [600.0, 600.0],
                y_length = [900.0, 900.0],
            ),
        ),
        river = Wflow.DomainRiver(;
            network = Wflow.NetworkRiver(;
                edges_at_node = Wflow.EdgesAtNode(; src = [[1]], dst = [[2]]),
            ),
            parameters = Wflow.RiverParameters(;
                flow_width = [30.0],
                flow_length = [800.0],
            ),
        ),
    )

    dt = 1000.0

    # update_river_cell_storage_and_depth!
    @test Wflow.compute_river_storage_change(land, river, domain, 1, dt) ≈ 290.0
    @test Wflow.compute_external_inflow(river, land, 1, 1, dt) |> collect ≈
          [-8.0e-4, 8.0e-4]

    river_h_expected = 1.0014614814814815
    land_h_expected = 0.001461481481481508
    river_storage_expected = 24035.07555555556
    @test Wflow.compute_water_depths(1289.2, 1, 1, river, domain) |> collect ≈
          [river_h_expected, land_h_expected, river_storage_expected]

    Wflow.update_river_cell_storage_and_depth!(land, river, domain, 1, dt)
    @test river.variables.h[1] ≈ river_h_expected
    @test land.variables.h[1] ≈ land_h_expected
    @test river.variables.storage[1] ≈ river_storage_expected

    # update_land_cell_storage_and_depth!
    @test Wflow.compute_land_storage_change(land, domain.land.network, 2, dt) ≈ 200.0

    Wflow.update_land_cell_storage_and_depth!(land, domain.land, 2, dt)
    @test land.variables.storage[2] ≈ 1450.0
    @test land.variables.h[2] ≈ 0.002685185185185185
end

@testitem "local inertial long channel MacDonald (1997)" begin
    using Wflow: Unit, to_SI
    using QuadGK: quadgk
    using Graphs: DiGraph, add_edge!, ne
    using Statistics: mean
    using Wflow: GRAVITATIONAL_ACCELERATION

    CM = Unit(; cm = 1)

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

    variables =
        Wflow.LocalInertialRiverFlowVariables(; n_cells = n, n_edges = _ne, h = h_init)

    boundary_conditions =
        Wflow.RiverFlowBC(; n, external_inflow = zeros(n), reservoir = nothing)

    sw_river = Wflow.LocalInertialRiverFlow(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain = nothing,
        allocation = Wflow.NoAllocationRiver(n),
    )

    # run until steady state is reached
    epsilon = 1e-12
    (; flow_length) = domain_river.parameters
    while true
        sw_river.boundary_conditions.inwater[1] = 20.0
        h0 = mean(sw_river.variables.h)
        dt = Wflow.stable_timestep(sw_river, flow_length)
        Wflow.local_inertial_river_update!(sw_river, domain, dt, true)
        d = abs(h0 - mean(sw_river.variables.h))
        if d <= epsilon
            break
        end
    end

    # test for mean absolute error [m]
    @test mean(abs.(sw_river.variables.h .- h_a)) ≈ to_SI(1.873574206931199, CM)
end

@testitem "unit: local_inertial_flow" begin
    # Case of general area
    q0 = 0.0004713562869434079
    zs0 = 206.10117949049967
    zs1 = 201.9003737619653
    hf = 0.0011733869840497846
    A = 0.04970535373017763
    R = 0.0011733219820725962
    length = 533.453125
    mannings_n_sq = 0.0008999999597668652
    froude_limit = true
    dt = 89.29563868855615

    @test Wflow.local_inertial_flow(
        q0,
        zs0,
        zs1,
        hf,
        A,
        R,
        length,
        mannings_n_sq,
        froude_limit,
        dt,
    ) ≈ 0.005331926324969742

    # Case of rectangular area
    theta = 1.0
    q0 = 0.0001769756305800402
    qd = 0.0
    qu = 0.0
    zs0 = 601.4761297394623
    zs1 = 601.4730243288751
    hf = 0.00310727852479431
    width = 620.6649135473787
    length = 926.602742473319
    mannings_n_sq = 0.1773345894316103
    froude_limit = true
    dt = 49.774905820268735

    @test Wflow.local_inertial_flow(
        theta,
        q0,
        qd,
        qu,
        zs0,
        zs1,
        hf,
        width,
        length,
        mannings_n_sq,
        froude_limit,
        dt,
    ) ≈ 0.00017992597962222483
end

@testitem "unit: local_inertial_flow" begin
    # Case of general area
    q0 = 0.0004713562869434079
    zs0 = 206.10117949049967
    zs1 = 201.9003737619653
    hf = 0.0011733869840497846
    A = 0.04970535373017763
    R = 0.0011733219820725962
    length = 533.453125
    mannings_n_sq = 0.0008999999597668652
    froude_limit = true
    dt = 89.29563868855615

    @test Wflow.local_inertial_flow(
        q0,
        zs0,
        zs1,
        hf,
        A,
        R,
        length,
        mannings_n_sq,
        froude_limit,
        dt,
    ) ≈ 0.005331926324969742

    # Case of rectangular area
    theta = 1.0
    q0 = 0.0001769756305800402
    qd = 0.0
    qu = 0.0
    zs0 = 601.4761297394623
    zs1 = 601.4730243288751
    hf = 0.00310727852479431
    width = 620.6649135473787
    length = 926.602742473319
    mannings_n_sq = 0.1773345894316103
    froude_limit = true
    dt = 49.774905820268735

    @test Wflow.local_inertial_flow(
        theta,
        q0,
        qd,
        qu,
        zs0,
        zs1,
        hf,
        width,
        length,
        mannings_n_sq,
        froude_limit,
        dt,
    ) ≈ 0.00017992597962222483
end
