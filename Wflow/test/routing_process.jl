@testitem "unit: kinematic_wave" begin
    alpha = 2.586
    dt = 600.0
    dx = 1061.375

    # Case !(q_in + q_prev + q_lat ≈ 0.0)
    q_in = 1.104e-6
    q_prev = 0.0
    q_lat = 1.142e-6

    q, crossarea = Wflow.kinematic_wave(q_in, q_prev, q_lat, alpha, dt, dx)
    @test q ≈ 1.09308660753423e-6
    @test crossarea ≈ 0.0006852061693892164

    # Case q_in + q_prev + q_lat ≈ 0.0
    q_in = 0.0
    q_prev = 0.0
    q_lat = 0.0
    @test Wflow.kinematic_wave(q_in, q_prev, q_lat, alpha, dt, dx) == (0.0, 0.0)
end

@testitem "unit: ssf_celerity" begin
    zi = 0.3
    theta_e = 0.274
    slope = 0.00586
    i = 1

    # Case kh_profile::KhExponential
    kh_profile = Wflow.KhExponential([0.0002795374658372667], [1.8001038115471601])
    @test Wflow.ssf_celerity(zi, slope, theta_e, kh_profile, i) ≈ 3.4838105601686665e-6

    # Case kh_profile::KhExponentialConstant
    kh_profile = Wflow.KhExponentialConstant(kh_profile, [0.2])
    @test Wflow.ssf_celerity(zi, slope, theta_e, kh_profile, i) ≈ 4.170921791220723e-6
end

@testitem "unit: kw_ssf_newton_raphson" begin
    ssf = 0.008738344907407408
    constant_term = 77.774
    celerity = 0.0001418287037037037
    dt = 86400.0
    dx = 1103.816

    @test Wflow.kw_ssf_newton_raphson(ssf, constant_term, celerity, dt, dx) ≈
          0.01090947420564454
end

@testitem "kinematic wave overland flow" begin
    using NCDatasets: NCDataset
    using Graphs: topological_sort_by_dfs

    dt_sec = 86400.0
    ldd_MISSING_VALUE = 255

    # read the staticmaps into memory
    nc = NCDataset(normpath(@__DIR__, "data/input/staticmaps-rhine.nc"))
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
    beta = Wflow.BETA_KINWAVE
    AlpPow = (2.0 / 3.0) * beta
    AlpTermR = (N ./ sqrt.(slope)) .^ beta
    P = Bw + (2.0 * waterlevel)
    alpha = AlpTermR .* P .^ AlpPow

    Q = zeros(n)
    Q = Wflow.kin_wave!(Q, graph, toposort, Qold, q, alpha, DCL, dt_sec)

    @test sum(Q) ≈ 2.957806043289641e6
    @test Q[toposort[1]] ≈ 0.007260052312634069
    @test Q[toposort[n - 100]] ≈ 3945.762718338739
    @test Q[sink] ≈ 4131.101474418251
end

@testitem "unit: kinematic_wave_ssf" begin
    using StaticArrays: SVector
    include("testing_utils.jl")

    ### Shared values
    n = 1
    N = 4
    act_thickl = [SVector(0.1, 0.3, 0.8, 0.8)]
    sumlayers = [SVector(0.0, 0.1, 0.4, 1.2, 2.0)]
    maxlayers = 4
    nlayers = [4]
    n_unsatlayers = [3]
    theta_s = [0.48642662167549133]
    theta_r = [0.11939866840839386]
    theta_fc = [0.28219206182657536]
    ssfin = 0.0

    ssf_prev = 0.30038365579798126
    zi_prev = 0.0005198340870375973
    q_net = 0.005618827458801466
    slope = 0.4522336721420288
    sy = 0.20423455984891598
    d = 2.0
    dt = 86400.0
    dx = 1117.0150713112287
    dw = 517.495693771673
    ssfmax = 0.0009215296489248933
    kh_profile = Wflow.KhExponential([0.002379589787235966], [1.0141291422769427])
    i = 1

    soil_model = init_sbm_soil_model(
        n,
        N;
        # Variables
        ustorelayerthickness = [SVector(0.1, 0.3, 0.11983408703759733, NaN)],
        ustorelayerdepth = [
            SVector(0.0001909439890049523, 0.01627933934181815, 0.019508197676020186, 0.0),
        ],
        n_unsatlayers,
        zi = [0.5198340870375974],
        # Parameters
        maxlayers,
        sumlayers,
        nlayers,
        theta_s,
        theta_r,
        theta_fc,
        act_thickl,
    )

    ssf, zi, exfilt, net_flux = Wflow.kinematic_wave_ssf(
        ssfin,
        ssf_prev,
        zi_prev,
        q_net,
        slope,
        sy,
        d,
        dt,
        dx,
        dw,
        ssfmax,
        kh_profile,
        soil_model,
        i,
    )
    @test ssf ≈ 0.23130576097772237
    @test zi ≈ 0.1656875455413981
    @test exfilt ≈ 0.0
    @test net_flux ≈ -3.904277181728481e-7

    # Case: ssfin + ssf_prev ≈ 0.0 && r <= 0
    ssf_prev = 0.0
    r = 0.0
    zi_prev = 0.5198340870375974
    ssfmax = 0.0009215296489248933
    kh_profile = Wflow.KhExponential([0.002379589787235966], [1.0141291422769427])
    ssf, zi, exfilt, sy_d = Wflow.kinematic_wave_ssf(
        ssfin,
        ssf_prev,
        zi_prev,
        q_net,
        slope,
        sy,
        d,
        dt,
        dx,
        dw,
        ssfmax,
        kh_profile,
        soil_model,
        i,
    )
    @test iszero(ssf)
    @test zi == d
    @test iszero(exfilt)
    @test sy_d ≈ 0.0

    # Case: !(ssfin + ssf_prev ≈ 0.0 && r <= 0)
    # Case: !(zi > d)
    ssf_prev = 0.30038365579798126
    ssf, zi, exfilt, sy_d = Wflow.kinematic_wave_ssf(
        ssfin,
        ssf_prev,
        zi_prev,
        q_net,
        slope,
        sy,
        d,
        dt,
        dx,
        dw,
        ssfmax,
        kh_profile,
        soil_model,
        i,
    )
    @test ssf ≈ 0.255794305836017
    @test zi ≈ 0.7029236021516849
    @test exfilt ≈ 0.0
    @test net_flux ≈ -3.904277181728481e-7

    soil = init_sbm_soil_model(
        n,
        N;
        # Variables
        ustorelayerthickness = [SVector(0.1, 0.3, 0.348312461531486, NaN)],
        ustorelayerdepth = [
            SVector(0.0001909439890049523, 0.01627933934181815, 0.058425012193036086, 0.0),
        ],
        n_unsatlayers,
        zi = [0.7588905603985703],
        # Parameters
        maxlayers,
        sumlayers,
        nlayers,
        theta_s,
        theta_r,
        theta_fc,
        act_thickl,
    )

    ssf_prev = 0.627032986563781
    zi_prev = 0.7588905603985703
    q_net = 773.9150244657355
    slope = 0.4522336721420288
    sy = 0.20423455984891598
    d = 2.0
    dt = 1.0
    dx = 1117.0150713112287
    dw = 517.495693771673
    ssfmax = 153.46698446681825
    kh_profile = Wflow.KhExponentialConstant(kh_profile, [0.2])
    i = 1

    ssf, zi, exfilt, net_flux = Wflow.kinematic_wave_ssf(
        ssfin,
        ssf_prev,
        zi_prev,
        q_net,
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

    @test ssf ≈ 0.6300110241913047
    @test zi ≈ 0.7521780183868452
    @test exfilt ≈ 0.0
    @test net_flux ≈ 0.00133774649062672
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

@testitem "unit: kinwave_river_update!" begin
    using Graphs: DiGraph, add_edge!
    n = 2
    model = Wflow.KinWaveRiverFlowModel(;
        timestepping = Wflow.TimeStepping(),
        boundary_conditions = Wflow.RiverFlowBC(;
            n,
            external_inflow = [-0.1],
            reservoir = Wflow.ReservoirModel(;
                boundary_conditions = Wflow.ReservoirBC(; n, external_inflow = [0.02]),
                parameters = Wflow.ReservoirParameters(;
                    id = [1],
                    storfunc = [Wflow.ReservoirProfileType.linear],
                    outflowfunc = [Wflow.ReservoirOutflowType.simple],
                    area = [2500.0],
                ),
                variables = Wflow.ReservoirVariables(;
                    waterlevel = [3.0],
                    storage = [7500.0],
                ),
            ),
        ),
        parameters = Wflow.RiverFlowParameters(;
            flow = Wflow.ManningFlowParameters(;
                slope = [0.01],
                mannings_n = [0.03],
                alpha_pow = 0.4,
                alpha_term = [0.5],
                alpha = [5.0],
            ),
            bankfull_depth = [10.0],
        ),
        variables = Wflow.FlowVariables(; n, q = [0.2]),
        allocation = Wflow.NoAllocationRiverModel(n),
    )
    graph = DiGraph(2)
    add_edge!(graph, 1, 2)
    domain = Wflow.DomainRiver(;
        network = Wflow.NetworkRiver(;
            graph,
            order_of_subdomains = [[1]],
            order_subdomain = [[1]],
            subdomain_indices = [[1]],
            upstream_nodes = [[]],
            reservoir_indices = [1],
        ),
        parameters = Wflow.RiverParameters(; flow_width = [30.0], flow_length = [800.0]),
    )
    dt = 1200.0

    Wflow.kinwave_river_update!(model, domain, dt)

    @test model.variables.q ≈ [0.1598124775930105]
    @test model.variables.h[1] ≈ 0.055464507410878765
    @test model.variables.storage[1] ≈ 1331.1481778610903
    @test model.variables.q_cumulative[1] ≈ 191.7749731116126
end

@testitem "unit: local_inertial_river_update!" begin
    dt = 86400.0
    n = 2
    river_flow_model = Wflow.LocalInertialRiverFlowModel(;
        timestepping = Wflow.TimeStepping(),
        boundary_conditions = Wflow.RiverFlowBC(;
            n,
            external_inflow = [-1.0, -1.0],
            inwater = [100.0, 100.0],
            reservoir = Wflow.ReservoirModel(;
                boundary_conditions = Wflow.ReservoirBC(;
                    n,
                    external_inflow = [-1.0],
                    inflow_overland = [3000.0],
                    inflow_subsurface = [5000.0],
                    precipitation = [2.3148148148148148e-8],
                    evaporation = [1.1574074074074074e-8],
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
        floodplain = Wflow.FloodPlainModel(;
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
        allocation = Wflow.AllocationRiverModel(; n),
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

    Wflow.update_river_channel_flow!(river_flow_model, domain.river, dt)

    @test river_flow_model.variables.zs_src[1] ≈ 1.0
    @test river_flow_model.variables.zs_dst[1] ≈ 2.0
    @test river_flow_model.variables.zs_max[1] ≈ 2.0
    @test river_flow_model.variables.hf[1] ≈ 1.5
    @test river_flow_model.variables.a[1] ≈ 150.0
    @test river_flow_model.variables.r[1] ≈ 1.4563106796116505
    @test river_flow_model.variables.q[1] ≈ -575.3037784510024
    @test river_flow_model.variables.q_cumulative[1] ≈ -575303.7784510023

    Wflow.update_floodplain_flow!(river_flow_model, domain.river, dt)

    @test river_flow_model.floodplain.variables.hf[1] ≈ 1.0
    @test river_flow_model.floodplain.variables.hf_index[1] == 1
    @test river_flow_model.floodplain.variables.a[1] ≈ 90.0
    @test river_flow_model.floodplain.variables.r[1] ≈ 0.2356020942408377
    @test river_flow_model.floodplain.variables.q[1] ≈ -281.84014086002725
    @test river_flow_model.floodplain.variables.q_cumulative[1] ≈ -281840.1408600272

    Wflow.update_bc_reservoir_model!(river_flow_model, domain, dt)

    @test river_flow_model.boundary_conditions.reservoir.variables.storage[1] ≈
          7.391913765967477e6
    @test river_flow_model.boundary_conditions.reservoir.variables.waterlevel[1] ≈
          2.428382753193495
    @test river_flow_model.boundary_conditions.reservoir.variables.outflow[1] ≈
          0.00018509186397934759
    @test river_flow_model.boundary_conditions.reservoir.boundary_conditions.inflow_cumulative[1] ≈
          7.141856080688971e6
    @test river_flow_model.boundary_conditions.reservoir.variables.outflow_cumulative[1] ≈
          0.1850918639793476
    @test river_flow_model.boundary_conditions.reservoir.variables.actevap_cumulative[1] ≈
          1.1574074074074073e-5
    @test river_flow_model.variables.q[1] ≈ 0.00018509186397934759
    @test river_flow_model.variables.q_cumulative[1] ≈ -575303.5933591384

    Wflow.update_water_depth_and_storage!(river_flow_model, domain.river, dt)

    @test river_flow_model.variables.storage[2] ≈ 99000.0
    @test river_flow_model.variables.h[2] ≈ 99.0

    Wflow.update_water_depth_and_storage!(
        river_flow_model.floodplain,
        river_flow_model,
        domain.river,
        dt,
    )

    @test river_flow_model.variables.h ≈ [1.0, -79.7]
    @test river_flow_model.variables.storage ≈ [0.0, -79700.0]
    @test river_flow_model.floodplain.variables.storage ≈ [0.0, 178700.0]
end

@testitem "unit: update_directional_flow!" begin
    n = 3
    overland_flow_model = Wflow.LocalInertialOverlandFlowModel(;
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

    Wflow.update_directional_flow!(overland_flow_model, domain, i, dt, is_x_direction)
    @test overland_flow_model.variables.qx_cumulative[1] ≈ 26493.90166029366
end

@testitem "unit: local_inertial_update_water_depth!" begin
    n = 2
    overland_flow_model = Wflow.LocalInertialOverlandFlowModel(;
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
    river_flow_model = Wflow.LocalInertialRiverFlowModel(;
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
        allocation = Wflow.NoAllocationRiverModel(1),
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
    @test Wflow.compute_river_storage_change(
        overland_flow_model,
        river_flow_model,
        domain,
        1,
        dt,
    ) ≈ 290.0
    @test Wflow.compute_external_inflow(river_flow_model, overland_flow_model, 1, 1, dt) |>
          collect ≈ [-0.2, 0.2]

    river_h_expected = 1.0014614814814815
    land_h_expected = 0.001461481481481508
    river_storage_expected = 24035.07555555556
    @test Wflow.compute_water_depths(1289.2, 1, 1, river_flow_model, domain) |> collect ≈
          [river_h_expected, land_h_expected, river_storage_expected]

    Wflow.update_river_cell_storage_and_depth!(
        overland_flow_model,
        river_flow_model,
        domain,
        1,
        dt,
    )
    @test river_flow_model.variables.h[1] ≈ 1.0010925925925926
    @test overland_flow_model.variables.h[1] ≈ 0.001092592592592645
    @test river_flow_model.variables.storage[1] ≈ 24026.222222222223

    # update_land_cell_storage_and_depth!
    @test Wflow.compute_land_storage_change(
        overland_flow_model,
        domain.land.network,
        2,
        dt,
    ) ≈ 200.0

    Wflow.update_land_cell_storage_and_depth!(overland_flow_model, domain.land, 2, dt)
    @test overland_flow_model.variables.storage[2] ≈ 1450.0
    @test overland_flow_model.variables.h[2] ≈ 0.002685185185185185
end

@testitem "unit: kinwave_river_update!" begin
    # Test river kinematic wave routing on a 2-node graph (1 → 2).
    # Node 1 has a simple reservoir (outflowfunc = simple) and a negative external inflow
    # (i.e. abstraction). The test verifies discharge, water depth, storage and averaged
    # discharge for the river, as well as waterlevel, storage, outflow and actual
    # evaporation for the reservoir.
    using Graphs: DiGraph, add_edge!

    dt = 86400.0
    n = 2

    river_flow_model = Wflow.KinWaveRiverFlowModel(;
        timestepping = Wflow.TimeStepping(; stable_timesteps = zeros(n)),
        boundary_conditions = Wflow.RiverFlowBC(;
            n,
            external_inflow = zeros(n),
            reservoir = Wflow.ReservoirModel(;
                boundary_conditions = Wflow.ReservoirBC(;
                    n,
                    external_inflow = [-0.1],
                    precipitation = [2.0833332436504178e-10],
                    evaporation = [5.324074170655674e-9],
                    inflow_overland = [0.0],
                    inflow_subsurface = [0.02735223635554672],
                ),
                parameters = Wflow.ReservoirParameters(;
                    id = [1],
                    storfunc = [Wflow.ReservoirProfileType.linear],
                    outflowfunc = [Wflow.ReservoirOutflowType.simple],
                    area = [1.498462875e6],
                    maxrelease = [24.007999420166016],
                    demand = [3.000999927520752],
                    targetminfrac = [0.07482631504535675],
                    targetfullfrac = [0.7536525130271912],
                    maxstorage = [6.2e7],
                ),
                variables = Wflow.ReservoirVariables(;
                    waterlevel = [29.656373296565086],
                    storage = [4.443897439204416e7],
                    outflow_obs = [NaN],
                ),
            ),
        ),
        parameters = Wflow.RiverFlowParameters(;
            flow = Wflow.ManningFlowParameters(;
                slope = [0.017522206529974937, 0.01738094352185726],
                mannings_n = [0.03, 0.03],
                alpha_pow = 0.4,
                alpha_term = [0.41038937516728013, 0.4113871693777356],
                alpha = [2.544585458995107, 2.5507721996678145],
            ),
            bankfull_depth = [1.0, 1.0],
        ),
        variables = Wflow.FlowVariables(; n, q = [0.5499295110293246, 3.0005238507869465]),
        allocation = Wflow.NoAllocationRiverModel(n),
    )
    graph = DiGraph(2)
    add_edge!(graph, 1, 2)
    domain = Wflow.DomainRiver(;
        network = Wflow.NetworkRiver(;
            graph,
            order_of_subdomains = [[1]],
            order_subdomain = [[1, 2]],
            subdomain_indices = [[1, 2]],
            upstream_nodes = [[], [1]],
            reservoir_indices = [1, 0],
        ),
        parameters = Wflow.RiverParameters(;
            flow_width = [94.73094177246094, 94.73094177246094],
            flow_length = [1059.8125, 951.96875],
        ),
    )
    dt = Wflow.stable_timestep(river_flow_model, domain.parameters.flow_length, 0.05)

    Wflow.kinwave_river_update!(river_flow_model, domain, dt)

    @test dt ≈ 994.6119029285007
    @test river_flow_model.variables.q ≈ [0.37903337592185243, 3.1969698861305855]
    @test river_flow_model.variables.h ≈ [0.01500830539624011, 0.05407828963124342]
    @test river_flow_model.variables.storage ≈ [1506.7893805755937, 4876.828625285123]
    @test river_flow_model.variables.q_cumulative ≈ [376.9911072990474, 3179.744302049454]

    (; reservoir) = river_flow_model.boundary_conditions
    @test reservoir.variables.waterlevel[1] ≈ 29.654579645252387
    @test reservoir.variables.storage[1] ≈ 4.443628667214138e7
    @test reservoir.variables.outflow[1] ≈ 3.0009999145314317
    @test reservoir.variables.actevap_cumulative[1] ≈ 5.295387542208319e-6
end

@testitem "unit: local inertial river flow with one reservoir" begin
    # Test local inertial river routing (no floodplain) on a 3-node graph (1 → 2 → 3) and a
    # simple reservoir (outflowfunc = simple) at node 2. Each sub-step of the local inertial
    # update is called and verified individually:
    #   1. update_river_channel_flow!  — edge discharge, water surface elevations
    #   2. update_floodplain_flow!     — no floodplain
    #   3. update_bc_reservoir_model!  — reservoir storage, outflow, evaporation
    #   4. update_water_depth_and_storage!
    model_dt = 86400.0
    n = 3
    river_flow_model = Wflow.LocalInertialRiverFlowModel(;
        timestepping = Wflow.TimeStepping(),
        boundary_conditions = Wflow.RiverFlowBC(;
            n,
            external_inflow = zeros(n),
            inwater = [0.012816561479797707, 0.01544689027914484, -0.0004760637654763434],
            reservoir = Wflow.ReservoirModel(;
                boundary_conditions = Wflow.ReservoirBC(;
                    n = 1,
                    external_inflow = [0.0],
                    inflow_overland = [0.0],
                    inflow_subsurface = [0.04279912663469156],
                    precipitation = [2.0833332436504186e-10],
                    evaporation = [5.324074170655674e-9],
                ),
                parameters = Wflow.ReservoirParameters(;
                    id = [2],
                    storfunc = [Wflow.ReservoirProfileType.linear],
                    outflowfunc = [Wflow.ReservoirOutflowType.simple],
                    area = [1.498462875e6],
                    maxrelease = [24.007999420166016],
                    demand = [3.000999927520752],
                    targetminfrac = [0.07482631504535675],
                    targetfullfrac = [0.7536525130271912],
                    maxstorage = [6.2e7],
                ),
                variables = Wflow.ReservoirVariables(;
                    waterlevel = [29.6558203236325],
                    storage = [4.443814578263413e7],
                ),
            ),
        ),
        parameters = Wflow.LocalInertialRiverFlowParameters(;
            n,
            ne = 2,
            active_n = [1, 3],
            active_e = [1],
            froude_limit = true,
            h_thresh = 0.001,
            zb = [315.1000061035156, 314.3999938964844, 278.8000183105469],
            zb_max = [315.1000061035156, 314.3999938964844],
            bankfull_storage = [57155.32165002823, 100397.03622722626, 90180.89622545242],
            bankfull_depth = [1.0, 1.0, 1.0],
            mannings_n_sq = [0.0008999999597668652, 0.0008999999597668652],
            mannings_n = [0.029999999329447746, 0.029999999329447746, 0.029999999329447746],
            flow_length_at_edge = [831.578125, 1005.890625],
            flow_width_at_edge = [94.73094177246094, 94.73094177246094],
        ),
        variables = Wflow.LocalInertialRiverFlowVariables(;
            n_cells = n,
            n_edges = 2,
            h = [0.04484241735240722, 0.0, 0.07939389691400389],
            storage = [2562.9827873416416, 0.0, 7159.812778536053],
            # set q reservoir edge at zero, it is computed and not a reservoir state
            q = [0.534560186001325, 0.0],
        ),
        floodplain = nothing,
        allocation = Wflow.AllocationRiverModel(; n),
    )
    domain = Wflow.Domain(;
        river = Wflow.DomainRiver(;
            network = Wflow.NetworkRiver(;
                nodes_at_edge = Wflow.NodesAtEdge(; src = [1, 2], dst = [2, 3]),
                edges_at_node = Wflow.EdgesAtNode(;
                    src = [[], [1], [2]],
                    dst = [[1], [2], []],
                ),
            ),
            parameters = Wflow.RiverParameters(;
                flow_width = [94.73094177246094, 94.73094177246094, 94.73094177246094],
                flow_length = [603.34375, 1059.8125, 951.96875],
            ),
        ),
        reservoir = Wflow.DomainReservoir(;
            network = Wflow.NetworkReservoir(; river_indices = [2]),
        ),
    )
    dt = Wflow.stable_timestep(river_flow_model, domain.river.parameters.flow_length)

    Wflow.update_river_channel_flow!(river_flow_model, domain.river, dt)

    @test dt ≈ 909.829412320351
    @test river_flow_model.variables.zs_src ≈ [315.14484852086804, 0.0]
    @test river_flow_model.variables.zs_dst ≈ [314.3999938964844, 0.0]
    @test river_flow_model.variables.zs_max ≈ [315.14484852086804, 0.0]
    @test river_flow_model.variables.hf ≈ [0.04484241735241312, 0.0]
    @test river_flow_model.variables.a ≈ [4.247964427147839, 0.0]
    @test river_flow_model.variables.r ≈ [0.04480000374545946, 0.0]
    @test river_flow_model.variables.q ≈ [0.534558444239785, 0.0]
    @test river_flow_model.variables.q_cumulative ≈ [486.35699517356477, 0.0]

    Wflow.update_floodplain_flow!(river_flow_model, domain.river, dt)

    Wflow.update_bc_reservoir_model!(
        river_flow_model.boundary_conditions.reservoir,
        river_flow_model,
        domain,
        dt,
    )

    # check river q and q_av has been set by reservoir model
    @test river_flow_model.variables.q[2] ≈ 3.0009999145276134
    @test river_flow_model.variables.q[2] ==
          river_flow_model.boundary_conditions.reservoir.variables.outflow[1]
    @test river_flow_model.variables.q_cumulative[2] ≈ 2730.397988608082
    @test river_flow_model.variables.q_cumulative[2] ==
          river_flow_model.boundary_conditions.reservoir.variables.outflow_cumulative[1]

    @test river_flow_model.boundary_conditions.reservoir.variables.storage[1] ≈
          4.443593370702217e7
    @test river_flow_model.boundary_conditions.reservoir.variables.waterlevel[1] ≈
          29.654344093791327
    @test river_flow_model.boundary_conditions.reservoir.variables.outflow[1] ≈
          3.0009999145276134
    @test river_flow_model.boundary_conditions.reservoir.boundary_conditions.inflow_cumulative[1] ≈
          525.2968994074305
    @test river_flow_model.boundary_conditions.reservoir.variables.actevap_cumulative[1] ≈
          4.843999273837613e-6

    Wflow.update_water_depth_and_storage!(river_flow_model, domain.river, dt)

    Wflow.update_water_depth_and_storage!(
        river_flow_model.floodplain,
        river_flow_model,
        domain.river,
        dt,
    )

    @test river_flow_model.variables.storage ≈ [2088.286676767209, 0.0, 9889.777630328164]
    @test river_flow_model.variables.h ≈ [0.03653704705843743, 0.0, 0.10966599406601261]
end

@testitem "unit: local inertial river flow with a floodplain" begin
    # Test local inertial river routing with a floodplain on a 3-node graph (1 → 2 → 3).
    # Each sub-step of the local inertial update is called and verified individually:
    #   1. update_river_channel_flow!  — edge discharge, water surface elevations
    #   2. update_floodplain_flow!     — floodplain discharge and geometry
    #   3. update_bc_reservoir_model!  — no reservoir
    #   4. update_water_depth_and_storage!
    n = 3
    river_flow_model = Wflow.LocalInertialRiverFlowModel(;
        timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.7),
        boundary_conditions = Wflow.RiverFlowBC(;
            n,
            external_inflow = zeros(n),
            inwater = [0.001214603164946035, 0.0069799723162954465, 0.00022016439899281862],
            reservoir = nothing,
        ),
        parameters = Wflow.LocalInertialRiverFlowParameters(;
            n,
            ne = 2,
            active_n = [1, 2, 3],
            active_e = [1, 2],
            froude_limit = true,
            h_thresh = 0.001,
            zb = [165.10784077644348, 165.10784077644348, 165.10784077644348],
            zb_max = [165.10784077644348, 165.10784077644348],
            bankfull_storage = [107921.00118967971, 200292.17449130033, 69688.46493603183],
            bankfull_depth = [1.592156171798706, 1.592156171798706, 1.592156171798706],
            mannings_n_sq = [0.0008999999597668652, 0.0008999999597668652],
            mannings_n = [0.029999999329447746, 0.029999999329447746, 0.029999999329447746],
            flow_length_at_edge = [648.828125, 568.34375],
            flow_width_at_edge = [149.17837524414062, 149.17837524414062],
        ),
        variables = Wflow.LocalInertialRiverFlowVariables(;
            n_cells = n,
            n_edges = 2,
            h = [1.8817912224982847, 1.8197068314233162, 1.7619620034455687],
            storage = [127553.31189184493, 228917.89427333255, 77120.8437153619],
            q = [137.1750567107639, 133.74797838974442],
        ),
        floodplain = Wflow.FloodPlainModel(;
            parameters = Wflow.FloodPlainParameters(;
                profile = Wflow.FloodPlainProfile(;
                    depth = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5],
                    storage = [
                        0.0 0.0 0.0
                        77602.0 281141.0 111313.0
                        189506.0 609512.0 256357.0
                        346960.0 981178.0 420515.0
                        526346.0 1.39783e6 602100.0
                        747343.0 1.87239e6 814606.0
                    ],
                    width = [
                        149.178 149.178 149.178
                        341.577 666.779 758.764
                        492.562 778.794 988.691
                        693.057 881.476 1118.98
                        789.594 988.161 1237.77
                        972.752 1125.52 1448.54
                    ],
                    a = [
                        0.0 0.0 0.0
                        170.788 333.389 379.382
                        417.07 722.786 873.727
                        763.598 1163.52 1433.22
                        1158.4 1657.6 2052.1
                        1644.77 2220.36 2776.38
                    ],
                    p = [
                        192.399 517.6 609.585
                        193.399 518.6 610.585
                        345.384 631.615 841.512
                        546.879 735.297 972.803
                        644.416 842.983 1092.59
                        828.573 981.337 1304.37
                    ],
                ),
                mannings_n = [0.072, 0.072, 0.072],
                mannings_n_sq = [0.005184, 0.005184],
                zb_max = [166.6999969482422, 166.6999969482422],
            ),
            variables = Wflow.FloodPlainVariables(;
                n,
                n_edges = 2,
                q = [3.306660222819796, 6.14041420989245],
                h = [0.2896350506995787, 0.22755065962461016, 0.16980583164686275],
                storage = [25320.207706612186, 99321.92021301284, 30370.81429688439],
            ),
        ),
        allocation = Wflow.AllocationRiverModel(; n),
    )
    domain = Wflow.Domain(;
        river = Wflow.DomainRiver(;
            network = Wflow.NetworkRiver(;
                nodes_at_edge = Wflow.NodesAtEdge(; src = [1, 2], dst = [2, 3]),
                edges_at_node = Wflow.EdgesAtNode(;
                    src = [[], [1], [2]],
                    dst = [[1], [2], []],
                ),
            ),
            parameters = Wflow.RiverParameters(;
                flow_width = [149.17837524414062, 149.17837524414062, 149.17837524414062],
                flow_length = [454.375, 843.28125, 293.40625],
            ),
        ),
        reservoir = Wflow.DomainReservoir(; network = Wflow.NetworkReservoir()),
    )
    dt = Wflow.stable_timestep(river_flow_model, domain.river.parameters.flow_length)

    Wflow.update_river_channel_flow!(river_flow_model, domain.river, dt)

    @test dt ≈ 49.40931052556788
    @test river_flow_model.variables.zs_src ≈ [166.98963199894177, 166.92754760786679]
    @test river_flow_model.variables.zs_dst ≈ [166.92754760786679, 166.86980277988906]
    @test river_flow_model.variables.zs_max ≈ [166.98963199894177, 166.92754760786679]
    @test river_flow_model.variables.hf ≈ [1.8817912224982933, 1.8197068314233036]
    @test river_flow_model.variables.a ≈ [280.7225571209805, 271.4609085323917]
    @test river_flow_model.variables.r ≈ [1.8354842671202385, 1.7763698223484754]
    @test river_flow_model.variables.q ≈ [137.1827776559179, 133.7538757670657]
    @test river_flow_model.variables.q_cumulative ≈ [6778.106459961183, 6608.686781773178]

    Wflow.update_floodplain_flow!(river_flow_model, domain.river, dt)

    @test river_flow_model.floodplain.variables.hf ≈
          [0.2896350506995873, 0.22755065962459753]
    @test river_flow_model.floodplain.variables.hf_index == [1, 2]
    @test river_flow_model.floodplain.variables.a ≈ [55.72538543569419, 117.7803635852996]
    @test river_flow_model.floodplain.variables.r ≈
          [0.28876507912737354, 0.22735103521877678]
    @test river_flow_model.floodplain.variables.q ≈ [3.3074651032168534, 6.14213116671929]
    @test river_flow_model.floodplain.variables.q_cumulative ≈
          [163.41957033732095, 303.47846610520196]

    Wflow.update_bc_reservoir_model!(
        river_flow_model.boundary_conditions.reservoir,
        river_flow_model,
        domain,
        dt,
    )

    Wflow.update_water_depth_and_storage!(river_flow_model, domain.river, dt)

    @test river_flow_model.variables.storage ≈
          [120775.26544458869, 229087.6588271402, 83729.54137530623]
    @test river_flow_model.variables.h ≈
          [1.7817948514048583, 1.8210563184054411, 1.9129493838748897]
    Wflow.update_water_depth_and_storage!(
        river_flow_model.floodplain,
        river_flow_model,
        domain.river,
        dt,
    )

    @test river_flow_model.variables.storage ≈
          [124521.72968438723, 228924.52867849727, 78479.82315286293]
    @test river_flow_model.variables.h ≈
          [1.8370663564507992, 1.8197595694254676, 1.7930102910029269]
    @test river_flow_model.floodplain.variables.storage ≈
          [21410.323896476315, 99344.9914658879, 35924.01098543289]
    @test river_flow_model.floodplain.variables.h ≈
          [0.2449101846520932, 0.2276033976267614, 0.2008541192042207]
end

@testitem "unit: update_directional_flow!" begin
    # Test local inertial overland flow routing in x-direction at edge 2 (i = 2) using 3
    # nodes.
    n = 3
    overland_flow_model = Wflow.LocalInertialOverlandFlowModel(;
        timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.7),
        boundary_conditions = Wflow.LocalInertialOverlandFlowBC(;
            n,
            runoff = [0.0, 0.0, 0.003001456821567986],
        ),
        parameters = Wflow.LocalInertialOverlandFlowParameters(;
            n,
            ywidth = [926.6857061478484, 869.7426481339323, 812.7995901200163],
            xwidth = [],
            zx_max = [257.3280029296875, 232.67100524902344, 232.67100524902344],
            theta = 1.0,
            h_thresh = 1e-3,
            zy_max = [],
            mannings_n_sq = [0.24167056670421605, 0.2883451232664811, 0.3928782408368683],
            z = [257.3280029296875, 227.5050048828125, 232.67100524902344],
            froude_limit = true,
        ),
        variables = Wflow.LocalInertialOverlandFlowVariables(;
            n,
            qx0 = [0.0, -3.6332616217117395, -0.7525806207906618, 0.0],
            qx = [0.0, -3.63341089804407, -0.7526187151790501, 0.0],
            h = [0.0, 1.3754708010382453, 0.11735139800699446],
            storage = [0.0, 783157.9568615163, 237954.47204911432],
        ),
    )
    domain = Wflow.Domain(;
        land = Wflow.DomainLand(;
            network = Wflow.NetworkLand(;
                edge_indices = Wflow.EdgeConnectivity(;
                    n = 3,
                    xu = [2, 3, 4],
                    xd = [4, 1, 2],
                ),
            ),
            parameters = Wflow.LandParameters(;
                x_length = [614.4202561305977, 614.4202561305977, 614.4202561305977],
                y_length = [926.6857061478484, 926.6857061478484, 926.6857061478484],
                river_location = [0, 0, 1],
            ),
        ),
    )
    i = 2
    dt = Wflow.stable_timestep(overland_flow_model, domain.land.parameters)
    is_x_direction = true

    @test dt ≈ 117.10556654947368
    overland_flow_model.variables.qx0 .= overland_flow_model.variables.qx
    Wflow.update_directional_flow!(overland_flow_model, domain, i, dt, is_x_direction)
    @test overland_flow_model.variables.qx[i] ≈ -3.633493490896127
end

@testitem "unit: local_inertial_update_water_depth!" begin
    # Test update of river and land water depth and storage using 3 nodes for local inertial
    # overland and river flow (subgrid channel) routing, with 2 land cells and 1 river cell:
    # | land | land | river |.
    n_land = 3
    overland_flow_model = Wflow.LocalInertialOverlandFlowModel(;
        timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.7),
        variables = Wflow.LocalInertialOverlandFlowVariables(;
            n = n_land,
            qx = [0.0, -3.63341089804407, -0.7526187151790501, 0.0],
            qy = [0.0, -0.7369647824685662, 0.0, 0.0],
            storage = [0.0, 783157.9568615163, 237954.47204911432],
            h = [0.0, 1.3754708010382453, 0.11735139800699446],
        ),
        boundary_conditions = Wflow.LocalInertialOverlandFlowBC(;
            n = n_land,
            runoff = [0.0, 0.0, 0.003001456821567986],
        ),
        parameters = Wflow.LocalInertialOverlandFlowParameters(;
            n = n_land,
            xwidth = [],
            ywidth = [],
            theta = 1.0,
            h_thresh = 1e-3,
            zx_max = [],
            zy_max = [],
            mannings_n_sq = [],
            z = [],
            froude_limit = true,
        ),
    )
    n_river = 1
    river_flow_model = Wflow.LocalInertialRiverFlowModel(;
        timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.7),
        boundary_conditions = Wflow.RiverFlowBC(;
            n = n_river,
            reservoir = nothing,
            external_inflow = zeros(n_river),
        ),
        parameters = Wflow.LocalInertialRiverFlowParameters(;
            n = n_river,
            ne = 2,
            active_n = [1],
            active_e = [1, 2],
            froude_limit = true,
            h_thresh = 1e-3,
            zb = [],
            zb_max = [],
            bankfull_storage = [171137.5821314017],
            bankfull_depth = [1.3683528900146484],
            mannings_n_sq = [],
            mannings_n = [],
            flow_length_at_edge = [],
            flow_width_at_edge = [],
        ),
        variables = Wflow.LocalInertialRiverFlowVariables(;
            n_cells = n_river,
            n_edges = 2,
            q = [56.685647296907476, 53.70963118023338],
            q_average = [55.860830141394764, 53.06943588871848],
            q_channel_average = [55.860830141394764, 53.06943588871848],
            h = [1.485704288021643],
            storage = [185814.5230442402],
        ),
        floodplain = nothing,
        allocation = Wflow.NoAllocationRiverModel(n_river),
    )
    domain = Wflow.Domain(;
        land = Wflow.DomainLand(;
            network = Wflow.NetworkLand(;
                river_indices = [0, 0, 1],
                edge_indices = Wflow.EdgeConnectivity(;
                    n = 3,
                    xd = [4, 1, 2],
                    yd = [4, 4, 4],
                ),
            ),
            parameters = Wflow.LandParameters(;
                x_length = [614.4202561305977, 614.4202561305977, 614.4202561305977],
                y_length = [926.6857061478484, 926.6857061478484, 926.6857061478484],
            ),
        ),
        river = Wflow.DomainRiver(;
            network = Wflow.NetworkRiver(;
                edges_at_node = Wflow.EdgesAtNode(; src = [[1]], dst = [[2]]),
            ),
            parameters = Wflow.RiverParameters(;
                flow_width = [113.88611602783203],
                flow_length = [1098.1875],
            ),
        ),
    )

    dt = Wflow.stable_timestep(river_flow_model, domain.river.parameters.flow_length)
    @test dt ≈ 201.394687315008

    # Test functions called by function update_river_cell_storage_and_depth!
    storage_change = Wflow.compute_river_storage_change(
        overland_flow_model,
        river_flow_model,
        domain,
        3,
        dt,
    )
    @test storage_change ≈ 19.782071832453088
    @test Wflow.compute_external_inflow(river_flow_model, overland_flow_model, 3, 1, dt) |>
          collect ≈ [0.0, 0.0]

    river_h_expected = 1.4857390315391559
    land_h_expected = 0.11738614152450744
    river_storage_expected = 185818.86835722585
    total_storage = overland_flow_model.variables.storage[3] + storage_change
    @test Wflow.compute_water_depths(total_storage, 1, 3, river_flow_model, domain) |>
          collect ≈ [river_h_expected, land_h_expected, river_storage_expected]

    # Test update_river_cell_storage_and_depth!
    Wflow.update_river_cell_storage_and_depth!(
        overland_flow_model,
        river_flow_model,
        domain,
        3,
        dt,
    )
    @test river_flow_model.variables.h[1] ≈ river_h_expected
    @test overland_flow_model.variables.h[3] ≈ land_h_expected
    @test river_flow_model.variables.storage[1] ≈ river_storage_expected

    # Test function called by function update_land_cell_storage_and_depth!
    @test Wflow.compute_land_storage_change(
        overland_flow_model,
        domain.land.network,
        2,
        dt,
    ) ≈ 880.1704436259577

    # Test update_land_cell_storage_and_depth!
    Wflow.update_land_cell_storage_and_depth!(overland_flow_model, domain.land, 2, dt)
    @test overland_flow_model.variables.storage[2] ≈ 784038.1273051423
    @test overland_flow_model.variables.h[2] ≈ 1.3770166561681556
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

    h_thresh = 1.0e-3
    froude_limit = true
    h_init = zeros(n - 1)
    push!(h_init, h_a[n])

    timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.7)
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

    sw_river = Wflow.LocalInertialRiverFlowModel(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain = nothing,
        allocation = Wflow.NoAllocationRiverModel(n),
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
    @test mean(abs.(sw_river.variables.h .- h_a)) ≈ 0.01873574206931199
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
