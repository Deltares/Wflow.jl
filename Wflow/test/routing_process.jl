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
    water_table_depth = 0.3
    theta_e = 0.274
    slope = 0.00586
    i = 1

    # Case kh_profile::KhExponential
    kh_profile = Wflow.KhExponential([0.0002795374658372667], [1.8001038115471601])
    @test Wflow.ssf_celerity(water_table_depth, slope, theta_e, kh_profile, i) ≈
          3.4838105601686665e-6

    # Case kh_profile::KhExponentialConstant
    kh_profile = Wflow.KhExponentialConstant(kh_profile, [0.2])
    @test Wflow.ssf_celerity(water_table_depth, slope, theta_e, kh_profile, i) ≈
          4.170921791220723e-6
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
    graph, ldd = Wflow.flowgraph(ldd, inds, Wflow.PCR_DIR)
    # a topological sort is used for visiting nodes in order from upstream to downstream
    toposort = topological_sort_by_dfs(graph)
    sink = toposort[end]
    @test ldd[sink] == Wflow.LDD_PIT  # the most downstream node must be a sink

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
    actual_layer_thickness = [SVector(0.1, 0.3, 0.8, 0.8)]
    cumulative_layer_depth = [SVector(0.0, 0.1, 0.4, 1.2, 2.0)]
    maximum_number_of_layers = 4
    number_of_layers = [4]
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
        unsaturated_layer_thickness = [SVector(0.1, 0.3, 0.11983408703759733, NaN)],
        unsaturated_layer_depth = [
            SVector(0.0001909439890049523, 0.01627933934181815, 0.019508197676020186, 0.0),
        ],
        n_unsatlayers,
        water_table_depth = [0.5198340870375974],
        # Parameters
        maximum_number_of_layers,
        cumulative_layer_depth,
        number_of_layers,
        theta_s,
        theta_r,
        theta_fc,
        actual_layer_thickness,
    )

    ssf, water_table_depth, exfilt, net_flux = Wflow.kinematic_wave_ssf(
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
    @test water_table_depth ≈ 0.1656875455413981
    @test exfilt ≈ 0.0
    @test net_flux ≈ -3.904277181728481e-7

    # Case: ssfin + ssf_prev ≈ 0.0 && r <= 0
    ssf_prev = 0.0
    hydraulic_radius = 0.0
    zi_prev = 0.5198340870375974
    ssfmax = 0.0009215296489248933
    kh_profile = Wflow.KhExponential([0.002379589787235966], [1.0141291422769427])
    ssf, water_table_depth, exfilt, sy_d = Wflow.kinematic_wave_ssf(
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
    @test water_table_depth == d
    @test iszero(exfilt)
    @test sy_d ≈ 0.0

    # Case: !(ssfin + ssf_prev ≈ 0.0 && r <= 0)
    # Case: !(zi > d)
    ssf_prev = 0.30038365579798126
    ssf, water_table_depth, exfilt, sy_d = Wflow.kinematic_wave_ssf(
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
    @test water_table_depth ≈ 0.7029236021516849
    @test exfilt ≈ 0.0
    @test net_flux ≈ -3.904277181728481e-7

    soil = init_sbm_soil_model(
        n,
        N;
        # Variables
        unsaturated_layer_thickness = [SVector(0.1, 0.3, 0.348312461531486, NaN)],
        unsaturated_layer_depth = [
            SVector(0.0001909439890049523, 0.01627933934181815, 0.058425012193036086, 0.0),
        ],
        n_unsatlayers,
        water_table_depth = [0.7588905603985703],
        # Parameters
        maximum_number_of_layers,
        cumulative_layer_depth,
        number_of_layers,
        theta_s,
        theta_r,
        theta_fc,
        actual_layer_thickness,
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

    ssf, water_table_depth, exfilt, net_flux = Wflow.kinematic_wave_ssf(
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
    @test water_table_depth ≈ 0.7521780183868452
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
    # Test river kinematic wave routing on a 2-node graph (1 → 2).
    # Node 1 has a simple reservoir (outflow_curve_type = simple) and a negative external inflow
    # (i.e. abstraction). The test verifies discharge, water depth, storage and averaged
    # discharge for the river, as well as waterlevel, storage, outflow and actual
    # evaporation for the reservoir.
    using Graphs: DiGraph, add_edge!

    n = 2
    bankfull_depth = [1.0, 1.0]
    flow_width = [94.73094177246094, 94.73094177246094]
    flow_length = [1059.8125, 951.96875]
    bankfull_storage = @. flow_length * flow_width * bankfull_depth

    river_flow_model = Wflow.RiverFlowModel(;
        routing_method = Wflow.KinematicWave(),
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
                    storage_curve_type = [Wflow.ReservoirProfileType.linear],
                    outflow_curve_type = [Wflow.ReservoirOutflowType.simple],
                    area = [1.498462875e6],
                    maximum_release = [24.007999420166016],
                    demand = [3.000999927520752],
                    target_minimum_fraction = [0.07482631504535675],
                    target_full_fraction = [0.7536525130271912],
                    maximum_storage = [6.2e7],
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
            bankfull_depth,
            bankfull_storage,
        ),
        variables = Wflow.RiverFlowVariables(;
            n,
            q = [0.5499295110293246, 3.0005238507869465],
        ),
        allocation = Wflow.NoAllocationRiverModel(n),
        floodplain = nothing,
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
        parameters = Wflow.RiverParameters(; flow_width, flow_length),
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
    # simple reservoir (outflow_curve_type = simple) at node 2. Each sub-step of the local inertial
    # update is called and verified individually:
    #   1. update_river_channel_flow!  — edge discharge, water surface elevations
    #   2. update_floodplain_flow!     — no floodplain
    #   3. update_bc_reservoir_model!  — reservoir storage, outflow, evaporation
    #   4. update_water_depth_and_storage!
    model_dt = 86400.0
    n = 3
    river_flow_model = Wflow.RiverFlowModel(;
        routing_method = Wflow.LocalInertial(),
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
                    storage_curve_type = [Wflow.ReservoirProfileType.linear],
                    outflow_curve_type = [Wflow.ReservoirOutflowType.simple],
                    area = [1.498462875e6],
                    maximum_release = [24.007999420166016],
                    demand = [3.000999927520752],
                    target_minimum_fraction = [0.07482631504535675],
                    target_full_fraction = [0.7536525130271912],
                    maximum_storage = [6.2e7],
                ),
                variables = Wflow.ReservoirVariables(;
                    waterlevel = [29.6558203236325],
                    storage = [4.443814578263413e7],
                ),
            ),
        ),
        parameters = Wflow.RiverFlowStaggeredParameters(;
            n,
            n_edges = 2,
            active_n = [1, 3],
            active_e = [1],
            froude_limit = true,
            h_thresh = 0.001,
            zb = [315.1000061035156, 314.3999938964844, 278.8000183105469],
            zb_at_edge = [315.1000061035156, 314.3999938964844],
            bankfull_storage = [57155.32165002823, 100397.03622722626, 90180.89622545242],
            bankfull_depth = [1.0, 1.0, 1.0],
            mannings_n_sq_at_edge = [0.0008999999597668652, 0.0008999999597668652],
            mannings_n_at_edge = [
                0.029999999329447746,
                0.029999999329447746,
                0.029999999329447746,
            ],
            flow_length_at_edge = [831.578125, 1005.890625],
            flow_width_at_edge = [94.73094177246094, 94.73094177246094],
        ),
        variables = Wflow.RiverFlowStaggeredVariables(;
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
    @test river_flow_model.variables.zs_at_edge ≈ [315.14484852086804, 0.0]
    @test river_flow_model.variables.water_depth_at_edge ≈ [0.04484241735241312, 0.0]
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

@testitem "unit: local inertial river flow with floodplain" begin
    # Test local inertial river routing with a floodplain on a 3-node graph (1 → 2 → 3).
    # Each sub-step of the local inertial update is called and verified individually:
    #   1. update_river_channel_flow!  — edge discharge, water surface elevations
    #   2. update_floodplain_flow!     — floodplain discharge and geometry
    #   3. update_bc_reservoir_model!  — no reservoir
    #   4. update_water_depth_and_storage!
    n = 3
    river_flow_model = Wflow.RiverFlowModel(;
        routing_method = Wflow.LocalInertial(),
        timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.7),
        boundary_conditions = Wflow.RiverFlowBC(;
            n,
            external_inflow = zeros(n),
            inwater = [0.001214603164946035, 0.0069799723162954465, 0.00022016439899281862],
            reservoir = nothing,
        ),
        parameters = Wflow.RiverFlowStaggeredParameters(;
            n,
            n_edges = 2,
            active_n = [1, 2, 3],
            active_e = [1, 2],
            froude_limit = true,
            h_thresh = 0.001,
            zb = [165.10784077644348, 165.10784077644348, 165.10784077644348],
            zb_at_edge = [165.10784077644348, 165.10784077644348],
            bankfull_storage = [107921.00118967971, 200292.17449130033, 69688.46493603183],
            bankfull_depth = [1.592156171798706, 1.592156171798706, 1.592156171798706],
            mannings_n = [0.029999999329447746, 0.029999999329447746, 0.029999999329447746],
            mannings_n_sq_at_edge = [0.0008999999597668652, 0.0008999999597668652],
            mannings_n_at_edge = [0.029999999329447746, 0.029999999329447746],
            flow_length_at_edge = [648.828125, 568.34375],
            flow_width_at_edge = [149.17837524414062, 149.17837524414062],
        ),
        variables = Wflow.RiverFlowStaggeredVariables(;
            n_cells = n,
            n_edges = 2,
            h = [1.8817912224982847, 1.8197068314233162, 1.7619620034455687],
            storage = [127553.31189184493, 228917.89427333255, 77120.8437153619],
            q = [137.1750567107639, 133.74797838974442],
        ),
        floodplain = Wflow.FloodPlainModel(;
            routing_method = Wflow.LocalInertial(),
            parameters = Wflow.FloodPlainStaggeredParameters(;
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
                        149.17837524414062 149.17837524414062 149.17837524414062;
                        341.5768913342503 666.778728923476 758.7636596016615;
                        492.562310866575 778.7935519733185 988.6905953775695;
                        693.0574965612104 881.4757828423199 1118.9809351368624;
                        789.5944979367263 988.1614230127849 1237.7718606880392;
                        972.7515818431912 1125.5153603853994 1448.5444669293854;
                    ],
                    flow_area = [
                        0.0 0.0 0.0;
                        170.78844566712516 333.389364461738 379.38182980083076;
                        417.0696011004127 722.7861404483972 873.7271274896154;
                        763.5983493810179 1163.5240318695571 1433.2175950580468;
                        1158.395598349381 1657.6047433759495 2052.103525402066;
                        1644.7713892709767 2220.3624235686493 2776.375758866759;
                    ],
                    wetted_perimeter = [
                        192.3985160901097 517.6003536793354 609.5852843575209;
                        193.3985160901097 518.6003536793354 610.5852843575209;
                        345.38393562243436 631.6151767291778 841.5122201334289;
                        546.8791213170698 735.2974075981792 972.8025598927218;
                        644.4161226925856 842.9830477686443 1092.5934854438985;
                        828.5732065990505 981.3369851412588 1304.3660916852448;
                    ],
                ),
                mannings_n = [0.072, 0.072, 0.072],
                mannings_n_at_edge = [0.072, 0.072],
                mannings_n_sq_at_edge = [0.005184, 0.005184],
                zb_at_edge = [166.6999969482422, 166.6999969482422],
            ),
            variables = Wflow.FloodPlainStaggeredVariables(;
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
    @test river_flow_model.variables.zs_at_edge ≈ [166.98963199894177, 166.92754760786679]
    @test river_flow_model.variables.water_depth_at_edge ≈
          [1.8817912224982933, 1.8197068314233036]
    @test river_flow_model.variables.q ≈ [137.1827776559179, 133.7538757670657]
    @test river_flow_model.variables.q_cumulative ≈ [6778.106459961183, 6608.686781773178]

    Wflow.update_floodplain_flow!(river_flow_model, domain.river, dt)

    @test river_flow_model.floodplain.variables.water_depth_at_edge ≈
          [0.2896350506995873, 0.22755065962459753]
    @test river_flow_model.floodplain.variables.q ≈ [3.3074672215578524, 6.1421232455587536]

    @test river_flow_model.floodplain.variables.q_cumulative ≈
          [163.41967500308917, 303.4780747261213]

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
          [124521.73491986065, 228924.54042985922, 78479.82701991481]
    @test river_flow_model.variables.h ≈
          [1.8370664336896243, 1.8197596628390198, 1.793010379352563]
    @test river_flow_model.floodplain.variables.storage ≈
          [21410.318556337137, 99344.98021057082, 35924.006727001935]
    @test river_flow_model.floodplain.variables.h ≈
          [0.2449102618909183, 0.22760349104031377, 0.20085420755385677]
end

@testitem "unit: river flow on a staggered grid using Manning's equation with floodplain" begin
    # Test river flow on a staggered grid using Manning's equation with a floodplain on a
    # 3-node graph (1 → 2 → 3).
    # Each sub-step of the update of Manning's flow on staggered grid  is called and
    # verified individually:
    #   1. update_river_channel_flow!  — edge discharge, water surface elevations
    #   2. update_floodplain_flow!     — floodplain discharge and geometry
    #   3. update_bc_reservoir_model!  — no reservoir
    #   4. update_water_depth_and_storage!
    n = 3
    river_flow_model = Wflow.RiverFlowModel(;
        routing_method = Wflow.ManningStaggered(),
        timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.7),
        boundary_conditions = Wflow.RiverFlowBC(;
            n,
            external_inflow = zeros(n),
            inwater = [0.5469041025637776, 1.2789538388011605, 0.1723191462506401],
            reservoir = nothing,
        ),
        parameters = Wflow.RiverFlowStaggeredParameters(;
            n,
            n_edges = 2,
            active_n = [1, 2, 3],
            active_e = [1, 2],
            h_thresh = 0.001,
            zb = [415.6000061035156, 415.6000061035156, 415.6000061035156],
            zb_at_edge = [415.6000061035156, 415.6000061035156],
            bankfull_storage = [36208.125, 30283.125, 25717.5],
            bankfull_depth = [1.0, 1.0, 1.0],
            mannings_n = [0.029999999329447746, 0.029999999329447746, 0.029999999329447746],
            mannings_n_at_edge = [0.029999999329447746, 0.029999999329447746],
            flow_length_at_edge = [1108.1875, 933.34375],
            flow_width_at_edge = [30.0, 30.0],
            slope_at_edge = [1.0e-5, 1.0e-5],
        ),
        variables = Wflow.RiverFlowStaggeredVariables(;
            n_cells = n,
            n_edges = 2,
            h = [1.0029472866230495, 2.432593354336019, 0.40220945912422745],
            storage = [36314.840722458204, 73666.52862352695, 10343.82176502732],
            q = [12.586919865724491, 12.586919865724491],
            water_depth_at_edge = [2.4325367293470777, 2.4325367293470777],
        ),
        floodplain = Wflow.FloodPlainModel(;
            routing_method = Wflow.ManningStaggered(),
            parameters = Wflow.FloodPlainStaggeredParameters(;
                profile = Wflow.FloodPlainProfile(;
                    depth = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5],
                    storage = [
                        0.0 0.0 0.0
                        37399.0 36242.0 34514.0;
                        79402.0 73635.0 69028.0;
                        150174.0 111029.0 103542.0;
                        233029.0 155325.0 138056.0;
                        341777.0 212852.0 172570.0;
                    ],
                    width = [
                        30.0 30.0 30.0;
                        61.97338304593237 71.8063277815615 80.52260134149898;
                        69.60260991144943 74.08680577054052 80.52260134149898;
                        117.27533530112372 74.08878707200793 80.52260134149898;
                        137.2979131065196 87.76372980001238 80.52260134149898;
                        180.20485733519754 113.9783295152003 80.52260134149898;
                    ],
                    flow_area = [
                        0.0 0.0 0.0;
                        30.986691522966186 35.90316389078075 40.26130067074949;
                        65.7879964786909 72.94656677605101 80.52260134149898;
                        124.42566412925277 109.99096031205498 120.78390201224848;
                        193.07462068251255 153.87282521206117 161.04520268299797;
                        283.1770493501113 210.86198996966132 201.30650335374747;
                    ],
                    wetted_perimeter = [
                        31.973383045932373 41.80632778156151 50.522601341498984;
                        32.97338304593237 42.80632778156151 51.522601341498984;
                        41.602609911449434 46.08680577054052 52.522601341498984;
                        90.27533530112372 47.08878707200793 53.522601341498984;
                        111.29791310651959 61.763729800012385 54.522601341498984;
                        155.20485733519754 88.9783295152003 55.522601341498984;
                    ],
                ),
                mannings_n = [0.072, 0.072, 0.072],
                mannings_n_at_edge = [0.072, 0.072],
                zb_at_edge = [416.6000061035156, 416.6000061035156],
                slope_at_edge = [1.0e-5, 1.0e-5],
            ),
            variables = Wflow.FloodPlainStaggeredVariables(;
                n,
                n_edges = 2,
                q = [2.0759155671111853, 3.279029180032925],
                water_depth_at_edge = [1.4325367293470777, 1.4325367293470777],
                h = [0.0029472866230494782, 1.432593354336019, 0.0],
                storage = [113.73542237265065, 62604.388160555245, 0.0],
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
                flow_width = [30.0, 30.0, 30.0],
                flow_length = [1206.9375, 1009.4375, 857.25],
            ),
        ),
        reservoir = Wflow.DomainReservoir(; network = Wflow.NetworkReservoir()),
    )
    dt = Wflow.stable_timestep(river_flow_model, domain.river.parameters.flow_length)

    Wflow.update_river_channel_flow!(river_flow_model, domain.river, dt)

    @test dt ≈ 2272.7523105527266
    @test river_flow_model.variables.zs_src ≈ [416.6029533901387, 418.03259945785163]
    @test river_flow_model.variables.zs_dst ≈ [418.03259945785163, 416.00221556263983]
    @test river_flow_model.variables.zs_at_edge ≈ [416.6029533901387, 418.03259945785163]
    @test river_flow_model.variables.water_depth_at_edge ≈
          [1.0029472866230495, 2.432593354336019]
    @test river_flow_model.variables.q ≈ [3.0436243193280137, 12.587380945649572]
    @test river_flow_model.variables.q_cumulative ≈ [6917.404204207213, 28607.999128032432]

    Wflow.update_floodplain_flow!(river_flow_model, domain.river, dt)

    @test river_flow_model.floodplain.variables.water_depth_at_edge ≈
          [0.0029472866230494782, 1.432593354336019]
    @test river_flow_model.floodplain.variables.q ≈ [8.50693948416018e-5, 3.279243909752732]

    @test river_flow_model.floodplain.variables.q_cumulative ≈
          [0.1933416636835727, 7452.909172756478]

    Wflow.update_bc_reservoir_model!(
        river_flow_model.boundary_conditions.reservoir,
        river_flow_model,
        domain,
        dt,
    )

    Wflow.update_water_depth_and_storage!(river_flow_model, domain.river, dt)

    @test river_flow_model.variables.storage ≈
          [30640.414081003582, 54882.67899192734, 39343.459630853366]
    @test river_flow_model.variables.h ≈
          [0.8462303441838974, 1.8123188736937599, 1.529832201063609]

    Wflow.update_water_depth_and_storage!(
        river_flow_model.floodplain,
        river_flow_model,
        domain.river,
        dt,
    )

    @test river_flow_model.variables.storage ≈
          [30753.95616171255, 63042.829748341144, 33570.77415623857]
    @test river_flow_model.variables.h ≈
          [0.849366161923948, 2.0817808514920815, 1.3053669352090433]
    @test river_flow_model.floodplain.variables.storage ≈
          [0.0, 46991.52157304865, 13225.594647371276]
    @test river_flow_model.floodplain.variables.h ≈
          [0.0, 1.0817808514920815, 0.3053669352090434]
end

@testitem "unit: kinematic river flow including 1D floodplain schematization" begin
    # Test river kinematic wave routing with floodplain on a 2-node graph (1 → 2).
    # Each sub-step of the update of kinematic wave routing with floodplain is called and
    # verified individually:
    #   1. river_channel_floodplain_exchange!
    #   2. kinwave_river_update!
    #   3. update_floodplain_model!
    using Graphs: DiGraph, add_edge!

    n = 2
    flow_length = [750.953125, 851.8125]
    flow_width = [229.91920471191406, 229.91920471191406]

    river_flow_model = Wflow.RiverFlowModel(;
        routing_method = Wflow.KinematicWave(),
        timestepping = Wflow.TimeStepping(; stable_timesteps = zeros(n)),
        boundary_conditions = Wflow.RiverFlowBC(;
            n,
            inwater = [0.02953203136486251, 0.0001849569866329713],
            reservoir = nothing,
        ),
        parameters = Wflow.RiverFlowParameters(;
            flow = Wflow.ManningFlowParameters(;
                slope = [1.0e-5, 1.0e-5],
                mannings_n = [0.03, 0.03],
                alpha_pow = 0.4,
                alpha_term = [3.8572052304976667, 3.8572052304976667],
                alpha = [34.0789466790827, 34.0789466790827],
            ),
            bankfull_depth = [2.1051321029663086, 2.1051321029663086],
            bankfull_storage = [363469.04651181493, 412286.02275520907],
        ),
        variables = Wflow.RiverFlowVariables(;
            n,
            q = [296.52948301601174, 192.8313119108856],
            qlat = [3.9326064945614956e-5, 2.1713344971219756e-7],
            h = [4.509741437854894, 3.4835130995322534],
            storage = [778645.3962305915, 682239.2566234164],
        ),
        allocation = Wflow.NoAllocationRiverModel(n),
        floodplain = Wflow.FloodPlainModel(;
            routing_method = Wflow.Manning(),
            parameters = Wflow.FloodPlainParameters(;
                profile = Wflow.FloodPlainProfile(;
                    depth = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5],

                    storage = [
                        0.0 0.0;
                        86329.2726379633 97924.02628183365;
                        172659.2726379633 207762.02628183365;
                        258989.2726379633 386716.02628183365;
                        369518.2726379633 605009.0262818336;
                        724648.2726379633 869843.0262818336;
                    ],
                    width = [
                        229.91920471191406 229.91920471191406;
                        229.91920471191406 229.91920471191406;
                        229.9211418821914 257.89243524836746;
                        229.9211418821914 420.1722796977034;
                        294.369904912507 512.5376770122533;
                        945.8113647239966 621.8128989654414;
                    ],
                    flow_area = [
                        0.0 0.0;
                        114.95960235595703 114.95960235595703;
                        229.9201732970527 243.90581998014076;
                        344.8807442381484 453.99195982899244;
                        492.06569669440194 710.260798335119;
                        964.9713790564002 1021.1672478178398;
                    ],
                    wetted_perimeter = [
                        0.0 0.0;
                        1.0 1.0;
                        2.00193717027733 29.9732305364534;
                        3.00193717027733 193.25307498578934;
                        68.45070020059296 286.61847230033925;
                        720.8921600120825 396.8936942535273;
                    ],
                ),
                mannings_n = [0.072, 0.072],
                slope = [1.0e-5, 1.0e-5],
            ),
            variables = Wflow.FloodPlainVariables(;
                n,
                q = [1.2861909826521447, 1.9846650910027395],
                h = [0.0, 0.0],
                storage = [0.0, 0.0],
            ),
        ),
    )
    graph = DiGraph(2)
    add_edge!(graph, 1, 2)
    domain = Wflow.DomainRiver(;
        network = Wflow.NetworkRiver(;
            graph,
            order = [1, 2],
            order_of_subdomains = [[1]],
            order_subdomain = [[1, 2]],
            subdomain_indices = [[1, 2]],
            upstream_nodes = [[], [1]],
            reservoir_indices = [0, 0],
        ),
        parameters = Wflow.RiverParameters(; flow_width, flow_length),
    )
    dt = Wflow.stable_timestep(river_flow_model, domain.parameters.flow_length, 0.05)

    Wflow.river_channel_floodplain_exchange!(river_flow_model, domain.parameters, dt)

    @test dt ≈ 1602.881460805217
    @test river_flow_model.variables.h ≈ [4.509741437854894, 3.4835130995322534]
    @test river_flow_model.variables.storage ≈ [778645.3962305915, 682239.2566234164]
    @test river_flow_model.floodplain.variables.h ≈ [2.0642836103410205, 1.1737631111525133]
    @test river_flow_model.floodplain.variables.storage ≈
          [58760.1445203583, 40074.01437791623]

    Wflow.kinwave_river_update!(river_flow_model, domain, dt)

    @test river_flow_model.variables.h ≈ [2.872002930358695, 3.1138609375883397]
    @test river_flow_model.variables.storage ≈ [495875.84798393055, 609843.6005807515]
    @test river_flow_model.variables.q ≈ [139.7837242183345, 159.9486203252721]
    @test river_flow_model.variables.q_cumulative ≈ [224056.7400718776, 256378.67820075116]

    Wflow.update_floodplain_model!(river_flow_model, domain, dt)

    @test river_flow_model.floodplain.variables.storage ≈
          [52745.30941770246, 41649.967280840756]
    @test river_flow_model.floodplain.variables.q ≈ [3.752513987924126, 2.7693140810933157]
    @test river_flow_model.floodplain.variables.q_cumulative ≈
          [6014.835102655834, 4438.882199731312]
end

@testitem "unit: update_directional_flow!" begin
    # Test local inertial overland flow routing in x-direction at edge 2 (i = 2) using 3
    # nodes.
    n = 3
    overland_flow_model = Wflow.OverlandFlowModel(;
        routing_method = Wflow.LocalInertial(),
        timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.7),
        boundary_conditions = Wflow.LocalInertialOverlandFlowBC(;
            n,
            runoff = [0.0, 0.0, 0.003001456821567986],
        ),
        parameters = Wflow.LocalInertialOverlandFlowParameters(;
            n,
            ywidth_at_edge = [926.6857061478484, 869.7426481339323, 812.7995901200163],
            xwidth_at_edge = [],
            zx_max_at_edge = [257.3280029296875, 232.67100524902344, 232.67100524902344],
            theta = 1.0,
            h_thresh = 1e-3,
            zy_max_at_edge = [],
            mannings_n_sq_at_edge = [
                0.24167056670421605,
                0.2883451232664811,
                0.3928782408368683,
            ],
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
                    ind_x_up = [2, 3, 4],
                    ind_x_down = [4, 1, 2],
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
    overland_flow_model = Wflow.OverlandFlowModel(;
        routing_method = Wflow.LocalInertial(),
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
            xwidth_at_edge = [],
            ywidth_at_edge = [],
            theta = 1.0,
            h_thresh = 1e-3,
            zx_max_at_edge = [],
            zy_max_at_edge = [],
            mannings_n_sq_at_edge = [],
            z = [],
            froude_limit = true,
        ),
    )
    n_river = 1
    river_flow_model = Wflow.RiverFlowModel(;
        routing_method = Wflow.LocalInertial(),
        timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.7),
        boundary_conditions = Wflow.RiverFlowBC(; n = n_river, reservoir = nothing),
        parameters = Wflow.RiverFlowStaggeredParameters(;
            n = n_river,
            n_edges = 2,
            active_n = [1],
            active_e = [1, 2],
            froude_limit = true,
            h_thresh = 1e-3,
            zb = [],
            zb_at_edge = [],
            bankfull_storage = [171137.5821314017],
            bankfull_depth = [1.3683528900146484],
            mannings_n_sq_at_edge = [],
            mannings_n_at_edge = [],
            flow_length_at_edge = [],
            flow_width_at_edge = [],
        ),
        variables = Wflow.RiverFlowStaggeredVariables(;
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
                    ind_x_down = [4, 1, 2],
                    ind_y_down = [4, 4, 4],
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

    # Test functions called by function update_river_and_land_storage_and_depth!
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

    # Test update_river_and_land_storage_and_depth!
    Wflow.update_river_and_land_storage_and_depth!(
        overland_flow_model,
        river_flow_model,
        domain,
        3,
        dt,
    )
    @test river_flow_model.variables.h[1] ≈ river_h_expected
    @test overland_flow_model.variables.h[3] ≈ land_h_expected
    @test river_flow_model.variables.storage[1] ≈ river_storage_expected

    # Test function called by function update_land_storage_and_depth!
    @test Wflow.compute_land_storage_change(
        overland_flow_model,
        domain.land.network,
        2,
        dt,
    ) ≈ 880.1704436259577

    # Test update_land_storage_and_depth!
    Wflow.update_land_storage_and_depth!(overland_flow_model, domain.land, 2, dt)
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
    n_edges = ne(graph)

    # determine z, width, length and manning's n at edges
    zb_at_edge = fill(0.0, n_edges)
    flow_width_at_edge = fill(0.0, n_edges)
    flow_length_at_edge = fill(0.0, n_edges)
    mannings_n_sq_at_edge = fill(0.0, n_edges)
    for i in 1:n_edges
        zb_at_edge[i] = max(zb[nodes_at_edge.src[i]], zb[nodes_at_edge.dst[i]])
        flow_width_at_edge[i] =
            min(width[nodes_at_edge.dst[i]], width[nodes_at_edge.src[i]])
        flow_length_at_edge[i] = 0.5 * (dl[nodes_at_edge.dst[i]] + dl[nodes_at_edge.src[i]])
        mannings_n_at_edge =
            (
                n_river[nodes_at_edge.dst[i]] * dl[nodes_at_edge.dst[i]] +
                n_river[nodes_at_edge.src[i]] * dl[nodes_at_edge.src[i]]
            ) / (dl[nodes_at_edge.dst[i]] + dl[nodes_at_edge.src[i]])
        mannings_n_sq_at_edge[i] = mannings_n_at_edge * mannings_n_at_edge
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
    parameters = Wflow.RiverFlowStaggeredParameters(;
        n,
        n_edges,
        active_n = collect(1:(n - 1)),
        active_e = collect(1:n_edges),
        h_thresh,
        zb_at_edge,
        mannings_n_sq_at_edge,
        flow_width_at_edge,
        flow_length_at_edge,
        bankfull_storage = fill(Wflow.MISSING_VALUE, n),
        bankfull_depth = fill(Wflow.MISSING_VALUE, n),
        zb,
        froude_limit,
    )

    variables = Wflow.RiverFlowStaggeredVariables(; n_cells = n, n_edges, h = h_init)

    boundary_conditions =
        Wflow.RiverFlowBC(; n, external_inflow = zeros(n), reservoir = nothing)

    sw_river = Wflow.RiverFlowModel(;
        routing_method = Wflow.LocalInertial(),
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
        Wflow.staggered_scheme_river_update!(sw_river, domain, dt, true)
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
    q_previous = 0.0004713562869434079
    zs0 = 206.10117949049967
    zs1 = 201.9003737619653
    water_depth_at_edge = 0.0011733869840497846
    A = 0.04970535373017763
    R = 0.0011733219820725962
    length = 533.453125
    mannings_n_sq_at_edge = 0.0008999999597668652
    froude_limit = true
    dt = 89.29563868855615

    @test Wflow.local_inertial_flow(
        q_previous,
        zs0,
        zs1,
        water_depth_at_edge,
        A,
        R,
        length,
        mannings_n_sq_at_edge,
        froude_limit,
        dt,
    ) ≈ 0.005331926324969742

    # Case of rectangular area
    theta = 1.0
    q_previous = 0.0001769756305800402
    qd = 0.0
    qu = 0.0
    zs0 = 601.4761297394623
    zs1 = 601.4730243288751
    water_depth_at_edge = 0.00310727852479431
    width = 620.6649135473787
    length = 926.602742473319
    mannings_n_sq_at_edge = 0.1773345894316103
    froude_limit = true
    dt = 49.774905820268735

    @test Wflow.local_inertial_flow(
        theta,
        q_previous,
        qd,
        qu,
        zs0,
        zs1,
        water_depth_at_edge,
        width,
        length,
        mannings_n_sq_at_edge,
        froude_limit,
        dt,
    ) ≈ 0.00017992597962222483
end

@testitem "unit: local_inertial_flow" begin
    # Case of general area
    q_previous = 0.0004713562869434079
    zs0 = 206.10117949049967
    zs1 = 201.9003737619653
    water_depth_at_edge = 0.0011733869840497846
    A = 0.04970535373017763
    R = 0.0011733219820725962
    length = 533.453125
    mannings_n_sq_at_edge = 0.0008999999597668652
    froude_limit = true
    dt = 89.29563868855615

    @test Wflow.local_inertial_flow(
        q_previous,
        zs0,
        zs1,
        water_depth_at_edge,
        A,
        R,
        length,
        mannings_n_sq_at_edge,
        froude_limit,
        dt,
    ) ≈ 0.005331926324969742

    # Case of rectangular area
    theta = 1.0
    q_previous = 0.0001769756305800402
    qd = 0.0
    qu = 0.0
    zs0 = 601.4761297394623
    zs1 = 601.4730243288751
    water_depth_at_edge = 0.00310727852479431
    width = 620.6649135473787
    length = 926.602742473319
    mannings_n_sq_at_edge = 0.1773345894316103
    froude_limit = true
    dt = 49.774905820268735

    @test Wflow.local_inertial_flow(
        theta,
        q_previous,
        qd,
        qu,
        zs0,
        zs1,
        water_depth_at_edge,
        width,
        length,
        mannings_n_sq_at_edge,
        froude_limit,
        dt,
    ) ≈ 0.00017992597962222483
end
