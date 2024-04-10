const dt_sec = 86400.0
const ldd_mv = 255

# read the staticmaps into memory
nc = NCDataset(staticmaps_rhine_path)
# helper function to get the axis order and directionality right
read_right(nc, var) = reverse(permutedims(Array(nc[var])); dims = 2)
ldd_2d = read_right(nc, "ldd")

inds, _ = Wflow.active_indices(ldd_2d, ldd_mv)
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
graph = Wflow.flowgraph(ldd, inds, Wflow.pcr_dir)
# a topological sort is used for visiting nodes in order from upstream to downstream
toposort = topological_sort_by_dfs(graph)
sink = toposort[end]
@test ldd[sink] == 5  # the most downstream node must be a sink

# calculate parameters of kinematic wave
const q = 0.000001
const beta = 0.6
const AlpPow = (2.0 / 3.0) * beta
AlpTermR = (N ./ sqrt.(slope)) .^ beta
P = Bw + (2.0 * waterlevel)
alpha = AlpTermR .* P .^ AlpPow

Q = zeros(n)
Q = Wflow.kin_wave!(Q, graph, toposort, Qold, q, alpha, beta, DCL, dt_sec)

@testset "flow rate" begin
    @test sum(Q) ≈ 2.957806043289641e6
    @test Q[toposort[1]] ≈ 0.007260052312634069f0
    @test Q[toposort[n-100]] ≈ 3945.762718338739f0
    @test Q[sink] ≈ 4131.101474418251
end

@testset "kinematic wave subsurface flow" begin
    @test all(
        isapprox.(
            Wflow.kinematic_wave_ssf(
                210128378079.0733,
                215395179156.82645,
                1540.34273559,
                1.238,
                18021.0,
                0.25,
                0.346,
                0.0017669756,
                1800.0,
                1.0,
                1697.05 * 1000.0,
                1200.0 * 1000.0,
                2443723.716252628,
                1.0,
                "exponential",
            ),
            (7.410313985168225e10, 1540.1496836278836, -0.0),
        ),
    )
end

@testset "accucapacity" begin
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
    flux, new_material = Wflow.accucapacityflux(material, network, capacity)
    @test new_material != material
    @test new_material == [0.0, 0.5, 0.5, 0.0, 3.5, 1.5]
    @test flux == Float64[0.5, 1.5, 1.5, 1, 1.5, 1.5]

    # example 2, accucapacityflux
    material = fill(10.0, 6)
    capacity = Float64[2, 30, 30, 2, 30, 2]
    flux, new_material = Wflow.accucapacityflux(material, network, capacity)
    @test new_material != material
    @test new_material == [8.0, 0.0, 0.0, 10.0, 0.0, 40.0]
    @test flux == Float64[2, 10, 10, 2, 30, 2]

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

@testset "local inertial long channel MacDonald (1997)" begin

    g = 9.80665
    L = 1000.0
    dx = 5.0
    n = Int(L / dx)

    # analytical solution MacDonald (1997) for channel with length L of 1000.0 m, Manning's
    # n of 0.03, constant inflow of 20.0 m3/s at upper boundary and channel width of 10.0 m
    # water depth profile h(x)
    h(x) = (4 / g)^(1.0 / 3.0) * (1.0 + 0.5 * exp(-16.0 * (x / L - 0.5)^2.0))
    # spatial derivative of h(x)
    h_acc(x) =
        -(4 / g)^(1.0 / 3.0) * 16.0 / L * (x / L - 0.5) * exp(-16 * (x / L - 0.5)^2.0)
    # solution for channel slope s(x)
    s(x) =
        (1.0 - 4.0 / (g * h(x)^(3.0))) * h_acc(x) +
        0.36 * (2 * h(x) + 10.0)^(4.0 / 3.0) / ((10.0 * h(x))^(10.0 / 3.0))

    h_a = h.([dx:dx:L;]) # water depth profile (analytical solution)
    # integrate slope to get elevation (bed level) z
    x = [dx:dx:L;]
    zb = first.([quadgk(s, xi, L, rtol = 1e-12) for xi in x])

    # initialize ShallowWaterRiver
    graph = DiGraph(n)
    for i = 1:n
        add_edge!(graph, i, i + 1)
    end

    dl = fill(dx, n)
    width = fill(10.0, n)
    n_river = fill(0.03, n)

    # for each link the src and dst node is required
    nodes_at_link = Wflow.adjacent_nodes_at_link(graph)
    _ne = ne(graph)

    # determine z, width, length and manning's n at links
    zb_max = fill(0.0, _ne)
    width_at_link = fill(0.0, _ne)
    length_at_link = fill(0.0, _ne)
    mannings_n_sq = fill(0.0, _ne)
    for i = 1:_ne
        zb_max[i] = max(zb[nodes_at_link.src[i]], zb[nodes_at_link.dst[i]])
        width_at_link[i] = min(width[nodes_at_link.dst[i]], width[nodes_at_link.src[i]])
        length_at_link[i] = 0.5 * (dl[nodes_at_link.dst[i]] + dl[nodes_at_link.src[i]])
        mannings_n =
            (
                n_river[nodes_at_link.dst[i]] * dl[nodes_at_link.dst[i]] +
                n_river[nodes_at_link.src[i]] * dl[nodes_at_link.src[i]]
            ) / (dl[nodes_at_link.dst[i]] + dl[nodes_at_link.src[i]])
        mannings_n_sq[i] = mannings_n * mannings_n
    end


    network = (
        nodes_at_link = nodes_at_link,
        links_at_node = Wflow.adjacent_links_at_node(graph, nodes_at_link),
    )

    alpha = 0.7
    dt = 1.0
    h_thresh = 1.0e-03
    froude_limit = true
    h_init = zeros(n - 1)
    push!(h_init, h_a[n])

    sw_river = Wflow.ShallowWaterRiver(
        n = n,
        ne = _ne,
        active_n = collect(1:n-1),
        active_e = collect(1:_ne),
        g = 9.80665,
        alpha = alpha,
        h_thresh = h_thresh,
        dt = dt,
        q0 = zeros(_ne),
        q = zeros(_ne),
        q_av = zeros(_ne),
        q_channel_av = zeros(_ne),
        zb_max = zb_max,
        mannings_n_sq = mannings_n_sq,
        mannings_n = n_river,
        h = h_init,
        eta_max = zeros(_ne),
        eta_src = zeros(_ne),
        eta_dst = zeros(_ne),
        hf = zeros(_ne),
        h_av = zeros(n),
        width = width,
        width_at_link = width_at_link,
        a = zeros(_ne),
        r = zeros(_ne),
        volume = fill(0.0, n),
        error = zeros(n),
        inflow = zeros(n),
        inflow_wb = zeros(n),
        inwater = zeros(n),
        dl = dl,
        dl_at_link = length_at_link,
        bankfull_volume = fill(Wflow.mv, n),
        bankfull_depth = fill(Wflow.mv, n),
        zb = zb,
        froude_limit = froude_limit,
        reservoir_index = Int[],
        lake_index = Int[],
        waterbody = zeros(n),
        reservoir = nothing,
        lake = nothing,
        floodplain = nothing,
    )

    # run until steady state is reached
    epsilon = 1.0e-12
    while true
        sw_river.inwater[1] = 20.0
        h0 = mean(sw_river.h)
        dt = Wflow.stable_timestep(sw_river)
        Wflow.shallowwater_river_update(sw_river, network, dt, 0.0, true)
        d = abs(h0 - mean(sw_river.h))
        if d <= epsilon
            break
        end
    end

    # test for mean absolute error [cm]
    @test mean(abs.(sw_river.h .- h_a)) * 100.0 ≈ 1.873574206931199

end
