const Δt_sec = 86400.0
const ldd_mv = 255

# read the staticmaps into memory
nc = NCDataset(staticmaps_rhine_path)
ldd_2d = Wflow.ncread(nc, "ldd")

inds, _ = Wflow.active_indices(ldd_2d, ldd_mv)
n = length(inds)

# take out only the active cells
ldd = ldd_2d[inds]
slope = Wflow.ncread(nc, "slope"; sel = inds)
N = Wflow.ncread(nc, "N"; sel = inds)
Qold = Wflow.ncread(nc, "Qold"; sel = inds)
Bw = Wflow.ncread(nc, "Bw"; sel = inds)
waterlevel = Wflow.ncread(nc, "waterlevel"; sel = inds)
DCL = Wflow.ncread(nc, "DCL"; sel = inds)
close(nc)

# create the directed acyclic graph from the drainage direction array
graph = Wflow.flowgraph(ldd, inds, Wflow.pcrdir)
# a topological sort is used for visiting nodes in order from upstream to downstream
toposort = topological_sort_by_dfs(graph)
sink = toposort[end]
@test ldd[sink] == 5  # the most downstream node must be a sink

# calculate parameters of kinematic wave
const q = 0.000001
const β = 0.6
const AlpPow = (2.0 / 3.0) * β
AlpTermR = (N ./ sqrt.(slope)) .^ β
P = Bw + (2.0 * waterlevel)
α = AlpTermR .* P .^ AlpPow

Q = zeros(n)
Q = Wflow.kin_wave!(Q, graph, toposort, Qold, q, α, β, DCL, Δt_sec)

@testset "flow rate" begin
    @test sum(Q) ≈ 2.957806043289641e6
    @test Q[toposort[1]] ≈ 0.00021453683608501235
    @test Q[toposort[n-100]] ≈ 4054.9507466731234
    @test Q[sink] ≈ 4131.101474418251
end

@testset "kinematic wave subsurface flow" begin
    @test all(isapprox.(
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
        ),
        (7.410313985168225e10, 1540.1496836278836, -0.0),
    ))
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
