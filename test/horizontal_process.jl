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
