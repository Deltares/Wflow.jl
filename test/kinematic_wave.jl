const Δt = 86400.0
const ldd_mv = 255

# read the staticmaps into memory
nc = NCDataset(staticmaps_rhine_path)
ldd_2d = nc["ldd"][:]
slope_2d = nc["slope"][:]
N_2d = nc["N"][:]
Qold_2d = nc["Qold"][:]
Bw_2d = nc["Bw"][:]
waterlevel_2d = nc["waterlevel"][:]
DCL_2d = nc["DCL"][:]
close(nc)

inds = Wflow.active_indices(ldd_2d, ldd_mv)
n = length(inds)

# take out only the active cells
ldd = ldd_2d[inds]
slope = slope_2d[inds]
N = N_2d[inds]
Qold = Qold_2d[inds]
Bw = Bw_2d[inds]
waterlevel = waterlevel_2d[inds]
DCL = DCL_2d[inds]

# create the directed acyclic graph from the drainage direction array
dag = Wflow.flowgraph(ldd, inds, Wflow.pcrdir)
# a topological sort is used for visiting nodes in order from upstream to downstream
toposort = topological_sort_by_dfs(dag)
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
Q = Wflow.kin_wave!(Q, dag, toposort, Qold, q, α, β, DCL, Δt)

@testset "flow rate" begin
    @test sum(Q) ≈ 2.957806043289641e6
    @test Q[toposort[1]] ≈ 0.00021453683608501235
    @test Q[toposort[n-100]] ≈ 4054.9507466731234
    @test Q[sink] ≈ 4131.101474418251
end
