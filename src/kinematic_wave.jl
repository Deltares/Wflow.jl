
"Map from PCRaster LDD value to a CartesianIndex"
const pcrdir = [
    CartesianIndex(1, -1),  # 1
    CartesianIndex(1, 0),  # 2
    CartesianIndex(1, 1),  # 3
    CartesianIndex(0, -1),  # 4
    CartesianIndex(0, 0),  # 5
    CartesianIndex(0, 1),  # 6
    CartesianIndex(-1, -1),  # 7
    CartesianIndex(-1, 0),  # 8
    CartesianIndex(-1, 1),  # 9
]

# in case the input data is permuted the other way
permute_indices(inds) = [CartesianIndex(i[2], i[1]) for i in inds]

"Get a list of indices that are active, based on a nodata value"
function active_indices(A, nodata)
    inds = CartesianIndices(size(A))
    filter(i -> !isequal(A[i], nodata), inds)
end

"Convert a gridded drainage direction to a directed graph"
function flowgraph(ldd::AbstractVector, inds::AbstractVector, pcrdir::AbstractVector)
    # prepare a directed graph to be filled
    n = length(inds)
    dag = DiGraph(n)

    # loop over ldd, adding the edge to the downstream node
    for (from_node, from_index) in enumerate(inds)
        ldd_val = ldd[from_node]
        # skip pits to prevent cycles
        ldd_val == 5 && continue
        to_index = from_index + pcrdir[ldd_val]
        # find the node id of the downstream cell
        to_node = searchsortedfirst(inds, to_index)
        add_edge!(dag, from_node, to_node)
    end
    @assert is_directed(dag)
    @assert !is_cyclic(dag)
    return dag
end

"Kinematic wave flow rate for a single cell and timestep"
function kinematic_wave(Qin, Qold, q, α, β, Δt, Δx)
    ϵ = 1.0e-12
    max_iters = 3000

    if Qin + Qold + q ≈ 0.0
        return 0.0
    else
        # common terms
        ab_pQ = α * β * pow(((Qold + Qin) / 2.0), (β - 1.0))
        Δtx = Δt / Δx
        C = Δtx * Qin + α * pow(Qold, β) + Δt * q

        Qkx = (Δtx * Qin + Qold * ab_pQ + Δt * q) / (Δtx + ab_pQ)
        if isnan(Qkx)
            Qkx = 0.0
        end
        Qkx = max(Qkx, 1.0e-30)
        count = 1

        while true
            fQkx = Δtx * Qkx + α * pow(Qkx, β) - C
            dfQkx = Δtx + α * β * pow(Qkx, (β - 1.0))
            Qkx = Qkx - fQkx / dfQkx
            Qkx = max(Qkx, 1.0e-30)
            if (abs(fQkx) <= ϵ) || (count >= max_iters)
                break
            end
            count += 1
        end

        return Qkx
    end
end

"Kinematic wave flow rate over the whole network for a single timestep"
function kin_wave!(Q, dag, toposort, Qold, q, α, β, DCL, Δt)
    for v in toposort
        upstream_nodes = inneighbors(dag, v)
        Qin = isempty(upstream_nodes) ? 0.0 : sum(Q[i] for i in upstream_nodes)
        Q[v] = kinematic_wave(Qin, Qold[v], q, α[v], β, Δt, DCL[v])
    end
    return Q
end
