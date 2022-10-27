function broadcast_to_ranks(obj, comm::MPI.Comm)
    return MPI.bcast(obj, root, comm)
end

broadcast_to_ranks(obj, comm::Nothing) = obj

function scatter_to_ranks(
    data,
    indices,
    nprocs::Int,
    comm::MPI.Comm,
    rank::Int;
    default_size = 1,
)

    reqs = Vector{MPI.Request}(undef, 0)

    size = length(indices[rank+1])
    remote_size = size == 0 ? default_size : length(indices[rank+1])
    A_local = zeros(eltype(data), remote_size)
    remote_buf = MPI.Buffer(A_local)

    if rank == root
        for sendrank = 0:nprocs-1
            # Get the indices on the root buffer to send to the remote buffer
            inds = indices[sendrank+1]
            data_on_root = data[inds]
            root_buf = MPI.Buffer(data_on_root)
            sendtag = sendrank + 1000
            sreq = MPI.Isend(root_buf, sendrank, sendtag, comm)
            push!(reqs, sreq)
        end
    end

    recievetag = rank + 1000
    rreq = MPI.Irecv!(remote_buf, root, recievetag, comm)
    push!(reqs, rreq)

    MPI.Waitall!(reqs)
    return A_local
end

scatter_to_ranks(A, indices, nprocs::Int, comm::Nothing, rank::Int) = A

function subgraphs!(subgraphs, sub_toposort, ranks, graph, local_indices)
    for r in ranks
        g = induced_subgraph(graph, local_indices[r])
        subgraphs[r] = g
        sub_toposort[r] = topological_sort_by_dfs(g[1])
    end
end

function set_sinks!(sinks, subgraphs, sub_toposort, ranks)
    for r in ranks
        sinks[r] = subgraphs[r][2][sub_toposort[r][end]]
    end
end

function connectivity(toposort, graph, graph_subbas, subgraphs, sinks)

    n = length(toposort)
    dst_nodes = fill(-1, n)
    neighbors = fill(-1, n)

    for v in toposort
        dst_node = outneighbors(graph, sinks[v])
        idx_ds = outneighbors(graph_subbas, v) # downstream rank

        if !isempty(idx_ds)
            idx_ds = only(idx_ds)
            neighbors[v] = idx_ds - 1 # MPI rank nr.
            g = subgraphs[idx_ds][1]
            d = findfirst(x -> x == only(dst_node), subgraphs[idx_ds][2])
            if !isnothing(d)
                add_vertex!(g)
                add_edge!(g, nv(g), d)
                dst_nodes[v] = nv(g)
            end
        end
    end
    return dst_nodes, neighbors
end

function set_subdomains(graph, subdomains, local_indices, toposort, nprocs, index_pit)

    sinks = Vector{Int}(undef, nprocs)
    subgraphs = Vector{Tuple{Graphs.DiGraph,Vector{Int}}}(undef, nprocs)
    sub_toposort = Vector{Vector{Int}}(undef, nprocs)
    ghost_node = fill(-1, nprocs)
    neighbor_rank = fill(-1, nprocs)
    order_subbas = Vector{Vector{Int}}()
    index = Vector{Int}()

    # assign basin ids
    n_pits = length(index_pit)
    basin = fill(0, length(toposort))
    basin[index_pit] = [1:n_pits;]
    basin_fill = fillnodata_upstream(graph, toposort, basin, 0)

    for i = 1:n_pits
        inds_basin = findall(x -> x == i, basin_fill)
        subdomains_basin = subdomains[inds_basin]
        ranks_basin = unique(subdomains_basin)
        g, _ = induced_subgraph(graph, inds_basin)
        subgraphs!(subgraphs, sub_toposort, ranks_basin, g, local_indices)
        set_sinks!(sinks, subgraphs, sub_toposort, ranks_basin)

        subbas = fill(0, length(inds_basin))
        subbas[sinks[ranks_basin]] = ranks_basin
        offset = minimum(ranks_basin) - 1
        graph_subbas = graph_from_nodes(g, subbas, subdomains_basin; offset = offset)
        toposort_subbas = topological_sort_by_dfs(graph_subbas)
        dist = Graphs.Experimental.Traversals.distances(
            Graph(graph_subbas),
            toposort_subbas[end],
        )
        max_dist = maximum([dist; 1])
        v_subbas = subbasins_order(graph_subbas, toposort_subbas[end], max_dist)
        for i in eachindex(v_subbas)
            v_subbas[i] .+= offset
        end
        append!(order_subbas, v_subbas)
        append!(index, 1:length(v_subbas))

        ghost_node[ranks_basin], neighbor_rank[ranks_basin] =
            connectivity(toposort_subbas, graph, graph_subbas, subgraphs, sinks)
    end
    subbas_order = Vector{Vector{Int}}(undef, maximum(index))
    for m = 1:maximum(index)
        subbas_order[m] = reduce(vcat, order_subbas[index.==m])
    end
    return subgraphs, sub_toposort, ghost_node, neighbor_rank, subbas_order
end
