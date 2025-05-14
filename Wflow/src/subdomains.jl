"""
    fillnodata_upstream(g, toposort, data, nodata)

Fill `nodata` upstream cells with the value from the first downstream valid cell, based on
directed acyclic graph `g`, topological order `toposort`, nodata value `nodata` and `data`
containing missing values. Returns filled data `data_out`.
"""
function fillnodata_upstream(g, toposort, data, nodata)
    data_out = copy(data)
    for v in reverse(toposort)  # down- to upstream
        idx_ds = outneighbors(g, v)
        if !isempty(idx_ds)
            if data_out[v] == nodata && data_out[only(idx_ds)] != nodata
                data_out[v] = data_out[only(idx_ds)]
            end
        end
    end
    return data_out
end

"""
    stream_order(g, toposort)

Return the Strahler streamorder based on directed acyclic graph `g` and topological order
`toposort`.
"""
function stream_order(g, toposort)
    n = length(toposort)
    strord = fill(1, n)
    for v in toposort
        inds_up = inneighbors(g, v)
        if !isempty(inds_up)
            sto_up = strord[inds_up]
            if length(findall(x -> x == maximum(sto_up), sto_up)) > 1
                strord[v] = (maximum(sto_up) + 1)
            else
                strord[v] = maximum(sto_up)
            end
        end
    end
    return strord
end

"""
    subbasins(g, streamorder, toposort, min_sto)

Return subbasins with a unique id starting at 1, based on a minimum streamorder `min_sto`,
directed acyclic graph `g` and topological order `toposort`.
"""
function subbasins(g, streamorder, toposort, min_sto)
    n = length(toposort)
    subbas = fill(Int(0), n)

    i = 1
    for v in toposort
        if streamorder[v] < min_sto
            continue
        end
        ds_nodes = outneighbors(g, v)
        if !isempty(ds_nodes)
            if streamorder[v] != streamorder[only(ds_nodes)]
                subbas[v] = i
                i += 1
            end
        else
            # also set pits (without a downstream node)
            subbas[v] = i
            i += 1
        end
    end
    return subbas
end

"""
    subbasins_order(g, outlet, max_dist)

Group subbasins (down- to upstream), starting at the `outlet` of the basin, per distance (1
to maximum distance `max_dist`) from the `outlet`, based on the directed acyclic graph `g`
of subbasins. Returns grouped subbasins `order`, ordered from `max_dist` (including
subbasins without upstream neighbor and distance < `max_dist`) to distance 0 (`outlet`).
"""
function subbasins_order(g, outlet, max_dist)
    order = Vector{Vector{Int}}(undef, max_dist + 1)
    order[1] = [outlet]
    for i in 1:max_dist
        v = Vector{Int}()
        for n in order[i]
            ups_nodes = inneighbors(g, n)
            if !isempty(ups_nodes)
                append!(v, ups_nodes)
            end
        end
        order[i + 1] = v
    end

    # move subbasins without upstream neighbor (headwater) to index [max_dist+1]
    for i in 1:max_dist
        for s in order[i]
            if isempty(inneighbors(g, s))
                append!(order[max_dist + 1], s)
                filter!(e -> e ≠ s, order[i])
            end
        end
    end
    return reverse(order)
end

"""
    graph_from_nodes(graph, subbas, subbas_fill)

Extract directed acyclic graph `g` representing the flow network at subbasin level from
`subbas` containing nodes with a unique id representing subbasin outlets, `subbas_fill` with
subbasin ids for the complete domain, and directed acyclic graph `graph` representing the
flow network for each subbasin cell.
"""
function graph_from_nodes(graph, subbas, subbas_fill)
    n = maximum(subbas)
    g = DiGraph(n)
    for i in 1:n
        idx = findall(x -> x == i, subbas)
        ds_idx = outneighbors(graph, only(idx))
        to_node = subbas_fill[ds_idx]
        if !isempty(to_node)
            add_edge!(g, i, only(to_node))
        end
    end
    return g
end

"""
    function kinwave_set_subdomains(graph, toposort, index_pit, streamorder, min_sto)

Setup subdomains for parallel execution (threading) of the kinematic wave calculation.
Subdomains are subbasins based on a minimum stream order `min_sto` (see also `subbasins(g,
streamorder, toposort, min_sto)`). Subbasins are extracted for each basin outlet in Vector
`index_pit`.

# Arguments
- `graph` directed acyclic graph of the kinematic wave domain
- `toposort` topological order of `graph`
- `index_pit` Vector with basin outlets (pits)
- `streamorder` stream order of the kinematic wave domain
- `min_sto` minimum `streamorder` value

# Output
- `subbas_order` grouped subbasin ids (`Vector{Vector{Int}}`) ordered upstream (first index)
  to downstream (last index)
- `indices_subbas` list of indices per subbasin id stored as `Vector{Vector{Int}}`
- `topo_subbas` topological order per subbasin id stored as `Vector{Vector{Int}}`
"""
function kinwave_set_subdomains(graph, toposort, index_pit, streamorder, min_sto)
    if nthreads() > 1
        # extract basins (per outlet/pit), assign unique basin id
        n_pits = length(index_pit)
        basin = fill(Int(0), length(toposort))
        basin[index_pit] = [1:n_pits;]
        basin_fill = fillnodata_upstream(graph, toposort, basin, 0)

        # pre-allocate the Vector with indices matching the topological order of the
        # complete domain (upstream neighbors are stored at these indices)
        index_toposort = fill(Int(0), length(toposort))
        for (i, j) in enumerate(toposort)
            index_toposort[j] = i
        end

        order_subbas = Vector{Vector{Int}}()
        indices_subbas = Vector{Vector{Int}}()
        topo_subbas = Vector{Vector{Int}}()
        index = Vector{Int}()
        total_subbas = 0
        for i in 1:n_pits
            # extract subbasins per basin, make a graph at the subbasin level, calculate the
            # maximum distance of this graph, and group and order the subbasin ids from
            # upstream to downstream
            basin = findall(x -> x == i, basin_fill)
            g, vmap = induced_subgraph(graph, basin)
            toposort_b = topological_sort_by_dfs(g)
            streamorder_subbas = streamorder[vmap]
            subbas = subbasins(g, streamorder_subbas, toposort_b, min_sto)
            subbas_fill = fillnodata_upstream(g, toposort_b, subbas, 0)
            n_subbas = max(length(subbas[subbas .> 0]), 1)
            if n_subbas > 1
                graph_subbas = graph_from_nodes(g, subbas, subbas_fill)
                toposort_subbas = topological_sort_by_dfs(graph_subbas)
                dist = Graphs.Experimental.Traversals.distances(
                    Graph(graph_subbas),
                    toposort_subbas[end],
                )
                max_dist = maximum([dist; 1])
                v_subbas = subbasins_order(graph_subbas, toposort_subbas[end], max_dist)
            else
                v_subbas = [[1]]
            end
            # subbasins need a unique id (in case of multiple basins/outlets in the
            # kinematic wave domain)
            for n in 1:length(v_subbas)
                v_subbas[n] .= v_subbas[n] .+ total_subbas
            end
            total_subbas += n_subbas
            append!(order_subbas, v_subbas)
            append!(index, 1:length(v_subbas))
            # in case of multiple subbasins calculate topological order per subbasin
            # (subgraph of the corresponding basin graph g), and the indices that match the
            # subbasin topological order
            if n_subbas > 1
                for s in 1:n_subbas
                    subbas_s = findall(x -> x == s, subbas_fill)
                    sg, _ = induced_subgraph(g, subbas_s)
                    toposort_sg = topological_sort_by_dfs(sg)
                    push!(topo_subbas, basin[subbas_s[toposort_sg]])
                    push!(indices_subbas, index_toposort[basin[subbas_s[toposort_sg]]])
                end
            else
                push!(topo_subbas, basin[toposort_b])
                push!(indices_subbas, index_toposort[basin[toposort_b]])
            end
        end
        # reduce the order of subbasin ids by merging groups of subbasins that have the same
        # index (multiple basins/outlets in the kinematic wave domain)
        subbas_order = Vector{Vector{Int}}(undef, maximum(index))
        for m in 1:maximum(index)
            subbas_order[m] = reduce(vcat, order_subbas[index .== m])
        end
    else
        subbas_order = [[1]]
        indices_subbas = [[1:length(toposort);]]
        topo_subbas = [toposort]
    end

    return subbas_order, indices_subbas, topo_subbas
end
