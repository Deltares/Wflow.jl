const KIN_WAVE_MIN_FLOW = 1e-30 # [m³ s⁻¹]
const KIN_WAVE_MIN_FLOW_QROOT = KIN_WAVE_MIN_FLOW^0.2

"Convert a gridded drainage direction to a directed graph. Also returns the possibly modified
drainage direction, where invalid values are set to 5 (pit)."
function flowgraph(ldd::AbstractVector, indices::AbstractVector, PCR_DIR::AbstractVector)
    # prepare a directed graph to be filled
    n = length(indices)
    graph = DiGraph(n)

    # loop over ldd, adding the edge to the downstream node
    for (from_node, from_index) in enumerate(indices)
        ldd_val = ldd[from_node]
        # skip pits to prevent cycles
        ldd_val == 5 && continue
        to_index = from_index + PCR_DIR[ldd_val]
        # find the node id of the downstream cell
        to_node = searchsortedfirst(indices, to_index)
        if to_node > length(indices) || indices[to_node] != to_index
            @warn "Invalid drainage direction value at node `$from_node` (LDD=`$ldd_val`), assign pit value at node"
            ldd[from_node] = 5
            continue
        end
        add_edge!(graph, from_node, to_node)
    end
    if is_cyclic(graph)
        error("""One or more cycles detected in flow graph.
            The provided local drainage direction map may be unsound.
            Verify that each active flow cell flows towards a pit.
            """)
    end
    return graph, ldd
end

"""
    accucapacitystate!(material, network, capacity, dt)

Transport of material downstream with a limited transport capacity over a directed graph.
Mutates the material input. The network is expected to hold a graph and order field, where
the graph implements the Graphs interface, and the order is a valid topological ordering
such as that returned by `Graphs.topological_sort_by_dfs`.
"""
function accucapacitystate!(material, network, capacity, dt)
    (; graph, order) = network
    for v in order
        downstream_nodes = outneighbors(graph, v)
        n = length(downstream_nodes)
        flux_val = min(material[v], capacity[v])
        material[v] -= flux_val * dt
        if n == 0
            # pit: material is transported out of the map if a capacity is set,
            # cannot add the material anywhere
        elseif n == 1
            material[only(downstream_nodes)] += flux_val
        else
            error("bifurcations not supported")
        end
    end
    return nothing
end

"""
    accucapacitystate!(material, network, capacity, dt) -> material

Non mutating version of `accucapacitystate!`.
"""
function accucapacitystate(material, network, capacity, dt)
    material = copy(material)
    accucapacitystate!(material, network, capacity, dt)
    return material
end

"""
    accucapacityflux!(flux, material, network, capacity)

Transport of material downstream with a limited transport capacity over a directed graph.
Updates the material input, and overwrites the flux input, not using existing values. The
network is expected to hold a graph and order field, where the graph implements the Graphs
interface, and the order is a valid topological ordering such as that returned by
`Graphs.topological_sort_by_dfs`.
"""
function accucapacityflux!(flux, material, network, capacity, dt)
    (; graph, order) = network
    for v in order
        downstream_nodes = outneighbors(graph, v)
        n = length(downstream_nodes)

        (n > 1) && error("bifurcations not supported")

        # pit: material is transported out of the map if a capacity is set,
        # cannot add the material anywhere
        to_pit = iszero(n)

        # Let [u] be the unit of material
        flux_val = min(material[v] / dt, capacity[v])
        material_update = flux_val * dt
        material[v] -= material_update
        flux[v] = flux_val
        if !to_pit
            material[only(downstream_nodes)] += material_update
        end
    end
    return nothing
end

"""
    accucapacityflux_rate!(flux, remaining, material, network, capacity)

Transport of material rate downstream with a limited transport capacity over a directed
graph. All of `material`, `capacity`, `flux` and `remaining` are rates [u s⁻¹]. Overwrites
`flux` with the transported rate and `remaining` with the deposited (untransported) rate.
Does not mutate `material`. The network is expected to hold a graph and order field, where
the graph implements the Graphs interface, and the order is a valid topological ordering
such as that returned by `Graphs.topological_sort_by_dfs`.
"""
function accucapacityflux_rate!(flux, remaining, material, network, capacity)
    (; graph, order) = network
    remaining .= material
    for v in order
        downstream_nodes = outneighbors(graph, v)
        n = length(downstream_nodes)

        (n > 1) && error("bifurcations not supported")

        # pit: material is transported out of the map if a capacity is set,
        # cannot add the material anywhere
        to_pit = iszero(n)

        flux_val = min(remaining[v], capacity[v])
        remaining[v] -= flux_val
        flux[v] = flux_val
        if !to_pit
            remaining[only(downstream_nodes)] += flux_val
        end
    end
    return nothing
end

"""
    accucapacityflux(material, network, capacity, dt) -> flux

Non mutating version of `accucapacityflux!`.
"""
function accucapacityflux(material, network, capacity, dt)
    flux = zero(material)
    accucapacityflux!(flux, material, network, capacity, dt)
    return flux
end

"""
    accucapacityflux_state(material, network, capacity) -> flux, material

Non mutating version of combined `accucapacityflux!` and `accucapacitystate!`.
"""
function accucapacityflux_state(material, network, capacity, dt)
    flux = zero(material)
    material = copy(material)
    accucapacityflux!(flux, material, network, capacity, dt)
    return flux, material
end

function flux_in!(flux_in, flux, network)
    (; upstream_nodes, order) = network
    for (i, v) in enumerate(order)
        flux_in[v] = sum_at(flux, upstream_nodes[i])
    end
    return nothing
end
