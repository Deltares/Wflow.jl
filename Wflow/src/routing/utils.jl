const KIN_WAVE_MIN_FLOW = 1e-30 # [m³ s⁻¹]

"Convert a gridded drainage direction to a directed graph"
function flowgraph(ldd::AbstractVector, indices::AbstractVector, PCR_DIR::AbstractVector)
    # prepare a directed graph to be filled
    n_cells = length(indices)
    graph = DiGraph(n_cells)

    # loop over ldd, adding the edge to the downstream node
    for (from_node, from_index) in enumerate(indices)
        ldd_val = ldd[from_node]
        # skip pits to prevent cycles
        ldd_val == 5 && continue
        to_index = from_index + PCR_DIR[ldd_val]
        # find the node id of the downstream cell
        to_node = searchsortedfirst(indices, to_index)
        add_edge!(graph, from_node, to_node)
    end
    if is_cyclic(graph)
        error("""One or more cycles detected in flow graph.
            The provided local drainage direction map may be unsound.
            Verify that each active flow cell flows towards a pit.
            """)
    end
    return graph
end

"""
    accucapacitystate!(material, network, capacity)

Transport of material downstream with a limited transport capacity over a directed graph.
Mutates the material input. The network is expected to hold a graph and order field, where
the graph implements the Graphs interface, and the order is a valid topological ordering
such as that returned by `Graphs.topological_sort_by_dfs`.
"""
function accucapacitystate!(material, network, capacity)
    (; graph, cell_order) = network
    for cell_idx in cell_order
        downstream_nodes = outneighbors(graph, cell_idx)
        n_downstream_cells = length(downstream_nodes)
        flux_val = min(material[cell_idx], capacity[cell_idx])
        material[cell_idx] -= flux_val
        if iszero(n_downstream_cells)
            # pit: material is transported out of the map if a capacity is set,
            # cannot add the material anywhere
        elseif isone(n_downstream_cells)
            material[only(downstream_nodes)] += flux_val
        else
            error("bifurcations not supported")
        end
    end
    return nothing
end

"""
    accucapacitystate!(material, network, capacity) -> material

Non mutating version of `accucapacitystate!`.
"""
function accucapacitystate(material, network, capacity)
    material = copy(material)
    accucapacitystate!(material, network, capacity)
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
function accucapacityflux!(flux, material, network, capacity)
    (; graph, cell_order) = network
    for cell_idx in cell_order
        downstream_nodes = outneighbors(graph, cell_idx)
        n_downstream_cells = length(downstream_nodes)
        flux_val = min(material[cell_idx], capacity[cell_idx])
        material[cell_idx] -= flux_val
        flux[cell_idx] = flux_val
        if iszero(n_downstream_cells)
            # pit: material is transported out of the map if a capacity is set,
            # cannot add the material anywhere
        elseif isone(n_downstream_cells)
            material[only(downstream_nodes)] += flux_val
        else
            error("bifurcations not supported")
        end
    end
    return nothing
end

"""
    accucapacityflux(material, network, capacity) -> flux

Non mutating version of `accucapacityflux!`.
"""
function accucapacityflux(material, network, capacity)
    flux = zero(material)
    accucapacityflux!(flux, material, network, capacity)
    return flux
end

"""
    accucapacityflux_state(material, network, capacity) -> flux, material

Non mutating version of combined `accucapacityflux!` and `accucapacitystate!`.
"""
function accucapacityflux_state(material, network, capacity)
    flux = zero(material)
    material = copy(material)
    accucapacityflux!(flux, material, network, capacity)
    return flux, material
end

function flux_in!(flux_in, flux, network)
    (; upstream_nodes, cell_order) = network
    for (i, cell_idx) in enumerate(cell_order)
        flux_in[cell_idx] = sum_at(flux, upstream_nodes[i])
    end
    return nothing
end
