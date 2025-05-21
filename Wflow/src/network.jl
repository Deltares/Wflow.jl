# maps the fields of struct `EdgeConnectivity` to the defined Wflow cartesian indices of
# const `neigbors`.
const DIRS = (:yd, :xd, :xu, :yu)

"""
Struct for storing 2D staggered grid edge connectivity in `x` and `y` directions. For
example used by the q-centered numerical scheme of Almeida et al. (2012), that uses
information on both neighboring cell edges (interfaces).

See also: de Almeida, G. A. M., P. D.Bates, J. Freer, and M. Souvignet (2012), Improving the
stability of a simple formulation of the shallow water equations for 2D flood modelling,
Water Resour. Res., 48, doi:10.1029/2011WR011570.

Edges without neigbors are handled by an extra index (at `n + 1`, with `n` edges). The
linear index `i` of the `EdgeConnectivity` fields represents the edge between node index `i`
and the neighboring nodes in the CartesianIndex(-1,0) and CartesianIndex(0,-1) directions.
The edges are defined as follows:
- `xu` is the edge between node `i` and node `xu` in the `CartesianIndex(1,0)` direction.
- `xd` is the edge between node `xd` in the `CartesianIndex(-1,0)` direction and the
  neighboring node (CartesianIndex(-2,0) direction).
- `yu` is the edge between node `i` and node `yu` in the `CartesianIndex(0,1)` direction.
- `yd` is the edge between node `yd` in the `CartesianIndex(0,-1)` direction and the
  neighboring node (`CartesianIndex(0,-2)` direction).
"""
@with_kw struct EdgeConnectivity
    xu::Vector{Int} = Int[]
    xd::Vector{Int} = Int[]
    yu::Vector{Int} = Int[]
    yd::Vector{Int} = Int[]
end

"Struct for storing source `src` node and destination `dst` node of an edge."
@kwdef struct NodesAtEdge{I <: DenseArray{Int}}
    src::I = Int[]
    dst::I = Int[]
end

function Adapt.adapt_structure(to, from::NodesAtEdge)
    return NodesAtEdge(adapt(to, from.src), adapt(to, from.dst))
end

"Struct for storing source `src` edge and destination `dst` edge of a node."
@kwdef struct EdgesAtNode{I <: DenseArray{Int}}
    src::I = Int[]
    dst::I = Int[]
end

function Adapt.adapt_structure(to, from::EdgesAtNode)
    return EdgesAtNode(adapt(to, from.src), adapt(to, from.dst))
end

"Struct for storing network information land domain."
@kwdef struct NetworkLand
    # dimension of 2D model
    modelsize::Tuple{Int, Int} = (Int(0), Int(0))
    # local drain direction using the 8 point pour algorithm
    local_drain_direction::Vector{UInt8} = UInt8[]
    # water allocation areas [-]
    allocation_area_indices::Vector{Vector{Int}} = Vector{Int}[]
    # directed acyclic graph
    graph::SimpleDiGraph{Int} = DiGraph{Int}(0)
    # Strahler stream order
    streamorder::Vector{Int} = Int[]
    # maps from the 1D internal land domain to the 2D model (external) domain
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    # traversion order of land domain
    order::Vector{Int} = Int[]
    # execution order of sub-domains for kinematic wave routing (land domain)
    order_of_subdomains::Vector{Vector{Int}} = Vector{Int}[]
    # traversion order per sub-domain
    order_subdomain::Vector{Vector{Int}} = Vector{Int}[]
    # maps from the 2D model (external) domain to the 1D internal land domain
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)
    # maps from the land domain to the river domain (zero value represents no river)
    river_indices::Vector{Int} = Int[]
    # maps from the land domain to the river domain excluding reservoir and lake locations
    river_inds_excl_waterbody::Vector{Int} = Int[]
    # 2D staggered grid edge indices
    edge_indices::EdgeConnectivity = EdgeConnectivity()
    # maps `order_subdomain` to traversion order of the complete domain
    subdomain_indices::Vector{Vector{Int}} = Vector{Int}[]
    # upstream nodes (directed graph)
    upstream_nodes::Vector{Vector{Int}} = Vector{Int}[]
end

"""
Initialize `NetworkLand` fields related to catchment (active indices model domain) and land
drainage network.
"""
function NetworkLand(dataset::NCDataset, config::Config, modelsettings::NamedTuple)
    lens = lens_input(config, "subbasin_location__count"; optional = false)
    subcatch_2d = ncread(dataset, config, lens; allow_missing = true)
    indices, reverse_indices = active_indices(subcatch_2d, missing)
    modelsize = size(subcatch_2d)
    graph, local_drain_direction =
        get_drainage_network(dataset, config, indices; do_pits = modelsettings.pits)
    order = topological_sort_by_dfs(graph)
    streamorder = stream_order(graph, order)

    network = NetworkLand(;
        modelsize,
        indices,
        reverse_indices,
        local_drain_direction,
        graph,
        order,
        streamorder,
    )
    return network
end

"""
Set subdomain fields of `NetworkLand` for the land domain. Subdomains are created when
nthreads > 1 to run the kinematic wave parallel, otherwise it is equal to the complete
domain.
"""
function network_subdomains(config::Config, network::NetworkLand)
    pit_inds = findall(x -> x == 5, network.local_drain_direction)
    min_streamorder = get(config.model, "land_streamorder__min_count", 5)
    order_of_subdomains, subdomain_inds, toposort_subdomain = kinwave_set_subdomains(
        network.graph,
        network.order,
        pit_inds,
        network.streamorder,
        min_streamorder,
    )
    @reset network.order_of_subdomains = order_of_subdomains
    @reset network.order_subdomain = toposort_subdomain
    @reset network.subdomain_indices = subdomain_inds

    return network
end

"Initialize `EdgeConnectivity`"
function EdgeConnectivity(network::NetworkLand)
    (; modelsize, indices, reverse_indices) = network
    n = Int(length(indices))
    edge_indices =
        EdgeConnectivity(; xu = zeros(n), xd = zeros(n), yu = zeros(n), yd = zeros(n))

    nrow, ncol = modelsize
    for (v, i) in enumerate(indices)
        for (m, neighbor) in enumerate(NEIGHBORS)
            j = i + neighbor
            dir = DIRS[m]
            if (1 <= j[1] <= nrow) && (1 <= j[2] <= ncol) && (reverse_indices[j] != 0)
                getfield(edge_indices, dir)[v] = reverse_indices[j]
            else
                getfield(edge_indices, dir)[v] = n + 1
            end
        end
    end
    return edge_indices
end

"Return local drain direction `ldd` and directed `graph` based on `ldd`"
function get_drainage_network(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    do_pits = false,
    logging = true,
)
    lens = lens_input(config, "basin__local_drain_direction"; optional = false)
    ldd_2d = ncread(dataset, config, lens; allow_missing = true, logging)
    ldd = convert(Array{UInt8}, ldd_2d[indices])
    if do_pits
        lens = lens_input(config, "basin_pit_location__mask"; optional = false)
        pits_2d = ncread(dataset, config, lens; type = Bool, fill = false)
        ldd = set_pit_ldd(pits_2d, ldd, indices)
    end
    graph = flowgraph(ldd, indices, PCR_DIR)
    return graph, ldd
end

@kwdef struct NetworkRiver{I <: DenseArray{Int}, I8 <: DenseArray{UInt8}}
    local_drain_direction::I8 = UInt8[]                           # local drain direction using the 8 point pour algorithm
    allocation_area_indices::Vector{Vector{Int}} = Vector{Int}[]  # water allocation areas [-]
    edges_at_node::EdgesAtNode = EdgesAtNode()                 # source and destination edge of a node
    graph::SimpleDiGraph{Int} = DiGraph{Int}(0)                   # directed graph
    streamorder::I = Int[]                                        # Strahler stream order
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]      # maps from the 1D internal river domain to the 2D model (external) domain
    lake_indices::I = Int[]                                       # maps lakes to the river domain (zero value represents no lake)
    land_indices::I = Int[]                                       # land domain indices masked by river domain (zero value represents no river)
    pit_indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]  # maps 1D pits (local drain direction) to the 2D model (external) domain
    nodes_at_edge::NodesAtEdge = NodesAtEdge()                 # source and destination node of an edge
    order::I = Int[]                                              # traversion order of river domain
    order_of_subdomains::Vector{Vector{Int}} = Vector{Int}[]      # execution order of sub-domains for kinematic wave river flow routing
    order_subdomain::Vector{Vector{Int}} = Vector{Int}[]          # traversion order per sub-domain
    reservoir_indices::I = Int[]                                  # maps reservoirs to the river domain (zero value represents no reservoir)
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)                         # maps from the 2D model (external) domain to the 1D internal river domain
    subdomain_indices::Vector{Vector{Int}} = Vector{Int}[]        # maps `order_subdomain` to traversion order of the complete domain
    upstream_nodes::Vector{Vector{Int}} = Vector{Int}[]           # upstream nodes (directed graph)
end

"""
Initialize `NetworkRiver` fields related to river location (active indices model domain) and
river drainage network.
"""
function NetworkRiver(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand;
    do_pits = false,
)
    logging = false
    lens = lens_input(config, "river_location__mask"; optional = false)
    river_location_2d = ncread(dataset, config, lens; type = Bool, fill = false, logging)
    indices, reverse_indices = active_indices(river_location_2d, 0)
    graph, local_drain_direction =
        get_drainage_network(dataset, config, indices; do_pits, logging)
    order = topological_sort_by_dfs(graph)
    river_location = river_location_2d[network.indices]
    land_indices =
        Vector{Int}(filter(i -> !isequal(river_location[i], 0), 1:length(network.indices)))
    streamorder = network.streamorder[land_indices]

    network = NetworkRiver(;
        local_drain_direction,
        graph,
        streamorder,
        indices,
        land_indices,
        order,
        reverse_indices,
    )

    return network
end

"""
Set subdomain fields of `NetworkRiver` for the river domain. Subdomains are created when
nthreads > 1 to run the kinematic wave parallel, otherwise it is equal to the complete
domain.
"""
function network_subdomains(config::Config, network::NetworkRiver)
    min_streamorder = get(config.model, "river_streamorder__min_count", 6)
    pit_inds = findall(x -> x == 5, network.local_drain_direction)
    order_of_subdomains, subdomain_inds, toposort_subdomain = kinwave_set_subdomains(
        network.graph,
        network.order,
        pit_inds,
        network.streamorder,
        min_streamorder,
    )
    @reset network.order_of_subdomains = order_of_subdomains
    @reset network.order_subdomain = toposort_subdomain
    @reset network.subdomain_indices = subdomain_inds
    return network
end

"Initialize `NodesAtEdge`"
function NodesAtEdge(network::NetworkRiver)
    index_pit = findall(x -> x == 5, network.local_drain_direction)
    add_vertex_edge_graph!(network.graph, index_pit)
    nodes_at_edge = NodesAtEdge(; adjacent_nodes_at_edge(network.graph)...)
    return nodes_at_edge, index_pit
end

"Initialize `EdgesAtNode`"
function EdgesAtNode(network::NetworkRiver)
    (; src, dst) = adjacent_edges_at_node(network.graph, network.nodes_at_edge)
    edges_at_node = EdgesAtNode(ragged_to_dense(src), ragged_to_dense(dst))
    return edges_at_node
end

"Struct for storing network information water body (reservoir or lake)."
@kwdef struct NetworkWaterBody
    # list of 2D indices representing water body area (coverage)
    indices_coverage::Vector{Vector{CartesianIndex{2}}} = Vector{CartesianIndex{2}}[]
    # list of 2D indices representing water body outlet
    indices_outlet::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    # maps from the 2D model (external) domain to a list of water bodies
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)
    # maps from the 1D river domain to a list of water bodies (zero value represents no water body)
    river_indices::Vector{Int} = Int[]
end

"Initialize `NetworkWaterBody`"
function NetworkWaterBody(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    waterbody_type::String,
)
    logging = false
    # allow waterbody only in river cells
    # note that these locations are only the waterbody outlet pixels
    lens = lens_input(config, "$(waterbody_type)_location__count"; optional = false)
    locs = ncread(dataset, config, lens; sel = indices, type = Int, fill = 0, logging)

    # this holds the same ids as locs, but covers the entire reservoir or lake
    lens = lens_input(config, "$(waterbody_type)_area__count"; optional = false)
    coverage_2d = ncread(dataset, config, lens; allow_missing = true, logging)
    # for each waterbody, a list of 2D indices, needed for getting the mean precipitation
    inds_coverage = Vector{CartesianIndex{2}}[]
    rev_inds = zeros(Int, size(coverage_2d))

    # construct a map from the rivers to the waterbody and
    # a map of the waterbody to the 2D model grid
    inds_map2river = fill(Int(0), length(indices))
    inds = CartesianIndex{2}[]
    counter = 0
    for (i, ind) in enumerate(indices)
        id = locs[i]
        if id > 0
            push!(inds, ind)
            counter += 1
            inds_map2river[i] = counter
            rev_inds[ind] = counter

            # get all indices related to this waterbody outlet
            # done in this loop to ensure that the order is equal to the order in the
            # waterbody model struct
            cov = findall(isequal(id), coverage_2d)
            push!(inds_coverage, cov)
        end
    end
    network = NetworkWaterBody(;
        indices_outlet = inds,
        indices_coverage = inds_coverage,
        reverse_indices = rev_inds,
        river_indices = findall(x -> x ≠ 0, inds_map2river),
    )
    return network, inds_map2river
end

"""
Struct for storing forward `indices` and reverse indices `reverse_indices` in the 2D
external model domain, and 1D land domain indices `land_indices` of `Drainage` cells
(boundary condition groundwater flow).
"""
@kwdef struct NetworkDrain
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)
    land_indices::Vector{Int} = Int[]
end

"Initialize `NetworkDrain`"
function NetworkDrain(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    surface_flow_width::Vector{Float},
)
    n_cells = length(indices)
    lens = lens_input_parameter(config, "land_drain_location__mask")
    drain_2d = ncread(dataset, config, lens; type = Bool, fill = false)
    drain = drain_2d[indices]

    # check if drain occurs where overland flow is not possible (surface_flow_width = 0.0)
    # and correct if this is the case
    false_drain =
        filter(i -> !isequal(drain[i], 0) && surface_flow_width[i] == 0.0, 1:n_cells)
    n_false_drain = length(false_drain)
    if n_false_drain > 0
        drain_2d[indices[false_drain]] .= 0
        drain[false_drain] .= 0
        @info "$n_false_drain drain locations are removed that occur where overland flow
         is not possible (overland flow width is zero)"
    end
    land_indices = filter(i -> !isequal(drain[i], 0), 1:n_cells)
    indices, reverse_indices = active_indices(drain_2d, 0)
    network = NetworkDrain(; indices, reverse_indices, land_indices)
    return network
end