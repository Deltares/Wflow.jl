# maps the fields of struct `EdgeConnectivity` to the defined Wflow cartesian indices of
# const `neigbors`.
const dirs = (:yd, :xd, :xu, :yu)

"""
Struct for storing forward `indices` and reverse indices `reverse_indices` of `Drainage`
cells (boundary condition groundwater flow) in the 2D external model domain.
"""
@kwdef struct NetworkDrain
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    reverse_indices::Matrix{Int64} = zeros(Int, 0, 0)
end

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

"Struct for storing network information land domain."
@kwdef struct NetworkLand
    # water allocation areas [-]
    allocation_area_indices::Vector{Vector{Int64}} = Vector{Int}[]
    # elevation [m]
    altitude::Vector{Float64} = Float64[]
    # grid cell area [m²]
    area::Vector{Float64} = Float64[]
    # flow fraction [-] to river
    frac_to_river::Vector{Float64} = Float64[]
    # directed acyclic graph
    graph::SimpleDiGraph{Int} = DiGraph(0)
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
    # slope [m m⁻¹]
    slope::Vector{Float64} = Float64[]
    # 2D staggered grid edge indices
    edge_indices::EdgeConnectivity = EdgeConnectivity()
    # maps `order_subdomain` to traversion order of the complete domain
    subdomain_indices::Vector{Vector{Int}} = Vector{Int}[]
    # upstream nodes (directed graph)
    upstream_nodes::Vector{Vector{Int}} = Vector{Int}[]
    # water body (reservoir and lake) location
    waterbody::Vector{Bool} = Bool[]
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

"Struct for storing source `src` node and destination `dst` node of an edge"
@kwdef struct NodesAtEdge
    src::Vector{Int} = Int[]
    dst::Vector{Int} = Int[]
end

"Struct for storing source `src` edge and destination `dst` edge of a node"
@kwdef struct EdgesAtNode
    src::Vector{Vector{Int}} = Vector{Int}[]
    dst::Vector{Vector{Int}} = Vector{Int}[]
end

"Struct for storing network information river domain."
@kwdef struct NetworkRiver
    # water allocation areas [-]
    allocation_area_indices::Vector{Vector{Int64}} = Vector{Int}[]
    # grid cell area [m²]
    cell_area::Vector{Float64} = Float64[]
    # source and destination edge of a node
    edges_at_node::EdgesAtNode = EdgesAtNode()
    # directed graph
    graph::SimpleDiGraph{Int} = DiGraph(0)
    # maps from the 1D internal river domain to the 2D model (external) domain
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    # maps lakes to the river domain (zero value represents no lake)
    lake_indices::Vector{Int} = Int[]
    # land domain indices masked by river domain (zero value represents no river)
    land_indices::Vector{Int} = Int[]
    # source and destination node of an edge
    nodes_at_edge::NodesAtEdge = NodesAtEdge()
    # traversion order of river domain
    order::Vector{Int} = Int[]
    # execution order of sub-domains for kinematic wave river flow routing
    order_of_subdomains::Vector{Vector{Int}} = Vector{Int}[]
    # traversion order per sub-domain
    order_subdomain::Vector{Vector{Int}} = Vector{Int}[]
    # maps reservoirs to the river domain (zero value represents no reservoir)
    reservoir_indices::Vector{Int} = Int[]
    # maps from the 2D model (external) domain to the 1D internal river domain
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)
    # maps `order_subdomain` to traversion order of the complete domain
    subdomain_indices::Vector{Vector{Int}} = Vector{Int}[]
    # upstream nodes (directed graph)
    upstream_nodes::Vector{Vector{Int}} = Vector{Int}[]
end

""" 
Struct for storing network information of different domains: `land`, `drain` (`Drainage`
boundary condition of groundwater flow), `river`, `reservoir` and `lake`.
"""
@kwdef struct Network
    drain::NetworkDrain = NetworkDrain()
    lake::NetworkWaterBody = NetworkWaterBody()
    land::NetworkLand = NetworkLand()
    reservoir::NetworkWaterBody = NetworkWaterBody()
    river::NetworkRiver = NetworkRiver()
end

function NetworkReservoir(
    config::Config,
    dataset::NCDataset,
    to_do::NamedTuple,
    data::NamedTuple,
    ::AbstractSbmModelType,
)::Tuple{NetworkWaterBody, Vector{Int}, Union{Nothing, SimpleReservoir}}
    (; river_data, pits) = data
    (; inds_river) = river_data
    n_river_cells = length(inds_river)
    if to_do.reservoirs
        reservoir, reservoir_network, inds_reservoir_map2river =
            SimpleReservoir(dataset, config, inds_river, n_river_cells, pits)
        network_reservoir = NetworkWaterBody(; reservoir_network...)
    else
        network_reservoir = NetworkWaterBody()
        inds_reservoir_map2river = zeros(Int, n_river_cells)
        reservoir = nothing
    end
    return network_reservoir, inds_reservoir_map2river, reservoir
end

function NetworkLake(
    config::Config,
    dataset::NCDataset,
    to_do::NamedTuple,
    data::NamedTuple,
    ::AbstractSbmModelType,
)
    (; inds_river) = data.river_data
    n_river_cells = length(inds_river)
    if to_do.lakes
        lake, lake_network, inds_lake_map2river, pits =
            Lake(dataset, config, inds_river, n_river_cells, pits)
        network_lake = NetworkWaterBody(; lake_network...)
    else
        network_lake = NetworkWaterBody()
        inds_lake_map2river = fill(0, n_river_cells)
        lake = nothing
    end
    network_lake, inds_lake_map2river, lake
end

function NetworkLand(
    config::Config,
    to_do::NamedTuple,
    data::NamedTuple,
    ::AbstractSbmModelType,
)::Tuple{NetworkLand, NamedTuple}
    (; sub_catchment_data, cell_data, river_data, pits, routing_types) = data
    (; ldd, indices, reverse_indices, land_slope) = sub_catchment_data

    graph = flowgraph(ldd, indices, pcr_dir)

    (;
        streamorder,
        toposort,
        order_of_subdomains,
        subdomain_inds,
        toposort_subdomain,
        min_streamorder_land,
    ) = get_kinwave_subdomains_land(config, graph, data, to_do)

    data = (; data..., min_streamorder_land, streamorder)

    frac_to_river =
        fraction_runoff_to_river(graph, ldd, river_data.inds_land_map2river, land_slope)

    network_land = NetworkLand(;
        graph,
        upstream_nodes = filter_upsteam_nodes(graph, pits[indices]),
        order_of_subdomains,
        order_subdomain = toposort_subdomain,
        subdomain_indices = subdomain_inds,
        order = toposort,
        indices,
        reverse_indices,
        area = cell_data.x_length .* cell_data.y_length,
        slope = land_slope,
        frac_to_river,
    )

    return network_land, data
end

function NetworkRiver(
    config::Config,
    data::NamedTuple,
    to_do::NamedTuple,
    ::AbstractSbmModelType,
)::Tuple{NetworkRiver, NamedTuple}
    (;
        sub_catchment_data,
        cell_data,
        reservoir_data,
        lake_data,
        river_data,
        land_hydrology,
    ) = data
    (; inds_land_map2river, ldd_river, inds_river, reverse_inds_river) = river_data

    if to_do.pits
        ldd_river = set_pit_ldd(pits_2d, ldd_river, inds_river)
    end

    allocation_area_inds, river_allocation_area_inds =
        get_allocation_area_inds(land_hydrology, data, to_do)

    graph_river = flowgraph(ldd_river, inds_river, pcr_dir)

    (;
        min_streamorder_river,
        order_of_river_subdomains,
        river_subdomain_inds,
        toposort_river,
    ) = get_kinwave_subdomains_river(config, graph_river, data)

    data = (; data..., min_streamorder_river)

    network_river = NetworkRiver(;
        graph = graph_river,
        indices = inds_river,
        reverse_indices = reverse_inds_river,
        reservoir_indices = reservoir_data.inds_reservoir_map2river,
        land_indices = inds_land_map2river,
        lake_indices = lake_data.inds_lake_map2river,
        allocation_area_indices = river_allocation_area_inds,
        cell_area = cell_data.x_length[inds_land_map2river] .*
                    cell_data.y_length[inds_land_map2river],
        order_of_subdomains = order_of_river_subdomains,
        subdomain_indices = river_subdomain_inds,
        order = toposort_river,
    )

    return network_river, data
end

function Network(
    config::Config,
    dataset::NCDataset,
    data::NamedTuple,
    to_do::NamedTuple,
    type::AbstractModelType,
)
    (; routing_types) = data
    (; river_routing, land_routing) = routing_types
    pits = zeros(Bool, size(data.sub_catchment_data.subcatch_2d))
    data = (; data..., pits)

    network_land, data = NetworkLand(config, to_do, data, type)
    network_reservoir, inds_reservoir_map2river, reservoir =
        NetworkReservoir(config, dataset, to_do, data, type)
    network_lake, inds_lake_map2river, lake =
        NetworkLake(config, dataset, to_do, data, type)

    reservoir_data = (; inds_reservoir_map2river)
    lake_data = (; inds_lake_map2river)
    data = (; data..., reservoir_data, lake_data, reservoir, lake)
    network_river, data = NetworkRiver(config, data, to_do, type)

    network_drain = NetworkDrain(config, dataset, data, to_do, type)
    network = Network(;
        land = network_land,
        reservoir = network_reservoir,
        lake = network_lake,
        river = network_river,
        drain = network_drain,
    )
    data = (; data..., network)

    if nthreads() > 1
        if river_routing == "kinematic-wave"
            @info "Parallel execution of kinematic wave" data.min_streamorder_land data.min_streamorder_river
        elseif land_routing == "kinematic-wave" || to_do.lateral_ssf
            @info "Parallel execution of kinematic wave" data.min_streamorder_land
        end
    end

    return network, data
end

NetworkReservoir(config, dataset, data, to_do, ::AbstractModelType) = NetworkWaterBody()
NetworkLake(config, dataset, data, to_do, ::AbstractModelType) = NetworkWaterBody()
NetworkDrain(config, dataset, data, to_do, ::AbstractModelType) = NetworkDrain()