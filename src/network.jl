"""
Struct for storing forward `indices` and reverse indices `reverse_indices` of `Drainage`
cells (boundary condition groundwater flow) in the 2D external model domain.
"""
@kwdef struct NetworkDrain
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    reverse_indices::Matrix{Int64} = zeros(Int, 0, 0)
end

"Struct for storing 2D staggered grid topology."
@with_kw struct Indices
    xu::Vector{Int} = Int[]     # index of neighbor cell in the (+1, 0) direction
    xd::Vector{Int} = Int[]     # index of neighbor cell in the (-1, 0) direction
    yu::Vector{Int} = Int[]     # index of neighbor cell in the (0, +1) direction
    yd::Vector{Int} = Int[]     # index of neighbor cell in the (0, -1) direction
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
    # staggered grid topology
    staggered_indices::Indices = Indices()
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