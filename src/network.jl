@kwdef struct NetworkDrain
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    reverse_indices::Matrix{Int64} = zeros(Int, 0, 0)
end

@kwdef struct NetworkLake
    river_indices::Vector{Int} = Int[]
end

# Stores edges in x and y direction between cells of a Vector with CartesianIndex(x, y), for
# staggered grid calculations.
@with_kw struct Indices
    xu::Vector{Int} = Int[]     # index of neighbor cell in the (+1, 0) direction
    xd::Vector{Int} = Int[]     # index of neighbor cell in the (-1, 0) direction
    yu::Vector{Int} = Int[]     # index of neighbor cell in the (0, +1) direction
    yd::Vector{Int} = Int[]     # index of neighbor cell in the (0, -1) direction
end

@kwdef struct NetworkLand
    allocation_area_indices::Vector{Vector{Int64}} = Vector{Int}[]
    altitude::Vector{Float64} = Float64[]
    area::Vector{Float64} = Float64[]
    frac_to_river::Vector{Float64} = Float64[]
    graph::SimpleDiGraph{Int} = DiGraph(0)
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    order::Vector{Int} = Int[]
    order_of_subdomains::Vector{Vector{Int}} = Vector{Int}[]
    order_subdomain::Vector{Vector{Int}} = Vector{Int}[]
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)
    river_indices::Vector{Int} = Int[]
    river_inds_excl_waterbody::Vector{Int} = Int[]
    slope::Vector{Float64} = Float64[]
    staggered_indices::Indices = Indices()
    subdomain_indices::Vector{Vector{Int}} = Vector{Int}[]
    upstream_nodes::Vector{Vector{Int}} = Vector{Int}[]
    waterbody::Vector{Bool} = Bool[]
end

@kwdef struct NetworkReservoir
    indices_coverage::Vector{Vector{CartesianIndex{2}}} = Vector{CartesianIndex{2}}[]
    indices_outlet::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)
    river_indices::Vector{Int} = Int[]
end

@kwdef struct NodesAtEdge
    src::Vector{Int} = Int[]
    dst::Vector{Int} = Int[]
end

@kwdef struct EdgesAtNode
    src::Vector{Vector{Int}} = Vector{Int}[]
    dst::Vector{Vector{Int}} = Vector{Int}[]
end

@kwdef struct NetworkRiver
    allocation_area_indices::Vector{Vector{Int64}} = Vector{Int}[]
    cell_area::Vector{Float64} = Float64[]
    edges_at_node::EdgesAtNode = EdgesAtNode()
    graph::SimpleDiGraph{Int} = DiGraph(0)
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    lake_indices::Vector{Int} = Int[]
    land_indices::Vector{Int} = Int[]
    nodes_at_edge::NodesAtEdge = NodesAtEdge()
    order::Vector{Int} = Int[]
    order_of_subdomains::Vector{Vector{Int}} = Vector{Int}[]
    order_subdomain::Vector{Vector{Int}} = Vector{Int}[]
    reservoir_indices::Vector{Int} = Int[]
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)
    subdomain_indices::Vector{Vector{Int}} = Vector{Int}[]
    upstream_nodes::Vector{Vector{Int}} = Vector{Int}[]
end

@kwdef struct Network
    drain::NetworkDrain = NetworkDrain()
    lake::NetworkLake = NetworkLake()
    land::NetworkLand = NetworkLand()
    reservoir::NetworkReservoir = NetworkReservoir()
    river::NetworkRiver = NetworkRiver()
    frac_to_river::Vector{Float64} = Float64[]
    index_river::Vector{Int} = Int[]
end

@kwdef struct Lateral{L, R, S}
    land::L = nothing
    river::R = nothing
    subsurface::S = nothing
end
