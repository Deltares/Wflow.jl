"Struct to store (shared) land parameters"
@with_kw struct LandParameters{T <: DenseArray{Float}, B <: DenseArray{Bool}}
    # cell length x direction [m]
    x_length::T = Float[]
    # cell length y direction [m]
    y_length::T = Float[]
    # cell area [m²]
    area::T = Float[]
    # flow width [m]
    flow_width::T = Float[]
    # suface flow width [m]
    surface_flow_width::T = Float[]
    # flow length [m]
    flow_length::T = Float[]
    # flow fraction to river [-] 
    flow_fraction_to_river::T = Float[]
    # slope [-]
    slope::T = Float[]
    # water body (reservoir and lake) location [-]
    waterbody_outlet::B = Bool[]
    # waterbody coverage [-]
    waterbody_coverage::B = Bool[]
    # river location [-]
    river_location::B = Bool[]
    # fraction of river [-]
    river_fraction::T = Float[]
    # fraction of open water (excluding rivers) [-]
    water_fraction::T = Float[]
end

function Adapt.adapt_structure(to, from::LandParameters)
    return LandParameters(
        adapt(to, from.x_length),
        adapt(to, from.y_length),
        adapt(to, from.area),
        adapt(to, from.flow_width),
        adapt(to, from.flow_length),
        adapt(to, from.flow_fraction_to_river),
        adapt(to, from.slope),
        adapt(to, from.waterbody_outlet),
        adapt(to, from.waterbody_coverage),
        adapt(to, from.river_location),
        adapt(to, from.river_fraction),
        adapt(to, from.water_fraction),
    )
end

"Struct to store (shared) river parameters"
@with_kw struct RiverParameters{T <: DenseArray{Float}, B <: DenseArray{Bool}}
    flow_width::T = Float[]          # river flow width [m]
    flow_length::T = Float[]         # river flow length [m]
    slope::T = Float[]               # slope [-]
    waterbody_outlet::B = Bool[]    # water body (reservoir and lake) location
    waterbody_coverage::B = Bool[]  # waterbody coverage [-]
    cell_area::T = Float[]          # grid cell area [m²]
end

@kwdef struct DomainLand
    network::NetworkLand = NetworkLand()
    parameters::LandParameters = LandParameters()
end

@kwdef struct DomainRiver{
    T <: DenseArray{Float},
    I <: DenseArray{Int},
    I8 <: DenseArray{UInt8},
    B <: DenseArray{Bool},
}
    network::NetworkRiver{I, I8} = ()
    parameters::RiverParameters{T, B} = RiverParameters()
end

function Adapt.adapt_structure(to, from::DomainRiver)
    return DomainRiver(adapt(to, from.network), adapt(to, from.parameters))
end

@kwdef struct DomainWaterBody
    network::NetworkWaterBody = NetworkWaterBody()
end

function Adapt.adapt_structure(to, from::DomainWaterBody)
    return DomainWaterBody(adapt(to, from.network))
end

@kwdef struct DomainDrain
    network::NetworkDrain = NetworkDrain()
end

function Adapt.adapt_structure(to, from::DomainDrain)
    return DomainDrain(adapt(to, from.network))
end

"""
Struct for storing information about different model domains 'land`, `river`, `reservoir`,
`lake` and `drain` (`Drainage` boundary condition of `GroundwaterFlow`). It holds network
information for each domain like active indices in the 2D model domain and connectivity
information for flow routing. The `land` and `river` domains contain shared parameters for
each domain that can used by different model components.
"""
@kwdef struct Domain
    land::DomainLand = DomainLand()
    river::DomainRiver = DomainRiver()
    reservoir::DomainWaterBody = DomainWaterBody()
    lake::DomainWaterBody = DomainWaterBody()
    drain::DomainDrain = DomainDrain()
end

function Adapt.adapt_structure(to, from::Domain)
    return Domain(
        adapt(to, from.land),
        adapt(to, from.river),
        adapt(to, from.reservoir),
        adapt(to, from.lake),
        adapt(to, from.drain),
    )
end

"Initialize `Domain` for model types `sbm` and `sbm_gwf`"
function Domain(
    dataset::NCDataset,
    config::Config,
    modelsettings::NamedTuple,
    routing_types::NamedTuple,
)
    network_land = NetworkLand(dataset, config, modelsettings)
    if routing_types.land == "kinematic-wave" ||
       routing_types.subsurface == "kinematic-wave"
        network_land = network_subdomains(config, network_land)
    end

    network_river =
        NetworkRiver(dataset, config, network_land; do_pits = modelsettings.pits)
    if routing_types.river == "kinematic-wave"
        network_river = network_subdomains(config, network_river)
    end

    pits = zeros(Bool, network_land.modelsize)
    nriv = Int(length(network_river.indices))
    if modelsettings.reservoirs
        network_reservoir, inds_reservoir_map2river =
            NetworkWaterBody(dataset, config, network_river.indices, "reservoir")
        pits[network_reservoir.indices_outlet] .= true
    else
        network_reservoir = NetworkWaterBody()
        inds_reservoir_map2river = fill(Int(0), nriv)
    end
    @reset network_river.reservoir_indices = inds_reservoir_map2river

    if modelsettings.lakes
        network_lake, inds_lake_map2river =
            NetworkWaterBody(dataset, config, network_river.indices, "lake")
        pits[network_lake.indices_outlet] .= true
    else
        network_lake = NetworkWaterBody()
        inds_lake_map2river = fill(Int(0), nriv)
    end
    @reset network_river.lake_indices = inds_lake_map2river

    if routing_types.river == "kinematic-wave"
        @reset network_river.upstream_nodes =
            filter_upsteam_nodes(network_river.graph, pits[network_river.indices])
    elseif routing_types.river == "local-inertial"
        nodes_at_edge, index_pit = NodesAtEdge(network_river)
        @reset network_river.nodes_at_edge = nodes_at_edge
        @reset network_river.pit_indices = network_river.indices[index_pit]
        @reset network_river.edges_at_node = EdgesAtNode(network_river)
    end

    if routing_types.land == "kinematic-wave" ||
       routing_types.subsurface == "kinematic-wave"
        @reset network_land.upstream_nodes =
            filter_upsteam_nodes(network_land.graph, pits[network_land.indices])
    end
    if routing_types.land == "local-inertial"
        @reset network_land.edge_indices = EdgeConnectivity(network_land)
        @reset network_land.river_indices =
            network_river.reverse_indices[network_land.indices]
    end

    domain_land = DomainLand(; network = network_land)
    domain_river = DomainRiver(; network = network_river)

    domain = Domain(;
        land = domain_land,
        river = domain_river,
        reservoir = DomainWaterBody(; network = network_reservoir),
        lake = DomainWaterBody(; network = network_lake),
    )

    land_params, river_params = initialize_shared_parameters(dataset, config, domain)
    @reset domain.land.parameters = land_params
    @reset domain.river.parameters = river_params

    if modelsettings.drains
        (; indices) = domain.land.network
        (; surface_flow_width) = domain.land.parameters
        @reset domain.drain.network =
            NetworkDrain(dataset, config, indices, surface_flow_width)
    end

    if modelsettings.water_demand
        allocation_area_indices, river_allocation_indices =
            get_allocation_area_indices(dataset, config, domain)
        @reset domain.land.network.allocation_area_indices = allocation_area_indices
        @reset domain.river.network.allocation_area_indices = river_allocation_indices

        mask = copy(domain.river.network.reverse_indices)
        mask_waterbody_coverage!(mask, config, domain)
        @reset domain.land.network.river_inds_excl_waterbody =
            mask[domain.land.network.indices]
    end

    if nthreads() > 1
        (; min_streamorder_land, min_streamorder_river) = modelsettings
        if routing_types.river == "kinematic-wave"
            @info "Parallel execution of kinematic wave" min_streamorder_land min_streamorder_river
        elseif routing_types.land == "kinematic-wave" ||
               routing_types.subsurface == "kinematic-wave"
            @info "Parallel execution of kinematic wave" min_streamorder_land
        end
    end

    return domain
end

"Initialize `Domain` for model type `sediment`"
function Domain(dataset::NCDataset, config::Config, modelsettings::NamedTuple)
    network_land = NetworkLand(dataset, config, modelsettings)
    network_river = NetworkRiver(dataset, config, network_land)

    domain_land = DomainLand(; network = network_land)
    domain_river = DomainRiver(; network = network_river)

    domain = Domain(; land = domain_land, river = domain_river)

    land_params, river_params = initialize_shared_parameters(dataset, config, domain)
    @reset domain.land.parameters = land_params
    @reset domain.river.parameters = river_params
end

"Initialize (shared) land parameters for model type `sediment`"
function LandParameters(dataset::NCDataset, config::Config, network::NetworkLand)
    x_length, y_length = get_cell_lengths(dataset, config, network)
    area = x_length .* y_length
    flow_width = map(get_flow_width, network.local_drain_direction, x_length, y_length)
    slope = get_landsurface_slope(dataset, config, network)
    waterbody_outlet = waterbody_mask(dataset, config, network)
    waterbody_coverage = waterbody_mask(dataset, config, network; region = "area")
    river_location = river_mask(dataset, config, network)

    land_parameters = LandParameters(;
        area,
        flow_width,
        slope,
        waterbody_outlet,
        waterbody_coverage,
        river_location,
    )
    return land_parameters
end

"Initialize (shared) land parameters for model types `sbm` and `sbm_gwf`"
function LandParameters(dataset::NCDataset, config::Config, domain::Domain)
    (; land_indices) = domain.river.network
    (; network) = domain.land
    x_length, y_length = get_cell_lengths(dataset, config, network)
    area = x_length .* y_length
    flow_width = map(get_flow_width, network.local_drain_direction, x_length, y_length)
    flow_length = map(get_flow_length, network.local_drain_direction, x_length, y_length)
    slope = get_landsurface_slope(dataset, config, network)
    river_location = river_mask(dataset, config, network)
    river_fraction = get_river_fraction(dataset, config, network, river_location, area)

    water_fraction = get_water_fraction(dataset, config, network, river_fraction)

    land_area = @. (Float(1.0) - river_fraction) * area
    surface_flow_width =
        map(get_surface_width, flow_width, flow_length, land_area, river_location)

    flow_fraction_to_river = get_flow_fraction_to_river(
        network.graph,
        network.local_drain_direction,
        land_indices,
        slope,
    )

    waterbody_outlet = waterbody_mask(dataset, config, network)
    waterbody_coverage = waterbody_mask(dataset, config, network; region = "area")

    land_parameters = LandParameters(;
        x_length,
        y_length,
        area,
        flow_width,
        surface_flow_width,
        flow_length,
        slope,
        river_location,
        flow_fraction_to_river,
        waterbody_outlet,
        waterbody_coverage,
        river_fraction,
        water_fraction,
    )
    return land_parameters
end

"Initialize (shared) river parameters"
function RiverParameters(dataset::NCDataset, config::Config, network::NetworkRiver)
    (; indices) = network
    lens = lens_input_parameter(config, "river__length"; optional = false)
    flow_length = ncread(dataset, config, lens; sel = indices, type = Float)
    minimum(flow_length) > 0 || error("river length must be positive on river cells")

    lens = lens_input_parameter(config, "river__width"; optional = false)
    flow_width = ncread(dataset, config, lens; sel = indices, type = Float)
    minimum(flow_width) > 0 || error("river width must be positive on river cells")

    lens = lens_input_parameter(config, "river__slope"; optional = false)
    slope = ncread(dataset, config, lens; sel = indices, type = Float)
    clamp!(slope, Float(0.00001), Float(Inf))

    river_parameters = RiverParameters(; flow_width, flow_length, slope)
    return river_parameters
end

"Initialize (shared) parameters for the river and land domains"
function initialize_shared_parameters(dataset::NCDataset, config::Config, domain::Domain)
    land_params = LandParameters(dataset, config, domain)
    river_params = RiverParameters(dataset, config, domain.river.network)

    @reset river_params.cell_area = land_params.area[domain.river.network.land_indices]
    @reset river_params.waterbody_coverage =
        land_params.waterbody_coverage[domain.river.network.land_indices]
    @reset river_params.waterbody_outlet =
        land_params.waterbody_outlet[domain.river.network.land_indices]

    return land_params, river_params
end

"Return open water fraction (excluding rivers)"
function get_water_fraction(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand,
    river_fraction::Vector{Float},
)
    lens = lens_input_parameter(config, "land~water-covered__area_fraction")
    water_fraction =
        ncread(dataset, config, lens; sel = network.indices, defaults = 0.0, type = Float)
    water_fraction = max.(water_fraction .- river_fraction, Float(0.0))
    return water_fraction
end

"Return river fraction"
function get_river_fraction(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand,
    river_location::Vector{Bool},
    area::Vector{Float},
)
    logging = false
    lens = lens_input_parameter(config, "river__width"; optional = false)
    river_width_2d = ncread(dataset, config, lens; type = Float, fill = 0, logging)
    river_width = river_width_2d[network.indices]

    lens = lens_input_parameter(config, "river__length"; optional = false)
    river_length_2d = ncread(dataset, config, lens; type = Float, fill = 0, logging)
    river_length = river_length_2d[network.indices]

    n = Int(length(river_location))
    river_fraction = fill(MISSING_VALUE, n)
    for i in 1:n
        river_fraction[i] = if river_location[i]
            min((river_length[i] * river_width[i]) / (area[i]), Float(1.0))
        else
            Float(0.0)
        end
    end
    return river_fraction
end

"Return cell lengths in x and y direction"
function get_cell_lengths(dataset::NCDataset, config::Config, network::NetworkLand)
    y_coords = read_y_axis(dataset)
    x_coords = read_x_axis(dataset)
    y = permutedims(repeat(y_coords; outer = (1, length(x_coords))))[network.indices]
    celllength = abs(mean(diff(x_coords)))

    cell_length_in_meter = get(config.model, "cell_length_in_meter__flag", false)::Bool
    x_length, y_length = cell_lengths(y, celllength, cell_length_in_meter)
    return x_length, y_length
end

"Return land surface slope"
function get_landsurface_slope(dataset::NCDataset, config::Config, network::NetworkLand)
    lens = lens_input_parameter(config, "land_surface__slope"; optional = false)
    slope = ncread(dataset, config, lens; sel = network.indices, type = Float)
    clamp!(slope, Float(0.00001), Float(Inf))
    return Vector{Float}(slope)
end

"Return river mask"
function river_mask(dataset::NCDataset, config::Config, network::NetworkLand)
    lens = lens_input(config, "river_location__mask"; optional = false)
    river_2d = ncread(dataset, config, lens; type = Bool, fill = false)
    river_location = river_2d[network.indices]
    return river_location
end

"Return waterbody (reservoir or lake) mask"
function waterbody_mask(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand;
    region = "location",
)
    do_reservoirs = get(config.model, "reservoir__flag", false)::Bool
    do_lakes = get(config.model, "lake__flag", false)::Bool
    waterbodies = fill(Int(0), length(network.indices))
    if do_reservoirs
        lens = lens_input(config, "reservoir_$(region)__count"; optional = false)
        reservoirs =
            ncread(dataset, config, lens; sel = network.indices, type = Float, fill = 0)
        waterbodies = waterbodies .+ reservoirs
    end
    if do_lakes
        lens = lens_input(config, "lake_$(region)__count"; optional = false)
        lakes = ncread(dataset, config, lens; sel = network.indices, type = Float, fill = 0)
        waterbodies = waterbodies .+ lakes
    end
    waterbodies = Vector{Bool}(waterbodies .> 0)
    return waterbodies
end

"Mask reservoir and lake coverage based on mask value `mask_value`"
function mask_waterbody_coverage!(
    mask::AbstractMatrix,
    config::Config,
    domain::Domain;
    mask_value = 0,
)
    do_reservoirs = get(config.model, "reservoir__flag", false)::Bool
    do_lakes = get(config.model, "lake__flag", false)::Bool
    if do_reservoirs
        inds_cov = collect(Iterators.flatten(domain.reservoir.network.indices_coverage))
        mask[inds_cov] .= mask_value
    end
    if do_lakes
        inds_cov = collect(Iterators.flatten(domain.lake.network.indices_coverage))
        mask[inds_cov] .= mask_value
    end
    return nothing
end

"Return indices of 1D land and river domains per allocation area number."
function get_allocation_area_indices(dataset::NCDataset, config::Config, domain::Domain)
    (; indices) = domain.land.network
    lens = lens_input_parameter(config, "land_water_allocation_area__count")
    logging = false
    areas = ncread(dataset, config, lens; sel = indices, defaults = 1, type = Int, logging)
    unique_areas = unique(areas)
    allocation_area_inds = Vector{Int}[]
    river_allocation_area_inds = Vector{Int}[]
    for a in unique_areas
        area_index = findall(x -> x == a, areas)
        push!(allocation_area_inds, area_index)
        area_river_index = findall(x -> x == a, areas[domain.river.network.land_indices])
        push!(river_allocation_area_inds, area_river_index)
    end
    return allocation_area_inds, river_allocation_area_inds
end