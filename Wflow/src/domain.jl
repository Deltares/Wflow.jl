"Struct to store (shared) land parameters"
@with_kw struct LandParameters
    # cell length x direction [m]
    x_length::Vector{Float64} = Float64[]
    # cell length y direction [m]
    y_length::Vector{Float64} = Float64[]
    # cell area [m²]
    area::Vector{Float64} = Float64[]
    # flow width [m]
    flow_width::Vector{Float64} = Float64[]
    # suface flow width [m]
    surface_flow_width::Vector{Float64} = Float64[]
    # flow length [m]
    flow_length::Vector{Float64} = Float64[]
    # flow fraction to river [-]
    flow_fraction_to_river::Vector{Float64} = Float64[]
    # slope [-]
    slope::Vector{Float64} = Float64[]
    # reservoir location [-]
    reservoir_outlet::Vector{Bool} = Bool[]
    # reservoir coverage [-]
    reservoir_coverage::Vector{Bool} = Bool[]
    # river location [-]
    river_location::Vector{Bool} = Bool[]
    # fraction of river [-]
    river_fraction::Vector{Float64} = Float64[]
    # fraction of open water (excluding rivers) [-]
    water_fraction::Vector{Float64} = Float64[]
end

"Struct to store (shared) river parameters"
@with_kw struct RiverParameters
    # river flow width [m]
    flow_width::Vector{Float64} = Float64[]
    # river flow length [m]
    flow_length::Vector{Float64} = Float64[]
    # slope [-]
    slope::Vector{Float64} = Float64[]
    # reservoir location
    reservoir_outlet::Vector{Bool} = Bool[]
    # reservoir coverage [-]
    reservoir_coverage::Vector{Bool} = Bool[]
    # grid cell area [m²]
    cell_area::Vector{Float64} = Float64[]
end

@kwdef struct DomainLand
    network::NetworkLand = NetworkLand()
    parameters::LandParameters = LandParameters()
end

@kwdef struct DomainRiver
    network::NetworkRiver = NetworkRiver()
    parameters::RiverParameters = RiverParameters()
end

@kwdef struct DomainReservoir
    network::NetworkReservoir = NetworkReservoir()
end

@kwdef struct DomainDrain
    network::NetworkDrain = NetworkDrain()
end

"""
Struct for storing information about different model domains 'land`, `river`, `reservoir`
and `drain` (`Drainage` boundary condition of `GroundwaterFlow`). It holds network
information for each domain like active indices in the 2D model domain and connectivity
information for flow routing. The `land` and `river` domains contain shared parameters for
each domain that can used by different model components.
"""
@kwdef struct Domain
    land::DomainLand = DomainLand()
    river::DomainRiver = DomainRiver()
    reservoir::DomainReservoir = DomainReservoir()
    drain::DomainDrain = DomainDrain()
end

"Initialize `Domain` for model types `sbm` and `sbm_gwf`"
function Domain(dataset::NCDataset, config::Config, ::Union{SbmModel, SbmGwfModel})
    (; land_routing, river_routing) = config.model

    network_land = NetworkLand(dataset, config)
    if land_routing == RoutingType.kinematic_wave ||
       subsurface_routing(config) == RoutingType.kinematic_wave
        network_land = network_subdomains(config, network_land)
    end

    network_river =
        NetworkRiver(dataset, config, network_land; do_pits = config.model.pit__flag)
    if river_routing == RoutingType.kinematic_wave
        network_river = network_subdomains(config, network_river)
    end

    pits = zeros(Bool, network_land.modelsize)
    nriv = length(network_river.indices)
    if config.model.reservoir__flag
        network_reservoir, inds_reservoir_map2river =
            NetworkReservoir(dataset, config, network_river)
        pits[network_reservoir.indices_outlet] .= true
    else
        network_reservoir = NetworkReservoir()
        inds_reservoir_map2river = fill(0, nriv)
    end
    @reset network_river.reservoir_indices = inds_reservoir_map2river

    if river_routing == RoutingType.kinematic_wave
        @reset network_river.upstream_nodes =
            filter_upsteam_nodes(network_river.graph, pits[network_river.indices])
    elseif river_routing == RoutingType.local_inertial
        nodes_at_edge, index_pit = NodesAtEdge(network_river)
        @reset network_river.nodes_at_edge = nodes_at_edge
        @reset network_river.pit_indices = network_river.indices[index_pit]
        @reset network_river.edges_at_node = EdgesAtNode(network_river)
    end

    if land_routing == RoutingType.kinematic_wave ||
       subsurface_routing(config) == RoutingType.kinematic_wave
        @reset network_land.upstream_nodes =
            filter_upsteam_nodes(network_land.graph, pits[network_land.indices])
    end
    if land_routing == RoutingType.local_inertial
        @reset network_land.edge_indices = EdgeConnectivity(network_land)
        @reset network_land.river_indices =
            network_river.reverse_indices[network_land.indices]
    end

    domain_land = DomainLand(; network = network_land)
    domain_river = DomainRiver(; network = network_river)

    domain = Domain(;
        land = domain_land,
        river = domain_river,
        reservoir = DomainReservoir(; network = network_reservoir),
    )

    land_params, river_params = initialize_shared_parameters(dataset, config, domain)
    @reset domain.land.parameters = land_params
    @reset domain.river.parameters = river_params

    if config.model.drain__flag
        (; indices) = domain.land.network
        (; surface_flow_width) = domain.land.parameters
        @reset domain.drain.network =
            NetworkDrain(dataset, config, indices, surface_flow_width)
    end

    if do_water_demand(config)
        allocation_area_indices, river_allocation_indices =
            get_allocation_area_indices(dataset, config, domain)
        @reset domain.land.network.allocation_area_indices = allocation_area_indices
        @reset domain.river.network.allocation_area_indices = river_allocation_indices

        mask = copy(domain.river.network.reverse_indices)
        mask_reservoir_coverage!(mask, config, domain)
        @reset domain.land.network.river_inds_excl_reservoir =
            mask[domain.land.network.indices]
    end

    if nthreads() > 1
        min_streamorder_land = config.model.land_streamorder__min_count
        min_streamorder_river = config.model.river_streamorder__min_count
        if river_routing == RoutingType.kinematic_wave
            @info "Parallel execution of kinematic wave." min_streamorder_land min_streamorder_river
        elseif land_routing == RoutingType.kinematic_wave ||
               subsurface_routing(config) == RoutingType.kinematic_wave
            @info "Parallel execution of kinematic wave." * min_streamorder_land
        end
    end

    return domain
end

"Initialize `Domain` for model type `sediment`"
function Domain(dataset::NCDataset, config::Config, ::SedimentModel)
    network_land = NetworkLand(dataset, config)
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
    reservoir_outlet = reservoir_mask(dataset, config, network)
    reservoir_coverage = reservoir_mask(dataset, config, network; region = "area")
    river_location = river_mask(dataset, config, network)

    land_parameters = LandParameters(;
        area,
        flow_width,
        slope,
        reservoir_outlet,
        reservoir_coverage,
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

    land_area = @. (1.0 - river_fraction) * area
    surface_flow_width =
        map(get_surface_width, flow_width, flow_length, land_area, river_location)

    flow_fraction_to_river = get_flow_fraction_to_river(
        network.graph,
        network.local_drain_direction,
        land_indices,
        slope,
    )

    reservoir_outlet = reservoir_mask(dataset, config, network)
    reservoir_coverage = reservoir_mask(dataset, config, network; region = "area")

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
        reservoir_outlet,
        reservoir_coverage,
        river_fraction,
        water_fraction,
    )
    return land_parameters
end

"Initialize (shared) river parameters"
function RiverParameters(dataset::NCDataset, config::Config, network::NetworkRiver)
    (; indices) = network
    flow_length = ncread(
        dataset,
        config,
        "river__length";
        optional = false,
        sel = indices,
        type = Float64,
    )
    minimum(flow_length) > 0 || error("river length must be positive on river cells")

    flow_width = ncread(
        dataset,
        config,
        "river__width";
        optional = false,
        sel = indices,
        type = Float64,
    )
    minimum(flow_width) > 0 || error("river width must be positive on river cells")

    slope = ncread(
        dataset,
        config,
        "river__slope";
        optional = false,
        sel = indices,
        type = Float64,
    )
    clamp!(slope, 0.00001, Inf)

    river_parameters = RiverParameters(; flow_width, flow_length, slope)
    return river_parameters
end

"Initialize (shared) parameters for the river and land domains"
function initialize_shared_parameters(dataset::NCDataset, config::Config, domain::Domain)
    land_params = LandParameters(dataset, config, domain)
    river_params = RiverParameters(dataset, config, domain.river.network)

    @reset river_params.cell_area = land_params.area[domain.river.network.land_indices]
    @reset river_params.reservoir_coverage =
        land_params.reservoir_coverage[domain.river.network.land_indices]
    @reset river_params.reservoir_outlet =
        land_params.reservoir_outlet[domain.river.network.land_indices]

    return land_params, river_params
end

"Return open water fraction (excluding rivers)"
function get_water_fraction(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand,
    river_fraction::Vector{Float64},
)
    water_fraction = ncread(
        dataset,
        config,
        "land_water_covered__area_fraction";
        sel = network.indices,
        defaults = 0.0,
        type = Float64,
    )
    water_fraction = max.(water_fraction .- river_fraction, 0.0)
    return water_fraction
end

"Return river fraction"
function get_river_fraction(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand,
    river_location::Vector{Bool},
    area::Vector{Float64},
)
    river_width_2d = ncread(
        dataset,
        config,
        "river__width";
        optional = false,
        type = Float64,
        fill = 0,
        logging = false,
    )
    river_width = river_width_2d[network.indices]

    river_length_2d = ncread(
        dataset,
        config,
        "river__length";
        optional = false,
        type = Float64,
        fill = 0,
        logging = false,
    )
    river_length = river_length_2d[network.indices]

    n = length(river_location)
    river_fraction = fill(MISSING_VALUE, n)
    for i in 1:n
        river_fraction[i] = if river_location[i]
            min((river_length[i] * river_width[i]) / (area[i]), 1.0)
        else
            0.0
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

    x_length, y_length =
        cell_lengths(y, celllength, config.model.cell_length_in_meter__flag)
    return x_length, y_length
end

"Return land surface slope"
function get_landsurface_slope(dataset::NCDataset, config::Config, network::NetworkLand)
    slope = ncread(
        dataset,
        config,
        "land_surface__slope";
        optional = false,
        sel = network.indices,
        type = Float64,
    )
    clamp!(slope, 0.00001, Inf)
    return slope
end

"Return river mask"
function river_mask(dataset::NCDataset, config::Config, network::NetworkLand)
    river_2d = ncread(
        dataset,
        config,
        "river_location__mask";
        optional = false,
        type = Bool,
        fill = false,
    )
    river_location = river_2d[network.indices]
    return river_location
end

"Return reservoir mask"
function reservoir_mask(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand;
    region::String = "location",
)
    reservoirs = fill(0, length(network.indices))
    if config.model.reservoir__flag
        reservoirs = ncread(
            dataset,
            config,
            "reservoir_$(region)__count";
            optional = false,
            sel = network.indices,
            type = Float64,
            fill = 0,
        )
    end
    reservoirs = Vector{Bool}(reservoirs .> 0)
    return reservoirs
end

"Mask reservoir coverage based on mask value `mask_value`"
function mask_reservoir_coverage!(mask::AbstractMatrix, config::Config, domain::Domain)
    if config.model.reservoir__flag
        inds_cov = collect(Iterators.flatten(domain.reservoir.network.indices_coverage))
        mask[inds_cov] .= 0
    end
    return nothing
end

"Return indices of 1D land and river domains per allocation area number."
function get_allocation_area_indices(dataset::NCDataset, config::Config, domain::Domain)
    (; indices) = domain.land.network
    areas = ncread(
        dataset,
        config,
        "land_water_allocation_area__count";
        sel = indices,
        defaults = 1,
        type = Int,
        logging = false,
    )
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
