abstract type AbstractFloodPlainModel end

"Struct for storing local inertial river flow model parameters"
@with_kw struct LocalInertialRiverFlowParameters
    n_cells::Int                      # number of cells [-]
    n_edges::Int                      # number of edges [-]
    active_n::Vector{Int}                   # active nodes [-]
    active_e::Vector{Int}                   # active edges [-]
    froude_limit::Bool                      # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
    h_thresh::Float64                       # depth threshold for calculating flow [m]
    zb::Vector{Float64}                     # river bed elevation [m]
    zb_max::Vector{Float64}                 # maximum channel bed elevation [m]
    bankfull_storage::Vector{Float64}       # bankfull storage [m³]
    bankfull_depth::Vector{Float64}         # bankfull depth [m]
    mannings_n_sq::Vector{Float64}          # Manning's roughness squared at edge [(s m-1/3)2]
    mannings_n::Vector{Float64}             # Manning's roughness [s m-1/3] at node
    flow_length_at_edge::Vector{Float64}    # flow (river) length at edge [m]
    flow_width_at_edge::Vector{Float64}     # flow (river) width at edge [m]
end

"Initialize local inertial river flow model parameters"
function LocalInertialRiverFlowParameters(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
)
    alpha = config.model.river_local_inertial_flow__alpha_coefficient # stability coefficient for model time step (0.2-0.7)
    waterdepth_threshold = config.model.river_water_flow_threshold__depth # depth threshold for flow at edge
    froude_limit = config.model.river_water_flow__froude_limit_flag # limit flow to subcritical according to Froude number
    floodplain_1d = config.model.floodplain_1d__flag

    @info "Local inertial approach is used for river flow." alpha waterdepth_threshold froude_limit floodplain_1d

    (; pit_indices, river_indices_2d, graph, local_drain_direction, nodes_at_edge) =
        domain.network
    (; flow_width, flow_length, reservoir_outlet) = domain.parameters

    riverlength_bc = ncread(
        dataset,
        config,
        "model_boundary_condition_river__length",
        Routing;
        sel = pit_indices,
    )
    bankfull_elevation_2d = ncread(dataset, config, "river_bank_water__elevation", Routing)
    bankfull_depth_2d = ncread(dataset, config, "river_bank_water__depth", Routing)
    bankfull_depth = bankfull_depth_2d[river_indices_2d]
    zb = bankfull_elevation_2d[river_indices_2d] - bankfull_depth # river bed elevation

    bankfull_storage = bankfull_depth .* flow_width .* flow_length
    mannings_n = ncread(
        dataset,
        config,
        "river_water_flow__manning_n_parameter",
        Routing;
        sel = river_indices_2d,
    )

    n_cells = length(river_indices_2d)
    index_pit = findall(x -> x == 5, local_drain_direction)
    # set ghost points for boundary condition (downstream river outlet): river width, bed
    # elevation, manning n is copied from the upstream cell.
    append!(flow_length, riverlength_bc)
    append!(zb, zb[index_pit])
    append!(flow_width, flow_width[index_pit])
    append!(mannings_n, mannings_n[index_pit])
    append!(bankfull_depth, bankfull_depth[index_pit])

    # determine z, width, length and manning's n at edges
    n_edges = ne(graph)
    zb_max = zeros(n_edges)
    width_at_edge = zeros(n_edges)
    length_at_edge = zeros(n_edges)
    mannings_n_sq = zeros(n_edges)
    for river_edge_idx in 1:n_edges
        src_node = nodes_at_edge.src[river_edge_idx]
        dst_node = nodes_at_edge.dst[river_edge_idx]
        zb_max[river_edge_idx] = max(zb[src_node], zb[dst_node])
        width_at_edge[river_edge_idx] = min(flow_width[src_node], flow_width[dst_node])
        length_at_edge[river_edge_idx] =
            0.5 * (flow_length[dst_node] + flow_length[src_node])
        mannings_n_i =
            (
                mannings_n[dst_node] * flow_length[dst_node] +
                mannings_n[src_node] * flow_length[src_node]
            ) / (flow_length[dst_node] + flow_length[src_node])
        mannings_n_sq[river_edge_idx] = mannings_n_i * mannings_n_i
    end
    active_index = findall(x -> x == 0, reservoir_outlet)

    parameters = LocalInertialRiverFlowParameters(;
        n_cells,
        n_edges,
        active_n = active_index,
        active_e = active_index,
        froude_limit,
        h_thresh = waterdepth_threshold,
        zb,
        zb_max,
        bankfull_storage,
        bankfull_depth,
        mannings_n,
        mannings_n_sq,
        flow_length_at_edge = length_at_edge,
        flow_width_at_edge = width_at_edge,
    )
    return parameters
end

"Struct for storing local inertial river flow model variables"
@with_kw struct LocalInertialRiverFlowVariables
    n_cells::Int
    n_edges::Int
    q::Vector{Float64} = zeros(n_edges)   # river discharge at edge (subgrid channel) [m³ s⁻¹]
    q0::Vector{Float64} = zeros(n_edges) # river discharge at edge (subgrid channel) at previous time step [m³ s⁻¹]
    q_av::Vector{Float64}                     # average river channel (+ floodplain) discharge at edge [m³ s⁻¹] (model timestep Δt)
    q_channel_av::Vector{Float64}             # average river channel discharge at edge [m³ s⁻¹] (for model timestep Δt)
    h::Vector{Float64}                        # water depth [m]
    zs_max::Vector{Float64} = zeros(n_edges)  # maximum water elevation at edge [m]
    zs_src::Vector{Float64} = zeros(n_edges)  # water elevation of source node of edge [m]
    zs_dst::Vector{Float64} = zeros(n_edges)  # water elevation of downstream node of edge [m]
    hf::Vector{Float64} = zeros(n_edges)       # water depth at edge [m]
    a::Vector{Float64} = zeros(n_edges)        # flow area at edge [m²]
    r::Vector{Float64} = zeros(n_edges)        # wetted perimeter at edge [m]
    storage::Vector{Float64} = zeros(n_cells)        # river storage [m³]
    error::Vector{Float64} = zeros(n_cells)          # error storage [m³]
end

"Initialize shallow water river flow model variables"
function LocalInertialRiverFlowVariables(
    dataset::NCDataset,
    config::Config,
    network::NetworkRiver,
)
    (; pit_indices, river_indices_2d, graph) = network

    riverdepth_bc = ncread(
        dataset,
        config,
        "model_boundary_condition_river_bank_water__depth",
        Routing;
        sel = pit_indices,
    )

    n_cells = length(river_indices_2d)
    n_edges = ne(graph)
    # set river depth h to zero (including reservoir locations)
    h = zeros(n_cells)
    q_av = zeros(n_edges)
    # set ghost points for boundary condition (downstream river outlet): river depth `h`
    append!(h, riverdepth_bc)
    variables = LocalInertialRiverFlowVariables(;
        n_cells,
        n_edges,
        q_av,
        q_channel_av = config.model.floodplain_1d__flag ? zeros(n_edges) : q_av,
        h,
    )
    return variables
end

"Shallow water river flow model using the local inertial method"
@with_kw struct LocalInertialRiverFlowModel{
    R <: RiverFlowBC,
    F <: Union{AbstractFloodPlainModel, Nothing},
    A <: AbstractAllocationModel,
} <: AbstractRiverFlowModel
    timestepping::TimeStepping
    boundary_conditions::R
    parameters::LocalInertialRiverFlowParameters
    variables::LocalInertialRiverFlowVariables
    floodplain::F                                       # Floodplain (1D) schematization
    allocation::A                                       # Water allocation
end

"Initialize shallow water river flow model `LocalInertialRiverFlowModel`"
function LocalInertialRiverFlowModel(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    reservoir_model::Union{ReservoirModel, Nothing},
)
    # The local inertial approach makes use of a staggered grid (Bates et al. (2010)),
    # with nodes and edges. This information is extracted from the directed graph of the
    # river. Discharge q is calculated at edges between nodes and mapped to the source
    # nodes for gridded output (index of edge is equal to source node index, e.g.:
    # Edge 1 => 5
    # Edge 2 => 1
    # Edge 3 => 2
    # Edge 4 => 9
    # ⋮ )

    # The following boundary conditions can be set at ghost nodes, downstream of river
    # outlets (pits): river length and river depth
    alpha_coefficient = config.model.river_local_inertial_flow__alpha_coefficient # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; alpha_coefficient)

    parameters = LocalInertialRiverFlowParameters(dataset, config, domain)
    variables = LocalInertialRiverFlowVariables(dataset, config, domain.network)
    boundary_conditions = RiverFlowBC(dataset, config, domain.network, reservoir_model)

    if config.model.floodplain_1d__flag
        zb_floodplain = parameters.zb .+ parameters.bankfull_depth
        floodplain = FloodPlainModel(dataset, config, domain, zb_floodplain)
    else
        floodplain = nothing
    end

    n_cells = length(domain.network.river_indices_2d)
    river_flow = LocalInertialRiverFlowModel(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain,
        allocation = do_water_demand(config) ? AllocationRiverModel(n_cells) :
                     NoAllocationRiverModel(n_cells),
    )
    return river_flow
end

"Return the upstream inflow for a reservoir in `LocalInertialRiverFlowModel`"
function get_inflow_reservoir(
    river_flow_model::LocalInertialRiverFlowModel,
    src_edge::Vector{Int},
)
    q_in = sum_at(river_flow_model.variables.q, src_edge)
    if !isnothing(river_flow_model.floodplain)
        q_in += sum_at(river_flow_model.floodplain.variables.q, src_edge)
    end
    return q_in
end

# For local inertial river routing, `to_river` is included, as reservoir cells are excluded
# (boundary condition).
get_inflow_reservoir(
    ::LocalInertialRiverFlowModel,
    overland_flow_model::KinWaveOverlandFlowModel,
    inds::Vector{Int},
) = overland_flow_model.variables.q_av[inds] .+ overland_flow_model.variables.to_river[inds]

get_inflow_reservoir(
    ::LocalInertialRiverFlowModel,
    subsurface_flow_model::LateralSSFModel,
    inds::Vector{Int},
) =
    (
        subsurface_flow_model.variables.q_av[inds] .+
        subsurface_flow_model.variables.to_river[inds]
    ) ./ tosecond(BASETIMESTEP)

"""
Update river channel flow for the local inertial river flow model.
"""
function update_river_channel_flow!(
    river_flow_model::LocalInertialRiverFlowModel,
    domain::DomainRiver,
    dt::Float64,
)
    (; nodes_at_edge) = domain.network
    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters

    river_v.q0 .= river_v.q
    if !isnothing(river_flow_model.floodplain)
        river_flow_model.floodplain.variables.q0 .= river_flow_model.floodplain.variables.q
    end
    @batch per = thread minbatch = 1000 for river_edge_idx in river_p.active_e
        cell_idx_src = nodes_at_edge.src[river_edge_idx]
        cell_idx_dst = nodes_at_edge.dst[river_edge_idx]
        river_v.zs_src[river_edge_idx] =
            river_p.zb[cell_idx_src] + river_v.h[cell_idx_src]
        river_v.zs_dst[river_edge_idx] =
            river_p.zb[cell_idx_dst] + river_v.h[cell_idx_dst]

        river_v.zs_max[river_edge_idx] =
            max(river_v.zs_src[river_edge_idx], river_v.zs_dst[river_edge_idx])
        river_v.hf[river_edge_idx] =
            (river_v.zs_max[river_edge_idx] - river_p.zb_max[river_edge_idx])

        river_v.a[river_edge_idx] =
            river_p.flow_width_at_edge[river_edge_idx] * river_v.hf[river_edge_idx] # flow area (rectangular channel)
        river_v.r[river_edge_idx] =
            river_v.a[river_edge_idx] /
            (river_p.flow_width_at_edge[river_edge_idx] + 2.0 * river_v.hf[river_edge_idx]) # hydraulic radius (rectangular channel)

        river_v.q[river_edge_idx] = ifelse(
            river_v.hf[river_edge_idx] > river_p.h_thresh,
            local_inertial_flow(
                river_v.q0[river_edge_idx],
                river_v.zs_src[river_edge_idx],
                river_v.zs_dst[river_edge_idx],
                river_v.hf[river_edge_idx],
                river_v.a[river_edge_idx],
                river_v.r[river_edge_idx],
                river_p.flow_length_at_edge[river_edge_idx],
                river_p.mannings_n_sq[river_edge_idx],
                river_p.froude_limit,
                dt,
            ),
            0.0,
        )

        # limit q in case water is not available
        river_v.q[river_edge_idx] = ifelse(
            river_v.h[cell_idx_src] <= 0.0,
            min(river_v.q[river_edge_idx], 0.0),
            river_v.q[river_edge_idx],
        )
        river_v.q[river_edge_idx] = ifelse(
            river_v.h[cell_idx_dst] <= 0.0,
            max(river_v.q[river_edge_idx], 0.0),
            river_v.q[river_edge_idx],
        )
        # average river discharge (here accumulated for model timestep Δt)
        river_v.q_av[river_edge_idx] += river_v.q[river_edge_idx] * dt
    end
    return nothing
end

"""
Update floodplain flow for the local inertial river flow model.
"""
function update_floodplain_flow!(
    river_flow_model::LocalInertialRiverFlowModel{R, F},
    domain::DomainRiver,
    dt::Float64,
) where {R, F <: AbstractFloodPlainModel}
    (; nodes_at_edge) = domain.network
    (; flow_width) = domain.parameters

    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters
    floodplain_p = river_flow_model.floodplain.parameters
    floodplain_v = river_flow_model.floodplain.variables

    @batch per = thread minbatch = 1000 for river_edge_idx in 1:(river_p.n_edges)
        floodplain_v.hf[river_edge_idx] =
            max(river_v.zs_max[river_edge_idx] - floodplain_p.zb_max[river_edge_idx], 0.0)
    end

    n_active_flood_edges = 0
    @inbounds for river_edge_idx in river_p.active_e
        @inbounds if river_v.hf[river_edge_idx] > river_p.h_thresh
            n_active_flood_edges += 1
            floodplain_v.hf_index[n_active_flood_edges] = river_edge_idx
        else
            floodplain_v.q[river_edge_idx] = 0.0
        end
    end

    get_area(river_edge_idx, lower_idx, upper_idx, cell_idx) = flow_area(
        floodplain_p.profile.width[upper_idx, cell_idx],
        floodplain_p.profile.a[lower_idx, cell_idx],
        floodplain_p.profile.depth[lower_idx],
        floodplain_v.hf[river_edge_idx],
    )

    get_wetted_perimeter(river_edge_idx, lower_idx, cell_idx) = wetted_perimeter(
        floodplain_p.profile.p[lower_idx, cell_idx],
        floodplain_p.profile.depth[lower_idx],
        floodplain_v.hf[river_edge_idx],
    )

    @batch per = thread minbatch = 1000 for flood_edge_idx in 1:n_active_flood_edges
        river_edge_idx = floodplain_v.hf_index[flood_edge_idx]
        cell_idx_src = nodes_at_edge.src[river_edge_idx]
        cell_idx_dst = nodes_at_edge.dst[river_edge_idx]

        depth_count = 0
        for depth_idx in eachindex(floodplain_p.profile.depth)
            depth_count +=
                1 *
                (floodplain_p.profile.depth[depth_idx] <= floodplain_v.hf[river_edge_idx])
        end
        lower_idx = max(depth_count, 1)
        upper_idx = ifelse(
            lower_idx == length(floodplain_p.profile.depth),
            lower_idx,
            lower_idx + 1,
        )

        a_src = get_area(river_edge_idx, lower_idx, upper_idx, cell_idx_src)
        a_src = max(
            a_src - (floodplain_v.hf[river_edge_idx] * flow_width[cell_idx_src]),
            0.0,
        )

        a_dst = get_area(river_edge_idx, lower_idx, upper_idx, cell_idx_dst)
        a_dst = max(
            a_dst - (floodplain_v.hf[river_edge_idx] * flow_width[cell_idx_dst]),
            0.0,
        )

        floodplain_v.a[river_edge_idx] = min(a_src, a_dst)

        floodplain_v.r[river_edge_idx] = if a_src < a_dst
            a_src / get_wetted_perimeter(river_edge_idx, lower_idx, cell_idx_src)
        else
            a_dst / get_wetted_perimeter(river_edge_idx, lower_idx, cell_idx_dst)
        end

        floodplain_v.q[river_edge_idx] = if floodplain_v.a[river_edge_idx] > 1.0e-05
            local_inertial_flow(
                floodplain_v.q0[river_edge_idx],
                river_v.zs_src[river_edge_idx],
                river_v.zs_dst[river_edge_idx],
                floodplain_v.hf[river_edge_idx],
                floodplain_v.a[river_edge_idx],
                floodplain_v.r[river_edge_idx],
                river_p.flow_length_at_edge[river_edge_idx],
                floodplain_p.mannings_n_sq[river_edge_idx],
                river_p.froude_limit,
                dt,
            )
        else
            0.0
        end

        # limit floodplain q in case water is not available
        if floodplain_v.h[cell_idx_src] <= 0.0
            floodplain_v.q[river_edge_idx] = min(floodplain_v.q[river_edge_idx], 0.0)
        end

        if floodplain_v.h[cell_idx_dst] <= 0.0
            floodplain_v.q[river_edge_idx] = max(floodplain_v.q[river_edge_idx], 0.0)
        end

        if floodplain_v.q[river_edge_idx] * river_v.q[river_edge_idx] < 0.0
            floodplain_v.q[river_edge_idx] = 0.0
        end

        # average floodplain discharge (here accumulated for model timestep Δt)
        floodplain_v.q_av[river_edge_idx] += floodplain_v.q[river_edge_idx] * dt
    end
    return nothing
end

update_floodplain_flow!(
    model::LocalInertialRiverFlowModel{R, F},
    domain::DomainRiver,
    dt::Float64,
) where {R, F <: Nothing} = nothing

"""
Update reservoir boundary conditions for the local inertial river flow model.
"""
function update_bc_reservoir_model!(
    reservoir_model::ReservoirModel,
    river_flow_model::LocalInertialRiverFlowModel,
    domain::Domain,
    dt::Float64,
    dt_forcing::Float64,
)
    (; edges_at_node) = domain.river.network
    inds_reservoir = domain.reservoir.network.river_cell_indices_containing_reservoir

    river_v = river_flow_model.variables
    res_bc = reservoir_model.boundary_conditions

    for (reservoir_idx, cell_idx) in enumerate(inds_reservoir)
        q_in = get_inflow_reservoir(river_flow_model, edges_at_node.src[cell_idx])
        # If external_inflow < 0, abstraction is limited
        if res_bc.external_inflow[reservoir_idx] < 0.0
            _abstraction = min(
                -res_bc.external_inflow[reservoir_idx],
                (reservoir_model.variables.storage[reservoir_idx] / dt) * 0.98,
            )
            res_bc.actual_external_abstraction_av[reservoir_idx] += _abstraction * dt
            _inflow = -_abstraction
        else
            _inflow = res_bc.external_inflow[reservoir_idx]
        end
        net_inflow =
            q_in +
            res_bc.inflow_overland[reservoir_idx] +
            res_bc.inflow_subsurface[reservoir_idx] +
            _inflow
        update_reservoir_model!(reservoir_model, reservoir_idx, net_inflow, dt, dt_forcing)
        river_v.q[cell_idx] = reservoir_model.variables.outflow[reservoir_idx]
        # average river discharge (here accumulated for model timestep Δt)
        river_v.q_av[cell_idx] += river_v.q[cell_idx] * dt
    end
    return nothing
end
update_bc_reservoir_model!(
    reservoir_model::Nothing,
    river_flow_model::LocalInertialRiverFlowModel,
    domain::Domain,
    dt::Float64,
    dt_forcing::Float64,
) = nothing

"""
Update floodplain water depth and storage.
"""
function update_water_depth_and_storage!(
    floodplain_model::AbstractFloodPlainModel,
    river_flow_model::LocalInertialRiverFlowModel,
    domain::DomainRiver,
    dt::Float64,
)
    (; edges_at_node) = domain.network
    (; flow_length, flow_width) = domain.parameters

    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters
    floodplain_v = floodplain_model.variables
    floodplain_p = floodplain_model.parameters

    @batch per = thread minbatch = 1000 for cell_idx in river_p.active_n
        q_src = sum_at(floodplain_v.q, edges_at_node.src[cell_idx])
        q_dst = sum_at(floodplain_v.q, edges_at_node.dst[cell_idx])
        floodplain_v.storage[cell_idx] =
            floodplain_v.storage[cell_idx] + (q_src - q_dst) * dt
        if floodplain_v.storage[cell_idx] < 0.0
            floodplain_v.error[cell_idx] += abs(floodplain_v.storage[cell_idx])
            floodplain_v.storage[cell_idx] = 0.0
        end
        storage_total =
            river_v.storage[cell_idx] + floodplain_v.storage[cell_idx]
        if storage_total > river_p.bankfull_storage[cell_idx]
            flood_storage = storage_total - river_p.bankfull_storage[cell_idx]
            h = flood_depth(
                floodplain_p.profile,
                flood_storage,
                flow_length[cell_idx],
                cell_idx,
            )
            river_v.h[cell_idx] = river_p.bankfull_depth[cell_idx] + h
            river_v.storage[cell_idx] =
                river_v.h[cell_idx] *
                flow_width[cell_idx] *
                flow_length[cell_idx]
            floodplain_v.storage[cell_idx] =
                max(storage_total - river_v.storage[cell_idx], 0.0)
            floodplain_v.h[cell_idx] =
                floodplain_v.storage[cell_idx] > 0.0 ? h : 0.0
        else
            river_v.h[cell_idx] =
                storage_total / (flow_length[cell_idx] * flow_width[cell_idx])
            river_v.storage[cell_idx] = storage_total
            floodplain_v.h[cell_idx] = 0.0
            floodplain_v.storage[cell_idx] = 0.0
        end
    end
    return nothing
end

"""
Update floodplain water depth and storage (no-op for Nothing floodplain).
"""
update_water_depth_and_storage!(
    ::Nothing,
    ::LocalInertialRiverFlowModel,
    ::DomainRiver,
    ::Float64,
) = nothing

"""
Update water depth and storage for river.
"""
function update_water_depth_and_storage!(
    river_flow_model::LocalInertialRiverFlowModel,
    domain::DomainRiver,
    dt::Float64,
)
    (; edges_at_node) = domain.network
    (; flow_length, flow_width) = domain.parameters
    (; inwater, abstraction, external_inflow, actual_external_abstraction_av) =
        river_flow_model.boundary_conditions

    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters

    @batch per = thread minbatch = 1000 for cell_idx in river_p.active_n
        q_src = sum_at(river_v.q, edges_at_node.src[cell_idx])
        q_dst = sum_at(river_v.q, edges_at_node.dst[cell_idx])
        # internal abstraction (water demand) is limited by river storage and negative
        # external inflow as part of water allocation computations.
        river_v.storage[cell_idx] =
            river_v.storage[cell_idx] +
            (q_src - q_dst + inwater[cell_idx] - abstraction[cell_idx]) * dt

        if river_v.storage[cell_idx] < 0.0
            river_v.error[cell_idx] =
                river_v.error[cell_idx] + abs(river_v.storage[cell_idx])
            river_v.storage[cell_idx] = 0.0 # set storage to zero
        end
        # limit negative external inflow
        if external_inflow[cell_idx] < 0.0
            _abstraction = min(
                -external_inflow[cell_idx],
                river_v.storage[cell_idx] / dt * 0.80,
            )
            actual_external_abstraction_av[cell_idx] += _abstraction * dt
            _inflow = -_abstraction
        else
            _inflow = external_inflow[cell_idx]
        end
        river_v.storage[cell_idx] += _inflow * dt # add external inflow
        river_v.h[cell_idx] =
            river_v.storage[cell_idx] /
            (flow_length[cell_idx] * flow_width[cell_idx])
    end
    return nothing
end

"Update local inertial river flow model `LocalInertialRiverFlowModel` for a single timestep"
function local_inertial_river_update!(
    river_flow_model::LocalInertialRiverFlowModel,
    domain::Domain,
    dt::Float64,
    dt_forcing::Float64,
    update_h::Bool,
)
    # Update river channel flow
    update_river_channel_flow!(river_flow_model, domain.river, dt)

    # Update floodplain flow if present
    update_floodplain_flow!(river_flow_model, domain.river, dt)

    # Handle reservoir boundary conditions
    update_bc_reservoir_model!(
        river_flow_model.boundary_conditions.reservoir,
        river_flow_model,
        domain,
        dt,
        dt_forcing,
    )

    # Update water depth and storage if requested
    if update_h
        update_water_depth_and_storage!(river_flow_model, domain.river, dt)
        update_water_depth_and_storage!(
            river_flow_model.floodplain,
            river_flow_model,
            domain.river,
            dt,
        )
    end

    return nothing
end

"""
Update local inertial river flow model `LocalInertialRiverFlow` for a single timestep `dt`. An adaptive
timestepping method is used (computing a sub timestep `dt_s`).
"""
function update_river_flow_model!(
    river_flow_model::LocalInertialRiverFlowModel,
    domain::Domain,
    clock::Clock;
    update_h = true,
)
    (; reservoir) = river_flow_model.boundary_conditions
    (; flow_length) = domain.river.parameters

    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)

    if !isnothing(river_flow_model.floodplain)
        river_flow_model.floodplain.variables.q_av .= 0.0
    end
    set_flow_vars!(river_flow_model)

    dt = tosecond(clock.dt)
    t = 0.0
    while t < dt
        dt_s = stable_timestep(river_flow_model, flow_length)
        dt_s = check_timestepsize(dt_s, t, dt)
        local_inertial_river_update!(river_flow_model, domain, dt_s, dt, update_h)
        t += dt_s
    end
    average_flow_vars!(river_flow_model, dt)
    average_reservoir_vars!(reservoir, dt)

    if !isnothing(river_flow_model.floodplain)
        river_flow_model.floodplain.variables.q_av ./= dt
        river_flow_model.variables.q_channel_av .= river_flow_model.variables.q_av
        river_flow_model.variables.q_av .=
            river_flow_model.variables.q_channel_av .+
            river_flow_model.floodplain.variables.q_av
    end

    return nothing
end

"Struct to store local inertial overland flow model variables"
@with_kw struct LocalInertialOverlandFlowVariables
    n_cells::Int
    # flow in y direction at edge at previous time step [m³ s⁻¹]
    qy0::Vector{Float64} = zeros(n_cells + 1)
    # flow in x direction at edge at previous time step [m³ s⁻¹]
    qx0::Vector{Float64} = zeros(n_cells + 1)
    # flow in x direction at edge [m³ s⁻¹]
    qx::Vector{Float64} = zeros(n_cells + 1)
    # average flow in x direction at edge [m³ s⁻¹] for model timestep Δt
    qx_av::Vector{Float64} = zeros(n_cells + 1)
    # flow in y direction at edge [m³ s⁻¹]
    qy::Vector{Float64} = zeros(n_cells + 1)
    # average flow in y direction at edge [m³ s⁻¹] for model timestep Δt
    qy_av::Vector{Float64} = zeros(n_cells + 1)
    # total storage of cell [m³] (including river storage for river cells)
    storage::Vector{Float64} = zeros(n_cells)
    # error storage [m³]
    error::Vector{Float64} = zeros(n_cells)
    # water depth of cell [m] (for river cells the reference is the river bed elevation `zb`)
    h::Vector{Float64} = zeros(n_cells)
end

"Struct to store local inertial overland flow model parameters"
@with_kw struct LocalInertialOverlandFlowParameters
    n_cells::Int                   # number of cells [-]
    xwidth::Vector{Float64}             # effective flow width x direction at edge (floodplain) [m]
    ywidth::Vector{Float64}             # effective flow width y direction at edge (floodplain) [m]
    theta::Float64                      # weighting factor (de Almeida et al., 2012) [-]
    h_thresh::Float64                   # depth threshold for calculating flow [m]
    zx_max::Vector{Float64}             # maximum cell elevation at edge [m] (x direction)
    zy_max::Vector{Float64}             # maximum cell elevation at edge [m] (y direction)
    mannings_n_sq::Vector{Float64}      # Manning's roughness squared at edge [(s m-1/3)2]
    z::Vector{Float64}                  # elevation [m] of cell
    froude_limit::Bool                  # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
end

"Initialize shallow water overland flow model parameters"
function LocalInertialOverlandFlowParameters(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
)
    froude_limit = config.model.land_surface_water_flow__froude_limit_flag # limit flow to subcritical according to Froude number
    alpha = config.model.land_local_inertial_flow__alpha_coefficient # stability coefficient for model time step (0.2-0.7)
    theta = config.model.land_local_inertial_flow__theta_coefficient # weighting factor
    waterdepth_threshold = config.model.land_surface_water_flow_threshold__depth # depth threshold for flow at edge

    (; edge_indices, land_indices_2d) = domain.land.network
    (; x_length, y_length) = domain.land.parameters

    @info "Local inertial approach is used for overland flow." alpha theta waterdepth_threshold froude_limit

    mannings_n = ncread(
        dataset,
        config,
        "land_surface_water_flow__manning_n_parameter",
        Routing;
        sel = land_indices_2d,
    )
    elevation_2d =
        ncread(dataset, config, "land_surface_water_flow__ground_elevation", Routing)
    elevation = elevation_2d[land_indices_2d]
    n_cells = length(domain.land.network.land_indices_2d)

    zx_max = zeros(n_cells)
    zy_max = zeros(n_cells)
    for cell_idx in 1:n_cells
        xu = edge_indices.xu[cell_idx]
        if xu <= n_cells
            zx_max[cell_idx] = max(elevation[cell_idx], elevation[xu])
        end
        yu = edge_indices.yu[cell_idx]
        if yu <= n_cells
            zy_max[cell_idx] = max(elevation[cell_idx], elevation[yu])
        end
    end

    # set the effective flow width for river cells in the x and y direction at cell edges.
    # for reservoir cells, h is set to zero (fixed) and not updated, and overland flow from
    # a downstream cell is not possible (effective flowwidth is zero).
    we_x = copy(x_length)
    we_y = copy(y_length)
    set_effective_flowwidth!(we_x, we_y, domain)
    parameters = LocalInertialOverlandFlowParameters(;
        n_cells,
        xwidth = we_x,
        ywidth = we_y,
        theta,
        h_thresh = waterdepth_threshold,
        zx_max,
        zy_max,
        mannings_n_sq = mannings_n .* mannings_n,
        z = elevation,
        froude_limit,
    )
    return parameters
end

"Struct to store local inertial overland flow model boundary conditions"
@with_kw struct LocalInertialOverlandFlowBC
    n_cells::Int
    runoff::Vector{Float64} = zeros(n_cells) # runoff from hydrological model [m³ s⁻¹]
end

"Local inertial overland flow model using the local inertial method"
@with_kw struct LocalInertialOverlandFlowModel <: AbstractOverlandFlowModel
    timestepping::TimeStepping
    boundary_conditions::LocalInertialOverlandFlowBC
    parameters::LocalInertialOverlandFlowParameters
    variables::LocalInertialOverlandFlowVariables
end

"Initialize local inertial overland flow model"
function LocalInertialOverlandFlowModel(dataset::NCDataset, config::Config, domain::Domain)
    alpha_coefficient = config.model.land_local_inertial_flow__alpha_coefficient # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; alpha_coefficient)

    n_cells = length(domain.land.network.land_indices_2d)
    boundary_conditions = LocalInertialOverlandFlowBC(; n_cells)
    parameters = LocalInertialOverlandFlowParameters(dataset, config, domain)
    variables = LocalInertialOverlandFlowVariables(; n_cells)

    overland_flow_model = LocalInertialOverlandFlowModel(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
    )

    return overland_flow_model
end

"""
    stable_timestep(river_flow_model::LocalInertialRiverFlowModel, flow_length::Vector{Float64})
    stable_timestep(overland_flow_model::LocalInertialOverlandFlowModel, parameters::LandParameters)

Compute a stable timestep size for the local inertial approach, based on Bates et al. (2010).

dt = α * (Δx / sqrt(g max(h))
"""
function stable_timestep(
    river_flow_model::LocalInertialRiverFlowModel,
    flow_length::Vector{Float64},
)
    dt_min = Inf
    (; alpha_coefficient) = river_flow_model.timestepping
    (; n_cells) = river_flow_model.parameters
    (; h) = river_flow_model.variables
    @batch per = thread reduction = ((min, dt_min),) for cell_idx in 1:n_cells
        @fastmath @inbounds dt =
            alpha_coefficient * flow_length[cell_idx] /
            sqrt(GRAVITATIONAL_ACCELERATION * h[cell_idx])
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? 60.0 : dt_min
    return dt_min
end

function stable_timestep(
    overland_flow_model::LocalInertialOverlandFlowModel,
    parameters::LandParameters,
)
    dt_min = Inf
    (; alpha_coefficient) = overland_flow_model.timestepping
    (; n_cells) = overland_flow_model.parameters
    (; x_length, y_length, river_location) = parameters
    (; h) = overland_flow_model.variables
    @batch per = thread reduction = ((min, dt_min),) for cell_idx in 1:n_cells
        @fastmath @inbounds dt = if river_location[cell_idx] == 0
            alpha_coefficient * min(x_length[cell_idx], y_length[cell_idx]) /
            sqrt(GRAVITATIONAL_ACCELERATION * h[cell_idx])
        else
            Inf
        end
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? 60.0 : dt_min
    return dt_min
end

"""
Update boundary condition `runoff` overland flow model `LocalInertialOverlandFlowModel` for a
single timestep.
"""
function update_bc_overland_flow_model!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    external_models::NamedTuple,
    domain::Domain,
    dt::Float64,
)
    (; soil, runoff, subsurface_flow) = external_models
    (; net_runoff) = soil.variables
    (; net_runoff_river) = runoff.variables
    (; area) = domain.land.parameters
    river_indices = domain.river.network.cell_indices_containing_river

    @. overland_flow_model.boundary_conditions.runoff =
        net_runoff / 1000.0 * area / dt + net_runoff_river * area * 0.001 / dt
    overland_flow_model.boundary_conditions.runoff[river_indices] .+=
        get_flux_to_river(subsurface_flow, river_indices)
    return nothing
end

"""
Update subsurface flow contribution to inflow of a reservoir model for a river flow model
`LocalInertialRiverFlowModel` for a single timestep.
"""
function update_inflow!(
    reservoir_model::ReservoirModel,
    river_flow_model::LocalInertialRiverFlowModel,
    subsurface_flow_model::AbstractSubsurfaceFlowModel,
    network::NetworkReservoir,
)
    (; cell_indices_containing_reservoir) = network
    (; inflow_subsurface) = reservoir_model.boundary_conditions
    inflow_subsurface .= get_inflow_reservoir(
        river_flow_model,
        subsurface_flow_model,
        cell_indices_containing_reservoir,
    )
    return nothing
end
update_inflow!(
    ::Nothing,
    ::LocalInertialRiverFlowModel,
    ::AbstractSubsurfaceFlowModel,
    ::NetworkReservoir,
) = nothing

"""
Helper function to set flow variables of the `LocalInertialOverlandFlowModel` model to zero. This
is done at the start of each simulation timestep, during the timestep the total (weighted)
sum is computed from values at each sub timestep.
"""
function set_flow_vars!(overland_flow_model::LocalInertialOverlandFlowModel)
    (; qx_av, qy_av) = overland_flow_model.variables
    qx_av .= 0.0
    qy_av .= 0.0
    return nothing
end

"""
Helper function to compute average flow variables of the `LocalInertialOverlandFlowModel` model.
This is done at the end of each simulation timestep.
"""
function average_flow_vars!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    dt::Float64,
)
    (; qx_av, qy_av) = overland_flow_model.variables
    qx_av ./= dt
    qy_av ./= dt
    return nothing
end

"""
Update combined river `LocalInertialRiverFlowModel` and overland flow `LocalInertialOverlandFlowModel`
models for a single timestep `dt`. An adaptive timestepping method is used (computing a sub
timestep `dt_s`).
"""
function update_overland_flow_model!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    river_flow_model::LocalInertialRiverFlowModel,
    domain::Domain,
    clock::Clock;
    update_h = false,
)
    (; reservoir) = river_flow_model.boundary_conditions
    (; flow_length) = domain.river.parameters
    (; parameters) = domain.land

    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)
    set_flow_vars!(river_flow_model)
    set_flow_vars!(overland_flow_model)

    dt = tosecond(clock.dt)
    t = 0.0
    while t < dt
        dt_river = stable_timestep(river_flow_model, flow_length)
        dt_land = stable_timestep(overland_flow_model, parameters)
        dt_s = min(dt_river, dt_land)
        dt_s = check_timestepsize(dt_s, t, dt)

        local_inertial_update_fluxes!(overland_flow_model, domain, dt_s)
        update_inflow_reservoir!(overland_flow_model, reservoir, domain)
        local_inertial_river_update!(river_flow_model, domain, dt_s, dt, update_h)
        local_inertial_update_water_depth!(
            overland_flow_model,
            river_flow_model,
            domain,
            dt_s,
        )

        t += dt_s
    end
    average_flow_vars!(river_flow_model, dt)
    average_flow_vars!(overland_flow_model, dt)
    average_reservoir_vars!(reservoir, dt)

    return nothing
end

"""
Update flow for a single direction in the local inertial overland flow model.
`is_x_direction`: true for x-direction (xu/xd), false for y-direction (yu/yd)
"""
@inline function update_directional_flow!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    domain::Domain,
    cell_idx::Int,
    dt::Float64,
    is_x_direction::Bool,
)
    indices = domain.land.network.edge_indices
    (; x_length, y_length) = domain.land.parameters
    land_v = overland_flow_model.variables
    land_p = overland_flow_model.parameters

    # Select direction-specific parameters based on the boolean flag
    if is_x_direction
        upstream_idx = indices.xu[cell_idx]
        downstream_idx = indices.xd[cell_idx]
        width = land_p.ywidth[cell_idx]
        z_max = land_p.zx_max[cell_idx]
        length_vec = x_length
        q_current = land_v.qx
        q_prev = land_v.qx0
        q_av = land_v.qx_av
    else
        upstream_idx = indices.yu[cell_idx]
        downstream_idx = indices.yd[cell_idx]
        width = land_p.xwidth[cell_idx]
        z_max = land_p.zy_max[cell_idx]
        length_vec = y_length
        q_current = land_v.qy
        q_prev = land_v.qy0
        q_av = land_v.qy_av
    end

    # the effective flow width is zero when the river width exceeds the cell width and
    # floodplain flow is not calculated.
    if upstream_idx <= land_p.n_cells && width != 0.0
        zs_current = land_p.z[cell_idx] + land_v.h[cell_idx]
        zs_upstream = land_p.z[upstream_idx] + land_v.h[upstream_idx]
        zs_max = max(zs_current, zs_upstream)
        hf = (zs_max - z_max)

        if hf > land_p.h_thresh
            length = 0.5 * (length_vec[cell_idx] + length_vec[upstream_idx]) # can be precalculated
            q_current[cell_idx] = local_inertial_flow(
                land_p.theta,
                q_prev[cell_idx],
                q_prev[downstream_idx],
                q_prev[upstream_idx],
                zs_current,
                zs_upstream,
                hf,
                width,
                length,
                land_p.mannings_n_sq[cell_idx],
                land_p.froude_limit,
                dt,
            )
            # limit q in case water is not available
            if land_v.h[cell_idx] <= 0.0
                q_current[cell_idx] = min(q_current[cell_idx], 0.0)
            end
            if land_v.h[upstream_idx] <= 0.0
                q_current[cell_idx] = max(q_current[cell_idx], 0.0)
            end
        else
            q_current[cell_idx] = 0.0
        end
        q_av[cell_idx] += q_current[cell_idx] * dt
    end
    return nothing
end

"""
Update fluxes for overland flow `LocalInertialOverlandFlowModel` model for a single timestep
`dt`.
"""
function local_inertial_update_fluxes!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    domain::Domain,
    dt::Float64,
)
    p = overland_flow_model.parameters
    v = overland_flow_model.variables

    v.qx0 .= v.qx
    v.qy0 .= v.qy

    @batch per = thread minbatch = 6000 for cell_idx in 1:(p.n_cells)
        # update qx (x-direction)
        update_directional_flow!(overland_flow_model, domain, cell_idx, dt, true)

        # update qy (y-direction)
        update_directional_flow!(overland_flow_model, domain, cell_idx, dt, false)
    end
    return nothing
end

"""
Update boundary condition inflow to a reservoir from land `inflow_reservoir` of combined
river `LocalInertialRiverFlowModel`and overland flow `LocalInertialOverlandFlowModel` models for a
single timestep.
"""
function update_inflow_reservoir!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    reservoir_model::Union{ReservoirModel, Nothing},
    domain::Domain,
)
    indices = domain.land.network.edge_indices
    reservoir_indices = domain.reservoir.network.cell_indices_containing_reservoir
    land_bc = overland_flow_model.boundary_conditions
    land_v = overland_flow_model.variables

    for (reservoir_idx, cell_idx) in enumerate(reservoir_indices)
        yd = indices.yd[cell_idx]
        xd = indices.xd[cell_idx]
        reservoir_model.boundary_conditions.inflow_overland[reservoir_idx] =
            land_bc.runoff[cell_idx] +
            (land_v.qx[xd] - land_v.qx[cell_idx] + land_v.qy[yd] - land_v.qy[cell_idx])
    end
    return nothing
end

"""
Compute storage change for a river cell from fluxes.
"""
@inline function compute_river_storage_change(
    overland_flow_model::LocalInertialOverlandFlowModel,
    river_flow_model::LocalInertialRiverFlowModel,
    domain::Domain,
    cell_idx::Int,
    dt::Float64,
)
    indices = domain.land.network.edge_indices
    inds_river = domain.land.network.river_cell_indices
    edges_at_node = domain.river.network.edges_at_node
    river_idx = inds_river[cell_idx]

    yd = indices.yd[cell_idx]
    xd = indices.xd[cell_idx]

    net_river_flow =
        sum_at(river_flow_model.variables.q, edges_at_node.src[river_idx]) -
        sum_at(river_flow_model.variables.q, edges_at_node.dst[river_idx])
    net_land_flow =
        overland_flow_model.variables.qx[xd] - overland_flow_model.variables.qx[cell_idx] +
        overland_flow_model.variables.qy[yd] - overland_flow_model.variables.qy[cell_idx]
    net_flow =
        net_river_flow +
        net_land_flow +
        overland_flow_model.boundary_conditions.runoff[cell_idx] -
        river_flow_model.boundary_conditions.abstraction[river_idx]
    storage_change = net_flow * dt

    return storage_change
end

"""
Compute external inflow for river cells, including negative inflow (abstraction).
Returns tuple: (inflow, abstraction_to_add)
"""
@inline function compute_external_inflow(
    river_flow_model::LocalInertialRiverFlowModel,
    overland_flow_model::LocalInertialOverlandFlowModel,
    cell_idx::Int,
    river_idx::Int,
    dt::Float64,
)
    if river_flow_model.boundary_conditions.external_inflow[river_idx] < 0.0
        available_volume =
            if overland_flow_model.variables.storage[cell_idx] >=
               river_flow_model.parameters.bankfull_storage[river_idx]
                river_flow_model.parameters.bankfull_storage[river_idx]
            else
                river_flow_model.variables.storage[river_idx]
            end
        _abstraction = min(
            -river_flow_model.boundary_conditions.external_inflow[river_idx],
            available_volume / dt * 0.80,
        )
        return (-_abstraction, _abstraction * dt)
    else
        return (river_flow_model.boundary_conditions.external_inflow[river_idx], 0.0)
    end
end

"""
Compute river and land water depths based on total storage and bankfull capacity.
Returns tuple: (river_h, land_h, river_storage)
"""
@inline function compute_water_depths(
    total_storage::Float64,
    river_idx::Int,
    cell_idx::Int,
    river::LocalInertialRiverFlowModel,
    domain::Domain,
)
    if total_storage >= river.parameters.bankfull_storage[river_idx]
        # Storage exceeds bankfull capacity - water spills onto floodplain
        river_h =
            river.parameters.bankfull_depth[river_idx] +
            (total_storage - river.parameters.bankfull_storage[river_idx]) / (
                domain.land.parameters.x_length[cell_idx] *
                domain.land.parameters.y_length[cell_idx]
            )
        land_h = river_h - river.parameters.bankfull_depth[river_idx]
        river_storage =
            river_h *
            domain.river.parameters.flow_length[river_idx] *
            domain.river.parameters.flow_width[river_idx]
        return (river_h, land_h, river_storage)
    else
        # Storage is within channel capacity
        river_h =
            total_storage / (
                domain.river.parameters.flow_length[river_idx] *
                domain.river.parameters.flow_width[river_idx]
            )
        return (river_h, 0.0, total_storage)
    end
end

"""
Compute storage change for a land cell from horizontal fluxes and runoff.
"""
@inline function compute_land_storage_change(
    overland_flow_model::LocalInertialOverlandFlowModel,
    network::NetworkLand,
    cell_idx::Int,
    dt::Float64,
)
    indices = network.edge_indices
    yd = indices.yd[cell_idx]
    xd = indices.xd[cell_idx]

    return (
        overland_flow_model.variables.qx[xd] - overland_flow_model.variables.qx[cell_idx] +
        overland_flow_model.variables.qy[yd] - overland_flow_model.variables.qy[cell_idx] +
        overland_flow_model.boundary_conditions.runoff[cell_idx]
    ) * dt
end

"""
Update storage and water depth for a single river cell.
"""
@inline function update_river_cell_storage_and_depth!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    river_flow_model::LocalInertialRiverFlowModel,
    domain::Domain,
    cell_idx::Int,
    dt::Float64,
)
    inds_river = domain.land.network.river_cell_indices
    river_idx = inds_river[cell_idx]

    # Compute and apply storage change from fluxes
    storage_change = compute_river_storage_change(
        overland_flow_model,
        river_flow_model,
        domain,
        cell_idx,
        dt,
    )
    overland_flow_model.variables.storage[cell_idx] += storage_change

    # Handle negative storage
    if overland_flow_model.variables.storage[cell_idx] < 0.0
        overland_flow_model.variables.error[cell_idx] +=
            abs(overland_flow_model.variables.storage[cell_idx])
        overland_flow_model.variables.storage[cell_idx] = 0.0 # set storage to zero
    end

    # Compute and apply external inflow
    inflow, abstraction_to_add = compute_external_inflow(
        river_flow_model,
        overland_flow_model,
        cell_idx,
        river_idx,
        dt,
    )
    overland_flow_model.variables.storage[cell_idx] += inflow * dt
    river_flow_model.boundary_conditions.actual_external_abstraction_av[river_idx] +=
        abstraction_to_add

    # Compute and apply water depths
    river_h, land_h, river_storage = compute_water_depths(
        overland_flow_model.variables.storage[cell_idx],
        river_idx,
        cell_idx,
        river_flow_model,
        domain,
    )
    river_flow_model.variables.h[river_idx] = river_h
    overland_flow_model.variables.h[cell_idx] = land_h
    river_flow_model.variables.storage[river_idx] = river_storage

    return nothing
end

"""
Update storage and water depth for a single land cell (non-river).
"""
@inline function update_land_cell_storage_and_depth!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    domain::DomainLand,
    cell_idx::Int,
    dt::Float64,
)
    # Compute and apply storage change
    storage_change =
        compute_land_storage_change(overland_flow_model, domain.network, cell_idx, dt)
    overland_flow_model.variables.storage[cell_idx] += storage_change

    # Handle negative storage
    if overland_flow_model.variables.storage[cell_idx] < 0.0
        overland_flow_model.variables.error[cell_idx] +=
            abs(overland_flow_model.variables.storage[cell_idx])
        overland_flow_model.variables.storage[cell_idx] = 0.0 # set storage to zero
    end

    # Update water depth
    overland_flow_model.variables.h[cell_idx] =
        overland_flow_model.variables.storage[cell_idx] /
        (domain.parameters.x_length[cell_idx] * domain.parameters.y_length[cell_idx])

    return nothing
end

"""
Update storage and water depth for combined river `LocalInertialRiverFlowModel` and overland flow
`LocalInertialOverlandFlowModel` models for a single timestep `dt`.
"""
function local_inertial_update_water_depth!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    river_flow_model::LocalInertialRiverFlowModel,
    domain::Domain,
    dt::Float64,
)
    (; river_location, reservoir_outlet) = domain.land.parameters

    @batch per = thread minbatch = 6000 for cell_idx in
                                            1:(overland_flow_model.parameters.n_cells)
        if river_location[cell_idx]
            # Process river cells (excluding reservoir outlets)
            if !reservoir_outlet[cell_idx]
                update_river_cell_storage_and_depth!(
                    overland_flow_model,
                    river_flow_model,
                    domain,
                    cell_idx,
                    dt,
                )
            end
        else
            # Process land cells (non-river)
            update_land_cell_storage_and_depth!(
                overland_flow_model,
                domain.land,
                cell_idx,
                dt,
            )
        end
    end
    return nothing
end

"""
    FloodPlainProfile

Floodplain `storage` is a function of `depth` (flood depth intervals). Based on the
cumulative floodplain `storage` a floodplain profile as a function of `flood_depth` is
derived with floodplain area `a` (cumulative) and wetted perimeter radius `p` (cumulative).
"""
@with_kw struct FloodPlainProfile
    depth::Vector{Float64}        # Flood depth [m]
    storage::Matrix{Float64}      # Flood storage (cumulative) [m³]
    width::Matrix{Float64}        # Flood width [m]
    a::Matrix{Float64}            # Flow area (cumulative) [m²]
    p::Matrix{Float64}            # Wetted perimeter (cumulative) [m]
end

"Initialize floodplain profile `FloodPlainProfile`"
function FloodPlainProfile(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    index_pit::Vector{Int},
)
    (; river_indices_2d) = domain.network
    (; flow_width, flow_length) = domain.parameters
    storage = ncread(
        dataset,
        config,
        "floodplain_water__sum_of_volume_per_depth",
        Routing;
        sel = river_indices_2d,
    )
    n_cells = length(river_indices_2d)

    # for convenience (interpolation) flood depth 0.0 m is added, with associated area (a),
    # storage, width (river width) and wetted perimeter (p).
    storage = vcat(zeros(n_cells)', storage)
    start_storage = storage
    flood_depths = Float64.(dataset["flood_depth"][:])
    pushfirst!(flood_depths, 0.0)
    n_depths = length(flood_depths)

    p = zeros(n_depths, n_cells)
    a = zeros(n_depths, n_cells)
    segment_storage = zeros(n_depths, n_cells)
    width = zeros(n_depths, n_cells)
    width[1, :] = flow_width[1:n_cells]

    # determine flow area (a), width and wetted perimeter (p) FloodPlainModel
    h = diff(flood_depths)
    incorrect_vol = 0
    riv_cells = 0
    error_vol = 0
    for cell_idx in 1:n_cells
        riv_cell = 0
        diff_storage = diff(storage[:, cell_idx])

        for flood_depth_idx in 1:(n_depths - 1)
            # assume rectangular shape of flood depth segment
            width[flood_depth_idx + 1, cell_idx] =
                diff_storage[flood_depth_idx] /
                (h[flood_depth_idx] * flow_length[cell_idx])
            # check provided flood storage (floodplain width should be constant or increasing
            # as a function of flood depth)
            if width[flood_depth_idx + 1, cell_idx] <
               width[flood_depth_idx, cell_idx]
                # raise warning only if difference is larger than rounding error of 0.01 m³
                if (
                    (
                        width[flood_depth_idx, cell_idx] -
                        width[flood_depth_idx + 1, cell_idx]
                    ) *
                    h[flood_depth_idx] *
                    flow_length[cell_idx]
                ) > 0.01
                    incorrect_vol += 1
                    riv_cell = 1
                    error_vol =
                        error_vol + (
                            (
                                width[flood_depth_idx, cell_idx] -
                                width[flood_depth_idx + 1, cell_idx]
                            ) *
                            h[flood_depth_idx] *
                            flow_length[cell_idx]
                        )
                end
                width[flood_depth_idx + 1, cell_idx] =
                    width[flood_depth_idx, cell_idx]
            end
            a[flood_depth_idx + 1, cell_idx] =
                width[flood_depth_idx + 1, cell_idx] * h[flood_depth_idx]
            p[flood_depth_idx + 1, cell_idx] =
                (
                    width[flood_depth_idx + 1, cell_idx] -
                    width[flood_depth_idx, cell_idx]
                ) + 2.0 * h[flood_depth_idx]
            segment_storage[flood_depth_idx + 1, cell_idx] =
                a[flood_depth_idx + 1, cell_idx] * flow_length[cell_idx]
            if flood_depth_idx == 1
                # for interpolation wetted perimeter at flood depth 0.0 is required
                p[flood_depth_idx, cell_idx] =
                    p[flood_depth_idx + 1, cell_idx] - 2.0 * h[flood_depth_idx]
            end
        end

        p[2:end, cell_idx] = cumsum(p[2:end, cell_idx])
        a[:, cell_idx] = cumsum(a[:, cell_idx])
        storage[:, cell_idx] = cumsum(segment_storage[:, cell_idx])

        riv_cells += riv_cell
    end

    if incorrect_vol > 0
        perc_riv_cells = round(100.0 * (riv_cells / n_cells); digits = 2)
        perc_error_vol = round(100.0 * (error_vol / sum(start_storage[end, :])); digits = 2)
        @warn string(
            "The provided storage of $incorrect_vol rectangular floodplain schematization",
            " segments for $riv_cells river cells ($perc_riv_cells % of total river cells)",
            " is not correct and has been increased with $perc_error_vol % of provided storage.",
        )
    end

    # set floodplain parameters for ghost points
    storage = hcat(storage, storage[:, index_pit])
    width = hcat(width, width[:, index_pit])
    a = hcat(a, a[:, index_pit])
    p = hcat(p, p[:, index_pit])

    # initialize floodplain profile parameters
    profile = FloodPlainProfile(; storage, width, depth = flood_depths, a, p)
    return profile
end

"Struct to store floodplain flow model parameters"
@with_kw struct FloodPlainParameters
    profile::FloodPlainProfile          # floodplain profile
    mannings_n::Vector{Float64}         # manning's roughness [s m-1/3]
    mannings_n_sq::Vector{Float64}      # manning's roughness squared at edge [(s m-1/3)2]
    zb_max::Vector{Float64}             # maximum bankfull elevation at edge [m]
end

"Initialize floodplain flow model parameters"
function FloodPlainParameters(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    zb_floodplain::Vector{Float64},
    index_pit::Vector{Int},
)
    (; river_indices_2d, nodes_at_edge, graph) = domain.network
    (; flow_length) = domain.parameters
    n_edges = ne(graph)
    profile = FloodPlainProfile(dataset, config, domain, index_pit)

    mannings_n = ncread(
        dataset,
        config,
        "floodplain_water_flow__manning_n_parameter",
        Routing;
        sel = river_indices_2d,
    )
    # manning roughness at edges
    append!(mannings_n, mannings_n[index_pit]) # copy to ghost nodes
    mannings_n_sq = zeros(n_edges)
    zb_max = zeros(n_edges)
    for river_edge_idx in 1:n_edges
        src_node = nodes_at_edge.src[river_edge_idx]
        dst_node = nodes_at_edge.dst[river_edge_idx]
        mannings_n_i =
            (
                mannings_n[dst_node] * flow_length[dst_node] +
                mannings_n[src_node] * flow_length[src_node]
            ) / (flow_length[dst_node] + flow_length[src_node])
        mannings_n_sq[river_edge_idx] = mannings_n_i * mannings_n_i
        zb_max[river_edge_idx] = max(zb_floodplain[src_node], zb_floodplain[dst_node])
    end
    parameters = FloodPlainParameters(profile, mannings_n, mannings_n_sq, zb_max)
    return parameters
end

"Struct to store floodplain flow model variables"
@with_kw struct FloodPlainVariables
    n_cells::Int
    n_edges::Int
    storage::Vector{Float64} = zeros(n_cells)    # storage [m³]
    h::Vector{Float64}                     # water depth [m]
    error::Vector{Float64} = zeros(n_cells)      # error storage [m³]
    a::Vector{Float64} = zeros(n_edges)    # flow area at edge [m²]
    r::Vector{Float64} = zeros(n_edges)    # hydraulic radius at edge [m]
    hf::Vector{Float64} = zeros(n_edges)   # water depth at edge [m]
    q0::Vector{Float64} = zeros(n_edges)   # discharge at edge at previous time step
    q::Vector{Float64} = zeros(n_edges)    # discharge at edge  [m³ s⁻¹]
    q_av::Vector{Float64} = zeros(n_edges) # average river discharge at edge  [m³ s⁻¹] for model timestep Δt
    hf_index::Vector{Int} = zeros(Int, n_edges) # edge index with `hf` [-] above depth threshold
end

"Initialize floodplain flow model variables"
function FloodPlainVariables(n_cells::Int, n_edges::Int, index_pit::Vector{Int})
    variables = FloodPlainVariables(;
        n_cells,
        n_edges,
        h = zeros(n_cells + length(index_pit)),
    )
    return variables
end

"Floodplain flow model"
@with_kw struct FloodPlainModel <: AbstractFloodPlainModel
    parameters::FloodPlainParameters
    variables::FloodPlainVariables
end

"Determine the initial floodplain storage"
function initialize_storage!(river, domain::Domain, n_cells::Int)
    (; flow_width, flow_length) = domain.river.parameters
    (; floodplain) = river
    (; profile) = floodplain.parameters
    for cell_idx in 1:n_cells
        lower_idx, upper_idx =
            interpolation_indices(floodplain.variables.h[cell_idx], profile.depth)
        a = flow_area(
            profile.width[upper_idx, cell_idx],
            profile.a[lower_idx, cell_idx],
            profile.depth[lower_idx],
            floodplain.variables.h[cell_idx],
        )
        a = max(
            a - (flow_width[cell_idx] * floodplain.variables.h[cell_idx]),
            0.0,
        )
        floodplain.variables.storage[cell_idx] = flow_length[cell_idx] * a
    end
    return nothing
end

"helper function to get interpolation indices"
function interpolation_indices(x, v::AbstractVector)
    lower_idx = 1
    for element_idx in eachindex(v)
        if v[element_idx] <= x
            lower_idx = element_idx
        end
    end
    if lower_idx == length(v)
        upper_idx = lower_idx
    else
        upper_idx = lower_idx + 1
    end
    return lower_idx, upper_idx
end

"""
    flow_area(width, area, depth, h)

Compute floodplain flow area based on flow depth `h` and floodplain `depth`, `area` and
`width` of a floodplain profile.
"""
function flow_area(width, area, depth, h)
    dh = h - depth  # depth at lower_idx
    area = area + (width * dh) # area at lower_idx, width at upper_idx
    return area
end

"""
    function wetted_perimeter(p, depth, h)

Compute floodplain wetted perimeter based on flow depth `h` and floodplain `depth` and
wetted perimeter `p` of a floodplain profile.
"""
function wetted_perimeter(p, depth, h)
    dh = h - depth # depth at lower_idx
    p += 2.0 * dh # p at lower_idx
    return p
end

"Compute flood depth by interpolating flood storage `flood_storage` using flood depth intervals."
function flood_depth(
    profile::FloodPlainProfile,
    flood_storage::Float64,
    flow_length::Float64,
    cell_idx::Int,
)
    lower_idx, upper_idx =
        interpolation_indices(flood_storage, @view profile.storage[:, cell_idx])
    ΔA = (flood_storage - profile.storage[lower_idx, cell_idx]) / flow_length
    dh = ΔA / profile.width[upper_idx, cell_idx]
    flood_depth = profile.depth[lower_idx] + dh
    return flood_depth
end

"Initialize floodplain geometry and `FloodPlainModel` variables and parameters"
function FloodPlainModel(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    zb_floodplain::Vector{Float64},
)
    (; river_indices_2d, local_drain_direction, graph) = domain.network
    n_cells = length(river_indices_2d)
    index_pit = findall(x -> x == 5, local_drain_direction)
    parameters = FloodPlainParameters(dataset, config, domain, zb_floodplain, index_pit)
    n_edges = ne(graph)
    variables = FloodPlainVariables(n_cells, n_edges, index_pit)

    floodplain = FloodPlainModel(; parameters, variables)
    return floodplain
end
