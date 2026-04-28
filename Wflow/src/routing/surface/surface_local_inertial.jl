abstract type AbstractFloodPlainModel end

"Struct for storing local inertial river flow model parameters"
@with_kw struct LocalInertialRiverFlowParameters
    n_river_cells::Int                      # number of cells [-]
    n_river_edges::Int                      # number of edges [-]
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

    n_river_cells = length(river_indices_2d)
    index_pit = findall(x -> x == 5, local_drain_direction)
    # set ghost points for boundary condition (downstream river outlet): river width, bed
    # elevation, manning n is copied from the upstream cell.
    append!(flow_length, riverlength_bc)
    append!(zb, zb[index_pit])
    append!(flow_width, flow_width[index_pit])
    append!(mannings_n, mannings_n[index_pit])
    append!(bankfull_depth, bankfull_depth[index_pit])

    # determine z, width, length and manning's n at edges
    n_river_edges = ne(graph)
    zb_max = zeros(n_river_edges)
    width_at_edge = zeros(n_river_edges)
    length_at_edge = zeros(n_river_edges)
    mannings_n_sq = zeros(n_river_edges)
    for river_edge_idx in 1:n_river_edges
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
        n_river_cells,
        n_river_edges,
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
    n_river_cells::Int
    n_river_edges::Int
    q::Vector{Float64} = zeros(n_river_edges)   # river discharge at edge (subgrid channel) [m³ s⁻¹]
    q0::Vector{Float64} = zeros(n_river_edges) # river discharge at edge (subgrid channel) at previous time step [m³ s⁻¹]
    q_av::Vector{Float64}                     # average river channel (+ floodplain) discharge at edge [m³ s⁻¹] (model timestep Δt)
    q_channel_av::Vector{Float64}             # average river channel discharge at edge [m³ s⁻¹] (for model timestep Δt)
    h::Vector{Float64}                        # water depth [m]
    zs_max::Vector{Float64} = zeros(n_river_edges)  # maximum water elevation at edge [m]
    zs_src::Vector{Float64} = zeros(n_river_edges)  # water elevation of source node of edge [m]
    zs_dst::Vector{Float64} = zeros(n_river_edges)  # water elevation of downstream node of edge [m]
    hf::Vector{Float64} = zeros(n_river_edges)       # water depth at edge [m]
    a::Vector{Float64} = zeros(n_river_edges)        # flow area at edge [m²]
    r::Vector{Float64} = zeros(n_river_edges)        # wetted perimeter at edge [m]
    storage::Vector{Float64} = zeros(n_river_cells)        # river storage [m³]
    error::Vector{Float64} = zeros(n_river_cells)          # error storage [m³]
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

    n_river_cells = length(river_indices_2d)
    n_river_edges = ne(graph)
    # set river depth h to zero (including reservoir locations)
    h = zeros(n_river_cells)
    q_av = zeros(n_river_edges)
    # set ghost points for boundary condition (downstream river outlet): river depth `h`
    append!(h, riverdepth_bc)
    variables = LocalInertialRiverFlowVariables(;
        n_river_cells,
        n_river_edges,
        q_av,
        q_channel_av = config.model.floodplain_1d__flag ? zeros(n_river_edges) : q_av,
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

"Initialize shallow water river flow model `LocalInertialRiverFlow`"
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

    n_river_cells = length(domain.network.river_indices_2d)
    river_flow = LocalInertialRiverFlowModel(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain,
        allocation = do_water_demand(config) ? AllocationRiverModel(n_river_cells) :
                     NoAllocationRiverModel(n_river_cells),
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

"Update local inertial river flow model `LocalInertialRiverFlow` for a single timestep"
function local_inertial_river_update!(
    river_flow_model::LocalInertialRiverFlowModel,
    domain::Domain,
    dt::Float64,
    dt_forcing::Float64,
    update_h::Bool,
)
    (; nodes_at_edge, edges_at_node) = domain.river.network
    (; flow_length, flow_width) = domain.river.parameters
    (; inwater, abstraction, external_inflow, actual_external_abstraction_av) =
        river_flow_model.boundary_conditions

    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters

    river_v.q0 .= river_v.q
    if !isnothing(river_flow_model.floodplain)
        river_flow_model.floodplain.variables.q0 .= river_flow_model.floodplain.variables.q
    end
    @batch per = thread minbatch = 1000 for river_edge_idx in river_p.active_e
        river_cell_idx_src = nodes_at_edge.src[river_edge_idx]
        river_cell_idx_dst = nodes_at_edge.dst[river_edge_idx]
        river_v.zs_src[river_edge_idx] =
            river_p.zb[river_cell_idx_src] + river_v.h[river_cell_idx_src]
        river_v.zs_dst[river_edge_idx] =
            river_p.zb[river_cell_idx_dst] + river_v.h[river_cell_idx_dst]

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
            river_v.h[river_cell_idx_src] <= 0.0,
            min(river_v.q[river_edge_idx], 0.0),
            river_v.q[river_edge_idx],
        )
        river_v.q[river_edge_idx] = ifelse(
            river_v.h[river_cell_idx_dst] <= 0.0,
            max(river_v.q[river_edge_idx], 0.0),
            river_v.q[river_edge_idx],
        )
        # average river discharge (here accumulated for model timestep Δt)
        river_v.q_av[river_edge_idx] += river_v.q[river_edge_idx] * dt
    end
    if !isnothing(river_flow_model.floodplain)
        floodplain_p = river_flow_model.floodplain.parameters
        floodplain_v = river_flow_model.floodplain.variables

        @batch per = thread minbatch = 1000 for river_edge_idx in 1:(river_p.n_river_edges)
            floodplain_v.hf[river_edge_idx] = max(
                river_v.zs_max[river_edge_idx] - floodplain_p.zb_max[river_edge_idx],
                0.0,
            )
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

        @batch per = thread minbatch = 1000 for flood_edge_idx in 1:n_active_flood_edges
            river_edge_idx = floodplain_v.hf_index[flood_edge_idx]
            river_cell_idx_src = nodes_at_edge.src[river_edge_idx]
            river_cell_idx_dst = nodes_at_edge.dst[river_edge_idx]

            i0 = 0
            for depth_idx in eachindex(floodplain_p.profile.depth)
                i0 +=
                    1 * (
                        floodplain_p.profile.depth[depth_idx] <=
                        floodplain_v.hf[river_edge_idx]
                    )
            end
            i1 = max(i0, 1)
            i2 = ifelse(i1 == length(floodplain_p.profile.depth), i1, i1 + 1)

            a_src = flow_area(
                floodplain_p.profile.width[i2, river_cell_idx_src],
                floodplain_p.profile.a[i1, river_cell_idx_src],
                floodplain_p.profile.depth[i1],
                floodplain_v.hf[river_edge_idx],
            )
            a_src = max(
                a_src - (floodplain_v.hf[river_edge_idx] * flow_width[river_cell_idx_src]),
                0.0,
            )

            a_dst = flow_area(
                floodplain_p.profile.width[i2, river_cell_idx_dst],
                floodplain_p.profile.a[i1, river_cell_idx_dst],
                floodplain_p.profile.depth[i1],
                floodplain_v.hf[river_edge_idx],
            )
            a_dst = max(
                a_dst - (floodplain_v.hf[river_edge_idx] * flow_width[river_cell_idx_dst]),
                0.0,
            )

            floodplain_v.a[river_edge_idx] = min(a_src, a_dst)

            floodplain_v.r[river_edge_idx] = ifelse(
                a_src < a_dst,
                a_src / wetted_perimeter(
                    floodplain_p.profile.p[i1, river_cell_idx_src],
                    floodplain_p.profile.depth[i1],
                    floodplain_v.hf[river_edge_idx],
                ),
                a_dst / wetted_perimeter(
                    floodplain_p.profile.p[i1, river_cell_idx_dst],
                    floodplain_p.profile.depth[i1],
                    floodplain_v.hf[river_edge_idx],
                ),
            )

            floodplain_v.q[river_edge_idx] = ifelse(
                floodplain_v.a[river_edge_idx] > 1.0e-05,
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
                ),
                0.0,
            )

            # limit floodplain q in case water is not available
            floodplain_v.q[river_edge_idx] = ifelse(
                floodplain_v.h[river_cell_idx_src] <= 0.0,
                min(floodplain_v.q[river_edge_idx], 0.0),
                floodplain_v.q[river_edge_idx],
            )
            floodplain_v.q[river_edge_idx] = ifelse(
                floodplain_v.h[river_cell_idx_dst] <= 0.0,
                max(floodplain_v.q[river_edge_idx], 0.0),
                floodplain_v.q[river_edge_idx],
            )

            floodplain_v.q[river_edge_idx] = ifelse(
                floodplain_v.q[river_edge_idx] * river_v.q[river_edge_idx] < 0.0,
                0.0,
                floodplain_v.q[river_edge_idx],
            )
            # average floodplain discharge (here accumulated for model timestep Δt)
            floodplain_v.q_av[river_edge_idx] += floodplain_v.q[river_edge_idx] * dt
        end
    end
    # For reservoir locations the local inertial solution is replaced by the reservoir
    # model. These locations are handled as boundary conditions in the local inertial model
    # (fixed h).
    (; reservoir) = river_flow_model.boundary_conditions
    inds_reservoir = domain.reservoir.network.river_cell_indices_containing_reservoir
    if !isnothing(reservoir)
        res_bc = reservoir.boundary_conditions
    end

    for (reservoir_idx, river_cell_idx) in enumerate(inds_reservoir)
        q_in = get_inflow_reservoir(river_flow_model, edges_at_node.src[river_cell_idx])
        # If external_inflow < 0, abstraction is limited
        if res_bc.external_inflow[reservoir_idx] < 0.0
            _abstraction = min(
                -res_bc.external_inflow[reservoir_idx],
                (reservoir.variables.storage[reservoir_idx] / dt) * 0.98,
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
        update_reservoir_model!(reservoir, reservoir_idx, net_inflow, dt, dt_forcing)
        river_v.q[river_cell_idx] = reservoir.variables.outflow[reservoir_idx]
        # average river discharge (here accumulated for model timestep Δt)
        river_v.q_av[river_cell_idx] += river_v.q[river_cell_idx] * dt
    end
    if update_h
        @batch per = thread minbatch = 1000 for river_cell_idx in river_p.active_n
            q_src = sum_at(river_v.q, edges_at_node.src[river_cell_idx])
            q_dst = sum_at(river_v.q, edges_at_node.dst[river_cell_idx])
            # internal abstraction (water demand) is limited by river storage and negative
            # external inflow as part of water allocation computations.
            river_v.storage[river_cell_idx] =
                river_v.storage[river_cell_idx] +
                (q_src - q_dst + inwater[river_cell_idx] - abstraction[river_cell_idx]) * dt

            if river_v.storage[river_cell_idx] < 0.0
                river_v.error[river_cell_idx] =
                    river_v.error[river_cell_idx] + abs(river_v.storage[river_cell_idx])
                river_v.storage[river_cell_idx] = 0.0 # set storage to zero
            end
            # limit negative external inflow
            if external_inflow[river_cell_idx] < 0.0
                _abstraction = min(
                    -external_inflow[river_cell_idx],
                    river_v.storage[river_cell_idx] / dt * 0.80,
                )
                actual_external_abstraction_av[river_cell_idx] += _abstraction * dt
                _inflow = -_abstraction
            else
                _inflow = external_inflow[river_cell_idx]
            end
            river_v.storage[river_cell_idx] += _inflow * dt # add external inflow
            river_v.h[river_cell_idx] =
                river_v.storage[river_cell_idx] /
                (flow_length[river_cell_idx] * flow_width[river_cell_idx])

            if !isnothing(river_flow_model.floodplain)
                floodplain_v = river_flow_model.floodplain.variables
                floodplain_p = river_flow_model.floodplain.parameters
                q_src = sum_at(floodplain_v.q, edges_at_node.src[river_cell_idx])
                q_dst = sum_at(floodplain_v.q, edges_at_node.dst[river_cell_idx])
                floodplain_v.storage[river_cell_idx] =
                    floodplain_v.storage[river_cell_idx] + (q_src - q_dst) * dt
                if floodplain_v.storage[river_cell_idx] < 0.0
                    floodplain_v.error[river_cell_idx] =
                        floodplain_v.error[river_cell_idx] +
                        abs(floodplain_v.storage[river_cell_idx])
                    floodplain_v.storage[river_cell_idx] = 0.0
                end
                storage_total =
                    river_v.storage[river_cell_idx] + floodplain_v.storage[river_cell_idx]
                if storage_total > river_p.bankfull_storage[river_cell_idx]
                    flood_storage = storage_total - river_p.bankfull_storage[river_cell_idx]
                    h = flood_depth(
                        floodplain_p.profile,
                        flood_storage,
                        flow_length[river_cell_idx],
                        river_cell_idx,
                    )
                    river_v.h[river_cell_idx] = river_p.bankfull_depth[river_cell_idx] + h
                    river_v.storage[river_cell_idx] =
                        river_v.h[river_cell_idx] *
                        flow_width[river_cell_idx] *
                        flow_length[river_cell_idx]
                    floodplain_v.storage[river_cell_idx] =
                        max(storage_total - river_v.storage[river_cell_idx], 0.0)
                    floodplain_v.h[river_cell_idx] =
                        floodplain_v.storage[river_cell_idx] > 0.0 ? h : 0.0
                else
                    river_v.h[river_cell_idx] =
                        storage_total /
                        (flow_length[river_cell_idx] * flow_width[river_cell_idx])
                    river_v.storage[river_cell_idx] = storage_total
                    floodplain_v.h[river_cell_idx] = 0.0
                    floodplain_v.storage[river_cell_idx] = 0.0
                end
            end
        end
    end
    return nothing
end

"""
Update local inertial river flow model `LocalInertialRiverFlowModel` for a single timestep `dt`. An adaptive
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
    n_land_cells::Int
    # flow in y direction at edge at previous time step [m³ s⁻¹]
    qy0::Vector{Float64} = zeros(n_land_cells + 1)
    # flow in x direction at edge at previous time step [m³ s⁻¹]
    qx0::Vector{Float64} = zeros(n_land_cells + 1)
    # flow in x direction at edge [m³ s⁻¹]
    qx::Vector{Float64} = zeros(n_land_cells + 1)
    # average flow in x direction at edge [m³ s⁻¹] for model timestep Δt
    qx_av::Vector{Float64} = zeros(n_land_cells + 1)
    # flow in y direction at edge [m³ s⁻¹]
    qy::Vector{Float64} = zeros(n_land_cells + 1)
    # average flow in y direction at edge [m³ s⁻¹] for model timestep Δt
    qy_av::Vector{Float64} = zeros(n_land_cells + 1)
    # total storage of cell [m³] (including river storage for river cells)
    storage::Vector{Float64} = zeros(n_land_cells)
    # error storage [m³]
    error::Vector{Float64} = zeros(n_land_cells)
    # water depth of cell [m] (for river cells the reference is the river bed elevation `zb`)
    h::Vector{Float64} = zeros(n_land_cells)
end

"Struct to store local inertial overland flow model parameters"
@with_kw struct LocalInertialOverlandFlowParameters
    n_land_cells::Int                   # number of cells [-]
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
    n_land_cells = length(domain.land.network.land_indices_2d)

    zx_max = zeros(n_land_cells)
    zy_max = zeros(n_land_cells)
    for land_cell_idx in 1:n_land_cells
        xu = edge_indices.xu[land_cell_idx]
        if xu <= n_land_cells
            zx_max[land_cell_idx] = max(elevation[land_cell_idx], elevation[xu])
        end
        yu = edge_indices.yu[land_cell_idx]
        if yu <= n_land_cells
            zy_max[land_cell_idx] = max(elevation[land_cell_idx], elevation[yu])
        end
    end

    # set the effective flow width for river cells in the x and y direction at cell edges.
    # for reservoir cells, h is set to zero (fixed) and not updated, and overland flow from
    # a downstream cell is not possible (effective flowwidth is zero).
    we_x = copy(x_length)
    we_y = copy(y_length)
    set_effective_flowwidth!(we_x, we_y, domain)
    parameters = LocalInertialOverlandFlowParameters(;
        n_land_cells,
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
    n_land_cells::Int
    runoff::Vector{Float64} = zeros(n_land_cells) # runoff from hydrological model [m³ s⁻¹]
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

    n_land_cells = length(domain.land.network.land_indices_2d)
    boundary_conditions = LocalInertialOverlandFlowBC(; n_land_cells)
    parameters = LocalInertialOverlandFlowParameters(dataset, config, domain)
    variables = LocalInertialOverlandFlowVariables(; n_land_cells)

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
    (; n_river_cells) = river_flow_model.parameters
    (; h) = river_flow_model.variables
    @batch per = thread reduction = ((min, dt_min),) for river_cell_idx in 1:n_river_cells
        @fastmath @inbounds dt =
            alpha_coefficient * flow_length[river_cell_idx] /
            sqrt(GRAVITATIONAL_ACCELERATION * h[river_cell_idx])
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
    (; n_land_cells) = overland_flow_model.parameters
    (; x_length, y_length, river_location) = parameters
    (; h) = overland_flow_model.variables
    @batch per = thread reduction = ((min, dt_min),) for land_cell_idx in 1:n_land_cells
        @fastmath @inbounds dt = if river_location[land_cell_idx] == 0
            alpha_coefficient * min(x_length[land_cell_idx], y_length[land_cell_idx]) /
            sqrt(GRAVITATIONAL_ACCELERATION * h[land_cell_idx])
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
    river_indices = domain.river.network.land_cell_indices_containing_river

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
    (; land_cell_indices_containing_reservoir) = network
    (; inflow_subsurface) = reservoir_model.boundary_conditions
    inflow_subsurface .= get_inflow_reservoir(
        river_flow_model,
        subsurface_flow_model,
        land_cell_indices_containing_reservoir,
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
Update fluxes for overland flow `LocalInertialOverlandFlowModel` model for a single timestep
`dt`.
"""
function local_inertial_update_fluxes!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    domain::Domain,
    dt::Float64,
)
    indices = domain.land.network.edge_indices
    (; x_length, y_length) = domain.land.parameters
    land_v = overland_flow_model.variables
    land_p = overland_flow_model.parameters

    land_v.qx0 .= land_v.qx
    land_v.qy0 .= land_v.qy

    # update qx
    @batch per = thread minbatch = 6000 for land_cell_idx in 1:(land_p.n_land_cells)
        yu = indices.yu[land_cell_idx]
        yd = indices.yd[land_cell_idx]
        xu = indices.xu[land_cell_idx]
        xd = indices.xd[land_cell_idx]

        # the effective flow width is zero when the river width exceeds the cell width (dy
        # for flow in x dir) and floodplain flow is not calculated.
        if xu <= land_p.n_land_cells && land_p.ywidth[land_cell_idx] != 0.0
            zs_x = land_p.z[land_cell_idx] + land_v.h[land_cell_idx]
            zs_xu = land_p.z[xu] + land_v.h[xu]
            zs_max = max(zs_x, zs_xu)
            hf = (zs_max - land_p.zx_max[land_cell_idx])

            if hf > land_p.h_thresh
                length = 0.5 * (x_length[land_cell_idx] + x_length[xu]) # can be precalculated
                land_v.qx[land_cell_idx] = local_inertial_flow(
                    land_p.theta,
                    land_v.qx0[land_cell_idx],
                    land_v.qx0[xd],
                    land_v.qx0[xu],
                    zs_x,
                    zs_xu,
                    hf,
                    land_p.ywidth[land_cell_idx],
                    length,
                    land_p.mannings_n_sq[land_cell_idx],
                    land_p.froude_limit,
                    dt,
                )
                # limit qx in case water is not available
                if land_v.h[land_cell_idx] <= 0.0
                    land_v.qx[land_cell_idx] = min(land_v.qx[land_cell_idx], 0.0)
                end
                if land_v.h[xu] <= 0.0
                    land_v.qx[land_cell_idx] = max(land_v.qx[land_cell_idx], 0.0)
                end
            else
                land_v.qx[land_cell_idx] = 0.0
            end
            land_v.qx_av[land_cell_idx] += land_v.qx[land_cell_idx] * dt
        end

        # update qy

        # the effective flow width is zero when the river width exceeds the cell width (dx
        # for flow in y dir) and floodplain flow is not calculated.
        if yu <= land_p.n_land_cells && land_p.xwidth[land_cell_idx] != 0.0
            zs_y = land_p.z[land_cell_idx] + land_v.h[land_cell_idx]
            zs_yu = land_p.z[yu] + land_v.h[yu]
            zs_max = max(zs_y, zs_yu)
            hf = (zs_max - land_p.zy_max[land_cell_idx])

            if hf > land_p.h_thresh
                length = 0.5 * (y_length[land_cell_idx] + y_length[yu]) # can be precalculated
                land_v.qy[land_cell_idx] = local_inertial_flow(
                    land_p.theta,
                    land_v.qy0[land_cell_idx],
                    land_v.qy0[yd],
                    land_v.qy0[yu],
                    zs_y,
                    zs_yu,
                    hf,
                    land_p.xwidth[land_cell_idx],
                    length,
                    land_p.mannings_n_sq[land_cell_idx],
                    land_p.froude_limit,
                    dt,
                )
                # limit qy in case water is not available
                if land_v.h[land_cell_idx] <= 0.0
                    land_v.qy[land_cell_idx] = min(land_v.qy[land_cell_idx], 0.0)
                end
                if land_v.h[yu] <= 0.0
                    land_v.qy[land_cell_idx] = max(land_v.qy[land_cell_idx], 0.0)
                end
            else
                land_v.qy[land_cell_idx] = 0.0
            end
            land_v.qy_av[land_cell_idx] += land_v.qy[land_cell_idx] * dt
        end
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
    reservoir_indices = domain.reservoir.network.land_cell_indices_containing_reservoir
    land_bc = overland_flow_model.boundary_conditions
    land_v = overland_flow_model.variables

    for (reservoir_idx, land_cell_idx) in enumerate(reservoir_indices)
        yd = indices.yd[land_cell_idx]
        xd = indices.xd[land_cell_idx]
        reservoir_model.boundary_conditions.inflow_overland[reservoir_idx] =
            land_bc.runoff[land_cell_idx] + (
                land_v.qx[xd] - land_v.qx[land_cell_idx] + land_v.qy[yd] -
                land_v.qy[land_cell_idx]
            )
    end
    return nothing
end

"""
Update storage and water depth for combined river `LocalInertialRiverFlowModel`and overland flow
`LocalInertialOverlandFlowModel` models for a single timestep `dt`.
"""
function local_inertial_update_water_depth!(
    overland_flow_model::LocalInertialOverlandFlowModel,
    river_flow_model::LocalInertialRiverFlowModel,
    domain::Domain,
    dt::Float64,
)
    indices = domain.land.network.edge_indices
    inds_river = domain.land.network.land_cell_river_indices
    (; edges_at_node) = domain.river.network
    (; river_location, reservoir_outlet, x_length, y_length) = domain.land.parameters
    (; flow_width, flow_length) = domain.river.parameters

    river_bc = river_flow_model.boundary_conditions
    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters
    land_bc = overland_flow_model.boundary_conditions
    land_v = overland_flow_model.variables
    land_p = overland_flow_model.parameters

    # change in storage and water depth based on horizontal fluxes for river and land cells
    @batch per = thread minbatch = 6000 for land_cell_idx in 1:(land_p.n_land_cells)
        yd = indices.yd[land_cell_idx]
        xd = indices.xd[land_cell_idx]

        if river_location[land_cell_idx]
            # reservoir locations are boundary points (update storage and h not required)
            reservoir_outlet[land_cell_idx] && continue
            # internal abstraction (water demand) is limited by river storage and negative
            # external inflow as part of water allocation computations.
            land_v.storage[land_cell_idx] +=
                (
                    sum_at(river_v.q, edges_at_node.src[inds_river[land_cell_idx]]) -
                    sum_at(river_v.q, edges_at_node.dst[inds_river[land_cell_idx]]) +
                    land_v.qx[xd] - land_v.qx[land_cell_idx] + land_v.qy[yd] -
                    land_v.qy[land_cell_idx] + land_bc.runoff[land_cell_idx] -
                    river_bc.abstraction[inds_river[land_cell_idx]]
                ) * dt
            if land_v.storage[land_cell_idx] < 0.0
                land_v.error[land_cell_idx] =
                    land_v.error[land_cell_idx] + abs(land_v.storage[land_cell_idx])
                land_v.storage[land_cell_idx] = 0.0 # set storage to zero
            end
            # limit negative external inflow
            if river_bc.external_inflow[inds_river[land_cell_idx]] < 0.0
                available_volume =
                    if land_v.storage[land_cell_idx] >=
                       river_p.bankfull_storage[inds_river[land_cell_idx]]
                        river_p.bankfull_depth[inds_river[land_cell_idx]]
                    else
                        river_v.storage[inds_river[land_cell_idx]]
                    end
                _abstraction = min(
                    -river_bc.external_inflow[inds_river[land_cell_idx]],
                    available_volume / dt * 0.80,
                )
                river_bc.actual_external_abstraction_av[inds_river[land_cell_idx]] +=
                    _abstraction * dt
                _inflow = -_abstraction
            else
                _inflow = river_bc.external_inflow[inds_river[land_cell_idx]]
            end
            land_v.storage[land_cell_idx] += _inflow * dt
            if land_v.storage[land_cell_idx] >=
               river_p.bankfull_storage[inds_river[land_cell_idx]]
                river_v.h[inds_river[land_cell_idx]] =
                    river_p.bankfull_depth[inds_river[land_cell_idx]] +
                    (
                        land_v.storage[land_cell_idx] -
                        river_p.bankfull_storage[inds_river[land_cell_idx]]
                    ) / (x_length[land_cell_idx] * y_length[land_cell_idx])
                land_v.h[land_cell_idx] =
                    river_v.h[inds_river[land_cell_idx]] -
                    river_p.bankfull_depth[inds_river[land_cell_idx]]
                river_v.storage[inds_river[land_cell_idx]] =
                    river_v.h[inds_river[land_cell_idx]] *
                    flow_length[inds_river[land_cell_idx]] *
                    flow_width[inds_river[land_cell_idx]]
            else
                river_v.h[inds_river[land_cell_idx]] =
                    land_v.storage[land_cell_idx] / (
                        flow_length[inds_river[land_cell_idx]] *
                        flow_width[inds_river[land_cell_idx]]
                    )
                land_v.h[land_cell_idx] = 0.0
                river_v.storage[inds_river[land_cell_idx]] = land_v.storage[land_cell_idx]
            end
        else
            land_v.storage[land_cell_idx] +=
                (
                    land_v.qx[xd] - land_v.qx[land_cell_idx] + land_v.qy[yd] -
                    land_v.qy[land_cell_idx] + land_bc.runoff[land_cell_idx]
                ) * dt
            if land_v.storage[land_cell_idx] < 0.0
                land_v.error[land_cell_idx] =
                    land_v.error[land_cell_idx] + abs(land_v.storage[land_cell_idx])
                land_v.storage[land_cell_idx] = 0.0 # set storage to zero
            end
            land_v.h[land_cell_idx] =
                land_v.storage[land_cell_idx] /
                (x_length[land_cell_idx] * y_length[land_cell_idx])
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
@with_kw struct FloodPlainProfile{N}
    depth::Vector{Float64}        # Flood depth [m]
    storage::Array{Float64, 2}    # Flood storage (cumulative) [m³]
    width::Array{Float64, 2}      # Flood width [m]
    a::Array{Float64, 2}          # Flow area (cumulative) [m²]
    p::Array{Float64, 2}          # Wetted perimeter (cumulative) [m]
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
    n_river_cells = length(river_indices_2d)

    # for convenience (interpolation) flood depth 0.0 m is added, with associated area (a),
    # storage, width (river width) and wetted perimeter (p).
    storage = vcat(zeros(n_river_cells)', storage)
    start_storage = storage
    flood_depths = Float64.(dataset["flood_depth"][:])
    pushfirst!(flood_depths, 0.0)
    n_depths = length(flood_depths)

    p = zeros(n_depths, n_river_cells)
    a = zeros(n_depths, n_river_cells)
    segment_storage = zeros(n_depths, n_river_cells)
    width = zeros(n_depths, n_river_cells)
    width[1, :] = flow_width[1:n_river_cells]

    # determine flow area (a), width and wetted perimeter (p) FloodPlainModel
    h = diff(flood_depths)
    incorrect_vol = 0
    riv_cells = 0
    error_vol = 0
    for river_cell_idx in 1:n_river_cells
        riv_cell = 0
        diff_storage = diff(storage[:, river_cell_idx])

        for flood_depth_idx in 1:(n_depths - 1)
            # assume rectangular shape of flood depth segment
            width[flood_depth_idx + 1, river_cell_idx] =
                diff_storage[flood_depth_idx] /
                (h[flood_depth_idx] * flow_length[river_cell_idx])
            # check provided flood storage (floodplain width should be constant or increasing
            # as a function of flood depth)
            if width[flood_depth_idx + 1, river_cell_idx] <
               width[flood_depth_idx, river_cell_idx]
                # raise warning only if difference is larger than rounding error of 0.01 m³
                if (
                    (
                        width[flood_depth_idx, river_cell_idx] -
                        width[flood_depth_idx + 1, river_cell_idx]
                    ) *
                    h[flood_depth_idx] *
                    flow_length[river_cell_idx]
                ) > 0.01
                    incorrect_vol += 1
                    riv_cell = 1
                    error_vol =
                        error_vol + (
                            (
                                width[flood_depth_idx, river_cell_idx] -
                                width[flood_depth_idx + 1, river_cell_idx]
                            ) *
                            h[flood_depth_idx] *
                            flow_length[river_cell_idx]
                        )
                end
                width[flood_depth_idx + 1, river_cell_idx] =
                    width[flood_depth_idx, river_cell_idx]
            end
            a[flood_depth_idx + 1, river_cell_idx] =
                width[flood_depth_idx + 1, river_cell_idx] * h[flood_depth_idx]
            p[flood_depth_idx + 1, river_cell_idx] =
                (
                    width[flood_depth_idx + 1, river_cell_idx] -
                    width[flood_depth_idx, river_cell_idx]
                ) + 2.0 * h[flood_depth_idx]
            segment_storage[flood_depth_idx + 1, river_cell_idx] =
                a[flood_depth_idx + 1, river_cell_idx] * flow_length[river_cell_idx]
            if flood_depth_idx == 1
                # for interpolation wetted perimeter at flood depth 0.0 is required
                p[flood_depth_idx, river_cell_idx] =
                    p[flood_depth_idx + 1, river_cell_idx] - 2.0 * h[flood_depth_idx]
            end
        end

        p[2:end, river_cell_idx] = cumsum(p[2:end, river_cell_idx])
        a[:, river_cell_idx] = cumsum(a[:, river_cell_idx])
        storage[:, river_cell_idx] = cumsum(segment_storage[:, river_cell_idx])

        riv_cells += riv_cell
    end

    if incorrect_vol > 0
        perc_riv_cells = round(100.0 * (riv_cells / n_river_cells); digits = 2)
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
    profile = FloodPlainProfile{n_depths}(; storage, width, depth = flood_depths, a, p)
    return profile
end

"Struct to store floodplain flow model parameters"
@with_kw struct FloodPlainParameters{P}
    profile::P                          # floodplain profile
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
    n_river_edges = ne(graph)
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
    mannings_n_sq = zeros(n_river_edges)
    zb_max = zeros(n_river_edges)
    for river_edge_idx in 1:n_river_edges
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
    n_river_cells::Int
    n_river_edges::Int
    storage::Vector{Float64} = zeros(n_river_cells)    # storage [m³]
    h::Vector{Float64}                     # water depth [m]
    error::Vector{Float64} = zeros(n_river_cells)      # error storage [m³]
    a::Vector{Float64} = zeros(n_river_edges)    # flow area at edge [m²]
    r::Vector{Float64} = zeros(n_river_edges)    # hydraulic radius at edge [m]
    hf::Vector{Float64} = zeros(n_river_edges)   # water depth at edge [m]
    q0::Vector{Float64} = zeros(n_river_edges)   # discharge at edge at previous time step
    q::Vector{Float64} = zeros(n_river_edges)    # discharge at edge  [m³ s⁻¹]
    q_av::Vector{Float64} = zeros(n_river_edges) # average river discharge at edge  [m³ s⁻¹] for model timestep Δt
    hf_index::Vector{Int} = zeros(Int, n_river_edges) # edge index with `hf` [-] above depth threshold
end

"Initialize floodplain flow model variables"
function FloodPlainVariables(n_river_cells::Int, n_river_edges::Int, index_pit::Vector{Int})
    variables = FloodPlainVariables(;
        n_river_cells,
        n_river_edges,
        h = zeros(n_river_cells + length(index_pit)),
    )
    return variables
end

"Floodplain flow model"
@with_kw struct FloodPlainModel{P} <: AbstractFloodPlainModel
    parameters::FloodPlainParameters{P}
    variables::FloodPlainVariables
end

"Determine the initial floodplain storage"
function initialize_storage!(river, domain::Domain, n_river_cells::Int)
    (; flow_width, flow_length) = domain.river.parameters
    (; floodplain) = river
    (; profile) = floodplain.parameters
    for river_cell_idx in 1:n_river_cells
        i1, i2 =
            interpolation_indices(floodplain.variables.h[river_cell_idx], profile.depth)
        a = flow_area(
            profile.width[i2, river_cell_idx],
            profile.a[i1, river_cell_idx],
            profile.depth[i1],
            floodplain.variables.h[river_cell_idx],
        )
        a = max(
            a - (flow_width[river_cell_idx] * floodplain.variables.h[river_cell_idx]),
            0.0,
        )
        floodplain.variables.storage[river_cell_idx] = flow_length[river_cell_idx] * a
    end
    return nothing
end

"helper function to get interpolation indices"
function interpolation_indices(x, v::AbstractVector)
    i1 = 1
    for idx in eachindex(v)
        if v[idx] <= x
            i1 = idx
        end
    end
    if i1 == length(v)
        i2 = i1
    else
        i2 = i1 + 1
    end
    return i1, i2
end

"""
    flow_area(width, area, depth, h)

Compute floodplain flow area based on flow depth `h` and floodplain `depth`, `area` and
`width` of a floodplain profile.
"""
function flow_area(width, area, depth, h)
    dh = h - depth  # depth at i1
    area = area + (width * dh) # area at i1, width at i2
    return area
end

"""
    function wetted_perimeter(p, depth, h)

Compute floodplain wetted perimeter based on flow depth `h` and floodplain `depth` and
wetted perimeter `p` of a floodplain profile.
"""
function wetted_perimeter(p, depth, h)
    dh = h - depth # depth at i1
    p += 2.0 * dh # p at i1
    return p
end

"Compute flood depth by interpolating flood storage `flood_storage` using flood depth intervals."
function flood_depth(
    profile::FloodPlainProfile,
    flood_storage::Float64,
    flow_length::Float64,
    river_cell_idx::Int,
)
    i1, i2 = interpolation_indices(flood_storage, @view profile.storage[:, river_cell_idx])
    ΔA = (flood_storage - profile.storage[i1, river_cell_idx]) / flow_length
    dh = ΔA / profile.width[i2, river_cell_idx]
    flood_depth = profile.depth[i1] + dh
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
    n_river_cells = length(river_indices_2d)
    index_pit = findall(x -> x == 5, local_drain_direction)
    parameters = FloodPlainParameters(dataset, config, domain, zb_floodplain, index_pit)
    n_edges = ne(graph)
    variables = FloodPlainVariables(n_river_cells, n_edges, index_pit)

    floodplain = FloodPlainModel(; parameters, variables)
    return floodplain
end
