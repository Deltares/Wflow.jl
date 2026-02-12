abstract type AbstractFloodPlain end

"Struct for storing local inertial river flow model parameters"
@with_kw struct LocalInertialRiverFlowParameters
    # number of cells [-]
    n::Int
    # number of edges [-]
    ne::Int
    # active nodes [-]
    active_n::Vector{Int}
    # active edges [-]
    active_e::Vector{Int}
    # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
    froude_limit::Bool
    # depth threshold for calculating flow [m]
    h_thresh::Float64
    # river bed elevation [m]
    zb::Vector{Float64}
    # maximum channel bed elevation [m]
    zb_max::Vector{Float64}
    # bankfull storage [m³]
    bankfull_storage::Vector{Float64}
    # bankfull depth [m]
    bankfull_depth::Vector{Float64}
    # Manning's roughness squared at edge [(s m-1/3)²]
    mannings_n_sq::Vector{Float64}
    # Manning's roughness [s m-1/3] at node
    mannings_n::Vector{Float64}
    # flow (river) length at edge [m]
    flow_length_at_edge::Vector{Float64}
    # flow (river) width at edge [m]
    flow_width_at_edge::Vector{Float64}
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

    (; pit_indices, indices, graph, local_drain_direction, nodes_at_edge) = domain.network
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
    bankfull_depth = bankfull_depth_2d[indices]
    zb = bankfull_elevation_2d[indices] - bankfull_depth # river bed elevation

    bankfull_storage = bankfull_depth .* flow_width .* flow_length
    mannings_n = ncread(
        dataset,
        config,
        "river_water_flow__manning_n_parameter",
        Routing;
        sel = indices,
    )

    n = length(indices)
    index_pit = findall(==(5), local_drain_direction)
    # set ghost points for boundary condition (downstream river outlet): river width, bed
    # elevation, manning n is copied from the upstream cell.
    append!(flow_length, riverlength_bc)
    append!(zb, zb[index_pit])
    append!(flow_width, flow_width[index_pit])
    append!(mannings_n, mannings_n[index_pit])
    append!(bankfull_depth, bankfull_depth[index_pit])

    # determine z, width, length and manning's n at edges
    n_edges = ne(graph)
    zb_max = fill(Float64(0), n_edges)
    width_at_edge = fill(Float64(0), n_edges)
    length_at_edge = fill(Float64(0), n_edges)
    mannings_n_sq = fill(Float64(0), n_edges)
    for i in 1:n_edges
        src_node = nodes_at_edge.src[i]
        dst_node = nodes_at_edge.dst[i]
        zb_max[i] = max(zb[src_node], zb[dst_node])
        width_at_edge[i] = min(flow_width[src_node], flow_width[dst_node])
        length_at_edge[i] = (flow_length[dst_node] + flow_length[src_node]) / 2
        mannings_n_i =
            (
                mannings_n[dst_node] * flow_length[dst_node] +
                mannings_n[src_node] * flow_length[src_node]
            ) / (flow_length[dst_node] + flow_length[src_node])
        mannings_n_sq[i] = mannings_n_i * mannings_n_i
    end
    active_index = findall(==(0), reservoir_outlet)

    parameters = LocalInertialRiverFlowParameters(;
        n,
        ne = n_edges,
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
    # river discharge at edge (subgrid channel) [m³ s⁻¹]
    q::Vector{Float64} = zeros(n_edges)
    # river discharge at edge (subgrid channel) at previous time step [m³ s⁻¹]
    q0::Vector{Float64} = zeros(n_edges)
    # average river channel (+ floodplain) discharge at edge [m³ s⁻¹] (model time step dt)
    q_av::AverageVector = AverageVector(; n = n_edges)
    # average river channel discharge at edge [m³ s⁻¹] (for model time step dt)
    q_channel_av::AverageVector = AverageVector(; n = n_edges)
    # water depth [m]
    h::Vector{Float64}
    # maximum water elevation at edge [m]
    zs_max::Vector{Float64} = zeros(n_edges)
    # water elevation of source node of edge [m]
    zs_src::Vector{Float64} = zeros(n_edges)
    # water elevation of downstream node of edge [m]
    zs_dst::Vector{Float64} = zeros(n_edges)
    # water depth at edge [m]
    hf::Vector{Float64} = zeros(n_edges)
    # flow area at edge [m²]
    a::Vector{Float64} = zeros(n_edges)
    # wetted perimeter at edge [m]
    r::Vector{Float64} = zeros(n_edges)
    # river storage [m³]
    storage::Vector{Float64} = zeros(n_cells)
    # error storage [m³]
    error::Vector{Float64} = zeros(n_cells)
end

"Initialize shallow water river flow model variables"
function LocalInertialRiverFlowVariables(
    dataset::NCDataset,
    config::Config,
    network::NetworkRiver,
)
    (; pit_indices, indices, graph) = network

    riverdepth_bc = ncread(
        dataset,
        config,
        "model_boundary_condition_river_bank_water__depth",
        Routing;
        sel = pit_indices,
    )

    n_cells = length(indices)
    n_edges = ne(graph)
    # set river depth h to zero (including reservoir locations)
    q_av = AverageVector(; n = n_edges)
    # set ghost points for boundary condition (downstream river outlet): river depth `h`
    h = zeros(n_cells)
    append!(h, riverdepth_bc)
    q_channel_av = config.model.floodplain_1d__flag ? AverageVector(; n = n_edges) : q_av
    variables = LocalInertialRiverFlowVariables(; n_cells, n_edges, q_av, q_channel_av, h)
    return variables
end

"Shallow water river flow model using the local inertial method"
@with_kw struct LocalInertialRiverFlow{
    R <: RiverFlowBC,
    F <: Union{AbstractFloodPlain, Nothing},
    A <: AbstractAllocationModel,
} <: AbstractRiverFlowModel
    timestepping::TimeStepping
    boundary_conditions::R
    parameters::LocalInertialRiverFlowParameters
    variables::LocalInertialRiverFlowVariables
    # Floodplain (1D) schematization
    floodplain::F
    # Water allocation
    allocation::A
end

"Initialize shallow water river flow model `LocalIntertialRiverFlow`"
function LocalInertialRiverFlow(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    reservoir::Union{Reservoir, Nothing},
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
    cfl = config.model.river_local_inertial_flow__alpha_coefficient # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; cfl)

    parameters = LocalInertialRiverFlowParameters(dataset, config, domain)
    variables = LocalInertialRiverFlowVariables(dataset, config, domain.network)
    boundary_conditions = RiverFlowBC(dataset, config, domain.network, reservoir)

    if config.model.floodplain_1d__flag
        zb_floodplain = parameters.zb .+ parameters.bankfull_depth
        floodplain = FloodPlain(dataset, config, domain, zb_floodplain)
    else
        floodplain = nothing
    end

    n = length(domain.network.indices)
    river_flow = LocalInertialRiverFlow(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain,
        allocation = do_water_demand(config) ? AllocationRiver(n) : NoAllocationRiver(n),
    )
    return river_flow
end

"Return the upstream inflow for a reservoir in `LocalInertialRiverFlow`"
function get_inflow_reservoir(model::LocalInertialRiverFlow, src_edge::Vector{Int})
    q_in = sum_at(model.variables.q, src_edge)
    if !isnothing(model.floodplain)
        q_in += sum_at(model.floodplain.variables.q, src_edge)
    end
    return q_in
end

# For local inertial river routing, `to_river` is included, as reservoir cells are excluded
# (boundary condition).
get_inflow_reservoir(
    ::LocalInertialRiverFlow,
    model::KinWaveOverlandFlow,
    inds::Vector{Int},
) = get_average(model.variables.q_av)[inds] .+ get_average(model.variables.to_river)[inds]

get_inflow_reservoir(::LocalInertialRiverFlow, model::LateralSSF, inds::Vector{Int}) =
    (model.variables.ssf[inds] .+ model.variables.to_river[inds])

"""
Update river channel flow for the local inertial river flow model.
"""
function update_river_channel_flow!(
    model::LocalInertialRiverFlow,
    domain::DomainRiver,
    dt::Float64,
)
    (; nodes_at_edge) = domain.network
    river_v = model.variables
    river_p = model.parameters

    # [m³ s⁻¹] = [m³ s⁻¹]
    river_v.q0 .= river_v.q
    if !isnothing(model.floodplain)
        # [m³ s⁻¹] = [m³ s⁻¹]
        model.floodplain.variables.q0 .= model.floodplain.variables.q
    end

    @batch per = thread minbatch = 1000 for j in eachindex(river_p.active_e)
        i = river_p.active_e[j]
        i_src = nodes_at_edge.src[i]
        i_dst = nodes_at_edge.dst[i]
        river_v.zs_src[i] = river_p.zb[i_src] + river_v.h[i_src]
        river_v.zs_dst[i] = river_p.zb[i_dst] + river_v.h[i_dst]

        river_v.zs_max[i] = max(river_v.zs_src[i], river_v.zs_dst[i])
        river_v.hf[i] = (river_v.zs_max[i] - river_p.zb_max[i])

        # [m²] = [m] * [m]
        river_v.a[i] = river_p.flow_width_at_edge[i] * river_v.hf[i] # flow area (rectangular channel)
        # [m] = [m²] / ([m] + [-] * [m])
        river_v.r[i] = river_v.a[i] / (river_p.flow_width_at_edge[i] + 2 * river_v.hf[i]) # hydraulic radius (rectangular channel)

        river_v.q[i] = ifelse(
            river_v.hf[i] > river_p.h_thresh,
            local_inertial_flow(
                river_v.q0[i],
                river_v.zs_src[i],
                river_v.zs_dst[i],
                river_v.hf[i],
                river_v.a[i],
                river_v.r[i],
                river_p.flow_length_at_edge[i],
                river_p.mannings_n_sq[i],
                river_p.froude_limit,
                dt,
            ),
            0.0,
        )

        # limit q in case water is not available
        river_v.q[i] = ifelse(river_v.h[i_src] <= 0.0, min(river_v.q[i], 0.0), river_v.q[i])
        river_v.q[i] = ifelse(river_v.h[i_dst] <= 0.0, max(river_v.q[i], 0.0), river_v.q[i])
        # average river discharge (here accumulated for model timestep dt)
        add_to_cumulative!(river_v.q_av, i, river_v.q[i], dt)
    end
    return nothing
end

"""
Update floodplain flow for the local inertial river flow model.
"""
function update_floodplain_flow!(
    model::LocalInertialRiverFlow{R, F},
    domain::DomainRiver,
    dt::Float64,
) where {R, F <: AbstractFloodPlain}
    (; nodes_at_edge) = domain.network
    (; flow_width) = domain.parameters

    river_v = model.variables
    river_p = model.parameters
    floodplain_p = model.floodplain.parameters
    floodplain_v = model.floodplain.variables

    @batch per = thread minbatch = 1000 for i in 1:length(floodplain_v.hf)
        floodplain_v.hf[i] = max(river_v.zs_max[i] - floodplain_p.zb_max[i], 0.0)
    end

    n = 0
    @inbounds for i in river_p.active_e
        @inbounds if river_v.hf[i] > river_p.h_thresh
            n += 1
            floodplain_v.hf_index[n] = i
        else
            floodplain_v.q[i] = 0.0
        end
    end

    get_area(i, i1, i2, idx) = flow_area(
        floodplain_p.profile.width[i2, idx],
        floodplain_p.profile.a[i1, idx],
        floodplain_p.profile.depth[i1],
        floodplain_v.hf[i],
    )

    get_wetted_perimeter(i, i1, idx) = wetted_perimeter(
        floodplain_p.profile.p[i1, idx],
        floodplain_p.profile.depth[i1],
        floodplain_v.hf[i],
    )

    @batch per = thread minbatch = 1000 for j in 1:n
        i = floodplain_v.hf_index[j]
        i_src = nodes_at_edge.src[i]
        i_dst = nodes_at_edge.dst[i]

        i0 = 0
        for k in eachindex(floodplain_p.profile.depth)
            i0 += 1 * (floodplain_p.profile.depth[k] <= floodplain_v.hf[i])
        end
        i1 = max(i0, 1)
        i2 = ifelse(i1 == length(floodplain_p.profile.depth), i1, i1 + 1)

        # [m²]
        a_src = get_area(i, i1, i2, i_src)
        a_src = max(a_src - (floodplain_v.hf[i] * flow_width[i_src]), 0.0)

        # [m²]
        a_dst = get_area(i, i1, i2, i_dst)
        a_dst = max(a_dst - (floodplain_v.hf[i] * flow_width[i_dst]), 0.0)

        # [m²]
        floodplain_v.a[i] = min(a_src, a_dst)

        # [m]
        floodplain_v.r[i] = if a_src < a_dst
            # [m²] / [m]
            a_src / get_wetted_perimeter(i, i1, i_src)
        else
            # [m²] / [m]
            a_dst / get_wetted_perimeter(i, i1, i_dst)
        end

        # [m³ s⁻¹]
        floodplain_v.q[i] = if floodplain_v.a[i] > 1.0e-05
            local_inertial_flow(
                floodplain_v.q0[i],
                river_v.zs_src[i],
                river_v.zs_dst[i],
                floodplain_v.hf[i],
                floodplain_v.a[i],
                floodplain_v.r[i],
                river_p.flow_length_at_edge[i],
                floodplain_p.mannings_n_sq[i],
                river_p.froude_limit,
                dt,
            )
        else
            0.0
        end

        # limit floodplain q in case water is not available
        if floodplain_v.h[i_src] <= 0.0
            floodplain_v.q[i] = min(floodplain_v.q[i], 0.0)
        end

        if floodplain_v.h[i_dst] <= 0.0
            floodplain_v.q[i] = max(floodplain_v.q[i], 0.0)
        end

        # If the directions are opposite
        if floodplain_v.q[i] * river_v.q[i] < 0.0
            floodplain_v.q[i] = 0.0
        end

        # average floodplain discharge (here accumulated for model timestep dt)
        add_to_cumulative!(floodplain_v.q_av, i, floodplain_v.q[i], dt)
    end
    return nothing
end

update_floodplain_flow!(
    model::LocalInertialRiverFlow{R, F},
    domain::DomainRiver,
    dt::Float64,
) where {R, F <: Nothing} = nothing

"""
Update reservoir boundary conditions for the local inertial river flow model.
"""
function update_boundary_conditions_reservoir!(
    model::LocalInertialRiverFlow,
    domain::Domain,
    dt::Float64,
)
    (; edges_at_node) = domain.river.network
    (; reservoir) = model.boundary_conditions
    inds_reservoir = domain.reservoir.network.river_indices
    isnothing(reservoir) && return nothing

    river_v = model.variables
    res_bc = reservoir.boundary_conditions

    for v in eachindex(inds_reservoir)
        i = inds_reservoir[v]

        q_in = get_inflow_reservoir(model, edges_at_node.src[i])
        # If external_inflow < 0, abstraction is limited
        if res_bc.external_inflow[v] < 0.0
            # [m³ s⁻¹] = min([m³ s⁻¹], [m³] / [s])
            abstraction = min(
                -res_bc.external_inflow[v],
                (reservoir.variables.storage[v] / dt) * 0.98,
            )
            add_to_cumulative!(res_bc.actual_external_abstraction_av, v, abstraction, dt)
            # [m³ s⁻¹] = [m³ s⁻¹]
            inflow = -abstraction
        else
            # [m³ s⁻¹] = [m³ s⁻¹]
            inflow = res_bc.external_inflow[v]
        end
        # [m³ s⁻¹] = ∑ [m³ s⁻¹]
        net_inflow = q_in + res_bc.inflow_overland[v] + res_bc.inflow_subsurface[v] + inflow
        update!(reservoir, v, net_inflow, dt)
        river_v.q[i] = reservoir.variables.outflow[v]
        # average river discharge (here accumulated for model timestep dt)
        add_to_cumulative!(river_v.q_av, i, river_v.q[i], dt)
    end
    return nothing
end

"""
Update floodplain water depth and storage.
"""
function update_water_depth_and_storage!(
    floodplain::AbstractFloodPlain,
    model::LocalInertialRiverFlow,
    domain::DomainRiver,
    dt::Float64,
)
    (; edges_at_node) = domain.network
    (; flow_length, flow_width) = domain.parameters

    river_v = model.variables
    river_p = model.parameters
    floodplain_v = floodplain.variables
    floodplain_p = floodplain.parameters

    @batch per = thread minbatch = 1000 for i in river_p.active_n
        # [m³ s⁻¹]
        q_src = sum_at(floodplain_v.q, edges_at_node.src[i])
        q_dst = sum_at(floodplain_v.q, edges_at_node.dst[i])
        # [m³] += ([m³ s⁻¹] - [m³ s⁻¹]) * [s]
        floodplain_v.storage[i] += (q_src - q_dst) * dt
        if floodplain_v.storage[i] < 0.0
            floodplain_v.error[i] += abs(floodplain_v.storage[i])
            floodplain_v.storage[i] = 0.0
        end
        # [m³] = [m³] + [m³]
        storage_total = river_v.storage[i] + floodplain_v.storage[i]
        if storage_total > river_p.bankfull_storage[i]
            # [m³] = [m³] - [m³]
            flood_storage = storage_total - river_p.bankfull_storage[i]
            # [m]
            h = flood_depth(floodplain_p.profile, flood_storage, flow_length[i], i)
            # [m] = [m] + [m]
            river_v.h[i] = river_p.bankfull_depth[i] + h
            # [m³] = [m] * [m] * [m]
            river_v.storage[i] = river_v.h[i] * flow_width[i] * flow_length[i]
            # [m³] = max([m³] - [m³], [m³])
            floodplain_v.storage[i] = max(storage_total - river_v.storage[i], 0.0)
            floodplain_v.h[i] = floodplain_v.storage[i] > 0.0 ? h : 0.0
        else
            # [m] = [m³] / ([m] * [m])
            river_v.h[i] = storage_total / (flow_length[i] * flow_width[i])
            river_v.storage[i] = storage_total
            floodplain_v.h[i] = 0.0
            floodplain_v.storage[i] = 0.0
        end
    end
    return nothing
end

"""
Update floodplain water depth and storage (no-op for Nothing floodplain).
"""
update_water_depth_and_storage!(
    floodplain::Nothing,
    model::LocalInertialRiverFlow,
    domain::DomainRiver,
    dt::Float64,
) = nothing

"""
Update water depth and storage for river.
"""
function update_water_depth_and_storage!(
    model::LocalInertialRiverFlow,
    domain::DomainRiver,
    dt::Float64,
)
    (; edges_at_node) = domain.network
    (; flow_length, flow_width) = domain.parameters
    (; inwater, abstraction, external_inflow, actual_external_abstraction_av) =
        model.boundary_conditions

    river_v = model.variables
    river_p = model.parameters

    @batch per = thread minbatch = 1000 for i in river_p.active_n
        # [m³ s⁻¹]
        q_src = sum_at(river_v.q, edges_at_node.src[i])
        q_dst = sum_at(river_v.q, edges_at_node.dst[i])
        # internal abstraction (water demand) is limited by river storage and negative
        # external inflow as part of water allocation computations.
        # [m³] = ([m³ s⁻¹] - [m³ s⁻¹] + [m³ s⁻¹] - [m³ s⁻¹]) * [s]
        river_v.storage[i] += (q_src - q_dst + inwater[i] - abstraction[i]) * dt

        if river_v.storage[i] < 0.0
            river_v.error[i] = river_v.error[i] + abs(river_v.storage[i])
            river_v.storage[i] = 0.0 # set storage to zero
        end
        # limit negative external inflow
        if external_inflow[i] < 0.0
            # [m³ s⁻¹] = min([m³ s⁻¹], [m³] / [s] * [-])
            abstraction = min(-external_inflow[i], river_v.storage[i] / dt * 0.80)
            # [m³] += [m³ s⁻¹] * [s]
            add_to_cumulative!(actual_external_abstraction_av, i, abstraction, dt)
            # [m³ s⁻¹] = [m³ s⁻¹]
            inflow = -abstraction
        else
            # [m³ s⁻¹] = [m³ s⁻¹]
            inflow = external_inflow[i]
        end
        # [m³] += [m³ s⁻¹] * [s]
        river_v.storage[i] += inflow * dt # add external inflow
        # [m] = [m³] / ([m] * [m])
        river_v.h[i] = river_v.storage[i] / (flow_length[i] * flow_width[i])
    end
    return nothing
end

"Update local inertial river flow model `LocalIntertialRiverFlow` for a single timestep"
function local_inertial_river_update!(
    model::LocalInertialRiverFlow,
    domain::Domain,
    dt::Float64,
    update_h::Bool,
)
    # Update river channel flow
    update_river_channel_flow!(model, domain.river, dt)

    # Update floodplain flow if present
    update_floodplain_flow!(model, domain.river, dt)

    # Handle reservoir boundary conditions
    update_boundary_conditions_reservoir!(model, domain, dt)

    # Update water depth and storage if requested
    if update_h
        update_water_depth_and_storage!(model, domain.river, dt)
        update_water_depth_and_storage!(model.floodplain, model, domain.river, dt)
    end

    return nothing
end

"""
Update local inertial river flow model `LocalInertialRiverFlow` for a single timestep `dt`. An adaptive
timestepping method is used (computing a sub timestep `dt_s`).
"""
function update!(
    model::LocalInertialRiverFlow,
    domain::Domain,
    clock::Clock,
    dt::Number;
    update_h = true,
)
    (; reservoir) = model.boundary_conditions
    (; flow_length) = domain.river.parameters

    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)

    if !isnothing(model.floodplain)
        zero!(model.floodplain.variables.q_av)
    end
    set_flow_vars!(model)

    t = 0.0
    while t < dt
        dt_s = stable_timestep(model, flow_length)
        dt_s = check_timestepsize(dt_s, t, dt)
        local_inertial_river_update!(model, domain, dt_s, update_h)
        t += dt_s
    end
    average_flow_vars!(model, dt)
    average_reservoir_vars!(reservoir, dt)

    if !isnothing(model.floodplain)
        average!(model.floodplain.variables.q_av, dt)
        set_average!(model.variables.q_channel_av, get_average(model.variables.q_av))
        set_average!(
            model.variables.q_av,
            get_average(model.variables.q_channel_av) .+
            get_average(model.floodplain.variables.q_av),
        )
    end

    return nothing
end

"Struct to store local inertial overland flow model variables"
@with_kw struct LocalInertialOverlandFlowVariables
    n::Int
    # flow in y direction at edge at previous time step [m³ s⁻¹]
    qy0::Vector{Float64} = zeros(n + 1)
    # flow in x direction at edge at previous time step [m³ s⁻¹]
    qx0::Vector{Float64} = zeros(n + 1)
    # flow in x direction at edge [m³ s⁻¹]
    qx::Vector{Float64} = zeros(n + 1)
    # average flow in x direction at edge [m³ s⁻¹] for model timestep Δt
    qx_av::AverageVector = AverageVector(; n = n + 1)
    # flow in y direction at edge [m³ s⁻¹]
    qy::Vector{Float64} = zeros(n + 1)
    # average flow in y direction at edge [m³ s⁻¹] for model timestep Δt
    qy_av::AverageVector = AverageVector(; n = n + 1)
    # total storage of cell [m³] (including river storage for river cells)
    storage::Vector{Float64} = zeros(n)
    # error storage [m³]
    error::Vector{Float64} = zeros(n)
    # water depth of cell [m] (for river cells the reference is the river bed elevation `zb`)
    h::Vector{Float64} = zeros(n)
end

"Struct to store local inertial overland flow model parameters"
@with_kw struct LocalInertialOverlandFlowParameters
    # number of cells [-]
    n::Int
    # effective flow width x direction at edge (floodplain) [m]
    xwidth::Vector{Float64}
    # effective flow width y direction at edge (floodplain) [m]
    ywidth::Vector{Float64}
    # acceleration due to gravity [m s⁻²]
    g::Float64 = 9.80665
    # weighting factor (de Almeida et al., 2012) [-]
    theta::Float64
    # depth threshold for calculating flow [m]
    h_thresh::Float64
    # maximum cell elevation at edge [m] (x direction)
    zx_max::Vector{Float64}
    # maximum cell elevation at edge [m] (y direction)
    zy_max::Vector{Float64}
    # Manning's roughness squared at edge [(s m-1/3)2]
    mannings_n_sq::Vector{Float64}
    # elevation [m] of cell
    z::Vector{Float64}
    # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
    froude_limit::Bool
end

"Initialize shallow water overland flow model parameters"
function LocalInertialOverlandFlowParameters(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
)
    # limit flow to subcritical according to Froude number
    froude_limit = config.model.land_surface_water_flow__froude_limit_flag
    # stability coefficient for model time step (0.2-0.7)
    alpha = config.model.land_local_inertial_flow__alpha_coefficient
    # weighting factor
    theta = config.model.land_local_inertial_flow__theta_coefficient
    # depth threshold for flow at edge
    waterdepth_threshold = config.model.land_surface_water_flow_threshold__depth

    (; edge_indices, indices) = domain.land.network
    (; x_length, y_length) = domain.land.parameters

    @info "Local inertial approach is used for overland flow." alpha theta waterdepth_threshold froude_limit

    mannings_n = ncread(
        dataset,
        config,
        "land_surface_water_flow__manning_n_parameter",
        Routing;
        sel = indices,
    )
    elevation_2d =
        ncread(dataset, config, "land_surface_water_flow__ground_elevation", Routing)
    elevation = elevation_2d[indices]
    n = length(domain.land.network.indices)

    zx_max = zeros(n)
    zy_max = zeros(n)
    for i in 1:n
        xu = edge_indices.xu[i]
        if xu <= n
            zx_max[i] = max(elevation[i], elevation[xu])
        end
        yu = edge_indices.yu[i]
        if yu <= n
            zy_max[i] = max(elevation[i], elevation[yu])
        end
    end

    # set the effective flow width for river cells in the x and y direction at cell edges.
    # for reservoir cells, h is set to zero (fixed) and not updated, and overland flow from
    # a downstream cell is not possible (effective flowwidth is zero).
    we_x = copy(x_length)
    we_y = copy(y_length)
    set_effective_flowwidth!(we_x, we_y, domain)
    parameters = LocalInertialOverlandFlowParameters(;
        n,
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
    n::Int
    # runoff from hydrological model [m³ s⁻¹]
    runoff::Vector{Float64} = zeros(n)
end

"Local inertial overland flow model using the local inertial method"
@with_kw struct LocalInertialOverlandFlow <: AbstractOverlandFlowModel
    timestepping::TimeStepping
    boundary_conditions::LocalInertialOverlandFlowBC
    parameters::LocalInertialOverlandFlowParameters
    variables::LocalInertialOverlandFlowVariables
end

"Initialize local inertial overland flow model"
function LocalInertialOverlandFlow(dataset::NCDataset, config::Config, domain::Domain)
    cfl = config.model.land_local_inertial_flow__alpha_coefficient # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; cfl)

    n = length(domain.land.network.indices)
    boundary_conditions = LocalInertialOverlandFlowBC(; n)
    parameters = LocalInertialOverlandFlowParameters(dataset, config, domain)
    variables = LocalInertialOverlandFlowVariables(; n)

    overland_flow = LocalInertialOverlandFlow(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
    )

    return overland_flow
end

"""
    stable_timestep(model::LocalInertialRiverFlow, flow_length::Vector{Float64})
    stable_timestep(model::LocalInertialOverlandFlow, parameters::LandParameters)

Compute a stable timestep size for the local inertial approach, based on Bates et al. (2010).

dt = cfl * (Δx / sqrt(g max(h))
"""
function stable_timestep(model::LocalInertialRiverFlow, flow_length::Vector{Float64})
    dt_min = Inf
    (; cfl) = model.timestepping
    (; n) = model.parameters
    (; h) = model.variables
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(n)
        @fastmath @inbounds dt =
            cfl * flow_length[i] / sqrt(GRAVITATIONAL_ACCELERATION * h[i])
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? 60.0 : dt_min
    return dt_min
end

function stable_timestep(model::LocalInertialOverlandFlow, parameters::LandParameters)
    dt_min = Inf
    (; cfl) = model.timestepping
    (; n) = model.parameters
    (; x_length, y_length, river_location) = parameters
    (; h) = model.variables
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(n)
        @fastmath @inbounds dt = if river_location[i] == 0
            cfl * min(x_length[i], y_length[i]) / sqrt(GRAVITATIONAL_ACCELERATION * h[i])
        else
            Inf
        end
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? 60.0 : dt_min
    return dt_min
end

"""
Update boundary condition `runoff` overland flow model `LocalInertialOverlandFlow` for a
single timestep.
"""
function update_boundary_conditions!(
    model::LocalInertialOverlandFlow,
    external_models::NamedTuple,
    domain::Domain,
    dt::Float64,
)
    (; soil, runoff, subsurface_flow) = external_models
    (; net_runoff) = soil.variables
    (; net_runoff_river) = runoff.variables
    (; area) = domain.land.parameters
    river_indices = domain.river.network.land_indices

    # [m³ s⁻¹] =  [m s⁻¹] * [m²]
    @. model.boundary_conditions.runoff = (net_runoff + net_runoff_river) * area
    # [m³ s⁻¹] += [m³ s⁻¹]
    model.boundary_conditions.runoff[river_indices] .+=
        get_flux_to_river(subsurface_flow, river_indices)
    return nothing
end

"""
Update subsurface flow contribution to inflow of a reservoir model for a river flow model
`LocalInertialRiverFlow` for a single timestep.
"""
function update_inflow!(
    model::Union{Reservoir, Nothing},
    river_flow::LocalInertialRiverFlow,
    subsurface_flow::AbstractSubsurfaceFlowModel,
    network::NetworkReservoir,
)
    (; land_indices) = network
    if !isnothing(model)
        (; inflow_subsurface) = model.boundary_conditions
        inflow_subsurface .= get_inflow_reservoir(river_flow, subsurface_flow, land_indices)
    end
    return nothing
end

"""
Helper function to set flow variables of the `LocalInertialOverlandFlow` model to zero. This
is done at the start of each simulation timestep, during the timestep the total (weighted)
sum is computed from values at each sub timestep.
"""
function set_flow_vars!(model::LocalInertialOverlandFlow)
    (; qx_av, qy_av) = model.variables
    zero!(qx_av)
    zero!(qy_av)
    return nothing
end

"""
Helper function to compute average flow variables of the `LocalInertialOverlandFlow` model.
This is done at the end of each simulation timestep.
"""
function average_flow_vars!(model::LocalInertialOverlandFlow, dt::Float64)
    (; qx_av, qy_av) = model.variables
    average!(qx_av, dt)
    average!(qy_av, dt)
    return nothing
end

"""
Update combined river `LocalInertialRiverFlow` and overland flow `LocalInertialOverlandFlow`
models for a single timestep `dt`. An adaptive timestepping method is used (computing a sub
timestep `dt_s`).
"""
function update!(
    land::LocalInertialOverlandFlow,
    river::LocalInertialRiverFlow,
    domain::Domain,
    clock::Clock,
    dt::Number;
    update_h = false,
)
    (; reservoir) = river.boundary_conditions
    (; flow_length) = domain.river.parameters
    (; parameters) = domain.land

    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)
    set_flow_vars!(river)
    set_flow_vars!(land)

    t = 0.0
    while t < dt
        dt_river = stable_timestep(river, flow_length)
        dt_land = stable_timestep(land, parameters)
        dt_s = min(dt_river, dt_land)
        dt_s = check_timestepsize(dt_s, t, dt)

        local_inertial_update_fluxes!(land, domain, dt_s)
        update_inflow_reservoir!(land, reservoir, domain)
        local_inertial_river_update!(river, domain, dt_s, dt, update_h)
        local_inertial_update_water_depth!(land, river, domain, dt_s)

        t += dt_s
    end

    average_flow_vars!(river)
    average_flow_vars!(land)
    average_reservoir_vars!(reservoir)

    return nothing
end

"""
Update flow for a single direction in the local inertial overland flow model.
`is_x_direction`: true for x-direction (xu/xd), false for y-direction (yu/yd)
"""
@inline function update_directional_flow!(
    land::LocalInertialOverlandFlow,
    domain::Domain,
    i::Int,
    dt::Float64,
    is_x_direction::Bool,
)
    indices = domain.land.network.edge_indices
    (; x_length, y_length) = domain.land.parameters
    land_v = land.variables
    land_p = land.parameters

    # Select direction-specific parameters based on the boolean flag
    if is_x_direction
        upstream_idx = indices.xu[i]
        downstream_idx = indices.xd[i]
        width = land_p.ywidth[i]
        z_max = land_p.zx_max[i]
        length_vec = x_length
        q_current = land_v.qx
        q_prev = land_v.qx0
        q_av = land_v.qx_av
    else
        upstream_idx = indices.yu[i]
        downstream_idx = indices.yd[i]
        width = land_p.xwidth[i]
        z_max = land_p.zy_max[i]
        length_vec = y_length
        q_current = land_v.qy
        q_prev = land_v.qy0
        q_av = land_v.qy_av
    end

    # the effective flow width is zero when the river width exceeds the cell width and
    # floodplain flow is not calculated.
    if upstream_idx <= land_p.n && width != 0.0
        # [m] = [m] + [m]
        zs_current = land_p.z[i] + land_v.h[i]
        zs_upstream = land_p.z[upstream_idx] + land_v.h[upstream_idx]
        # [m] = max([m], [m])
        zs_max = max(zs_current, zs_upstream)
        # [m] = [m] - [m]
        hf = zs_max - z_max

        if hf > land_p.h_thresh'
            # [m] = [-] * ([m] + [m])
            length = 0.5 * (length_vec[i] + length_vec[upstream_idx]) # can be precalculated
            q_current[i] = local_inertial_flow(
                land_p.theta,
                q_prev[i],
                q_prev[downstream_idx],
                q_prev[upstream_idx],
                zs_current,
                zs_upstream,
                hf,
                width,
                length,
                land_p.mannings_n_sq[i],
                land_p.froude_limit,
                dt,
            )
            # limit q in case water is not available
            if land_v.h[i] <= 0.0
                q_current[i] = min(q_current[i], 0.0)
            end
            if land_v.h[upstream_idx] <= 0.0
                q_current[i] = max(q_current[i], 0.0)
            end
        else
            q_current[i] = 0.0
        end
        add_to_cumulative!(q_av, i, q_current[i], dt)
    end
    return nothing
end

"""
Update fluxes for overland flow `LocalInertialOverlandFlow` model for a single timestep
`dt`.
"""
function local_inertial_update_fluxes!(
    land::LocalInertialOverlandFlow,
    domain::Domain,
    dt::Float64,
)
    land_v = land.variables
    land_p = land.parameters

    # [m³ s⁻¹] = [m³ s⁻¹]
    land_v.qx0 .= land_v.qx
    land_v.qy0 .= land_v.qy

    @batch per = thread minbatch = 6000 for i in 1:(land_p.n)
        # update qx (x-direction)
        update_directional_flow!(land, domain, i, dt, true)

        # update qy (y-direction)
        update_directional_flow!(land, domain, i, dt, false)
    end
    return nothing
end

"""
Update boundary condition inflow to a reservoir from land `inflow_reservoir` of combined
river `LocalInertialRiverFlow`and overland flow `LocalInertialOverlandFlow` models for a
single timestep.
"""
function update_inflow_reservoir!(
    land::LocalInertialOverlandFlow,
    reservoir::Union{Reservoir, Nothing},
    domain::Domain,
)
    indices = domain.land.network.edge_indices
    reservoir_indices = domain.reservoir.network.land_indices
    land_bc = land.boundary_conditions
    land_v = land.variables

    for (i, j) in enumerate(reservoir_indices)
        yd = indices.yd[j]
        xd = indices.xd[j]
        reservoir.boundary_conditions.inflow_overland[i] =
            land_bc.runoff[j] +
            (land_v.qx[xd] - land_v.qx[j] + land_v.qy[yd] - land_v.qy[j])
    end
    return nothing
end

"""
Compute storage change for a river cell from fluxes.
"""
@inline function compute_river_storage_change(
    land::LocalInertialOverlandFlow,
    river::LocalInertialRiverFlow,
    domain::Domain,
    i::Int,
    dt::Float64,
)
    indices = domain.land.network.edge_indices
    inds_river = domain.land.network.river_indices
    edges_at_node = domain.river.network.edges_at_node
    river_idx = inds_river[i]

    yd = indices.yd[i]
    xd = indices.xd[i]

    return (
        sum_at(river.variables.q, edges_at_node.src[river_idx]) -
        sum_at(river.variables.q, edges_at_node.dst[river_idx]) + land.variables.qx[xd] -
        land.variables.qx[i] + land.variables.qy[yd] - land.variables.qy[i] +
        land.boundary_conditions.runoff[i] -
        river.boundary_conditions.abstraction[river_idx]
    ) * dt
end

"""
Compute external inflow for river cells, including negative inflow (abstraction).
Returns tuple: (inflow, abstraction_to_add)
"""
@inline function compute_external_inflow(
    river::LocalInertialRiverFlow,
    land::LocalInertialOverlandFlow,
    i::Int,
    river_idx::Int,
    dt::Float64,
)
    if river.boundary_conditions.external_inflow[river_idx] < 0.0
        # [m]
        available_volume =
            if land.variables.storage[i] >= river.parameters.bankfull_storage[river_idx]
                # [m]
                river.parameters.bankfull_depth[river_idx]
            else
                # [m]
                river.variables.storage[river_idx]
            end
        # [m³ s⁻¹] = min([m³ s⁻¹], [m³] / [s] * [-])
        abstraction = min(
            -river.boundary_conditions.external_inflow[river_idx],
            available_volume / dt * 0.80,
        )
        return (-abstraction, abstraction)
    else
        return (river.boundary_conditions.external_inflow[river_idx], 0.0)
    end
end

"""
Compute river and land water depths based on total storage and bankfull capacity.
Returns tuple: (river_h, land_h, river_storage)
"""
@inline function compute_water_depths(
    total_storage::Float64,
    river_idx::Int,
    i::Int,
    river::LocalInertialRiverFlow,
    domain::Domain,
)
    if total_storage >= river.parameters.bankfull_storage[river_idx]
        # Storage exceeds bankfull capacity - water spills onto floodplain
        # [m] = [m] + ([m³] - [m³]) / ([m] * [m])
        river_h =
            river.parameters.bankfull_depth[river_idx] +
            (total_storage - river.parameters.bankfull_storage[river_idx]) /
            (domain.land.parameters.x_length[i] * domain.land.parameters.y_length[i])
        land_h = river_h - river.parameters.bankfull_depth[river_idx]
        # [m³] = [m] * [m] * [m]
        river_storage =
            river_h *
            domain.river.parameters.flow_length[river_idx] *
            domain.river.parameters.flow_width[river_idx]
        return (river_h, land_h, river_storage)
    else
        # Storage is within channel capacity
        # [m] = [m³] / ([m] * [m])
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
    land::LocalInertialOverlandFlow,
    network::NetworkLand,
    i::Int,
    dt::Float64,
)
    indices = network.edge_indices
    yd = indices.yd[i]
    xd = indices.xd[i]
    # [m³] = (∑ [m³ s⁻¹]) * [s]
    return (
        land.variables.qx[xd] - land.variables.qx[i] + land.variables.qy[yd] -
        land.variables.qy[i] + land.boundary_conditions.runoff[i]
    ) * dt
end

"""
Update storage and water depth for a single river cell.
"""
@inline function update_river_cell_storage_and_depth!(
    land::LocalInertialOverlandFlow,
    river::LocalInertialRiverFlow,
    domain::Domain,
    i::Int,
    dt::Float64,
)
    inds_river = domain.land.network.river_indices
    river_idx = inds_river[i]

    # Compute and apply storage change from fluxes
    storage_change = compute_river_storage_change(land, river, domain, i, dt)
    land.variables.storage[i] += storage_change

    # Handle negative storage
    if land.variables.storage[i] < 0.0
        land.variables.error[i] += abs(land.variables.storage[i])
        land.variables.storage[i] = 0.0 # set storage to zero
    end

    # Compute and apply external inflow
    # [m³ s⁻¹], [m³ s⁻¹]
    inflow, abstraction_to_add = compute_external_inflow(river, land, i, river_idx, dt)
    # [m³] += [m³ s⁻¹] * [s]
    land.variables.storage[i] += inflow * dt
    add_to_cumulative!(
        river.boundary_conditions.actual_external_abstraction_av,
        river_idx,
        abstraction_to_add,
        dt,
    )

    # Compute and apply water depths
    river_h, land_h, river_storage =
        compute_water_depths(land.variables.storage[i], river_idx, i, river, domain)
    river.variables.h[river_idx] = river_h
    land.variables.h[i] = land_h
    river.variables.storage[river_idx] = river_storage

    return nothing
end

"""
Update storage and water depth for a single land cell (non-river).
"""
@inline function update_land_cell_storage_and_depth!(
    land::LocalInertialOverlandFlow,
    domain::DomainLand,
    i::Int,
    dt::Float64,
)
    # Compute and apply storage change
    storage_change = compute_land_storage_change(land, domain.network, i, dt)
    land.variables.storage[i] += storage_change

    # Handle negative storage
    if land.variables.storage[i] < 0.0
        land.variables.error[i] += abs(land.variables.storage[i])
        land.variables.storage[i] = 0.0 # set storage to zero
    end

    # Update water depth
    # [m] = [m³] / ([m] * [m])
    land.variables.h[i] =
        land.variables.storage[i] /
        (domain.parameters.x_length[i] * domain.parameters.y_length[i])

    return nothing
end

"""
Update storage and water depth for combined river `LocalInertialRiverFlow` and overland flow
`LocalInertialOverlandFlow` models for a single timestep `dt`.
"""
function local_inertial_update_water_depth!(
    land::LocalInertialOverlandFlow,
    river::LocalInertialRiverFlow,
    domain::Domain,
    dt::Float64,
)
    (; river_location, reservoir_outlet) = domain.land.parameters

    @batch per = thread minbatch = 6000 for i in 1:(land.parameters.n)
        if river_location[i]
            # Process river cells (excluding reservoir outlets)
            if !reservoir_outlet[i]
                update_river_cell_storage_and_depth!(land, river, domain, i, dt)
            end
        else
            # Process land cells (non-river)
            update_land_cell_storage_and_depth!(land, domain.land, i, dt)
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
    (; indices) = domain.network
    (; flow_width, flow_length) = domain.parameters
    storage = ncread(
        dataset,
        config,
        "floodplain_water__sum_of_volume_per_depth",
        Routing;
        sel = indices,
    )
    n = length(indices)

    # for convenience (interpolation) flood depth 0.0 m is added, with associated area (a),
    # storage, width (river width) and wetted perimeter (p).
    storage = vcat(fill(Float64(0), n)', storage)
    start_storage = storage
    flood_depths = Float64.(dataset["flood_depth"][:])
    pushfirst!(flood_depths, 0.0)
    n_depths = length(flood_depths)

    p = zeros(n_depths, n)
    a = zeros(n_depths, n)
    segment_storage = zeros(n_depths, n)
    width = zeros(n_depths, n)
    width[1, :] = flow_width[1:n]

    # determine flow area (a), width and wetted perimeter (p) FloodPlain
    h = diff(flood_depths)
    incorrect_vol = 0
    riv_cells = 0
    error_vol = 0
    for i in 1:n
        riv_cell = 0
        diff_storage = diff(storage[:, i])

        for j in 1:(n_depths - 1)
            # assume rectangular shape of flood depth segment
            width[j + 1, i] = diff_storage[j] / (h[j] * flow_length[i])
            # check provided flood storage (floodplain width should be constant or increasing
            # as a function of flood depth)
            if width[j + 1, i] < width[j, i]
                # raise warning only if difference is larger than rounding error of 0.01 m³
                if ((width[j, i] - width[j + 1, i]) * h[j] * flow_length[i]) > 0.01
                    incorrect_vol += 1
                    riv_cell = 1
                    error_vol =
                        error_vol +
                        ((width[j, i] - width[j + 1, i]) * h[j] * flow_length[i])
                end
                width[j + 1, i] = width[j, i]
            end
            a[j + 1, i] = width[j + 1, i] * h[j]
            p[j + 1, i] = (width[j + 1, i] - width[j, i]) + 2.0 * h[j]
            segment_storage[j + 1, i] = a[j + 1, i] * flow_length[i]
            if j == 1
                # for interpolation wetted perimeter at flood depth 0.0 is required
                p[j, i] = p[j + 1, i] - 2.0 * h[j]
            end
        end

        p[2:end, i] = cumsum(p[2:end, i])
        a[:, i] = cumsum(a[:, i])
        storage[:, i] = cumsum(segment_storage[:, i])

        riv_cells += riv_cell
    end

    if incorrect_vol > 0
        perc_riv_cells = round(100.0 * (riv_cells / n); digits = 2)
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
    (; indices, nodes_at_edge, graph) = domain.network
    (; flow_length) = domain.parameters
    n_edges = ne(graph)
    profile = FloodPlainProfile(dataset, config, domain, index_pit)

    mannings_n = ncread(
        dataset,
        config,
        "floodplain_water_flow__manning_n_parameter",
        Routing;
        sel = indices,
    )
    # manning roughness at edges
    append!(mannings_n, mannings_n[index_pit]) # copy to ghost nodes
    mannings_n_sq = fill(Float64(0), n_edges)
    zb_max = fill(Float64(0), n_edges)
    for i in 1:n_edges
        src_node = nodes_at_edge.src[i]
        dst_node = nodes_at_edge.dst[i]
        mannings_n_i =
            (
                mannings_n[dst_node] * flow_length[dst_node] +
                mannings_n[src_node] * flow_length[src_node]
            ) / (flow_length[dst_node] + flow_length[src_node])
        mannings_n_sq[i] = mannings_n_i * mannings_n_i
        zb_max[i] = max(zb_floodplain[src_node], zb_floodplain[dst_node])
    end
    parameters = FloodPlainParameters(profile, mannings_n, mannings_n_sq, zb_max)
    return parameters
end

"Struct to store floodplain flow model variables"
@with_kw struct FloodPlainVariables
    n::Int
    n_edges::Int
    # storage [m³]
    storage::Vector{Float64} = zeros(n)
    # water depth [m]
    h::Vector{Float64}
    # error storage [m³]
    error::Vector{Float64} = zeros(n)
    # flow area at egde [m²]
    a::Vector{Float64} = zeros(n_edges)
    # hydraulic radius at edge [m]
    r::Vector{Float64} = zeros(n_edges)
    # water depth at edge [m]
    hf::Vector{Float64} = zeros(n_edges)
    # discharge at edge at previous time step
    q0::Vector{Float64} = zeros(n_edges)
    # discharge at edge  [m³ s⁻¹]
    q::Vector{Float64} = zeros(n_edges)
    # average river discharge at edge  [m³ s⁻¹] for model timestep dt
    q_av::AverageVector = AverageVector(; n = n_edges)
    # edge index with `hf` [-] above depth threshold
    hf_index::Vector{Int} = zeros(Int, n_edges)
end

"Floodplain flow model"
@with_kw struct FloodPlain <: AbstractFloodPlain
    parameters::FloodPlainParameters
    variables::FloodPlainVariables
end

"Determine the initial floodplain storage"
function initialize_storage!(river, domain::Domain, nriv::Int)
    (; flow_width, flow_length) = domain.river.parameters
    (; floodplain) = river
    (; profile) = floodplain.parameters
    for i in 1:nriv
        i1, i2 = interpolation_indices(floodplain.variables.h[i], profile.depth)
        a = flow_area(
            profile.width[i2, i],
            profile.a[i1, i],
            profile.depth[i1],
            floodplain.variables.h[i],
        )
        a = max(a - (flow_width[i] * floodplain.variables.h[i]), 0.0)
        floodplain.variables.storage[i] = flow_length[i] * a
    end
    return nothing
end

"helper function to get interpolation indices"
function interpolation_indices(x, v::AbstractVector)
    i1 = 1
    for i in eachindex(v)
        if v[i] <= x
            i1 = i
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
    i::Int,
)
    i1, i2 = interpolation_indices(flood_storage, @view profile.storage[:, i])
    # [m²] = ([m³] - [m³]) / [m]
    ΔA = (flood_storage - profile.storage[i1, i]) / flow_length
    # [m] = [m²] / [m]
    dh = ΔA / profile.width[i2, i]
    # [m] = [m] + [m]
    flood_depth = profile.depth[i1] + dh
    return flood_depth
end

"Initialize floodplain geometry and `FloodPlain` variables and parameters"
function FloodPlain(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    zb_floodplain::Vector{Float64},
)
    (; indices, local_drain_direction, graph) = domain.network
    n = length(indices)
    index_pit = findall(x -> x == 5, local_drain_direction)
    parameters = FloodPlainParameters(dataset, config, domain, zb_floodplain, index_pit)
    h = zeros(n + length(index_pit))
    n_edges = ne(graph)
    variables = FloodPlainVariables(; n, n_edges, h)

    floodplain = FloodPlain(; parameters, variables)
    return floodplain
end
