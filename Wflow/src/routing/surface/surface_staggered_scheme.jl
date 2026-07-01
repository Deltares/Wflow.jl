"Struct for storing river flow model parameters on a staggered grid"
@with_kw struct RiverFlowStaggeredParameters <: AbstractRiverFlowParameters
    # number of cells [-]
    n::Int
    # number of edges [-]
    n_edges::Int
    # active nodes [-]
    active_n::Vector{Int}
    # active edges [-]
    active_e::Vector{Int}
    # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
    froude_limit::Bool = false
    # depth threshold for calculating flow [m]
    h_thresh::Float64 = 0.0
    # river bed elevation [m]
    zb::Vector{Float64} = Float64[]
    # bankfull storage [m³]
    bankfull_storage::Vector{Float64} = Float64[]
    # bankfull depth [m]
    bankfull_depth::Vector{Float64} = Float64[]
    # maximum channel bed elevation at edge [m]
    zb_max_at_edge::Vector{Float64} = Float64[]
    # Manning's roughness [s m-1/3]
    mannings_n::Vector{Float64} = Float64[]
    # Manning's roughness at edge [s m-1/3]
    mannings_n_at_edge::Vector{Float64} = Float64[]
    # Manning's roughness squared at edge [(s m-1/3)2]
    mannings_n_sq_at_edge::Vector{Float64} = Float64[]
    # slope at edge [-]
    slope_at_edge::Vector{Float64} = Float64[]
    # flow (river) length at edge [m]
    flow_length_at_edge::Vector{Float64} = Float64[]
    # flow (river) width at edge [m]
    flow_width_at_edge::Vector{Float64} = Float64[]
end

"Get river parameters for a river flow model on a staggered grid"
function get_river_parameters(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    index_pit::Vector{Int},
)
    (; pit_indices, indices) = domain.network
    (; flow_length, flow_width) = domain.parameters

    bankfull_elevation =
        ncread(dataset, config, "river_bank_water__elevation", Routing; sel = indices)
    bankfull_depth =
        ncread(dataset, config, "river_bank_water__depth", Routing; sel = indices)
    mannings_n = ncread(
        dataset,
        config,
        "river_water_flow__manning_n_parameter",
        Routing;
        sel = indices,
    )
    riverlength_bc = ncread(
        dataset,
        config,
        "model_boundary_condition_river__length",
        Routing;
        sel = pit_indices,
    )
    append!(flow_length, riverlength_bc)

    # set ghost points for boundary condition (downstream river outlet): flow width,
    # mannings n, bankfull depth and elevation is copied from the upstream cell.
    append!(flow_width, flow_width[index_pit])
    append!(mannings_n, mannings_n[index_pit])
    append!(bankfull_depth, bankfull_depth[index_pit])
    append!(bankfull_elevation, bankfull_elevation[index_pit])

    return mannings_n, bankfull_depth, bankfull_elevation, flow_length, flow_width
end

"""
Generate info log message for river flow routing on a staggered grid, solving the local
inertial or Manning's equation.
"""
function log_message_staggered_flow(config::Config)
    (; river_routing) = config.model
    waterdepth_threshold = config.model.river_water_flow_threshold__depth # depth threshold for flow at edge
    floodplain_1d = config.model.floodplain_1d__flag

    if river_routing == RoutingType.local_inertial
        alpha = config.model.river_local_inertial_flow__alpha_coefficient # stability coefficient for model time step (0.2-0.7)
        froude_limit = config.model.river_water_flow__froude_limit_flag # limit flow to subcritical according to Froude number
        @info "Local inertial approach is used for river flow." alpha waterdepth_threshold froude_limit floodplain_1d
    elseif river_routing == RoutingType.manning_staggered
        alpha = config.model.river_staggered_manning_flow__alpha_coefficient
        froude_limit = false
        @info "Manning's equation on a staggered grid is used for river flow." alpha floodplain_1d
    end
    return waterdepth_threshold, froude_limit
end

"Initialize river flow model parameters on a staggered grid."
function RiverFlowStaggeredParameters(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
)
    (; river_routing) = config.model
    (; graph, indices, local_drain_direction, nodes_at_edge) = domain.network
    (; reservoir_outlet) = domain.parameters

    waterdepth_threshold, froude_limit = log_message_staggered_flow(config)

    n = length(indices)
    n_edges = ne(graph)
    active_index = findall(x -> x == 0, reservoir_outlet)
    index_pit = findall(x -> x == LDD_PIT, local_drain_direction)

    mannings_n, bankfull_depth, bankfull_elevation, flow_length, flow_width =
        get_river_parameters(dataset, config, domain, index_pit)
    zb = bankfull_elevation - bankfull_depth # river bed elevation
    bankfull_storage = bankfull_depth .* flow_width .* flow_length

    # determine parameters at edges
    zb_max_at_edge = compute_value_at_edge(zb, nodes_at_edge, n_edges, maximum)
    flow_width_at_edge = compute_value_at_edge(flow_width, nodes_at_edge, n_edges, minimum)
    flow_length_at_edge = compute_value_at_edge(flow_length, nodes_at_edge, n_edges, mean)
    mannings_n_at_edge =
        compute_mannings_n_at_edge(mannings_n, flow_length, nodes_at_edge, n_edges)
    if river_routing == RoutingType.local_inertial
        mannings_n_sq_at_edge = mannings_n_at_edge .* mannings_n_at_edge
        slope_at_edge = []
    elseif river_routing == RoutingType.manning_staggered
        mannings_n_sq_at_edge = []
        slope_at_edge =
            compute_slope_at_edge(zb, flow_length_at_edge, nodes_at_edge, n_edges)
    end

    parameters = RiverFlowStaggeredParameters(;
        n,
        n_edges,
        active_n = active_index,
        active_e = active_index,
        froude_limit,
        h_thresh = waterdepth_threshold,
        zb,
        bankfull_storage,
        bankfull_depth,
        zb_max_at_edge,
        mannings_n,
        mannings_n_at_edge,
        mannings_n_sq_at_edge,
        slope_at_edge,
        flow_length_at_edge,
        flow_width_at_edge,
    )
    return parameters
end

"Struct for storing river flow model variables on a staggered grid."
@with_kw struct RiverFlowStaggeredVariables <: AbstractRiverFlowVariables
    n_cells::Int
    n_edges::Int
    # river discharge at edge (subgrid channel) [m³ s⁻¹]
    q::Vector{Float64} = zeros(n_edges)
    # river discharge at edge (subgrid channel) at previous time step [m³ s⁻¹]
    q_previous::Vector{Float64} = zeros(n_edges)
    # cumulative river channel (+ floodplain) discharge at edge [m³] (model time step dt)
    q_cumulative::Vector{Float64} = zeros(n_edges)
    # average river channel (+ floodplain) discharge at edge [m³ s⁻¹] (model time step dt)
    q_average::Vector{Float64} = zeros(n_edges)
    # cumulative river channel discharge at edge [m³] (for model time step dt)
    q_channel_cumulative::Vector{Float64} = q_cumulative
    # average river channel discharge at edge [m³ s⁻¹] (for model time step dt)
    q_channel_average::Vector{Float64} = q_average
    # water depth [m]
    h::Vector{Float64}
    # maximum water elevation at edge [m]
    zs_max_at_edge::Vector{Float64} = zeros(n_edges)
    # water elevation of source node of edge [m]
    zs_src::Vector{Float64} = zeros(n_edges)
    # water elevation of downstream node of edge [m]
    zs_dst::Vector{Float64} = zeros(n_edges)
    # water depth at edge [m]
    water_depth_at_edge::Vector{Float64} = zeros(n_edges)
    # river storage [m³]
    storage::Vector{Float64} = zeros(n_cells)
    # error storage [m³]
    error::Vector{Float64} = zeros(n_cells)
end

"Initialize river flow model variables on a staggered grid."
function RiverFlowStaggeredVariables(
    dataset::NCDataset,
    config::Config,
    network::NetworkRiver,
)
    (; pit_indices, indices, graph) = network
    (; river_routing) = config.model

    riverdepth_bc = ncread(
        dataset,
        config,
        "model_boundary_condition_river_bank_water__depth",
        Routing;
        sel = pit_indices,
    )

    n_cells = length(indices)
    n_edges = ne(graph)
    # set ghost points for boundary condition (downstream river outlet): river depth `h`
    h = zeros(n_cells)
    append!(h, riverdepth_bc)

    variables = if config.model.floodplain_1d__flag
        # When floodplain is enabled, q_channel fields must be separate from q fields
        # because q_average will hold channel + floodplain combined discharge.
        RiverFlowStaggeredVariables(;
            n_cells,
            n_edges,
            h,
            q_channel_cumulative = zeros(n_edges),
            q_channel_average = zeros(n_edges),
        )
    else
        # When floodplain is disabled, q_channel aliases q (they are identical)
        RiverFlowStaggeredVariables(; n_cells, n_edges, h)
    end
    return variables
end

"Initialize river flow model on a staggered grid."
function init_staggered_river_flow(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    reservoir_model::Union{ReservoirModel, Nothing},
)
    # This river flow model makes use of a staggered grid (Bates et al. (2010)), with nodes
    # and edges. This information is extracted from the directed graph of the river.
    # Discharge q is calculated at edges between nodes and mapped to the source nodes for
    # gridded output (index of edge is equal to source node index, e.g.:
    # Edge 1 => 5
    # Edge 2 => 1
    # Edge 3 => 2
    # Edge 4 => 9
    # ⋮ )

    # The following boundary conditions can be set at ghost nodes, downstream of river
    # outlets (pits): river length and river depth
    alpha_coefficient = config.model.river_local_inertial_flow__alpha_coefficient # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; alpha_coefficient)

    parameters = RiverFlowStaggeredParameters(dataset, config, domain)
    variables = RiverFlowStaggeredVariables(dataset, config, domain.network)
    boundary_conditions = RiverFlowBC(dataset, config, domain.network, reservoir_model)

    if config.model.floodplain_1d__flag
        zb_floodplain = parameters.zb .+ parameters.bankfull_depth
        floodplain = FloodPlainModel(dataset, config, domain, zb_floodplain)
    else
        floodplain = nothing
    end

    (; river_routing) = config.model
    routing_method = if river_routing == RoutingType.local_inertial
        LocalInertial()
    elseif river_routing == RoutingType.manning_staggered
        ManningStaggered()
    end

    n = length(domain.network.indices)
    river_flow = RiverFlowModel(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain,
        allocation = do_water_demand(config) ? AllocationRiverModel(n) :
                     NoAllocationRiverModel(n),
        routing_method,
    )
    return river_flow
end

"""
Return the upstream inflow for a reservoir in a river flow model on a staggered grid.
"""
function get_inflow_reservoir(
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
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
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    overland_flow_model::OverlandFlowModel{<:KinematicWave},
    inds::Vector{Int},
) =
    overland_flow_model.variables.q_average[inds] .+
    overland_flow_model.variables.to_river_average[inds]

get_inflow_reservoir(
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    subsurface_flow_model::LateralSSFModel,
    inds::Vector{Int},
) = (
    subsurface_flow_model.variables.q_average[inds] .+
    subsurface_flow_model.variables.to_river_average[inds]
)

"""
Update river channel flow for the local inertial river flow model.
"""
function update_river_channel_flow!(
    river_flow_model::RiverFlowModel{<:LocalInertial},
    domain::DomainRiver,
    dt::Float64,
)
    (; nodes_at_edge) = domain.network
    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters

    river_v.q_previous .= river_v.q
    if !isnothing(river_flow_model.floodplain)
        river_flow_model.floodplain.variables.q_previous .=
            river_flow_model.floodplain.variables.q
    end

    @batch per = thread minbatch = 1000 for j in eachindex(river_p.active_e)
        i = river_p.active_e[j]
        i_src = nodes_at_edge.src[i]
        i_dst = nodes_at_edge.dst[i]
        river_v.zs_src[i] = river_p.zb[i_src] + river_v.h[i_src]
        river_v.zs_dst[i] = river_p.zb[i_dst] + river_v.h[i_dst]

        river_v.zs_max_at_edge[i] = max(river_v.zs_src[i], river_v.zs_dst[i])
        river_v.water_depth_at_edge[i] =
            (river_v.zs_max_at_edge[i] - river_p.zb_max_at_edge[i])

        # rectangular channel
        flow_area_at_edge = river_p.flow_width_at_edge[i] * river_v.water_depth_at_edge[i]
        hydraulic_radius_at_edge =
            flow_area_at_edge / wetted_perimeter_channel(
                river_v.water_depth_at_edge[i],
                river_p.flow_width_at_edge[i],
            )

        river_v.q[i] = ifelse(
            river_v.water_depth_at_edge[i] > river_p.h_thresh,
            local_inertial_flow(
                river_v.q_previous[i],
                river_v.zs_src[i],
                river_v.zs_dst[i],
                river_v.water_depth_at_edge[i],
                flow_area_at_edge,
                hydraulic_radius_at_edge,
                river_p.flow_length_at_edge[i],
                river_p.mannings_n_sq_at_edge[i],
                river_p.froude_limit,
                dt,
            ),
            0.0,
        )

        # limit q in case water is not available
        river_v.q[i] = ifelse(river_v.h[i_src] <= 0.0, min(river_v.q[i], 0.0), river_v.q[i])
        river_v.q[i] = ifelse(river_v.h[i_dst] <= 0.0, max(river_v.q[i], 0.0), river_v.q[i])
        # average river discharge (here accumulated for model timestep dt)
        river_v.q_cumulative[i] += river_v.q[i] * dt
    end
    return nothing
end

"""
Update river channel flow using Manning's equation on a staggered grid.
"""
function update_river_channel_flow!(
    river_flow_model::RiverFlowModel{<:ManningStaggered},
    domain::DomainRiver,
    dt::Float64,
)
    (; nodes_at_edge) = domain.network
    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters
    (; inwater, abstraction) = river_flow_model.boundary_conditions

    @batch per = thread minbatch = 1000 for j in eachindex(river_p.active_e)
        i = river_p.active_e[j]
        i_src = nodes_at_edge.src[i]
        i_dst = nodes_at_edge.dst[i]
        river_v.zs_src[i] = river_p.zb[i_src] + river_v.h[i_src]
        river_v.zs_dst[i] = river_p.zb[i_dst] + river_v.h[i_dst]

        river_v.zs_max_at_edge[i] = max(river_v.zs_src[i], river_v.zs_dst[i])
        river_v.water_depth_at_edge[i] =
            (river_v.zs_max_at_edge[i] - river_p.zb_max_at_edge[i])

        # rectangular channel
        flow_area_at_edge = river_p.flow_width_at_edge[i] * river_v.water_depth_at_edge[i] # flow area (rectangular channel)
        hydraulic_radius_at_edge =
            flow_area_at_edge / wetted_perimeter_channel(
                river_v.water_depth_at_edge[i],
                river_p.flow_width_at_edge[i],
            )

        river_v.q[i] = ifelse(
            river_v.water_depth_at_edge[i] > river_p.h_thresh,
            manning_flow(
                river_p.mannings_n_at_edge[i],
                hydraulic_radius_at_edge,
                river_p.slope_at_edge[i],
                flow_area_at_edge,
            ),
            0.0,
        )

        # limit q in case water is not available
        river_v.q[i] = min(river_v.q[i], river_v.storage[i_src]/dt)
        # average river discharge (here accumulated for model timestep Δt)
        river_v.q_cumulative[i] += river_v.q[i] * dt
    end
    return nothing
end

"""
Update floodplain flow for the local inertial river flow model.
"""
function update_floodplain_flow!(
    river_flow_model::RiverFlowModel{T, F},
    domain::DomainRiver,
    dt::Float64,
) where {T <: LocalInertial, F <: FloodPlainModel{<:LocalInertial}}
    (; nodes_at_edge) = domain.network
    (; flow_width) = domain.parameters

    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters
    floodplain_p = river_flow_model.floodplain.parameters
    (; profile) = floodplain_p
    floodplain_v = river_flow_model.floodplain.variables

    @batch per = thread minbatch = 1000 for j in eachindex(river_p.active_e)
        i = river_p.active_e[j]
        i_src = nodes_at_edge.src[i]
        i_dst = nodes_at_edge.dst[i]

        floodplain_v.water_depth_at_edge[i] =
            max(river_v.zs_max_at_edge[i] - floodplain_p.zb_max_at_edge[i], 0.0)

        i1, i2 = interpolation_indices(
            floodplain_v.water_depth_at_edge[i],
            @view profile.depth[:]
        )
        a_src = compute_floodplain_flow_area(
            profile,
            floodplain_v.water_depth_at_edge[i],
            i_src,
            i1,
            i2,
        )
        a_dst = compute_floodplain_flow_area(
            profile,
            floodplain_v.water_depth_at_edge[i],
            i_dst,
            i1,
            i2,
        )
        flow_area_at_edge = min(a_src, a_dst)

        hydraulic_radius_at_edge = ifelse(
            a_src < a_dst,
            a_src / compute_wetted_perimeter(
                profile,
                floodplain_v.water_depth_at_edge[i],
                i_src,
                i1,
            ),
            a_dst / compute_wetted_perimeter(
                profile,
                floodplain_v.water_depth_at_edge[i],
                i_dst,
                i1,
            ),
        )
        hydraulic_radius_at_edge = ifelse(
            floodplain_v.water_depth_at_edge[i] > profile.depth[i1],
            hydraulic_radius_at_edge,
            0.0,
        )

        floodplain_v.q[i] = ifelse(
            flow_area_at_edge > 1.0e-05,
            local_inertial_flow(
                floodplain_v.q_previous[i],
                river_v.zs_src[i],
                river_v.zs_dst[i],
                floodplain_v.water_depth_at_edge[i],
                flow_area_at_edge,
                hydraulic_radius_at_edge,
                river_p.flow_length_at_edge[i],
                floodplain_p.mannings_n_sq_at_edge[i],
                river_p.froude_limit,
                dt,
            ),
            0.0,
        )

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
        floodplain_v.q_cumulative[i] += floodplain_v.q[i] * dt
    end
    return nothing
end

"""
Update floodplain flow for a river flow model using Manning's equation on a staggered grid.
"""
function update_floodplain_flow!(
    river_flow_model::RiverFlowModel{T, F},
    domain::DomainRiver,
    dt::Float64,
) where {T <: ManningStaggered, F <: FloodPlainModel{<:ManningStaggered}}
    (; nodes_at_edge) = domain.network
    (; flow_width) = domain.parameters

    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters
    floodplain_p = river_flow_model.floodplain.parameters
    (; profile) = floodplain_p
    floodplain_v = river_flow_model.floodplain.variables

    @batch per = thread minbatch = 1000 for j in eachindex(river_p.active_e)
        i = river_p.active_e[j]
        i_src = nodes_at_edge.src[i]
        i_dst = nodes_at_edge.dst[i]

        floodplain_v.water_depth_at_edge[i] =
            max(river_v.zs_max_at_edge[i] - floodplain_p.zb_max_at_edge[i], 0.0)

        i1, i2 = interpolation_indices(
            floodplain_v.water_depth_at_edge[i],
            @view profile.depth[:]
        )
        a_src = compute_floodplain_flow_area(
            profile,
            floodplain_v.water_depth_at_edge[i],
            i_src,
            i1,
            i2,
        )
        a_dst = compute_floodplain_flow_area(
            profile,
            floodplain_v.water_depth_at_edge[i],
            i_dst,
            i1,
            i2,
        )
        flow_area_at_edge = min(a_src, a_dst)

        hydraulic_radius_at_edge = ifelse(
            a_src < a_dst,
            a_src / compute_wetted_perimeter(
                profile,
                floodplain_v.water_depth_at_edge[i],
                i_src,
                i1,
            ),
            a_dst / compute_wetted_perimeter(
                profile,
                floodplain_v.water_depth_at_edge[i],
                i_dst,
                i1,
            ),
        )
        hydraulic_radius_at_edge = ifelse(
            floodplain_v.water_depth_at_edge[i] > profile.depth[i1],
            hydraulic_radius_at_edge,
            0.0,
        )

        floodplain_v.q[i] = ifelse(
            flow_area_at_edge > 1.0e-05,
            manning_flow(
                floodplain_p.mannings_n_at_edge[i],
                hydraulic_radius_at_edge,
                floodplain_p.slope_at_edge[i],
                flow_area_at_edge,
            ),
            0.0,
        )

        # limit floodplain q in case water is not available
        floodplain_v.q[i] = min(floodplain_v.q[i], floodplain_v.storage[i_src]/dt)

        # average floodplain discharge (here accumulated for model timestep Δt)
        floodplain_v.q_cumulative[i] += floodplain_v.q[i] * dt
    end
    return nothing
end

update_floodplain_flow!(
    model::RiverFlowModel{T, F},
    domain::DomainRiver,
    dt::Float64,
) where {T <: AbstractStaggeredRoutingMethod, F <: Nothing} = nothing

"""
Update reservoir boundary conditions for a river flow model on a staggered grid.
"""
function update_bc_reservoir_model!(
    reservoir_model::ReservoirModel,
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    domain::Domain,
    dt::Float64,
)
    (; edges_at_node) = domain.river.network
    inds_reservoir = domain.reservoir.network.river_indices

    river_v = river_flow_model.variables
    res_bc = reservoir_model.boundary_conditions

    for v in eachindex(inds_reservoir)
        i = inds_reservoir[v]

        q_in = get_inflow_reservoir(river_flow_model, edges_at_node.src[i])
        # If external_inflow < 0, abstraction is limited
        if res_bc.external_inflow[v] < 0.0
            abstraction = min(
                -res_bc.external_inflow[v],
                (reservoir_model.variables.storage[v] / dt) * 0.98,
            )
            res_bc.actual_external_abstraction_cumulative[v] += abstraction * dt
            inflow = -abstraction
        else
            inflow = res_bc.external_inflow[v]
        end
        net_inflow = q_in + res_bc.inflow_overland[v] + res_bc.inflow_subsurface[v] + inflow
        update_reservoir_model!(reservoir_model, v, net_inflow, dt)
        river_v.q[i] = reservoir_model.variables.outflow[v]
        # average river discharge (here accumulated for model timestep Δt)
        river_v.q_cumulative[i] += river_v.q[i] * dt
    end
    return nothing
end

update_bc_reservoir_model!(
    reservoir_model::Nothing,
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    domain::Domain,
    dt::Float64,
) = nothing

"""
Update river and floodplain water depth and storage for river and floodplain models on a
staggered grid.
"""
function update_water_depth_and_storage!(
    floodplain_model::AbstractFloodPlainModel,
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    domain::DomainRiver,
    dt::Float64,
)
    (; edges_at_node) = domain.network
    (; flow_length, flow_width) = domain.parameters

    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters
    floodplain_v = floodplain_model.variables
    floodplain_p = floodplain_model.parameters

    @batch per = thread minbatch = 1000 for i in river_p.active_n
        q_src = sum_at(floodplain_v.q, edges_at_node.src[i])
        q_dst = sum_at(floodplain_v.q, edges_at_node.dst[i])
        floodplain_v.storage[i] += (q_src - q_dst) * dt
        if floodplain_v.storage[i] < 0.0
            floodplain_v.error[i] += abs(floodplain_v.storage[i])
            floodplain_v.storage[i] = 0.0
        end
        storage_total = river_v.storage[i] + floodplain_v.storage[i]
        if storage_total > river_p.bankfull_storage[i]
            flood_storage = storage_total - river_p.bankfull_storage[i]
            h = compute_flood_depth(floodplain_p.profile, flood_storage, flow_length[i], i)
            river_v.h[i] = river_p.bankfull_depth[i] + h
            river_v.storage[i] = river_v.h[i] * flow_width[i] * flow_length[i]
            floodplain_v.storage[i] = max(storage_total - river_v.storage[i], 0.0)
            floodplain_v.h[i] = floodplain_v.storage[i] > 0.0 ? h : 0.0
        else
            river_v.h[i] = storage_total / (flow_length[i] * flow_width[i])
            river_v.storage[i] = storage_total
            floodplain_v.h[i] = 0.0
            floodplain_v.storage[i] = 0.0
        end
    end
    return nothing
end
update_water_depth_and_storage!(
    floodplain_model::Nothing,
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    domain::DomainRiver,
    dt::Float64,
) = nothing

"""
Update water depth and storage for river flow model on a staggered grid.
"""
function update_water_depth_and_storage!(
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    domain::DomainRiver,
    dt::Float64,
)
    (; edges_at_node) = domain.network
    (; flow_length, flow_width) = domain.parameters
    (; inwater, abstraction, external_inflow, actual_external_abstraction_cumulative) =
        river_flow_model.boundary_conditions

    river_v = river_flow_model.variables
    river_p = river_flow_model.parameters

    @batch per = thread minbatch = 1000 for i in river_p.active_n
        q_src = sum_at(river_v.q, edges_at_node.src[i])
        q_dst = sum_at(river_v.q, edges_at_node.dst[i])
        # internal abstraction (water demand) is limited by river storage and negative
        # external inflow as part of water allocation computations.
        river_v.storage[i] += (q_src - q_dst + inwater[i] - abstraction[i]) * dt

        if river_v.storage[i] < 0.0
            river_v.error[i] = river_v.error[i] + abs(river_v.storage[i])
            river_v.storage[i] = 0.0 # set storage to zero
        end
        # limit negative external inflow
        if external_inflow[i] < 0.0
            _abstraction = min(-external_inflow[i], river_v.storage[i] / dt * 0.80)
            actual_external_abstraction_cumulative[i] += _abstraction * dt
            inflow = -_abstraction
        else
            inflow = external_inflow[i]
        end
        river_v.storage[i] += inflow * dt # add external inflow
        river_v.h[i] = river_v.storage[i] / (flow_length[i] * flow_width[i])
    end
    return nothing
end

"Update river flow model on a staggered grid for a single timestep"
function staggered_scheme_river_update!(
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    domain::Domain,
    dt::Float64,
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
Update river flow model on a staggered grid for a single timestep `dt`. An adaptive
timestepping method is used (computing a sub timestep `dt_s`).
"""
function update_river_flow_model!(
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    domain::Domain,
    clock::Clock,
    dt::Float64;
    update_h = true,
)
    (; reservoir) = river_flow_model.boundary_conditions
    (; flow_length) = domain.river.parameters

    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)

    if !isnothing(river_flow_model.floodplain)
        river_flow_model.floodplain.variables.q_cumulative .= 0.0
    end
    set_flow_vars!(river_flow_model)

    t = 0.0
    while t < dt
        dt_s = stable_timestep(river_flow_model, flow_length)
        dt_s = check_timestepsize(dt_s, t, dt)
        staggered_scheme_river_update!(river_flow_model, domain, dt_s, update_h)
        t += dt_s
    end
    average_flow_vars!(river_flow_model, dt)
    average_reservoir_vars!(reservoir, dt)

    if !isnothing(river_flow_model.floodplain)
        v_floodplain = river_flow_model.floodplain.variables
        v_river = river_flow_model.variables

        @. v_floodplain.q_average = v_floodplain.q_cumulative / dt
        @. v_river.q_channel_average = v_river.q_average
        @. v_river.q_average = v_river.q_channel_average + v_floodplain.q_average
    end

    return nothing
end

"Struct to store local inertial overland flow model variables"
@with_kw struct LocalInertialOverlandFlowVariables <: AbstractOverlandFlowVariables
    n::Int
    # flow in y direction at edge at previous time step [m³ s⁻¹]
    qy0::Vector{Float64} = zeros(n + 1)
    # flow in x direction at edge at previous time step [m³ s⁻¹]
    qx0::Vector{Float64} = zeros(n + 1)
    # flow in x direction at edge [m³ s⁻¹]
    qx::Vector{Float64} = zeros(n + 1)
    # cumulative flow in x direction at edge [m³] for model timestep dt
    qx_cumulative::Vector{Float64} = zeros(n + 1)
    # average flow in x direction at edge [m³ s⁻¹] for model timestep dt
    qx_average::Vector{Float64} = zeros(n + 1)
    # flow in y direction at edge [m³ s⁻¹]
    qy::Vector{Float64} = zeros(n + 1)
    # cumulative flow in y direction at edge [m³] for model timestep dt
    qy_cumulative::Vector{Float64} = zeros(n + 1)
    # average flow in y direction at edge [m³ s⁻¹] for model timestep dt
    qy_average::Vector{Float64} = zeros(n + 1)
    # total storage of cell [m³] (including river storage for river cells)
    storage::Vector{Float64} = zeros(n)
    # error storage [m³]
    error::Vector{Float64} = zeros(n)
    # water depth of cell [m] (for river cells the reference is the river bed elevation `zb`)
    h::Vector{Float64} = zeros(n)
end

"Struct to store local inertial overland flow model parameters"
@with_kw struct LocalInertialOverlandFlowParameters <: AbstractOverlandFlowParameters
    # number of cells [-]
    n::Int
    # effective flow width x direction at edge (floodplain) [m]
    xwidth_at_edge::Vector{Float64}
    # effective flow width y direction at edge (floodplain) [m]
    ywidth_at_edge::Vector{Float64}
    # acceleration due to gravity [m s⁻²]
    g::Float64 = 9.80665
    # weighting factor (de Almeida et al., 2012) [-]
    theta::Float64
    # depth threshold for calculating flow [m]
    h_thresh::Float64
    # maximum cell elevation at edge [m] (x direction)
    zx_max_at_edge::Vector{Float64}
    # maximum cell elevation at edge [m] (y direction)
    zy_max_at_edge::Vector{Float64}
    # Manning's roughness squared at edge [(s m-1/3)2]
    mannings_n_sq_at_edge::Vector{Float64}
    # elevation [m] of cell
    z::Vector{Float64}
    # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
    froude_limit::Bool
end

"Initialize local inertial overland flow model parameters"
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
    elevation = ncread(
        dataset,
        config,
        "land_surface_water_flow__ground_elevation",
        Routing;
        sel = indices,
    )
    n = length(domain.land.network.indices)

    zx_max_at_edge = zeros(n)
    zy_max_at_edge = zeros(n)
    for i in 1:n
        ind_x_up = edge_indices.ind_x_up[i]
        if ind_x_up <= n
            zx_max_at_edge[i] = max(elevation[i], elevation[ind_x_up])
        end
        ind_y_up = edge_indices.ind_y_up[i]
        if ind_y_up <= n
            zy_max_at_edge[i] = max(elevation[i], elevation[ind_y_up])
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
        xwidth_at_edge = we_x,
        ywidth_at_edge = we_y,
        theta,
        h_thresh = waterdepth_threshold,
        zx_max_at_edge,
        zy_max_at_edge,
        mannings_n_sq_at_edge = mannings_n .* mannings_n,
        z = elevation,
        froude_limit,
    )
    return parameters
end

"Struct to store local inertial overland flow model boundary conditions"
@with_kw struct LocalInertialOverlandFlowBC <: AbstractOverlandFlowBC
    n::Int
    # runoff from hydrological model [m³ s⁻¹]
    runoff::Vector{Float64} = zeros(n)
end

"Initialize local inertial overland flow model"
function init_local_inertial_overland_flow(
    dataset::NCDataset,
    config::Config,
    domain::Domain,
)
    alpha_coefficient = config.model.land_local_inertial_flow__alpha_coefficient # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; alpha_coefficient)

    n = length(domain.land.network.indices)
    boundary_conditions = LocalInertialOverlandFlowBC(; n)
    parameters = LocalInertialOverlandFlowParameters(dataset, config, domain)
    variables = LocalInertialOverlandFlowVariables(; n)
    routing_method = LocalInertial()

    overland_flow_model = OverlandFlowModel(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        routing_method,
    )

    return overland_flow_model
end

"""
    stable_timestep(river_flow_model::RiverFlowModel{<:LocalInertial}, flow_length::Vector{Float64})
    stable_timestep(overland_flow_model::OverlandFlowModel{<:LocalInertial}, parameters::LandParameters)

Compute stable timestep size for the local inertial approach, based on Bates et al. (2010).

dt = α * (Δx / sqrt(g max(h))
"""
function stable_timestep(
    river_flow_model::RiverFlowModel{<:LocalInertial},
    flow_length::Vector{Float64},
)
    dt_min = Inf
    dt_min_default = 60.0
    (; alpha_coefficient) = river_flow_model.timestepping
    (; n) = river_flow_model.parameters
    (; h) = river_flow_model.variables
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(n)
        @fastmath @inbounds dt =
            alpha_coefficient * flow_length[i] / sqrt(GRAVITATIONAL_ACCELERATION * h[i])
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? dt_min_default : dt_min
    return dt_min
end

function stable_timestep(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    parameters::LandParameters,
)
    dt_min = Inf
    dt_min_default = 60.0
    (; alpha_coefficient) = overland_flow_model.timestepping
    (; n) = overland_flow_model.parameters
    (; x_length, y_length, river_location) = parameters
    (; h) = overland_flow_model.variables
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(n)
        @fastmath @inbounds dt = if river_location[i] == 0
            alpha_coefficient * min(x_length[i], y_length[i]) /
            sqrt(GRAVITATIONAL_ACCELERATION * h[i])
        else
            Inf
        end
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? dt_min_default : dt_min
    return dt_min
end

"""
    stable_timestep(river_flow_model::RiverFlowModel{<:ManningStaggered}, _)

Compute stable timestep size for a river flow model using the Manning's equation on staggered grid.
"""
function stable_timestep(river_flow_model::RiverFlowModel{<:ManningStaggered}, _)
    dt_min = Inf
    dt_min_default = 60.0
    (; alpha_coefficient) = river_flow_model.timestepping
    (;
        active_e,
        mannings_n_at_edge,
        flow_length_at_edge,
        flow_width_at_edge,
        slope_at_edge,
    ) = river_flow_model.parameters
    (; water_depth_at_edge) = river_flow_model.variables

    @batch per = thread reduction = ((min, dt_min),) for i in active_e
        @fastmath @inbounds h_r =
            (flow_width_at_edge[i] * water_depth_at_edge[i]) /
            wetted_perimeter_channel(water_depth_at_edge[i], flow_width_at_edge[i])
        # estimate wave celerity for wide rectangular channel
        celerity =
            cbrt(h_r^2) * sqrt(slope_at_edge[i]) / mannings_n_at_edge[i] / BETA_KINWAVE
        dt = alpha_coefficient * flow_length_at_edge[i] / celerity
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? dt_min_default : dt_min
    return dt_min
end

"""
Update boundary condition `runoff` local inertial overland flow model for a single timestep.
"""
function update_bc_overland_flow_model!(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    external_models::NamedTuple,
    domain::Domain,
    dt::Float64,
)
    (; soil, runoff, subsurface_flow) = external_models
    (; net_runoff) = soil.variables
    (; net_runoff_river) = runoff.variables
    (; area) = domain.land.parameters
    river_indices = domain.river.network.land_indices

    @. overland_flow_model.boundary_conditions.runoff =
        (net_runoff + net_runoff_river) * area
    overland_flow_model.boundary_conditions.runoff[river_indices] .+=
        get_flux_to_river(subsurface_flow, river_indices)
    return nothing
end

"""
Update subsurface flow contribution to inflow of a reservoir model for a river flow model on
a staggered grid for a single timestep.
"""
function update_inflow!(
    reservoir_model::ReservoirModel,
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    subsurface_flow_model::AbstractSubsurfaceFlowModel,
    network::NetworkReservoir,
)
    (; land_indices) = network
    (; inflow_subsurface) = reservoir_model.boundary_conditions
    inflow_subsurface .=
        get_inflow_reservoir(river_flow_model, subsurface_flow_model, land_indices)
    return nothing
end
update_inflow!(
    ::Nothing,
    ::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    ::AbstractSubsurfaceFlowModel,
    ::NetworkReservoir,
) = nothing

"""
Helper function to set flow variables of the local inertial overland flow model to zero.
This is done at the start of each simulation timestep, during the timestep the total
(weighted) sum is computed from values at each sub timestep.
"""
function set_flow_vars!(overland_flow_model::OverlandFlowModel{<:LocalInertial})
    (; qx_cumulative, qy_cumulative) = overland_flow_model.variables
    qx_cumulative .= 0.0
    qy_cumulative .= 0.0
    return nothing
end

"""
Helper function to compute average flow variables of the local inertial overland flow model.
This is done at the end of each simulation timestep.
"""
function average_flow_vars!(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    dt::Float64,
)
    (; qx_average, qx_cumulative, qy_average, qy_cumulative) = overland_flow_model.variables

    @. qx_average .= qx_cumulative / dt
    @. qy_average .= qy_cumulative / dt
    return nothing
end

"""
Update local inertial river and overland flow model for a single timestep `dt`. An adaptive
timestepping method is used (computing a sub timestep `dt_s`).
"""
function update_overland_flow_model!(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    river_flow_model::RiverFlowModel{<:LocalInertial},
    domain::Domain,
    clock::Clock,
    dt::Float64;
    update_h = false,
)
    (; reservoir) = river_flow_model.boundary_conditions
    (; flow_length) = domain.river.parameters
    land_parameters = domain.land.parameters

    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)
    set_flow_vars!(river_flow_model)
    set_flow_vars!(overland_flow_model)

    t = 0.0
    while t < dt
        dt_river = stable_timestep(river_flow_model, flow_length)
        dt_land = stable_timestep(overland_flow_model, land_parameters)
        dt_s = min(dt_river, dt_land)
        dt_s = check_timestepsize(dt_s, t, dt)

        local_inertial_update_fluxes!(overland_flow_model, domain, dt_s)
        update_inflow_reservoir!(overland_flow_model, reservoir, domain)
        staggered_scheme_river_update!(river_flow_model, domain, dt_s, update_h)
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
Update flow for the local inertial overland flow model at edge `i` in a single direction.
`is_x_direction`: true for x-direction (ind_x_up/ind_x_down), false for y-direction
(ind_y_up/ind_y_down)
"""
@inline function update_directional_flow!(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    domain::Domain,
    i::Int,
    dt::Float64,
    is_x_direction::Bool,
)
    indices = domain.land.network.edge_indices
    (; x_length, y_length) = domain.land.parameters
    land_v = overland_flow_model.variables
    land_p = overland_flow_model.parameters

    # Select direction-specific parameters based on the boolean flag
    if is_x_direction
        upstream_idx = indices.ind_x_up[i]
        downstream_idx = indices.ind_x_down[i]
        width_at_edge = land_p.ywidth_at_edge[i]
        z_max_at_edge = land_p.zx_max_at_edge[i]
        length_vec = x_length
        q_current = land_v.qx
        q_prev = land_v.qx0
        q_cumulative = land_v.qx_cumulative
    else
        upstream_idx = indices.ind_y_up[i]
        downstream_idx = indices.ind_y_down[i]
        width_at_edge = land_p.xwidth_at_edge[i]
        z_max_at_edge = land_p.zy_max_at_edge[i]
        length_vec = y_length
        q_current = land_v.qy
        q_prev = land_v.qy0
        q_cumulative = land_v.qy_cumulative
    end

    # the effective flow width is zero when the river width exceeds the cell width and
    # floodplain flow is not calculated.
    if upstream_idx <= land_p.n && width_at_edge != 0.0
        zs_current = land_p.z[i] + land_v.h[i]
        zs_upstream = land_p.z[upstream_idx] + land_v.h[upstream_idx]
        zs_max_at_edge = max(zs_current, zs_upstream)
        water_depth_at_edge = (zs_max_at_edge - z_max_at_edge)

        if water_depth_at_edge > land_p.h_thresh
            length_at_edge = 0.5 * (length_vec[i] + length_vec[upstream_idx]) # can be precalculated
            q_current[i] = local_inertial_flow(
                land_p.theta,
                q_prev[i],
                q_prev[downstream_idx],
                q_prev[upstream_idx],
                zs_current,
                zs_upstream,
                water_depth_at_edge,
                width_at_edge,
                length_at_edge,
                land_p.mannings_n_sq_at_edge[i],
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
        q_cumulative[i] += q_current[i] * dt
    end
    return nothing
end

"""
Update fluxes for local inertial overland flow model for a single timestep `dt`.
"""
function local_inertial_update_fluxes!(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    domain::Domain,
    dt::Float64,
)
    land_v = overland_flow_model.variables
    land_p = overland_flow_model.parameters

    land_v.qx0 .= land_v.qx
    land_v.qy0 .= land_v.qy

    @batch per = thread minbatch = 6000 for i in 1:(land_p.n)
        # update qx (x-direction)
        update_directional_flow!(overland_flow_model, domain, i, dt, true)

        # update qy (y-direction)
        update_directional_flow!(overland_flow_model, domain, i, dt, false)
    end
    return nothing
end

"""
Update boundary condition inflow to a reservoir from land of the local inertial river and
overland flow model for a single timestep.
"""
function update_inflow_reservoir!(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    reservoir_model::Union{ReservoirModel, Nothing},
    domain::Domain,
)
    indices = domain.land.network.edge_indices
    reservoir_indices = domain.reservoir.network.land_indices
    land_bc = overland_flow_model.boundary_conditions
    land_v = overland_flow_model.variables

    for (i, j) in enumerate(reservoir_indices)
        ind_y_down = indices.ind_y_down[j]
        ind_x_down = indices.ind_x_down[j]
        reservoir_model.boundary_conditions.inflow_overland[i] =
            land_bc.runoff[j] +
            (land_v.qx[ind_x_down] - land_v.qx[j] + land_v.qy[ind_y_down] - land_v.qy[j])
    end
    return nothing
end

"""
Compute storage change at node index `i` containing a river, from river and overland flow
fluxes of the local inertial river and overland flow model.
"""
@inline function compute_river_storage_change(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    river_flow_model::RiverFlowModel{<:LocalInertial},
    domain::Domain,
    i::Int,
    dt::Float64,
)
    indices = domain.land.network.edge_indices
    inds_river = domain.land.network.river_indices
    edges_at_node = domain.river.network.edges_at_node
    river_idx = inds_river[i]

    ind_y_down = indices.ind_y_down[i]
    ind_x_down = indices.ind_x_down[i]

    net_river_flow =
        sum_at(river_flow_model.variables.q, edges_at_node.src[river_idx]) -
        sum_at(river_flow_model.variables.q, edges_at_node.dst[river_idx])
    net_land_flow =
        overland_flow_model.variables.qx[ind_x_down] - overland_flow_model.variables.qx[i] +
        overland_flow_model.variables.qy[ind_y_down] - overland_flow_model.variables.qy[i]
    net_flow =
        net_river_flow + net_land_flow + overland_flow_model.boundary_conditions.runoff[i] -
        river_flow_model.boundary_conditions.abstraction[river_idx]
    storage_change = net_flow * dt

    return storage_change
end

"""
Compute external inflow at river node index `river_idx`, including negative inflow
(abstraction), of the local inertial river and overland flow model.
Returns tuple: (inflow, abstraction_to_add)
"""
@inline function compute_external_inflow(
    river_flow_model::RiverFlowModel{<:LocalInertial},
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    i::Int,
    river_idx::Int,
    dt::Float64,
)
    if river_flow_model.boundary_conditions.external_inflow[river_idx] < 0.0
        available_volume =
            if overland_flow_model.variables.storage[i] >=
               river_flow_model.parameters.bankfull_storage[river_idx]
                river_flow_model.parameters.bankfull_storage[river_idx]
            else
                river_flow_model.variables.storage[river_idx]
            end
        _abstraction = min(
            -river_flow_model.boundary_conditions.external_inflow[river_idx],
            available_volume / dt * 0.80,
        )
        return (-_abstraction, _abstraction)
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
    i::Int,
    river::RiverFlowModel{<:LocalInertial},
    domain::Domain,
)
    if total_storage >= river.parameters.bankfull_storage[river_idx]
        # Storage exceeds bankfull capacity - water spills onto floodplain
        river_h =
            river.parameters.bankfull_depth[river_idx] +
            (total_storage - river.parameters.bankfull_storage[river_idx]) /
            (domain.land.parameters.x_length[i] * domain.land.parameters.y_length[i])
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
Compute storage change for local inertial overland flow model at node index `i` (non-river)
from horizontal fluxes and runoff.
"""
@inline function compute_land_storage_change(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    network::NetworkLand,
    i::Int,
    dt::Float64,
)
    indices = network.edge_indices
    ind_y_down = indices.ind_y_down[i]
    ind_x_down = indices.ind_x_down[i]

    return (
        overland_flow_model.variables.qx[ind_x_down] - overland_flow_model.variables.qx[i] +
        overland_flow_model.variables.qy[ind_y_down] - overland_flow_model.variables.qy[i] +
        overland_flow_model.boundary_conditions.runoff[i]
    ) * dt
end

"""
Update storage and water depth for the local inertial overland flow model at node index `i`
containing a river.
"""
@inline function update_river_and_land_storage_and_depth!(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    river_flow_model::RiverFlowModel{<:LocalInertial},
    domain::Domain,
    i::Int,
    dt::Float64,
)
    inds_river = domain.land.network.river_indices
    river_idx = inds_river[i]

    # Compute and apply storage change from fluxes
    storage_change =
        compute_river_storage_change(overland_flow_model, river_flow_model, domain, i, dt)
    overland_flow_model.variables.storage[i] += storage_change

    # Handle negative storage
    if overland_flow_model.variables.storage[i] < 0.0
        overland_flow_model.variables.error[i] +=
            abs(overland_flow_model.variables.storage[i])
        overland_flow_model.variables.storage[i] = 0.0 # set storage to zero
    end

    # Compute and apply external inflow
    inflow, abstraction_to_add =
        compute_external_inflow(river_flow_model, overland_flow_model, i, river_idx, dt)
    overland_flow_model.variables.storage[i] += inflow * dt
    river_flow_model.boundary_conditions.actual_external_abstraction_cumulative[river_idx] +=
        abstraction_to_add * dt

    # Compute and apply water depths
    river_h, land_h, river_storage = compute_water_depths(
        overland_flow_model.variables.storage[i],
        river_idx,
        i,
        river_flow_model,
        domain,
    )
    river_flow_model.variables.h[river_idx] = river_h
    overland_flow_model.variables.h[i] = land_h
    river_flow_model.variables.storage[river_idx] = river_storage

    return nothing
end

"""
Update storage and water depth for local inertial overland flow model at node index `i`
(non-river).
"""
@inline function update_land_storage_and_depth!(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    domain::DomainLand,
    i::Int,
    dt::Float64,
)
    # Compute and apply storage change
    storage_change = compute_land_storage_change(overland_flow_model, domain.network, i, dt)
    overland_flow_model.variables.storage[i] += storage_change

    # Handle negative storage
    if overland_flow_model.variables.storage[i] < 0.0
        overland_flow_model.variables.error[i] +=
            abs(overland_flow_model.variables.storage[i])
        overland_flow_model.variables.storage[i] = 0.0 # set storage to zero
    end

    # Update water depth
    overland_flow_model.variables.h[i] =
        overland_flow_model.variables.storage[i] /
        (domain.parameters.x_length[i] * domain.parameters.y_length[i])

    return nothing
end

"""
Update storage and water depth for local inertial river and overland flow model for a single
timestep `dt`.
"""
function local_inertial_update_water_depth!(
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    river_flow_model::RiverFlowModel{<:LocalInertial},
    domain::Domain,
    dt::Float64,
)
    (; river_location, reservoir_outlet) = domain.land.parameters

    @batch per = thread minbatch = 6000 for i in 1:(overland_flow_model.parameters.n)
        if river_location[i]
            # Process river cells (excluding reservoir outlets)
            if !reservoir_outlet[i]
                update_river_and_land_storage_and_depth!(
                    overland_flow_model,
                    river_flow_model,
                    domain,
                    i,
                    dt,
                )
            end
        else
            # Process land cells (non-river)
            update_land_storage_and_depth!(overland_flow_model, domain.land, i, dt)
        end
    end
    return nothing
end
