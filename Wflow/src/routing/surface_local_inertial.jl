abstract type AbstractFloodPlain end

"Struct for storing local inertial river flow model parameters"
@with_kw struct LocalInertialRiverFlowParameters
    n::Int                                  # number of cells [-]
    ne::Int                                 # number of edges [-]
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

    (; pit_indices, indices, graph, local_drain_direction, nodes_at_edge) = domain.network
    (; flow_width, flow_length, reservoir_outlet) = domain.parameters

    riverlength_bc = ncread(
        dataset,
        config,
        "model_boundary_condition_river__length";
        sel = pit_indices,
        defaults = 1.0e04,
        type = Float64,
    )
    bankfull_elevation_2d = ncread(
        dataset,
        config,
        "river_bank_water__elevation";
        optional = false,
        type = Float64,
        fill = 0,
    )
    bankfull_depth_2d = ncread(
        dataset,
        config,
        "river_bank_water__depth";
        optional = false,
        type = Float64,
        fill = 0,
    )
    bankfull_depth = bankfull_depth_2d[indices]
    zb = bankfull_elevation_2d[indices] - bankfull_depth # river bed elevation

    bankfull_storage = bankfull_depth .* flow_width .* flow_length
    mannings_n = ncread(
        dataset,
        config,
        "river_water_flow__manning_n_parameter";
        sel = indices,
        defaults = 0.036,
        type = Float64,
    )

    n = length(indices)
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
    zb_max = fill(Float64(0), n_edges)
    width_at_edge = fill(Float64(0), n_edges)
    length_at_edge = fill(Float64(0), n_edges)
    mannings_n_sq = fill(Float64(0), n_edges)
    for i in 1:n_edges
        src_node = nodes_at_edge.src[i]
        dst_node = nodes_at_edge.dst[i]
        zb_max[i] = max(zb[src_node], zb[dst_node])
        width_at_edge[i] = min(flow_width[src_node], flow_width[dst_node])
        length_at_edge[i] = 0.5 * (flow_length[dst_node] + flow_length[src_node])
        mannings_n_i =
            (
                mannings_n[dst_node] * flow_length[dst_node] +
                mannings_n[src_node] * flow_length[src_node]
            ) / (flow_length[dst_node] + flow_length[src_node])
        mannings_n_sq[i] = mannings_n_i * mannings_n_i
    end
    active_index = findall(x -> x == 0, reservoir_outlet)

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
    n::Int
    n_edges::Int
    q::Vector{Float64} = zeros(n_edges)                        # river discharge at edge (subgrid channel) [m³ s⁻¹]
    q0::Vector{Float64} = zeros(n_edges)                       # river discharge at edge (subgrid channel) at previous time step [m³ s⁻¹]
    q_av::Vector{Float64}                     # average river channel (+ floodplain) discharge at edge [m³ s⁻¹] (model timestep Δt)
    q_channel_av::Vector{Float64}             # average river channel discharge at edge [m³ s⁻¹] (for model timestep Δt)
    h::Vector{Float64}                        # water depth [m]
    zs_max::Vector{Float64} = zeros(n_edges)  # maximum water elevation at edge [m]
    zs_src::Vector{Float64} = zeros(n_edges)  # water elevation of source node of edge [m]
    zs_dst::Vector{Float64} = zeros(n_edges)  # water elevation of downstream node of edge [m]
    hf::Vector{Float64} = zeros(n_edges)       # water depth at edge [m]
    a::Vector{Float64} = zeros(n_edges)        # flow area at edge [m²]
    r::Vector{Float64} = zeros(n_edges)        # wetted perimeter at edge [m]
    storage::Vector{Float64} = zeros(n)        # river storage [m³]
    error::Vector{Float64} = zeros(n)          # error storage [m³]
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
        "model_boundary_condition_river_bank_water__depth";
        sel = pit_indices,
        defaults = 0.0,
        type = Float64,
    )

    n = length(indices)
    n_edges = ne(graph)
    # set river depth h to zero (including reservoir locations)
    h = zeros(n)
    q_av = zeros(n_edges)
    # set ghost points for boundary condition (downstream river outlet): river depth `h`
    append!(h, riverdepth_bc)
    variables = LocalInertialRiverFlowVariables(;
        n,
        n_edges,
        q_av,
        q_channel_av = config.model.floodplain_1d__flag ? zeros(n_edges) : q_av,
        h,
    )
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
    floodplain::F                                       # Floodplain (1D) schematization
    allocation::A                                       # Water allocation
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
) = model.variables.q_av[inds] .+ model.variables.to_river[inds]

get_inflow_reservoir(::LocalInertialRiverFlow, model::LateralSSF, inds::Vector{Int}) =
    (model.variables.ssf[inds] .+ model.variables.to_river[inds]) ./ tosecond(BASETIMESTEP)

"Update local inertial river flow model `LocalIntertialRiverFlow` for a single timestep"
function local_inertial_river_update!(
    model::LocalInertialRiverFlow,
    domain::Domain,
    dt::Float64,
    dt_forcing::Float64,
    update_h::Bool,
)
    (; nodes_at_edge, edges_at_node) = domain.river.network
    (; flow_length, flow_width) = domain.river.parameters
    (; inwater, abstraction, external_inflow, actual_external_abstraction_av) =
        model.boundary_conditions

    river_v = model.variables
    river_p = model.parameters

    river_v.q0 .= river_v.q
    if !isnothing(model.floodplain)
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

        river_v.a[i] = river_p.flow_width_at_edge[i] * river_v.hf[i] # flow area (rectangular channel)
        river_v.r[i] = river_v.a[i] / (river_p.flow_width_at_edge[i] + 2.0 * river_v.hf[i]) # hydraulic radius (rectangular channel)

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
        # average river discharge (here accumulated for model timestep Δt)
        river_v.q_av[i] += river_v.q[i] * dt
    end
    if !isnothing(model.floodplain)
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

            a_src = flow_area(
                floodplain_p.profile.width[i2, i_src],
                floodplain_p.profile.a[i1, i_src],
                floodplain_p.profile.depth[i1],
                floodplain_v.hf[i],
            )
            a_src = max(a_src - (floodplain_v.hf[i] * flow_width[i_src]), 0.0)

            a_dst = flow_area(
                floodplain_p.profile.width[i2, i_dst],
                floodplain_p.profile.a[i1, i_dst],
                floodplain_p.profile.depth[i1],
                floodplain_v.hf[i],
            )
            a_dst = max(a_dst - (floodplain_v.hf[i] * flow_width[i_dst]), 0.0)

            floodplain_v.a[i] = min(a_src, a_dst)

            floodplain_v.r[i] = ifelse(
                a_src < a_dst,
                a_src / wetted_perimeter(
                    floodplain_p.profile.p[i1, i_src],
                    floodplain_p.profile.depth[i1],
                    floodplain_v.hf[i],
                ),
                a_dst / wetted_perimeter(
                    floodplain_p.profile.p[i1, i_dst],
                    floodplain_p.profile.depth[i1],
                    floodplain_v.hf[i],
                ),
            )

            floodplain_v.q[i] = ifelse(
                floodplain_v.a[i] > 1.0e-05,
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
                ),
                0.0,
            )

            # limit floodplain q in case water is not available
            floodplain_v.q[i] = ifelse(
                floodplain_v.h[i_src] <= 0.0,
                min(floodplain_v.q[i], 0.0),
                floodplain_v.q[i],
            )
            floodplain_v.q[i] = ifelse(
                floodplain_v.h[i_dst] <= 0.0,
                max(floodplain_v.q[i], 0.0),
                floodplain_v.q[i],
            )

            floodplain_v.q[i] =
                ifelse(floodplain_v.q[i] * river_v.q[i] < 0.0, 0.0, floodplain_v.q[i])
            # average floodplain discharge (here accumulated for model timestep Δt)
            floodplain_v.q_av[i] += floodplain_v.q[i] * dt
        end
    end
    # For reservoir locations the local inertial solution is replaced by the reservoir
    # model. These locations are handled as boundary conditions in the local inertial model
    # (fixed h).
    (; reservoir) = model.boundary_conditions
    inds_reservoir = domain.reservoir.network.river_indices
    if !isnothing(reservoir)
        res_bc = reservoir.boundary_conditions
    end
    for v in eachindex(inds_reservoir)
        i = inds_reservoir[v]

        q_in = get_inflow_reservoir(model, edges_at_node.src[i])
        # If external_inflow < 0, abstraction is limited
        if res_bc.external_inflow[v] < 0.0
            _abstraction = min(
                -res_bc.external_inflow[v],
                (reservoir.variables.storage[v] / dt) * 0.98,
            )
            res_bc.actual_external_abstraction_av[v] += _abstraction * dt
            _inflow = -_abstraction
        else
            _inflow = res_bc.external_inflow[v]
        end
        net_inflow =
            q_in + res_bc.inflow_overland[v] + res_bc.inflow_subsurface[v] + _inflow
        update!(reservoir, v, net_inflow, dt, dt_forcing)
        river_v.q[i] = reservoir.variables.outflow[v]
        # average river discharge (here accumulated for model timestep Δt)
        river_v.q_av[i] += river_v.q[i] * dt
    end
    if update_h
        @batch per = thread minbatch = 1000 for i in river_p.active_n
            q_src = sum_at(river_v.q, edges_at_node.src[i])
            q_dst = sum_at(river_v.q, edges_at_node.dst[i])
            # internal abstraction (water demand) is limited by river storage and negative
            # external inflow as part of water allocation computations.
            river_v.storage[i] =
                river_v.storage[i] + (q_src - q_dst + inwater[i] - abstraction[i]) * dt

            if river_v.storage[i] < 0.0
                river_v.error[i] = river_v.error[i] + abs(river_v.storage[i])
                river_v.storage[i] = 0.0 # set storage to zero
            end
            # limit negative external inflow
            if external_inflow[i] < 0.0
                _abstraction = min(-external_inflow[i], river_v.storage[i] / dt * 0.80)
                actual_external_abstraction_av[i] += _abstraction * dt
                _inflow = -_abstraction
            else
                _inflow = external_inflow[i]
            end
            river_v.storage[i] += _inflow * dt # add external inflow
            river_v.h[i] = river_v.storage[i] / (flow_length[i] * flow_width[i])

            if !isnothing(model.floodplain)
                floodplain_v = model.floodplain.variables
                floodplain_p = model.floodplain.parameters
                q_src = sum_at(floodplain_v.q, edges_at_node.src[i])
                q_dst = sum_at(floodplain_v.q, edges_at_node.dst[i])
                floodplain_v.storage[i] = floodplain_v.storage[i] + (q_src - q_dst) * dt
                if floodplain_v.storage[i] < 0.0
                    floodplain_v.error[i] =
                        floodplain_v.error[i] + abs(floodplain_v.storage[i])
                    floodplain_v.storage[i] = 0.0
                end
                storage_total = river_v.storage[i] + floodplain_v.storage[i]
                if storage_total > river_p.bankfull_storage[i]
                    flood_storage = storage_total - river_p.bankfull_storage[i]
                    h = flood_depth(floodplain_p.profile, flood_storage, flow_length[i], i)
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
        end
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
    clock::Clock;
    update_h = true,
)
    (; reservoir) = model.boundary_conditions
    (; flow_length) = domain.river.parameters

    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)

    if !isnothing(model.floodplain)
        model.floodplain.variables.q_av .= 0.0
    end
    set_flow_vars!(model)

    dt = tosecond(clock.dt)
    t = 0.0
    while t < dt
        dt_s = stable_timestep(model, flow_length)
        dt_s = check_timestepsize(dt_s, t, dt)
        local_inertial_river_update!(model, domain, dt_s, dt, update_h)
        t += dt_s
    end
    average_flow_vars!(model, dt)
    average_reservoir_vars!(reservoir, dt)

    if !isnothing(model.floodplain)
        model.floodplain.variables.q_av ./= dt
        model.variables.q_channel_av .= model.variables.q_av
        model.variables.q_av .=
            model.variables.q_channel_av .+ model.floodplain.variables.q_av
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
    qx_av::Vector{Float64} = zeros(n + 1)
    # flow in y direction at edge [m³ s⁻¹]
    qy::Vector{Float64} = zeros(n + 1)
    # average flow in y direction at edge [m³ s⁻¹] for model timestep Δt
    qy_av::Vector{Float64} = zeros(n + 1)
    # total storage of cell [m³] (including river storage for river cells)
    storage::Vector{Float64} = zeros(n)
    # error storage [m³]
    error::Vector{Float64} = zeros(n)
    # water depth of cell [m] (for river cells the reference is the river bed elevation `zb`)
    h::Vector{Float64} = zeros(n)
end

"Struct to store local inertial overland flow model parameters"
@with_kw struct LocalInertialOverlandFlowParameters
    n::Int                              # number of cells [-]
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

    (; edge_indices, indices) = domain.land.network
    (; x_length, y_length) = domain.land.parameters

    @info "Local inertial approach is used for overland flow." alpha theta waterdepth_threshold froude_limit

    mannings_n = ncread(
        dataset,
        config,
        "land_surface_water_flow__manning_n_parameter";
        sel = indices,
        defaults = 0.072,
        type = Float64,
    )
    elevation_2d = ncread(
        dataset,
        config,
        "land_surface_water_flow__ground_elevation";
        optional = false,
        type = Float64,
        fill = 0,
    )
    elevation = elevation_2d[indices]
    n = length(domain.land.network.indices)

    zx_max = fill(Float64(0), n)
    zy_max = fill(Float64(0), n)
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
    runoff::Vector{Float64} = zeros(n) # runoff from hydrological model [m³ s⁻¹]
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

    @. model.boundary_conditions.runoff =
        net_runoff / 1000.0 * area / dt + net_runoff_river * area * 0.001 / dt
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
    qx_av .= 0.0
    qy_av .= 0.0
    return nothing
end

"""
Helper function to compute average flow variables of the `LocalInertialOverlandFlow` model.
This is done at the end of each simulation timestep.
"""
function average_flow_vars!(model::LocalInertialOverlandFlow, dt::Float64)
    (; qx_av, qy_av) = model.variables
    qx_av ./= dt
    qy_av ./= dt
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
    clock::Clock;
    update_h = false,
)
    (; reservoir) = river.boundary_conditions
    (; flow_length) = domain.river.parameters
    (; parameters) = domain.land

    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)
    set_flow_vars!(river)
    set_flow_vars!(land)

    dt = tosecond(clock.dt)
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
    average_flow_vars!(river, dt)
    average_flow_vars!(land, dt)
    average_reservoir_vars!(reservoir, dt)

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
    indices = domain.land.network.edge_indices
    (; x_length, y_length) = domain.land.parameters
    land_v = land.variables
    land_p = land.parameters

    land_v.qx0 .= land_v.qx
    land_v.qy0 .= land_v.qy

    # update qx
    @batch per = thread minbatch = 6000 for i in 1:(land_p.n)
        yu = indices.yu[i]
        yd = indices.yd[i]
        xu = indices.xu[i]
        xd = indices.xd[i]

        # the effective flow width is zero when the river width exceeds the cell width (dy
        # for flow in x dir) and floodplain flow is not calculated.
        if xu <= land_p.n && land_p.ywidth[i] != 0.0
            zs_x = land_p.z[i] + land_v.h[i]
            zs_xu = land_p.z[xu] + land_v.h[xu]
            zs_max = max(zs_x, zs_xu)
            hf = (zs_max - land_p.zx_max[i])

            if hf > land_p.h_thresh
                length = 0.5 * (x_length[i] + x_length[xu]) # can be precalculated
                land_v.qx[i] = local_inertial_flow(
                    land_p.theta,
                    land_v.qx0[i],
                    land_v.qx0[xd],
                    land_v.qx0[xu],
                    zs_x,
                    zs_xu,
                    hf,
                    land_p.ywidth[i],
                    length,
                    land_p.mannings_n_sq[i],
                    land_p.froude_limit,
                    dt,
                )
                # limit qx in case water is not available
                if land_v.h[i] <= 0.0
                    land_v.qx[i] = min(land_v.qx[i], 0.0)
                end
                if land_v.h[xu] <= 0.0
                    land_v.qx[i] = max(land_v.qx[i], 0.0)
                end
            else
                land_v.qx[i] = 0.0
            end
            land_v.qx_av[i] += land_v.qx[i] * dt
        end

        # update qy

        # the effective flow width is zero when the river width exceeds the cell width (dx
        # for flow in y dir) and floodplain flow is not calculated.
        if yu <= land_p.n && land_p.xwidth[i] != 0.0
            zs_y = land_p.z[i] + land_v.h[i]
            zs_yu = land_p.z[yu] + land_v.h[yu]
            zs_max = max(zs_y, zs_yu)
            hf = (zs_max - land_p.zy_max[i])

            if hf > land_p.h_thresh
                length = 0.5 * (y_length[i] + y_length[yu]) # can be precalculated
                land_v.qy[i] = local_inertial_flow(
                    land_p.theta,
                    land_v.qy0[i],
                    land_v.qy0[yd],
                    land_v.qy0[yu],
                    zs_y,
                    zs_yu,
                    hf,
                    land_p.xwidth[i],
                    length,
                    land_p.mannings_n_sq[i],
                    land_p.froude_limit,
                    dt,
                )
                # limit qy in case water is not available
                if land_v.h[i] <= 0.0
                    land_v.qy[i] = min(land_v.qy[i], 0.0)
                end
                if land_v.h[yu] <= 0.0
                    land_v.qy[i] = max(land_v.qy[i], 0.0)
                end
            else
                land_v.qy[i] = 0.0
            end
            land_v.qy_av[i] += land_v.qy[i] * dt
        end
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
Update storage and water depth for combined river `LocalInertialRiverFlow`and overland flow
`LocalInertialOverlandFlow` models for a single timestep `dt`.
"""
function local_inertial_update_water_depth!(
    land::LocalInertialOverlandFlow,
    river::LocalInertialRiverFlow,
    domain::Domain,
    dt::Float64,
)
    indices = domain.land.network.edge_indices
    inds_river = domain.land.network.river_indices
    (; edges_at_node) = domain.river.network
    (; river_location, reservoir_outlet, x_length, y_length) = domain.land.parameters
    (; flow_width, flow_length) = domain.river.parameters

    river_bc = river.boundary_conditions
    river_v = river.variables
    river_p = river.parameters
    land_bc = land.boundary_conditions
    land_v = land.variables
    land_p = land.parameters

    # change in storage and water depth based on horizontal fluxes for river and land cells
    @batch per = thread minbatch = 6000 for i in 1:(land_p.n)
        yd = indices.yd[i]
        xd = indices.xd[i]

        if river_location[i]
            # reservoir locations are boundary points (update storage and h not required)
            reservoir_outlet[i] && continue
            # internal abstraction (water demand) is limited by river storage and negative
            # external inflow as part of water allocation computations.
            land_v.storage[i] +=
                (
                    sum_at(river_v.q, edges_at_node.src[inds_river[i]]) -
                    sum_at(river_v.q, edges_at_node.dst[inds_river[i]]) + land_v.qx[xd] -
                    land_v.qx[i] + land_v.qy[yd] - land_v.qy[i] + land_bc.runoff[i] -
                    river_bc.abstraction[inds_river[i]]
                ) * dt
            if land_v.storage[i] < 0.0
                land_v.error[i] = land_v.error[i] + abs(land_v.storage[i])
                land_v.storage[i] = 0.0 # set storage to zero
            end
            # limit negative external inflow
            if river_bc.external_inflow[inds_river[i]] < 0.0
                available_volume =
                    if land_v.storage[i] >= river_p.bankfull_storage[inds_river[i]]
                        river_p.bankfull_depth[inds_river[i]]
                    else
                        river_v.storage[inds_river[i]]
                    end
                _abstraction = min(
                    -river_bc.external_inflow[inds_river[i]],
                    available_volume / dt * 0.80,
                )
                river_bc.actual_external_abstraction_av[inds_river[i]] += _abstraction * dt
                _inflow = -_abstraction
            else
                _inflow = river_bc.external_inflow[inds_river[i]]
            end
            land_v.storage[i] += _inflow * dt
            if land_v.storage[i] >= river_p.bankfull_storage[inds_river[i]]
                river_v.h[inds_river[i]] =
                    river_p.bankfull_depth[inds_river[i]] +
                    (land_v.storage[i] - river_p.bankfull_storage[inds_river[i]]) /
                    (x_length[i] * y_length[i])
                land_v.h[i] =
                    river_v.h[inds_river[i]] - river_p.bankfull_depth[inds_river[i]]
                river_v.storage[inds_river[i]] =
                    river_v.h[inds_river[i]] *
                    flow_length[inds_river[i]] *
                    flow_width[inds_river[i]]
            else
                river_v.h[inds_river[i]] =
                    land_v.storage[i] /
                    (flow_length[inds_river[i]] * flow_width[inds_river[i]])
                land_v.h[i] = 0.0
                river_v.storage[inds_river[i]] = land_v.storage[i]
            end
        else
            land_v.storage[i] +=
                (
                    land_v.qx[xd] - land_v.qx[i] + land_v.qy[yd] - land_v.qy[i] +
                    land_bc.runoff[i]
                ) * dt
            if land_v.storage[i] < 0.0
                land_v.error[i] = land_v.error[i] + abs(land_v.storage[i])
                land_v.storage[i] = 0.0 # set storage to zero
            end
            land_v.h[i] = land_v.storage[i] / (x_length[i] * y_length[i])
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
    (; indices) = domain.network
    (; flow_width, flow_length) = domain.parameters
    storage = ncread(
        dataset,
        config,
        "floodplain_water__sum_of_volume_per_depth";
        optional = false,
        sel = indices,
        type = Float64,
        dimname = :flood_depth,
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
    (; indices, nodes_at_edge, graph) = domain.network
    (; flow_length) = domain.parameters
    n_edges = ne(graph)
    profile = FloodPlainProfile(dataset, config, domain, index_pit)

    mannings_n = ncread(
        dataset,
        config,
        "floodplain_water_flow__manning_n_parameter";
        sel = indices,
        defaults = 0.072,
        type = Float64,
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
    storage::Vector{Float64} = zeros(n)    # storage [m³]
    h::Vector{Float64}                     # water depth [m]
    error::Vector{Float64} = zeros(n)      # error storage [m³]
    a::Vector{Float64} = zeros(n_edges)    # flow area at egde [m²]
    r::Vector{Float64} = zeros(n_edges)    # hydraulic radius at edge [m]
    hf::Vector{Float64} = zeros(n_edges)   # water depth at edge [m]
    q0::Vector{Float64} = zeros(n_edges)   # discharge at edge at previous time step
    q::Vector{Float64} = zeros(n_edges)    # discharge at edge  [m³ s⁻¹]
    q_av::Vector{Float64} = zeros(n_edges) # average river discharge at edge  [m³ s⁻¹] for model timestep Δt
    hf_index::Vector{Int} = zeros(Int, n_edges) # edge index with `hf` [-] above depth threshold
end

"Initialize floodplain flow model variables"
function FloodPlainVariables(n::Int, n_edges::Int, index_pit::Vector{Int})
    variables = FloodPlainVariables(; n, n_edges, h = zeros(n + length(index_pit)))
    return variables
end

"Floodplain flow model"
@with_kw struct FloodPlain{P} <: AbstractFloodPlain
    parameters::FloodPlainParameters{P}
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
    ΔA = (flood_storage - profile.storage[i1, i]) / flow_length
    dh = ΔA / profile.width[i2, i]
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
    n_edges = ne(graph)
    variables = FloodPlainVariables(n, n_edges, index_pit)

    floodplain = FloodPlain(; parameters, variables)
    return floodplain
end
