"Struct for storing local inertial river flow model parameters"
@with_kw struct LocalInertialRiverFlowParameters
    n::Int                                  # number of cells [-]
    ne::Int                                 # number of edges [-]
    active_n::Vector{Int}                   # active nodes [-]
    active_e::Vector{Int}                   # active edges [-]
    g::Float64                              # acceleration due to gravity [m s⁻²]
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
    alpha = get(config.model, "river_local_inertial_flow__alpha_coefficient", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    waterdepth_threshold =
        get(config.model, "river_water_flow_threshold__depth", 1.0e-03)::Float64 # depth threshold for flow at edge
    froude_limit = get(config.model, "river_water_flow__froude_limit_flag", true)::Bool # limit flow to subcritical according to Froude number
    floodplain_1d = get(config.model, "floodplain_1d__flag", false)::Bool

    @info "Local inertial approach is used for river flow." alpha waterdepth_threshold froude_limit floodplain_1d

    (; pit_indices, indices, graph, local_drain_direction, nodes_at_edge) = domain.network
    (; flow_width, flow_length, reservoir_outlet) = domain.parameters

    lens = lens_input_parameter(config, "model_boundary_condition~river__length")
    riverlength_bc =
        ncread(dataset, config, lens; sel = pit_indices, defaults = 1.0e04, type = Float64)
    lens = lens_input_parameter(config, "river_bank_water__elevation"; optional = false)
    bankfull_elevation_2d = ncread(dataset, config, lens; type = Float64, fill = 0)
    lens = lens_input_parameter(config, "river_bank_water__depth"; optional = false)
    bankfull_depth_2d = ncread(dataset, config, lens; type = Float64, fill = 0)
    bankfull_depth = bankfull_depth_2d[indices]
    zb = bankfull_elevation_2d[indices] - bankfull_depth # river bed elevation

    bankfull_storage = bankfull_depth .* flow_width .* flow_length
    lens = lens_input_parameter(config, "river_water_flow__manning_n_parameter")
    mannings_n =
        ncread(dataset, config, lens; sel = indices, defaults = 0.036, type = Float64)

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
        g = 9.80665,
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
    q::Vector{Float64}                        # river discharge at edge (subgrid channel) [m³ s⁻¹]
    q0::Vector{Float64}                       # river discharge at edge (subgrid channel) at previous time step [m³ s⁻¹]
    q_av::Vector{Float64}                     # average river channel (+ floodplain) discharge at edge [m³ s⁻¹] (model timestep Δt)
    q_channel_av::Vector{Float64}             # average river channel discharge at edge [m³ s⁻¹] (for model timestep Δt)
    h::Vector{Float64}                        # water depth [m]
    zs_max::Vector{Float64}                   # maximum water elevation at edge [m]
    zs_src::Vector{Float64}                   # water elevation of source node of edge [m]
    zs_dst::Vector{Float64}                   # water elevation of downstream node of edge [m]
    hf::Vector{Float64}                       # water depth at edge [m]
    h_av::Vector{Float64}                     # average water depth for model timestep Δt [m]
    a::Vector{Float64}                        # flow area at edge [m²]
    r::Vector{Float64}                        # wetted perimeter at edge [m]
    storage::Vector{Float64}                  # river storage [m³]
    storage_av::Vector{Float64}               # average river storage for model timestep Δt [m³]
    error::Vector{Float64}                    # error storage [m³]
end

"Initialize shallow water river flow model variables"
function LocalInertialRiverFlowVariables(
    dataset::NCDataset,
    config::Config,
    network::NetworkRiver,
)
    (; pit_indices, indices, graph) = network
    floodplain_1d = get(config.model, "floodplain_1d", false)::Bool

    lens = lens_input_parameter(config, "model_boundary_condition~river_bank_water__depth")
    riverdepth_bc =
        ncread(dataset, config, lens; sel = pit_indices, defaults = 0.0, type = Float64)

    n = length(indices)
    n_edges = ne(graph)
    # set river depth h to zero (including reservoir locations)
    h = zeros(n)
    q_av = zeros(n_edges)
    # set ghost points for boundary condition (downstream river outlet): river depth `h`
    append!(h, riverdepth_bc)
    variables = LocalInertialRiverFlowVariables(;
        q = zeros(n_edges),
        q0 = zeros(n_edges),
        q_av = q_av,
        q_channel_av = floodplain_1d ? zeros(n_edges) : q_av,
        h = h,
        zs_max = zeros(n_edges),
        zs_src = zeros(n_edges),
        zs_dst = zeros(n_edges),
        hf = zeros(n_edges),
        h_av = zeros(n),
        a = zeros(n_edges),
        r = zeros(n_edges),
        storage = zeros(n),
        storage_av = zeros(n),
        error = zeros(n),
    )
    return variables
end

"Shallow water river flow model using the local inertial method"
@with_kw struct LocalInertialRiverFlow{R, F, A} <: AbstractRiverFlowModel
    timestepping::TimeStepping
    boundary_conditions::RiverFlowBC{R}
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
    cfl = get(config.model, "river_local_inertial_flow__alpha_coefficient", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; cfl)

    parameters = LocalInertialRiverFlowParameters(dataset, config, domain)
    variables = LocalInertialRiverFlowVariables(dataset, config, domain.network)

    n = length(domain.network.indices)
    boundary_conditions = RiverFlowBC(n, reservoir)

    floodplain_1d = get(config.model, "floodplain_1d__flag", false)::Bool
    if floodplain_1d
        zb_floodplain = parameters.zb .+ parameters.bankfull_depth
        floodplain = FloodPlain(dataset, config, domain, zb_floodplain)
    else
        floodplain = nothing
    end

    do_water_demand = haskey(config.model, "water_demand")
    river_flow = LocalInertialRiverFlow(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain,
        allocation = do_water_demand ? AllocationRiver(n) : NoAllocationRiver(),
    )
    return river_flow
end

"Return the upstream inflow for a reservoir in `LocalInertialRiverFlow`"
function get_inflow_reservoir(model::LocalInertialRiverFlow, src_edge::Vector{Int})
    q_in = sum_at(model.variables.q, src_edge)
    if !isnothing(model.floodplain)
        q_in = q_in + sum_at(model.floodplain.variables.q, src_edge)
    end
    return q_in
end

# For local inertial river routing, `to_river` is included, as reservoir cells are excluded
# (boundary condition).
get_inflow_reservoir(::LocalInertialRiverFlow, model::KinWaveOverlandFlow) =
    model.variables.q_av .+ model.variables.to_river
get_inflow_reservoir(::LocalInertialRiverFlow, model::LateralSSF) =
    (model.variables.ssf .+ model.variables.to_river) ./ tosecond(BASETIMESTEP)

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
    (; inwater, abstraction, inflow, actual_external_abstraction_av) =
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
                river_p.g,
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
                    river_p.g,
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
    (; reservoir, inflow_reservoir) = model.boundary_conditions
    inds_reservoir = domain.reservoir.network.river_indices
    for v in eachindex(inds_reservoir)
        i = inds_reservoir[v]

        q_in = get_inflow_reservoir(model, edges_at_node.src[i])
        update!(reservoir, v, q_in + inflow_reservoir[i], dt, dt_forcing)
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
            if inflow[i] < 0.0
                _inflow = max(-(0.80 * river_v.storage[i] / dt), inflow[i])
                actual_external_abstraction_av[i] += _inflow * dt
            else
                _inflow = inflow[i]
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
                # average variables (here accumulated for model timestep Δt)
                floodplain_v.storage_av[i] += floodplain_v.storage[i] * dt
                floodplain_v.h_av[i] += floodplain_v.h_av[i] * dt
            end
            # average variables (here accumulated for model timestep Δt)
            river_v.storage_av[i] += river_v.storage[i] * dt
            river_v.h_av[i] += river_v.h[i] * dt
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
    (; reservoir, actual_external_abstraction_av) = model.boundary_conditions
    (; flow_length) = domain.river.parameters

    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)

    if !isnothing(model.floodplain)
        set_flow_vars!(model.floodplain.variables)
    end
    set_flow_vars!(model.variables, actual_external_abstraction_av)

    dt = tosecond(clock.dt)
    t = 0.0
    while t < dt
        dt_s = stable_timestep(model, flow_length)
        dt_s = check_timestepsize(dt_s, t, dt)
        local_inertial_river_update!(model, domain, dt_s, dt, update_h)
        t = t + dt_s
    end
    average_flow_vars!(model.variables, actual_external_abstraction_av, dt)
    average_reservoir_vars!(reservoir, dt)

    if !isnothing(model.floodplain)
        average_flow_vars!(model.floodplain.variables, dt)
        model.variables.q_channel_av .= model.variables.q_av
        model.variables.q_av .=
            model.variables.q_channel_av .+ model.floodplain.variables.q_av
    end

    return nothing
end

"Struct to store local inertial overland flow model variables"
@with_kw struct LocalInertialOverlandFlowVariables
    qy0::Vector{Float64}              # flow in y direction at edge at previous time step [m³ s⁻¹]
    qx0::Vector{Float64}              # flow in x direction at edge at previous time step [m³ s⁻¹]
    qx::Vector{Float64}               # flow in x direction at egde [m³ s⁻¹]
    qy::Vector{Float64}               # flow in y direction at edge [m³ s⁻¹]
    storage::Vector{Float64}          # total storage of cell [m³] (including river storage for river cells)
    storage_av::Vector{Float64}       # average total storage of cell [m³] (including river storage for river cells) (model timestep Δt)
    error::Vector{Float64}            # error storage [m³]
    h::Vector{Float64}                # water depth of cell [m] (for river cells the reference is the river bed elevation `zb`)
    h_av::Vector{Float64}             # average water depth [m] (for river cells the reference is the river bed elevation `zb`) (model timestep Δt)
end

"Initialize local inertial overland flow model variables"
function LocalInertialOverlandFlowVariables(n::Int)
    variables = LocalInertialOverlandFlowVariables(;
        qx0 = zeros(n + 1),
        qy0 = zeros(n + 1),
        qx = zeros(n + 1),
        qy = zeros(n + 1),
        storage = zeros(n),
        storage_av = zeros(n),
        error = zeros(n),
        h = zeros(n),
        h_av = zeros(n),
    )
    return variables
end

"Struct to store local inertial overland flow model parameters"
@with_kw struct LocalInertialOverlandFlowParameters
    n::Int                              # number of cells [-]
    xwidth::Vector{Float64}             # effective flow width x direction at edge (floodplain) [m]
    ywidth::Vector{Float64}             # effective flow width y direction at edge (floodplain) [m]
    g::Float64                          # acceleration due to gravity [m s⁻²]
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
    froude_limit =
        get(config.model, "land_surface_water_flow__froude_limit_flag", true)::Bool # limit flow to subcritical according to Froude number
    alpha = get(config.model, "land_local_inertial_flow__alpha_coefficient", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    theta = get(config.model, "land_local_inertial_flow__theta_coefficient", 0.8)::Float64 # weighting factor
    waterdepth_threshold =
        get(config.model, "land_surface_water_flow_threshold__depth", 1.0e-03)::Float64 # depth threshold for flow at edge

    (; edge_indices, indices) = domain.land.network
    (; x_length, y_length) = domain.land.parameters

    @info "Local inertial approach is used for overland flow." alpha theta waterdepth_threshold froude_limit

    lens = lens_input_parameter(config, "land_surface_water_flow__manning_n_parameter")
    mannings_n =
        ncread(dataset, config, lens; sel = indices, defaults = 0.072, type = Float64)
    lens = lens_input_parameter(
        config,
        "land_surface_water_flow__ground_elevation";
        optional = false,
    )
    elevation_2d = ncread(dataset, config, lens; type = Float64, fill = 0)
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
        g = 9.80665,
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
    runoff::Vector{Float64}           # runoff from hydrological model [m³ s⁻¹]
    inflow_reservoir::Vector{Float64} # inflow to reservoir from hydrological model [m³ s⁻¹]
end

"Struct to store shallow water overland flow model boundary conditions"
function LocalInertialOverlandFlowBC(n::Int)
    bc = LocalInertialOverlandFlowBC(; runoff = zeros(n), inflow_reservoir = zeros(n))
    return bc
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
    cfl = get(config.model, "land_local_inertial_flow__alpha_coefficient", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; cfl)

    n = length(domain.land.network.indices)
    boundary_conditions = LocalInertialOverlandFlowBC(n)
    parameters = LocalInertialOverlandFlowParameters(dataset, config, domain)
    variables = LocalInertialOverlandFlowVariables(n)

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
    (; n, g) = model.parameters
    (; h) = model.variables
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(n)
        @fastmath @inbounds dt = cfl * flow_length[i] / sqrt(g * h[i])
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? 60.0 : dt_min
    return dt_min
end

function stable_timestep(model::LocalInertialOverlandFlow, parameters::LandParameters)
    dt_min = Inf
    (; cfl) = model.timestepping
    (; n, g) = model.parameters
    (; x_length, y_length, river_location) = parameters
    (; h) = model.variables
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(n)
        @fastmath @inbounds dt = if river_location[i] == 0
            cfl * min(x_length[i], y_length[i]) / sqrt(g * h[i])
        else
            Inf
        end
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? 60.0 : dt_min
    return dt_min
end

"""
Update boundary conditions `runoff` and inflow to a reservoir from land `inflow_reservoir`
for overland flow model `LocalInertialOverlandFlow` for a single timestep.
"""
function update_boundary_conditions!(
    model::LocalInertialOverlandFlow,
    external_models::NamedTuple,
    domain::Domain,
    dt::Float64,
)
    (; river_flow, soil, subsurface_flow, runoff) = external_models
    (; inflow_reservoir) = model.boundary_conditions
    (; reservoir) = river_flow.boundary_conditions
    (; net_runoff) = soil.variables
    (; net_runoff_river) = runoff.variables

    (; area) = domain.land.parameters
    (; network) = domain.river

    model.boundary_conditions.runoff .=
        net_runoff ./ 1000.0 .* area ./ dt .+ get_flux_to_river(subsurface_flow) .+
        net_runoff_river .* area .* 0.001 ./ dt

    if !isnothing(reservoir)
        inflow_subsurface = get_inflow_reservoir(river_flow, subsurface_flow)

        @. inflow_reservoir[network.land_indices] = inflow_subsurface[network.land_indices]
    end
    return nothing
end

"""
Helper function to set storage and water depth variables of the `LocalInertialOverlandFlow`
model to zero. This is done at the start of each simulation timestep, during the timestep
the total (weighted) sum is computed from values at each sub timestep.
"""
function set_flow_vars!(variables::LocalInertialOverlandFlowVariables)
    variables.h_av .= 0.0
    variables.storage_av .= 0.0
    return nothing
end

"""
Helper function to compute average flow variables of the `LocalInertialOverlandFlow` model.
This is done at the end of each simulation timestep.
"""
function average_flow_vars!(variables::LocalInertialOverlandFlowVariables, dt::Float64)
    variables.h_av ./= dt
    variables.storage_av ./= dt
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
    (; reservoir, actual_external_abstraction_av) = river.boundary_conditions
    (; flow_length) = domain.river.parameters
    (; parameters) = domain.land

    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)
    set_flow_vars!(river.variables, actual_external_abstraction_av)
    set_flow_vars!(land.variables)

    dt = tosecond(clock.dt)
    t = 0.0
    while t < dt
        dt_river = stable_timestep(river, flow_length)
        dt_land = stable_timestep(land, parameters)
        dt_s = min(dt_river, dt_land)
        dt_s = check_timestepsize(dt_s, t, dt)

        local_inertial_update_fluxes!(land, domain, dt_s)
        update_inflow_reservoir!(land, river, domain)
        local_inertial_river_update!(river, domain, dt_s, dt, update_h)
        local_inertial_update_water_depth!(land, river, domain, dt_s)

        t = t + dt_s
    end
    average_flow_vars!(river.variables, actual_external_abstraction_av, dt)
    average_flow_vars!(land.variables, dt)
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
                    land_p.g,
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
                    land_p.g,
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
    river::LocalInertialRiverFlow,
    domain::Domain,
)
    indices = domain.land.network.edge_indices
    inds_river = domain.land.network.river_indices
    (; river_location, reservoir_outlet) = domain.land.parameters

    river_bc = river.boundary_conditions
    land_bc = land.boundary_conditions
    land_v = land.variables
    land_p = land.parameters

    @batch per = thread minbatch = 6000 for i in 1:(land_p.n)
        yd = indices.yd[i]
        xd = indices.xd[i]
        if river_location[i] && reservoir_outlet[i]
            river_bc.inflow_reservoir[inds_river[i]] =
                land_bc.inflow_reservoir[i] +
                land_bc.runoff[i] +
                (land_v.qx[xd] - land_v.qx[i] + land_v.qy[yd] - land_v.qy[i])
        end
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
            if river_bc.inflow[inds_river[i]] < 0.0
                available_volume =
                    if land_v.storage[i] >= river_p.bankfull_storage[inds_river[i]]
                        river_p.bankfull_depth[inds_river[i]]
                    else
                        river_v.storage[inds_river[i]]
                    end
                _inflow =
                    max(-(0.80 * available_volume / dt), river_bc.inflow[inds_river[i]])
                river_bc.actual_external_abstraction_av[inds_river[i]] += _inflow * dt
            else
                _inflow = river_bc.inflow[inds_river[i]]
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
            # average variables (here accumulated for model timestep Δt)
            river_v.h_av[inds_river[i]] += river_v.h[inds_river[i]] * dt
            river_v.storage_av[inds_river[i]] += river_v.storage[inds_river[i]] * dt
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
        # average variables (here accumulated for model timestep Δt)
        land_v.h_av[i] += land_v.h[i] * dt
        land_v.storage_av[i] += land_v.storage[i] * dt
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
    lens = lens_input_parameter(config, "floodplain_water__sum_of_volume-per-depth")
    storage =
        ncread(dataset, config, lens; sel = indices, type = Float64, dimname = :flood_depth)
    n = length(indices)

    # for convenience (interpolation) flood depth 0.0 m is added, with associated area (a),
    # storage, width (river width) and wetted perimeter (p).
    storage = vcat(fill(Float64(0), n)', storage)
    start_storage = storage
    flood_depths = Float64.(dataset["flood_depth"][:])
    pushfirst!(flood_depths, 0.0)
    n_depths = length(flood_depths)

    p = zeros(Float64, n_depths, n)
    a = zeros(Float64, n_depths, n)
    segment_storage = zeros(Float64, n_depths, n)
    width = zeros(Float64, n_depths, n)
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

    lens = lens_input_parameter(config, "floodplain_water_flow__manning_n_parameter")
    mannings_n =
        ncread(dataset, config, lens; sel = indices, defaults = 0.072, type = Float64)
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
    storage::Vector{Float64}        # storage [m³]
    storage_av::Vector{Float64}     # average storage for model timestep Δt [m³]
    h::Vector{Float64}              # water depth [m]
    h_av::Vector{Float64}           # average water depth [m] for model timestep Δt
    error::Vector{Float64}          # error storage [m³]
    a::Vector{Float64}              # flow area at egde [m²]
    r::Vector{Float64}              # hydraulic radius at edge [m]
    hf::Vector{Float64}             # water depth at edge [m]
    q0::Vector{Float64}             # discharge at edge at previous time step
    q::Vector{Float64}              # discharge at edge  [m³ s⁻¹]
    q_av::Vector{Float64}           # average river discharge at edge  [m³ s⁻¹] for model timestep Δt
    hf_index::Vector{Int}           # edge index with `hf` [-] above depth threshold
end

"Initialize floodplain flow model variables"
function FloodPlainVariables(n::Int, n_edges::Int, index_pit::Vector{Int})
    variables = FloodPlainVariables(;
        storage = zeros(n),
        storage_av = zeros(n),
        error = zeros(n),
        h = zeros(n + length(index_pit)),
        h_av = zeros(n),
        a = zeros(n_edges),
        r = zeros(n_edges),
        hf = zeros(n_edges),
        q = zeros(n_edges),
        q_av = zeros(n_edges),
        q0 = zeros(n_edges),
        hf_index = zeros(Int, n_edges),
    )
    return variables
end

"Floodplain flow model"
@with_kw struct FloodPlain{P}
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
    p = p + (2.0 * dh) # p at i1
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
