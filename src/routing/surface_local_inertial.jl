"Struct for storing local inertial river flow model parameters"
@get_units @grid_loc @with_kw struct LocalInertialRiverFlowParameters{T}
    n::Int                                              # number of cells [-]
    ne::Int                                             # number of edges [-]
    active_n::Vector{Int} | "-"                         # active nodes [-]
    active_e::Vector{Int} | "-" | "edge"                # active edges [-]
    g::T                                                # acceleration due to gravity [m s⁻²]
    froude_limit::Bool                                  # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
    h_thresh::T                                         # depth threshold for calculating flow [m]
    zb::Vector{T} | "m"                                 # river bed elevation   
    zb_max::Vector{T} | "m"                             # maximum channel bed elevation
    bankfull_volume::Vector{T} | "m3"                   # bankfull volume
    bankfull_depth::Vector{T} | "m"                     # bankfull depth
    mannings_n_sq::Vector{T} | "(s m-1/3)2" | "edge"    # Manning's roughness squared at edge
    mannings_n::Vector{T} | "s m-1/3"                   # Manning's roughness at node
    flow_length::Vector{T} | "m"                        # flow (river) length
    flow_length_at_edge::Vector{T} | "m" | "edge"       # flow (river) length at edge
    flow_width::Vector{T} | "m"                         # flow (river) width
    flow_width_at_edge::Vector{T} | "m" | "edge"        # flow (river) width at edge
    waterbody::Vector{Bool} | "-"                       # water body cells (reservoir or lake)
end

"Initialize local inertial river flow model parameters"
function LocalInertialRiverFlowParameters(
    dataset,
    config,
    indices;
    river_length,
    river_width,
    waterbody,
    n_edges,
    nodes_at_edge,
    index_pit,
    inds_pit,
)
    cfl = get(config.model, "inertial_flow_alpha", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    h_thresh = get(config.model, "h_thresh", 1.0e-03)::Float64 # depth threshold for flow at edge
    froude_limit = get(config.model, "froude_limit", true)::Bool # limit flow to subcritical according to Froude number
    floodplain_1d = get(config.model, "floodplain_1d", false)::Bool

    @info "Local inertial approach is used for river flow." cfl h_thresh froude_limit floodplain_1d
    @warn string(
        "Providing the boundary condition `riverlength_bc` as part of the `[model]` setting ",
        "in the TOML file has been deprecated as of Wflow v0.8.0.\n The boundary condition should ",
        "be provided as part of the file `$(config.input.path_static)`.",
    )

    riverlength_bc = ncread(
        dataset,
        config,
        "routing.river_flow.riverlength_bc";
        sel = inds_pit,
        defaults = 1.0e04,
        type = Float64,
    )
    bankfull_elevation_2d = ncread(
        dataset,
        config,
        "routing.river_flow.bankfull_elevation";
        optional = false,
        type = Float64,
        fill = 0,
    )
    bankfull_depth_2d = ncread(
        dataset,
        config,
        "routing.river_flow.bankfull_depth";
        optional = false,
        type = Float64,
        fill = 0,
    )
    bankfull_depth = bankfull_depth_2d[indices]
    zb = bankfull_elevation_2d[indices] - bankfull_depth # river bed elevation

    bankfull_volume = bankfull_depth .* river_width .* river_length
    mannings_n = ncread(
        dataset,
        config,
        "routing.river_flow.mannings_n";
        sel = indices,
        defaults = 0.036,
        type = Float64,
    )

    n = length(indices)

    # set ghost points for boundary condition (downstream river outlet): river width, bed
    # elevation, manning n is copied from the upstream cell.
    append!(river_length, riverlength_bc)
    append!(zb, zb[index_pit])
    append!(river_width, river_width[index_pit])
    append!(mannings_n, mannings_n[index_pit])
    append!(bankfull_depth, bankfull_depth[index_pit])

    # determine z, width, length and manning's n at edges
    zb_max = fill(Float64(0), n_edges)
    width_at_edge = fill(Float64(0), n_edges)
    length_at_edge = fill(Float64(0), n_edges)
    mannings_n_sq = fill(Float64(0), n_edges)
    for i in 1:n_edges
        src_node = nodes_at_edge.src[i]
        dst_node = nodes_at_edge.dst[i]
        zb_max[i] = max(zb[src_node], zb[dst_node])
        width_at_edge[i] = min(river_width[src_node], river_width[dst_node])
        length_at_edge[i] = 0.5 * (river_length[dst_node] + river_length[src_node])
        mannings_n_i =
            (
                mannings_n[dst_node] * river_length[dst_node] +
                mannings_n[src_node] * river_length[src_node]
            ) / (river_length[dst_node] + river_length[src_node])
        mannings_n_sq[i] = mannings_n_i * mannings_n_i
    end
    active_index = findall(x -> x == 0, waterbody)

    parameters = LocalInertialRiverFlowParameters(;
        n,
        ne = n_edges,
        active_n = active_index,
        active_e = active_index,
        g = 9.80665,
        froude_limit,
        h_thresh,
        zb,
        zb_max,
        bankfull_volume,
        bankfull_depth,
        mannings_n,
        mannings_n_sq,
        flow_length = river_length,
        flow_length_at_edge = length_at_edge,
        flow_width = river_width,
        flow_width_at_edge = width_at_edge,
        waterbody,
    )
    return parameters
end

"Struct for storing local inertial river flow model variables"
@get_units @grid_loc @with_kw struct LocalInertialRiverFlowVariables{T}
    q::Vector{T} | "m3 s-1" | "edge"                    # river discharge (subgrid channel)
    q0::Vector{T} | "m3 s-1" | "edge"                   # river discharge (subgrid channel) at previous time step
    q_av::Vector{T} | "m3 s-1" | "edge"                 # average river channel (+ floodplain) discharge [m³ s⁻¹]
    q_channel_av::Vector{T} | "m3 s-1"                  # average river channel discharge [m³ s⁻¹]
    h::Vector{T} | "m"                                  # water depth
    zs_max::Vector{T} | "m" | "edge"                    # maximum water elevation at edge
    zs_src::Vector{T} | "m"                             # water elevation of source node of edge
    zs_dst::Vector{T} | "m"                             # water elevation of downstream node of edge
    hf::Vector{T} | "m" | "edge"                        # water depth at edge
    h_av::Vector{T} | "m"                               # average water depth
    a::Vector{T} | "m2" | "edge"                        # flow area at edge
    r::Vector{T} | "m" | "edge"                         # wetted perimeter at edge
    volume::Vector{T} | "m3"                            # river volume
    error::Vector{T} | "m3"                             # error volume    
end

"Initialize shallow water river flow model variables"
function LocalInertialRiverFlowVariables(dataset, config, indices, n_edges, inds_pit)
    floodplain_1d = get(config.model, "floodplain_1d", false)::Bool
    riverdepth_bc = ncread(
        dataset,
        config,
        "routing.river_flow.riverdepth_bc";
        sel = inds_pit,
        defaults = 0.0,
        type = Float64,
    )

    n = length(indices)
    # set river depth h to zero (including reservoir and lake locations)
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
        volume = zeros(n),
        error = zeros(n),
    )
    return variables
end

"Shallow water river flow model using the local inertial method"
@with_kw struct LocalInertialRiverFlow{T, R, L, F, A} <: AbstractRiverFlowModel
    timestepping::TimeStepping{T}
    boundary_conditions::RiverFlowBC{T, R, L}
    parameters::LocalInertialRiverFlowParameters{T}
    variables::LocalInertialRiverFlowVariables{T}
    floodplain::F                                       # Floodplain (1D) schematization
    allocation::A                                       # Water allocation
end

"Initialize shallow water river flow model `LocalIntertialRiverFlow`"
function LocalInertialRiverFlow(
    dataset,
    config,
    indices;
    graph_river,
    ldd_river,
    river_length,
    river_width,
    reservoir,
    lake,
    waterbody,
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
    cfl = get(config.model, "inertial_flow_alpha", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; cfl)

    index_pit = findall(x -> x == 5, ldd_river)
    inds_pit = indices[index_pit]

    add_vertex_edge_graph!(graph_river, index_pit)
    nodes_at_edge = adjacent_nodes_at_edge(graph_river)
    n_edges = ne(graph_river)

    parameters = LocalInertialRiverFlowParameters(
        dataset,
        config,
        indices;
        river_length,
        river_width,
        waterbody,
        n_edges,
        nodes_at_edge,
        index_pit,
        inds_pit,
    )
    variables = LocalInertialRiverFlowVariables(dataset, config, indices, n_edges, inds_pit)

    n = length(indices)
    boundary_conditions = RiverFlowBC(n, reservoir, lake)

    floodplain_1d = get(config.model, "floodplain_1d", false)::Bool
    if floodplain_1d
        zb_floodplain = parameters.zb .+ parameters.bankfull_depth
        floodplain = FloodPlain(
            dataset,
            config,
            indices;
            river_width,
            river_length,
            zb_floodplain,
            index_pit,
            n_edges,
            nodes_at_edge,
        )
    else
        floodplain = nothing
    end

    do_water_demand = haskey(config.model, "water_demand")
    sw_river = LocalInertialRiverFlow(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain,
        allocation = do_water_demand ? AllocationRiver(n) : NoAllocationRiver{Float64}(),
    )
    return sw_river, nodes_at_edge
end

"Return the upstream inflow for a waterbody in `LocalInertialRiverFlow`"
function get_inflow_waterbody(model::LocalInertialRiverFlow, src_edge)
    q_in = sum_at(model.variables.q, src_edge)
    if !isnothing(model.floodplain)
        q_in = q_in + sum_at(model.floodplain.variables.q, src_edge)
    end
    return q_in
end

# For local inertial river routing, `to_river` is included, as water body cells are excluded
# (boundary condition).
get_inflow_waterbody(::LocalInertialRiverFlow, model::KinWaveOverlandFlow) =
    model.variables.q_av .+ model.variables.to_river
get_inflow_waterbody(::LocalInertialRiverFlow, model::LateralSSF) =
    (model.variables.ssf .+ model.variables.to_river) ./ tosecond(BASETIMESTEP)

"Update local inertial river flow model `LocalIntertialRiverFlow` for a single timestep"
function local_inertial_river_update!(
    model::LocalInertialRiverFlow,
    network,
    dt,
    dt_forcing,
    doy,
    update_h,
)
    (; nodes_at_edge, edges_at_node) = network.river
    (; inwater, abstraction, inflow) = model.boundary_conditions
    river_v = model.variables
    river_p = model.parameters

    river_v.q0 .= river_v.q
    if !isnothing(model.floodplain)
        model.floodplain.variables.q0 .= model.floodplain.variables.q
    end
    @tturbo for j in eachindex(river_p.active_e)
        i = river_p.active_e[j]
        i_src = nodes_at_edge.src[i]
        i_dst = nodes_at_edge.dst[i]
        river_v.zs_src[i] = river_p.zb[i_src] + river_v.h[i_src]
        river_v.zs_dst[i] = river_p.zb[i_dst] + river_v.h[i_dst]

        river_v.zs_max[i] = max(river_v.zs_src[i], river_v.zs_dst[i])
        river_v.hf[i] = (river_v.zs_max[i] - river_p.zb_max[i])

        river_v.a[i] = river_p.flow_width_at_edge[i] * river_v.hf[i] # flow area (rectangular channel)
        river_v.r[i] = river_v.a[i] / (river_p.flow_width_at_edge[i] + 2.0 * river_v.hf[i]) # hydraulic radius (rectangular channel)

        river_v.q[i] = IfElse.ifelse(
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
        river_v.q[i] =
            IfElse.ifelse(river_v.h[i_src] <= 0.0, min(river_v.q[i], 0.0), river_v.q[i])
        river_v.q[i] =
            IfElse.ifelse(river_v.h[i_dst] <= 0.0, max(river_v.q[i], 0.0), river_v.q[i])

        river_v.q_av[i] += river_v.q[i] * dt
    end
    if !isnothing(model.floodplain)
        floodplain_p = model.floodplain.parameters
        floodplain_v = model.floodplain.variables

        @tturbo @. floodplain_v.hf = max(river_v.zs_max - floodplain_p.zb_max, 0.0)

        n = 0
        @inbounds for i in river_p.active_e
            @inbounds if river_v.hf[i] > river_p.h_thresh
                n += 1
                floodplain_v.hf_index[n] = i
            else
                floodplain_v.q[i] = 0.0
            end
        end

        @tturbo for j in 1:n
            i = floodplain_v.hf_index[j]
            i_src = nodes_at_edge.src[i]
            i_dst = nodes_at_edge.dst[i]

            i0 = 0
            for k in eachindex(floodplain_p.profile.depth)
                i0 += 1 * (floodplain_p.profile.depth[k] <= floodplain_v.hf[i])
            end
            i1 = max(i0, 1)
            i2 = IfElse.ifelse(i1 == length(floodplain_p.profile.depth), i1, i1 + 1)

            a_src = flow_area(
                floodplain_p.profile.width[i2, i_src],
                floodplain_p.profile.a[i1, i_src],
                floodplain_p.profile.depth[i1],
                floodplain_v.hf[i],
            )
            a_src = max(a_src - (floodplain_v.hf[i] * river_p.flow_width[i_src]), 0.0)

            a_dst = flow_area(
                floodplain_p.profile.width[i2, i_dst],
                floodplain_p.profile.a[i1, i_dst],
                floodplain_p.profile.depth[i1],
                floodplain_v.hf[i],
            )
            a_dst = max(a_dst - (floodplain_v.hf[i] * river_p.flow_width[i_dst]), 0.0)

            floodplain_v.a[i] = min(a_src, a_dst)

            floodplain_v.r[i] = IfElse.ifelse(
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

            floodplain_v.q[i] = IfElse.ifelse(
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
            floodplain_v.q[i] = IfElse.ifelse(
                floodplain_v.h[i_src] <= 0.0,
                min(floodplain_v.q[i], 0.0),
                floodplain_v.q[i],
            )
            floodplain_v.q[i] = IfElse.ifelse(
                floodplain_v.h[i_dst] <= 0.0,
                max(floodplain_v.q[i], 0.0),
                floodplain_v.q[i],
            )

            floodplain_v.q[i] = IfElse.ifelse(
                floodplain_v.q[i] * river_v.q[i] < 0.0,
                0.0,
                floodplain_v.q[i],
            )
            floodplain_v.q_av[i] += floodplain_v.q[i] * dt
        end
    end
    # For reservoir and lake locations the local inertial solution is replaced by the
    # reservoir or lake model. These locations are handled as boundary conditions in the
    # local inertial model (fixed h).
    (; reservoir, inflow_waterbody) = model.boundary_conditions
    inds_reservoir = network.reservoir.river_indices
    for v in eachindex(inds_reservoir)
        i = inds_reservoir[v]

        q_in = get_inflow_waterbody(model, edges_at_node.src[i])
        update!(reservoir, v, q_in + inflow_waterbody[i], dt, dt_forcing)
        river_v.q[i] = reservoir.variables.outflow[v]
        river_v.q_av[i] += river_v.q[i] * dt
    end
    (; lake, inflow_waterbody) = model.boundary_conditions
    inds_lake = network.lake.river_indices
    for v in eachindex(inds_lake)
        i = inds_lake[v]

        q_in = get_inflow_waterbody(model, edges_at_node.src[i])
        update!(lake, v, q_in + inflow_waterbody[i], doy, dt, dt_forcing)
        river_v.q[i] = max(lake.variables.outflow[v], 0.0)
        river_v.q_av[i] += river_v.q[i] * dt
    end
    if update_h
        @batch per = thread minbatch = 2000 for i in river_p.active_n
            q_src = sum_at(river_v.q, edges_at_node.src[i])
            q_dst = sum_at(river_v.q, edges_at_node.dst[i])
            river_v.volume[i] =
                river_v.volume[i] + (q_src - q_dst + inwater[i] - abstraction[i]) * dt

            if river_v.volume[i] < 0.0
                river_v.error[i] = river_v.error[i] + abs(river_v.volume[i])
                river_v.volume[i] = 0.0 # set volume to zero
            end
            river_v.volume[i] = max(river_v.volume[i] + inflow[i] * dt, 0.0) # add external inflow

            if !isnothing(model.floodplain)
                floodplain_v = model.floodplain.variables
                floodplain_p = model.floodplain.parameters
                q_src = sum_at(floodplain_v.q, edges_at_node.src[i])
                q_dst = sum_at(floodplain_v.q, edges_at_node.dst[i])
                floodplain_v.volume[i] = floodplain_v.volume[i] + (q_src - q_dst) * dt
                # TODO check following approach:
                # if floodplain volume negative, extract from river volume first
                if floodplain_v.volume[i] < 0.0
                    floodplain_v.error[i] =
                        floodplain_v.error[i] + abs(floodplain_v.volume[i])
                    floodplain_v.volume[i] = 0.0
                end
                volume_total = river_v.volume[i] + floodplain_v.volume[i]
                if volume_total > river_p.bankfull_volume[i]
                    flood_volume = volume_total - river_p.bankfull_volume[i]
                    h = flood_depth(
                        floodplain_p.profile,
                        flood_volume,
                        river_p.flow_length[i],
                        i,
                    )
                    river_v.h[i] = river_p.bankfull_depth[i] + h
                    river_v.volume[i] =
                        river_v.h[i] * river_p.flow_width[i] * river_p.flow_length[i]
                    floodplain_v.volume[i] = max(volume_total - river_v.volume[i], 0.0)
                    floodplain_v.h[i] = floodplain_v.volume[i] > 0.0 ? h : 0.0
                else
                    river_v.h[i] =
                        volume_total / (river_p.flow_length[i] * river_p.flow_width[i])
                    river_v.volume[i] = volume_total
                    floodplain_v.h[i] = 0.0
                    floodplain_v.volume[i] = 0.0
                end
                floodplain_v.h_av[i] += floodplain_v.h[i] * dt
            else
                river_v.h[i] =
                    river_v.volume[i] / (river_p.flow_length[i] * river_p.flow_width[i])
            end
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
    model::LocalInertialRiverFlow{T},
    network,
    doy,
    dt;
    update_h = true,
) where {T}
    (; reservoir, lake) = model.boundary_conditions

    set_waterbody_vars!(reservoir)
    set_waterbody_vars!(lake)

    if !isnothing(model.floodplain)
        model.floodplain.variables.q_av .= 0.0
        model.floodplain.variables.h_av .= 0.0
    end
    model.variables.q_av .= 0.0
    model.variables.h_av .= 0.0

    t = T(0.0)
    while t < dt
        dt_s = stable_timestep(model)
        if t + dt_s > dt
            dt_s = dt - t
        end
        local_inertial_river_update!(model, network, dt_s, dt, doy, update_h)
        t = t + dt_s
    end
    model.variables.q_av ./= dt
    model.variables.h_av ./= dt

    average_waterbody_vars!(reservoir, dt)
    average_waterbody_vars!(lake, dt)

    if !isnothing(model.floodplain)
        model.floodplain.variables.q_av ./= dt
        model.floodplain.variables.h_av ./= dt
        model.variables.q_channel_av .= model.variables.q_av
        model.variables.q_av .=
            model.variables.q_channel_av .+ model.floodplain.variables.q_av
    end

    return nothing
end

"Struct to store local inertial overland flow model variables"
@get_units @grid_loc @with_kw struct LocalInertialOverlandFlowVariables{T}
    qy0::Vector{T} | "m3 s-1" | "edge"      # flow in y direction at previous time step
    qx0::Vector{T} | "m3 s-1" | "edge"      # flow in x direction at previous time step
    qx::Vector{T} | "m3 s-1" | "edge"       # flow in x direction
    qy::Vector{T} | "m3 s-1" | "edge"       # flow in y direction
    volume::Vector{T} | "m3"                # total volume of cell (including river volume for river cells)
    error::Vector{T} | "m3"                 # error volume
    h::Vector{T} | "m"                      # water depth of cell (for river cells the reference is the river bed elevation `zb`)
    h_av::Vector{T} | "m"                   # average water depth (for river cells the reference is the river bed elevation `zb`)
end

"Initialize local inertial overland flow model variables"
function LocalInertialOverlandFlowVariables(n)
    variables = LocalInertialOverlandFlowVariables(;
        qx0 = zeros(n + 1),
        qy0 = zeros(n + 1),
        qx = zeros(n + 1),
        qy = zeros(n + 1),
        volume = zeros(n),
        error = zeros(n),
        h = zeros(n),
        h_av = zeros(n),
    )
    return variables
end

"Struct to store local inertial overland flow model parameters"
@get_units @grid_loc @with_kw struct LocalInertialOverlandFlowParameters{T}
    n::Int                                              # number of cells [-]
    x_length::Vector{T} | "m"                           # cell length x direction [m]
    y_length::Vector{T} | "m"                           # cell length y direction [m]
    xwidth::Vector{T} | "m" | "edge"                    # effective flow width x direction (floodplain) [m]
    ywidth::Vector{T} | "m" | "edge"                    # effective flow width y direction (floodplain) [m]
    g::T                                                # acceleration due to gravity [m s⁻²]
    theta::T                                            # weighting factor (de Almeida et al., 2012) [-]
    h_thresh::T                                         # depth threshold for calculating flow [m]
    zx_max::Vector{T} | "m" | "edge"                    # maximum cell elevation (x direction)
    zy_max::Vector{T} | "m" | "edge"                    # maximum cell elevation (y direction)
    mannings_n_sq::Vector{T} | "(s m-1/3)2" | "edge"    # Manning's roughness squared
    z::Vector{T} | "m"                                  # elevation of cell
    froude_limit::Bool                                  # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
    rivercells::Vector{Bool} | "-"                      # river cells
end

"Initialize shallow water overland flow model parameters"
function LocalInertialOverlandFlowParameters(
    dataset,
    config,
    indices;
    modelsize_2d,
    reverse_indices, # maps from the 2D external domain to the 1D internal domain (Int for linear indexing).
    x_length,
    y_length,
    river_width,
    graph_river,
    ldd_river,
    inds_river,
    river_location,
    waterbody,
)
    froude_limit = get(config.model, "froude_limit", true)::Bool # limit flow to subcritical according to Froude number
    cfl = get(config.model, "inertial_flow_alpha", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    theta = get(config.model, "inertial_flow_theta", 0.8)::Float64 # weighting factor
    h_thresh = get(config.model, "h_thresh", 1.0e-03)::Float64 # depth threshold for flow at edge

    @info "Local inertial approach is used for overlandflow." cfl theta h_thresh froude_limit

    mannings_n = ncread(
        dataset,
        config,
        "routing.overland_flow.mannings_n";
        sel = indices,
        defaults = 0.072,
        type = Float64,
    )
    elevation_2d = ncread(
        dataset,
        config,
        "routing.overland_flow.elevation";
        optional = false,
        type = Float64,
        fill = 0,
    )
    elevation = elevation_2d[indices]
    n = length(indices)

    # initialize edge connectivity of 2D staggered grid
    edge_indices =
        EdgeConnectivity(; xu = zeros(n), xd = zeros(n), yu = zeros(n), yd = zeros(n))

    nrow, ncol = modelsize_2d
    for (v, i) in enumerate(indices)
        for (m, neighbor) in enumerate(NEIGHBORS)
            j = i + neighbor
            dir = DIRS[m]
            if (1 <= j[1] <= nrow) && (1 <= j[2] <= ncol) && (reverse_indices[j] != 0)
                getfield(edge_indices, dir)[v] = reverse_indices[j]
            else
                getfield(edge_indices, dir)[v] = n + 1
            end
        end
    end

    # determine z at edges in x and y direction
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
    # for waterbody cells (reservoir or lake), h is set to zero (fixed) and not updated, and
    # overland flow from a downstream cell is not possible (effective flowwidth is zero).
    we_x = copy(x_length)
    we_y = copy(y_length)
    set_effective_flowwidth!(
        we_x,
        we_y,
        edge_indices,
        graph_river,
        river_width,
        ldd_river,
        waterbody,
        reverse_indices[inds_river],
    )
    parameters = LocalInertialOverlandFlowParameters(;
        n,
        x_length,
        y_length,
        xwidth = we_x,
        ywidth = we_y,
        g = 9.80665,
        theta,
        h_thresh,
        zx_max,
        zy_max,
        mannings_n_sq = mannings_n .* mannings_n,
        z = elevation,
        froude_limit,
        rivercells = river_location,
    )
    return parameters, edge_indices
end

"Struct to store local inertial overland flow model boundary conditions"
@get_units @grid_loc @with_kw struct LocalInertialOverlandFlowBC{T}
    runoff::Vector{T} | "m3 s-1"               # runoff from hydrological model
    inflow_waterbody::Vector{T} | "m3 s-1"     # inflow to water body from hydrological model
end

"Struct to store shallow water overland flow model boundary conditions"
function LocalInertialOverlandFlowBC(n)
    bc = LocalInertialOverlandFlowBC(; runoff = zeros(n), inflow_waterbody = zeros(n))
    return bc
end

"Local inertial overland flow model using the local inertial method"
@with_kw struct LocalInertialOverlandFlow{T} <: AbstractOverlandFlowModel
    timestepping::TimeStepping{T}
    boundary_conditions::LocalInertialOverlandFlowBC{T}
    parameters::LocalInertialOverlandFlowParameters{T}
    variables::LocalInertialOverlandFlowVariables{T}
end

"Initialize local inertial overland flow model"
function LocalInertialOverlandFlow(
    dataset,
    config,
    indices;
    modelsize_2d,
    reverse_indices, # maps from the 2D external domain to the 1D internal domain (Int for linear indexing).
    x_length,
    y_length,
    river_width,
    graph_river,
    ldd_river,
    inds_river,
    river_location,
    waterbody,
)
    cfl = get(config.model, "inertial_flow_alpha", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; cfl)

    n = length(indices)
    boundary_conditions = LocalInertialOverlandFlowBC(n)
    parameters, edge_indices = LocalInertialOverlandFlowParameters(
        dataset,
        config,
        indices;
        modelsize_2d,
        reverse_indices, # maps from the 2D external domain to the 1D internal domain (Int for linear indexing).
        x_length,
        y_length,
        river_width,
        graph_river,
        ldd_river,
        inds_river,
        river_location,
        waterbody,
    )
    variables = LocalInertialOverlandFlowVariables(n)

    sw_land = LocalInertialOverlandFlow{Float64}(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
    )

    return sw_land, edge_indices
end

"""
    stable_timestep(model::LocalInertialRiverFlow)
    stable_timestep(model::LocalInertialOverlandFlow)

Compute a stable timestep size for the local inertial approach, based on Bates et al. (2010).

dt = cfl * (Δx / sqrt(g max(h))
"""
function stable_timestep(model::LocalInertialRiverFlow{T})::T where {T}
    dt_min = T(Inf)
    (; cfl) = model.timestepping
    (; n, flow_length, g) = model.parameters
    (; h) = model.variables
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(n)
        @fastmath @inbounds dt = cfl * flow_length[i] / sqrt(g * h[i])
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? T(60.0) : dt_min
    return dt_min
end

function stable_timestep(model::LocalInertialOverlandFlow{T})::T where {T}
    dt_min = T(Inf)
    (; cfl) = model.timestepping
    (; n, g, x_length, y_length, rivercells) = model.parameters
    (; h) = model.variables
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(n)
        @fastmath @inbounds dt = if rivercells[i] == 0
            cfl * min(x_length[i], y_length[i]) / sqrt(g * h[i])
        else
            T(Inf)
        end
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? T(60.0) : dt_min
    return dt_min
end

"""
Update boundary conditions `runoff` and inflow to a waterbody from land `inflow_waterbody` for
overland flow model `LocalInertialOverlandFlow` for a single timestep.
"""
function update_boundary_conditions!(
    model::LocalInertialOverlandFlow,
    external_models::NamedTuple,
    network,
    dt,
)
    (; river_flow, soil, subsurface_flow, runoff) = external_models
    (; inflow_waterbody) = model.boundary_conditions
    (; reservoir, lake) = river_flow.boundary_conditions
    (; net_runoff) = soil.variables
    (; net_runoff_river) = runoff.variables

    model.boundary_conditions.runoff .=
        net_runoff ./ 1000.0 .* network.land.area ./ dt .+
        get_flux_to_river(subsurface_flow) .+
        net_runoff_river .* network.land.area .* 0.001 ./ dt

    if !isnothing(reservoir) || !isnothing(lake)
        inflow_land = get_inflow_waterbody(river_flow, model)
        inflow_subsurface = get_inflow_waterbody(river_flow, subsurface_flow)

        @. inflow_waterbody[network.river_indices] =
            inflow_land[network.river_indices] + inflow_subsurface[network.river_indices]
    end
    return nothing
end

"""
Update combined river `LocalInertialRiverFlow` and overland flow `LocalInertialOverlandFlow` models for a
single timestep `dt`. An adaptive timestepping method is used (computing a sub timestep
`dt_s`).
"""
function update!(
    land::LocalInertialOverlandFlow{T},
    river::LocalInertialRiverFlow{T},
    network,
    doy,
    dt;
    update_h = false,
) where {T}
    (; reservoir, lake) = river.boundary_conditions

    if !isnothing(reservoir)
        reservoir.boundary_conditions.inflow .= 0.0
        reservoir.variables.totaloutflow .= 0.0
        reservoir.variables.actevap .= 0.0
    end
    if !isnothing(lake)
        lake.boundary_conditions.inflow .= 0.0
        lake.variables.totaloutflow .= 0.0
        lake.variables.actevap .= 0.0
    end
    river.variables.q_av .= 0.0
    river.variables.h_av .= 0.0
    land.variables.h_av .= 0.0

    t = T(0.0)
    while t < dt
        dt_river = stable_timestep(river)
        dt_land = stable_timestep(land)
        dt_s = min(dt_river, dt_land)
        if t + dt_s > dt
            dt_s = dt - t
        end
        local_inertial_river_update!(river, network, dt_s, dt, doy, update_h)
        local_inertial_update!(land, river, network, dt_s)
        t = t + dt_s
    end
    river.variables.q_av ./= dt
    river.variables.h_av ./= dt
    land.variables.h_av ./= dt

    return nothing
end

"""
Update combined river `LocalInertialRiverFlow`and overland flow `LocalInertialOverlandFlow` models for a
single timestep `dt`.
"""
function local_inertial_update!(
    land::LocalInertialOverlandFlow{T},
    river::LocalInertialRiverFlow{T},
    network,
    dt,
) where {T}
    indices = network.land.edge_indices
    inds_river = network.land.river_indices

    (; edges_at_node) = network.river

    river_bc = river.boundary_conditions
    river_v = river.variables
    river_p = river.parameters
    land_bc = land.boundary_conditions
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
        if xu <= land_p.n && land_p.ywidth[i] != T(0.0)
            zs_x = land_p.z[i] + land_v.h[i]
            zs_xu = land_p.z[xu] + land_v.h[xu]
            zs_max = max(zs_x, zs_xu)
            hf = (zs_max - land_p.zx_max[i])

            if hf > land_p.h_thresh
                length = T(0.5) * (land_p.x_length[i] + land_p.x_length[xu]) # can be precalculated
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
                if land_v.h[i] <= T(0.0)
                    land_v.qx[i] = min(land_v.qx[i], T(0.0))
                end
                if land_v.h[xu] <= T(0.0)
                    land_v.qx[i] = max(land_v.qx[i], T(0.0))
                end
            else
                land_v.qx[i] = T(0.0)
            end
        end

        # update qy

        # the effective flow width is zero when the river width exceeds the cell width (dx
        # for flow in y dir) and floodplain flow is not calculated.
        if yu <= land_p.n && land_p.xwidth[i] != T(0.0)
            zs_y = land_p.z[i] + land_v.h[i]
            zs_yu = land_p.z[yu] + land_v.h[yu]
            zs_max = max(zs_y, zs_yu)
            hf = (zs_max - land_p.zy_max[i])

            if hf > land_p.h_thresh
                length = T(0.5) * (land_p.y_length[i] + land_p.y_length[yu]) # can be precalculated
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
                if land_v.h[i] <= T(0.0)
                    land_v.qy[i] = min(land_v.qy[i], T(0.0))
                end
                if land_v.h[yu] <= T(0.0)
                    land_v.qy[i] = max(land_v.qy[i], T(0.0))
                end
            else
                land_v.qy[i] = T(0.0)
            end
        end
    end

    # change in volume and water levels based on horizontal fluxes for river and land cells
    @batch per = thread minbatch = 6000 for i in 1:(land_p.n)
        yd = indices.yd[i]
        xd = indices.xd[i]

        if land_p.rivercells[i]
            if river_p.waterbody[inds_river[i]]
                # for reservoir or lake set inflow from land part, these are boundary points
                # and update of volume and h is not required
                river_bc.inflow_waterbody[inds_river[i]] =
                    land_bc.inflow_waterbody[i] +
                    land_bc.runoff[i] +
                    (land_v.qx[xd] - land_v.qx[i] + land_v.qy[yd] - land_v.qy[i])
            else
                land_v.volume[i] +=
                    (
                        sum_at(river_v.q, edges_at_node.src[inds_river[i]]) -
                        sum_at(river_v.q, edges_at_node.dst[inds_river[i]]) +
                        land_v.qx[xd] - land_v.qx[i] + land_v.qy[yd] - land_v.qy[i] +
                        river_bc.inflow[inds_river[i]] +
                        land_bc.runoff[i] - river_bc.abstraction[inds_river[i]]
                    ) * dt
                if land_v.volume[i] < T(0.0)
                    land_v.error[i] = land_v.error[i] + abs(land_v.volume[i])
                    land_v.volume[i] = T(0.0) # set volume to zero
                end
                if land_v.volume[i] >= river_p.bankfull_volume[inds_river[i]]
                    river_v.h[inds_river[i]] =
                        river_p.bankfull_depth[inds_river[i]] +
                        (land_v.volume[i] - river_p.bankfull_volume[inds_river[i]]) /
                        (land_p.x_length[i] * land_p.y_length[i])
                    land_v.h[i] =
                        river_v.h[inds_river[i]] - river_p.bankfull_depth[inds_river[i]]
                    river_v.volume[inds_river[i]] =
                        river_v.h[inds_river[i]] *
                        river_p.flow_length[inds_river[i]] *
                        river_p.flow_width[inds_river[i]]
                else
                    river_v.h[inds_river[i]] =
                        land_v.volume[i] / (
                            river_p.flow_length[inds_river[i]] *
                            river_p.flow_width[inds_river[i]]
                        )
                    land_v.h[i] = T(0.0)
                    river_v.volume[inds_river[i]] = land_v.volume[i]
                end
                river_v.h_av[inds_river[i]] += river_v.h[inds_river[i]] * dt
            end
        else
            land_v.volume[i] +=
                (
                    land_v.qx[xd] - land_v.qx[i] + land_v.qy[yd] - land_v.qy[i] +
                    land_bc.runoff[i]
                ) * dt
            if land_v.volume[i] < T(0.0)
                land_v.error[i] = land_v.error[i] + abs(land_v.volume[i])
                land_v.volume[i] = T(0.0) # set volume to zero
            end
            land_v.h[i] = land_v.volume[i] / (land_p.x_length[i] * land_p.y_length[i])
        end
        land_v.h_av[i] += land_v.h[i] * dt
    end
    return nothing
end

"""
    FloodPlainProfile

Floodplain `volume` is a function of `depth` (flood depth intervals). Based on the
cumulative floodplain `volume` a floodplain profile as a function of `flood_depth` is
derived with floodplain area `a` (cumulative) and wetted perimeter radius `p` (cumulative).
"""
@get_units @grid_loc @with_kw struct FloodPlainProfile{T, N}
    depth::Vector{T} | "m"                     # Flood depth
    volume::Array{T, 2} | "m3"                 # Flood volume (cumulative)
    width::Array{T, 2} | "m"                   # Flood width
    a::Array{T, 2} | "m2"                      # Flow area (cumulative)
    p::Array{T, 2} | "m"                       # Wetted perimeter (cumulative)
end

"Initialize floodplain profile `FloodPlainProfile`"
function FloodPlainProfile(dataset, config, indices; river_width, river_length, index_pit)
    volume = ncread(
        dataset,
        config,
        "routing.river_flow.floodplain.volume";
        sel = indices,
        type = Float64,
        dimname = :flood_depth,
    )
    n = length(indices)

    # for convenience (interpolation) flood depth 0.0 m is added, with associated area (a),
    # volume, width (river width) and wetted perimeter (p).
    volume = vcat(fill(Float64(0), n)', volume)
    start_volume = volume
    flood_depths = Float64.(dataset["flood_depth"][:])
    pushfirst!(flood_depths, 0.0)
    n_depths = length(flood_depths)

    p = zeros(Float64, n_depths, n)
    a = zeros(Float64, n_depths, n)
    segment_volume = zeros(Float64, n_depths, n)
    width = zeros(Float64, n_depths, n)
    width[1, :] = river_width[1:n]

    # determine flow area (a), width and wetted perimeter (p) FloodPlain
    h = diff(flood_depths)
    incorrect_vol = 0
    riv_cells = 0
    error_vol = 0
    for i in 1:n
        riv_cell = 0
        diff_volume = diff(volume[:, i])

        for j in 1:(n_depths - 1)
            # assume rectangular shape of flood depth segment
            width[j + 1, i] = diff_volume[j] / (h[j] * river_length[i])
            # check provided flood volume (floodplain width should be constant or increasing
            # as a function of flood depth)
            if width[j + 1, i] < width[j, i]
                # raise warning only if difference is larger than rounding error of 0.01 m³
                if ((width[j, i] - width[j + 1, i]) * h[j] * river_length[i]) > 0.01
                    incorrect_vol += 1
                    riv_cell = 1
                    error_vol =
                        error_vol +
                        ((width[j, i] - width[j + 1, i]) * h[j] * river_length[i])
                end
                width[j + 1, i] = width[j, i]
            end
            a[j + 1, i] = width[j + 1, i] * h[j]
            p[j + 1, i] = (width[j + 1, i] - width[j, i]) + 2.0 * h[j]
            segment_volume[j + 1, i] = a[j + 1, i] * river_length[i]
            if j == 1
                # for interpolation wetted perimeter at flood depth 0.0 is required
                p[j, i] = p[j + 1, i] - 2.0 * h[j]
            end
        end

        p[2:end, i] = cumsum(p[2:end, i])
        a[:, i] = cumsum(a[:, i])
        volume[:, i] = cumsum(segment_volume[:, i])

        riv_cells += riv_cell
    end

    if incorrect_vol > 0
        perc_riv_cells = round(100.0 * (riv_cells / n); digits = 2)
        perc_error_vol = round(100.0 * (error_vol / sum(start_volume[end, :])); digits = 2)
        @warn string(
            "The provided volume of $incorrect_vol rectangular floodplain schematization",
            " segments for $riv_cells river cells ($perc_riv_cells % of total river cells)",
            " is not correct and has been increased with $perc_error_vol % of provided volume.",
        )
    end

    # set floodplain parameters for ghost points
    volume = hcat(volume, volume[:, index_pit])
    width = hcat(width, width[:, index_pit])
    a = hcat(a, a[:, index_pit])
    p = hcat(p, p[:, index_pit])

    # initialize floodplain profile parameters
    profile =
        FloodPlainProfile{Float64, n_depths}(; volume, width, depth = flood_depths, a, p)
    return profile
end

"Struct to store floodplain flow model parameters"
@get_units @grid_loc @with_kw struct FloodPlainParameters{T, P}
    profile::P                                          # floodplain profile
    mannings_n::Vector{T} | "s m-1/3"                   # manning's roughness
    mannings_n_sq::Vector{T} | "(s m-1/3)2" | "edge"    # manning's roughness squared
    zb_max::Vector{T} | "m" | "edge"                    # maximum bankfull elevation (edge)
end

"Initialize floodplain flow model parameters"
function FloodPlainParameters(
    dataset,
    config,
    indices;
    river_width,
    river_length,
    zb_floodplain,
    nodes_at_edge,
    n_edges,
    index_pit,
)
    profile =
        FloodPlainProfile(dataset, config, indices; river_width, river_length, index_pit)

    mannings_n = ncread(
        dataset,
        config,
        "routing.river_flow.floodplain.mannings_n";
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
                mannings_n[dst_node] * river_length[dst_node] +
                mannings_n[src_node] * river_length[src_node]
            ) / (river_length[dst_node] + river_length[src_node])
        mannings_n_sq[i] = mannings_n_i * mannings_n_i
        zb_max[i] = max(zb_floodplain[src_node], zb_floodplain[dst_node])
    end
    parameters = FloodPlainParameters(profile, mannings_n, mannings_n_sq, zb_max)
    return parameters
end

"Struct to store floodplain flow model variables"
@get_units @grid_loc @with_kw struct FloodPlainVariables{T}
    volume::Vector{T} | "m3"                            # volume
    h::Vector{T} | "m"                                  # water depth
    h_av::Vector{T} | "m"                               # average water depth
    error::Vector{T} | "m3"                             # error volume
    a::Vector{T} | "m2" | "edge"                        # flow area
    r::Vector{T} | "m" | "edge"                         # hydraulic radius
    hf::Vector{T} | "m" | "edge"                        # water depth at edge
    q0::Vector{T} | "m3 s-1" | "edge"                   # discharge at previous time step
    q::Vector{T} | "m3 s-1" | "edge"                    # discharge
    q_av::Vector{T} | "m" | "edge"                      # average river discharge
    hf_index::Vector{Int} | "-" | "edge"                # index with `hf` above depth threshold
end

"Initialize floodplain flow model variables"
function FloodPlainVariables(n, n_edges, index_pit)
    variables = FloodPlainVariables(;
        volume = zeros(n),
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
@with_kw struct FloodPlain{T, P}
    parameters::FloodPlainParameters{T, P}
    variables::FloodPlainVariables{T}
end

"Determine the initial floodplain volume"
function initialize_volume!(river, nriv::Int)
    (; flow_width, flow_length) = river.parameters
    (; floodplain) = river
    profile = floodplain.parameters
    river = for i in 1:nriv
        i1, i2 = interpolation_indices(floodplain.variables.h[i], profile.depth)
        a = flow_area(
            profile.width[i2, i],
            profile.a[i1, i],
            profile.depth[i1],
            floodplain.variables.h[i],
        )
        a = max(a - (flow_width[i] * floodplain.h[i]), 0.0)
        floodplain.variables.volume[i] = flow_length[i] * a
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

"Compute flood depth by interpolating flood volume `flood_volume` using flood depth intervals."
function flood_depth(
    profile::FloodPlainProfile{T},
    flood_volume,
    flow_length,
    i::Int,
)::T where {T}
    i1, i2 = interpolation_indices(flood_volume, @view profile.volume[:, i])
    ΔA = (flood_volume - profile.volume[i1, i]) / flow_length
    dh = ΔA / profile.width[i2, i]
    flood_depth = profile.depth[i1] + dh
    return flood_depth
end

"Initialize floodplain geometry and `FloodPlain` variables and parameters"
function FloodPlain(
    dataset,
    config,
    indices;
    river_width,
    river_length,
    zb_floodplain,
    index_pit,
    n_edges,
    nodes_at_edge,
)
    n = length(indices)
    parameters = FloodPlainParameters(
        dataset,
        config,
        indices;
        river_width,
        river_length,
        zb_floodplain,
        nodes_at_edge,
        n_edges,
        index_pit,
    )
    variables = FloodPlainVariables(n, n_edges, index_pit)

    floodplain = FloodPlain(; parameters, variables)
    return floodplain
end