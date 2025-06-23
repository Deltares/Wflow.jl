"Struct for storing local inertial river flow model parameters"
@with_kw struct LocalInertialRiverFlowParameters{
    T <: DenseArray{Float},
    I <: DenseArray{Int},
}
    n::Int                                  # number of cells [-]
    ne::Int                                 # number of edges [-]
    active_n::I                             # active nodes [-]
    active_e::I                             # active edges [-]
    g::Float                                # acceleration due to gravity [m s⁻²]
    froude_limit::Bool                      # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
    h_thresh::Float                         # depth threshold for calculating flow [m]
    zb::T                                   # river bed elevation [m]
    zb_max::T                               # maximum channel bed elevation [m]
    bankfull_storage::T                     # bankfull storage [m³]
    bankfull_depth::T                       # bankfull depth [m]
    mannings_n_sq::T                        # Manning's roughness squared at edge [(s m-1/3)2]
    mannings_n::T                           # Manning's roughness [s m-1/3] at node
    flow_length_at_edge::T                  # flow (river) length at edge [m]
    flow_width_at_edge::T                   # flow (river) width at edge [m]
end

function Adapt.adapt_structure(to, from::LocalInertialRiverFlowParameters)
    return LocalInertialRiverFlowParameters(
        adapt(to, from.n),
        adapt(to, from.ne),
        adapt(to, from.active_n),
        adapt(to, from.active_e),
        adapt(to, from.g),
        adapt(to, from.froude_limit),
        adapt(to, from.h_thresh),
        adapt(to, from.zb),
        adapt(to, from.zb_max),
        adapt(to, from.bankfull_storage),
        adapt(to, from.bankfull_depth),
        adapt(to, from.mannings_n_sq),
        adapt(to, from.mannings_n),
        adapt(to, from.flow_length_at_edge),
        adapt(to, from.flow_width_at_edge),
    )
end

"Initialize local inertial river flow model parameters"
function LocalInertialRiverFlowParameters(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
)
    alpha =
        get(config.model, "river_local_inertial_flow__alpha_coefficient", Float(0.7))::Float # stability coefficient for model time step (0.2-0.7)
    waterdepth_threshold =
        get(config.model, "river_water_flow_threshold__depth", Float(1.0e-03))::Float # depth threshold for flow at edge
    froude_limit = get(config.model, "river_water_flow__froude_limit_flag", true)::Bool # limit flow to subcritical according to Froude number
    floodplain_1d = get(config.model, "floodplain_1d__flag", false)::Bool

    @info "Local inertial approach is used for river flow." alpha waterdepth_threshold froude_limit floodplain_1d

    (; pit_indices, indices, graph, local_drain_direction, nodes_at_edge) = domain.network
    (; flow_width, flow_length, waterbody_outlet) = domain.parameters

    lens = lens_input_parameter(config, "model_boundary_condition~river__length")
    riverlength_bc =
        ncread(dataset, config, lens; sel = pit_indices, defaults = 1.0e04, type = Float)
    lens = lens_input_parameter(config, "river_bank_water__elevation"; optional = false)
    bankfull_elevation_2d = ncread(dataset, config, lens; type = Float, fill = 0)
    lens = lens_input_parameter(config, "river_bank_water__depth"; optional = false)
    bankfull_depth_2d = ncread(dataset, config, lens; type = Float, fill = 0)
    bankfull_depth = bankfull_depth_2d[indices]
    zb = bankfull_elevation_2d[indices] - bankfull_depth # river bed elevation

    bankfull_storage = bankfull_depth .* flow_width .* flow_length
    lens = lens_input_parameter(config, "river_water_flow__manning_n_parameter")
    mannings_n =
        ncread(dataset, config, lens; sel = indices, defaults = 0.036, type = Float)

    n = Int(length(indices))
    index_pit = findall(x -> x == 5, local_drain_direction)
    # set ghost points for boundary condition (downstream river outlet): river width, bed
    # elevation, manning n is copied from the upstream cell.
    append!(flow_length, riverlength_bc)
    append!(zb, zb[index_pit])
    append!(flow_width, flow_width[index_pit])
    append!(mannings_n, mannings_n[index_pit])
    append!(bankfull_depth, bankfull_depth[index_pit])

    # determine z, width, length and manning's n at edges
    n_edges = Int(ne(graph))
    zb_max = fill(Float(0), n_edges)
    width_at_edge = fill(Float(0), n_edges)
    length_at_edge = fill(Float(0), n_edges)
    mannings_n_sq = fill(Float(0), n_edges)
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
    active_index = Vector{Int}(findall(x -> x == 0, waterbody_outlet))

    parameters = LocalInertialRiverFlowParameters(;
        n,
        ne = n_edges,
        active_n = active_index,
        active_e = active_index,
        g = Float(9.80665),
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
@with_kw struct LocalInertialRiverFlowVariables{T <: DenseArray{Float}}
    q::T                                      # river discharge at edge (subgrid channel) [m³ s⁻¹]
    q0::T                                     # river discharge at edge (subgrid channel) at previous time step [m³ s⁻¹]
    q_av::T                                   # average river channel (+ floodplain) discharge at edge [m³ s⁻¹] (model timestep Δt)
    q_channel_av::T                           # average river channel discharge at edge [m³ s⁻¹] (for model timestep Δt)
    h::T                                      # water depth [m]
    zs_max::T                                 # maximum water elevation at edge [m]
    zs_src::T                                 # water elevation of source node of edge [m]
    zs_dst::T                                 # water elevation of downstream node of edge [m]
    hf::T                                     # water depth at edge [m]
    h_av::T                                   # average water depth for model timestep Δt [m]
    a::T                                      # flow area at edge [m²]
    r::T                                      # wetted perimeter at edge [m]
    storage::T                                # river storage [m³]
    storage_av::T                             # average river storage for model timestep Δt [m³]
    error::T                                  # error storage [m³]
end

function Adapt.adapt_structure(to, from::LocalInertialRiverFlowVariables)
    return LocalInertialRiverFlowVariables(
        adapt(to, from.q),
        adapt(to, from.q0),
        adapt(to, from.q_av),
        adapt(to, from.q_channel_av),
        adapt(to, from.h),
        adapt(to, from.zs_max),
        adapt(to, from.zs_src),
        adapt(to, from.zs_dst),
        adapt(to, from.hf),
        adapt(to, from.h_av),
        adapt(to, from.a),
        adapt(to, from.r),
        adapt(to, from.storage),
        adapt(to, from.storage_av),
        adapt(to, from.error),
    )
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
        ncread(dataset, config, lens; sel = pit_indices, defaults = 0.0, type = Float)

    n = length(indices)
    n_edges = ne(graph)
    # set river depth h to zero (including reservoir and lake locations)
    h = zeros(Float, n)
    q_av = zeros(Float, n_edges)
    # set ghost points for boundary condition (downstream river outlet): river depth `h`
    append!(h, riverdepth_bc)
    variables = LocalInertialRiverFlowVariables(;
        q = zeros(Float, n_edges),
        q0 = zeros(Float, n_edges),
        q_av = q_av,
        q_channel_av = floodplain_1d ? zeros(Float, n_edges) : q_av,
        h = h,
        zs_max = zeros(Float, n_edges),
        zs_src = zeros(Float, n_edges),
        zs_dst = zeros(Float, n_edges),
        hf = zeros(Float, n_edges),
        h_av = zeros(Float, n),
        a = zeros(Float, n_edges),
        r = zeros(Float, n_edges),
        storage = zeros(Float, n),
        storage_av = zeros(Float, n),
        error = zeros(Float, n),
    )
    return variables
end

"Shallow water river flow model using the local inertial method"
@with_kw struct LocalInertialRiverFlow{
    T <: DenseArray{Float},
    I <: DenseArray{Int},
    R,
    L,
    F,
    A,
} <: AbstractRiverFlowModel
    timestepping::TimeStepping{T}
    boundary_conditions::RiverFlowBC{T, R, L}
    parameters::LocalInertialRiverFlowParameters{T, I}
    variables::LocalInertialRiverFlowVariables{T}
    floodplain::F                                       # Floodplain (1D) schematization
    allocation::A                                       # Water allocation
end

function Adapt.adapt_structure(to, from::LocalInertialRiverFlow)
    return LocalInertialRiverFlow(
        adapt(to, from.timestepping),
        adapt(to, from.boundary_conditions),
        adapt(to, from.parameters),
        adapt(to, from.variables),
        adapt(to, from.floodplain),
        adapt(to, from.allocation),
    )
end

"Initialize shallow water river flow model `LocalIntertialRiverFlow`"
function LocalInertialRiverFlow(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    reservoir::Union{SimpleReservoir, Nothing},
    lake::Union{Lake, Nothing},
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
    cfl =
        get(config.model, "river_local_inertial_flow__alpha_coefficient", Float(0.7))::Float # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; cfl)

    parameters = LocalInertialRiverFlowParameters(dataset, config, domain)
    variables = LocalInertialRiverFlowVariables(dataset, config, domain.network)

    n = Int(length(domain.network.indices))
    boundary_conditions = RiverFlowBC(n, reservoir, lake)

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

"Return the upstream inflow for a waterbody in `LocalInertialRiverFlow`"
function get_inflow_waterbody(
    model::LocalInertialRiverFlow,
    src_edge::DenseArray{Int},
)::Float
    q_in = sum_at(model.variables.q, src_edge)
    if !isnothing(model.floodplain)
        q_in += sum_at(model.floodplain.variables.q, src_edge)
    end
    return q_in
end

# For local inertial river routing, `to_river` is included, as water body cells are excluded
# (boundary condition).
get_inflow_waterbody(::LocalInertialRiverFlow, model::KinWaveOverlandFlow) =
    model.variables.q_av .+ model.variables.to_river
get_inflow_waterbody(::LocalInertialRiverFlow, model::LateralSSF) =
    (model.variables.ssf .+ model.variables.to_river) ./ Float(tosecond(BASETIMESTEP))

"Update local inertial river flow model `LocalIntertialRiverFlow` for a single timestep"
function local_inertial_river_update!(
    model::LocalInertialRiverFlow,
    # domain::Domain,
    dt::Float,
    dt_forcing::Float,
    doy::Int,
    update_h::Bool,
    nodes_at_edge::NodesAtEdge,
    edges_at_node::EdgesAtNode,
    flow_length::DenseArray{Float},
    flow_width::DenseArray{Float},
    inds_lake::DenseArray{Int},
    inds_reservoir::DenseArray{Int},
)
    # (; nodes_at_edge, edges_at_node) = domain.river.network
    # (; flow_length, flow_width) = domain.river.parameters
    (; inwater, abstraction, inflow) = model.boundary_conditions

    river_v = model.variables
    river_p = model.parameters

    river_v.q0 .= river_v.q
    if !isnothing(model.floodplain)
        model.floodplain.variables.q0 .= model.floodplain.variables.q
    end
    @inbounds AK.foreachindex(
        river_p.active_e;
        scheduler = :polyester,  # Use Polyester.jl on cpu
        min_elems = 1000,  # Same arg as `minbatch = 1000` in original Wflow.jl code
    ) do j
        i = river_p.active_e[j]
        i_src = nodes_at_edge.src[i]
        i_dst = nodes_at_edge.dst[i]
        river_v.zs_src[i] = river_p.zb[i_src] + river_v.h[i_src]
        river_v.zs_dst[i] = river_p.zb[i_dst] + river_v.h[i_dst]

        river_v.zs_max[i] = max(river_v.zs_src[i], river_v.zs_dst[i])
        river_v.hf[i] = (river_v.zs_max[i] - river_p.zb_max[i])

        river_v.a[i] = river_p.flow_width_at_edge[i] * river_v.hf[i] # flow area (rectangular channel)
        river_v.r[i] =
            river_v.a[i] / (river_p.flow_width_at_edge[i] + Float(2.0) * river_v.hf[i]) # hydraulic radius (rectangular channel)

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
            Float(0.0),
        )

        # limit q in case water is not available
        river_v.q[i] = ifelse(
            river_v.h[i_src] <= Float(0.0),
            min(river_v.q[i], Float(0.0)),
            river_v.q[i],
        )
        river_v.q[i] = ifelse(
            river_v.h[i_dst] <= Float(0.0),
            max(river_v.q[i], Float(0.0)),
            river_v.q[i],
        )
        # average river discharge (here accumulated for model timestep Δt)
        river_v.q_av[i] += river_v.q[i] * dt
    end

    if !isnothing(model.floodplain)
        floodplain_p = model.floodplain.parameters
        floodplain_v = model.floodplain.variables

        @inbounds AK.foreachindex(
            floodplain_v.hf;
            scheduler = :polyester,
            min_elems = 1000,
        ) do i
            floodplain_v.hf[i] = max(river_v.zs_max[i] - floodplain_p.zb_max[i], 0.0f0)
        end

        @inbounds AK.foreachindex(
            river_p.active_e;
            scheduler = :polyester,
            min_elems = 1000,
        ) do i
            if river_v.hf[i] <= river_p.h_thresh
                floodplain_v.q[i] = 0.0f0
            else
                i_src = nodes_at_edge.src[i]
                i_dst = nodes_at_edge.dst[i]

                i0 = Int(0)
                for k in eachindex(floodplain_p.profile.depth)
                    i0 += 1 * (floodplain_p.profile.depth[k] <= floodplain_v.hf[i])
                end
                i1 = max(i0, Int(1))
                i2 = ifelse(i1 == length(floodplain_p.profile.depth), i1, i1 + 1)

                a_src = flow_area(
                    floodplain_p.profile.width[i2, i_src],
                    floodplain_p.profile.a[i1, i_src],
                    floodplain_p.profile.depth[i1],
                    floodplain_v.hf[i],
                )
                a_src = max(a_src - (floodplain_v.hf[i] * flow_width[i_src]), 0.0f0)

                a_dst = flow_area(
                    floodplain_p.profile.width[i2, i_dst],
                    floodplain_p.profile.a[i1, i_dst],
                    floodplain_p.profile.depth[i1],
                    floodplain_v.hf[i],
                )
                a_dst = max(a_dst - (floodplain_v.hf[i] * flow_width[i_dst]), 0.0f0)

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
                    floodplain_v.a[i] > Float(1.0e-05),
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
                    0.0f0,
                )

                # limit floodplain q in case water is not available
                floodplain_v.q[i] = ifelse(
                    floodplain_v.h[i_src] <= 0.0f0,
                    min(floodplain_v.q[i], 0.0f0),
                    floodplain_v.q[i],
                )
                floodplain_v.q[i] = ifelse(
                    floodplain_v.h[i_dst] <= 0.0f0,
                    max(floodplain_v.q[i], 0.0f0),
                    floodplain_v.q[i],
                )

                floodplain_v.q[i] = ifelse(
                    floodplain_v.q[i] * river_v.q[i] < 0.0f0,
                    0.0f0,
                    floodplain_v.q[i],
                )
                # average floodplain discharge (here accumulated for model timestep Δt)
                floodplain_v.q_av[i] += floodplain_v.q[i] * dt
            end
        end
    end
    # For reservoir and lake locations the local inertial solution is replaced by the
    # reservoir or lake model. These locations are handled as boundary conditions in the
    # local inertial model (fixed h).
    (; reservoir, inflow_waterbody) = model.boundary_conditions
    # inds_reservoir = domain.reservoir.network.river_indices

    if !isnothing(reservoir)  # should be declared before GPU kernel, otherwise it can't compile.
        AK.foreachindex(
            inds_reservoir;
            scheduler = :polyester,  # Use Polyester.jl on cpu
            min_elems = 1000,  # Same arg as `minbatch = 1000` in original Wflow.jl code
        ) do v
            i = inds_reservoir[v]
            # q_in_reservoir = get_inflow_waterbody(model, edges_at_node.src[i, :])
            q_in_reservoir = 0.0f0 # get_inflow_waterbody depends on !isnothing(floodplain). Not GPU compat.
            update!(reservoir, Int(v), q_in_reservoir + inflow_waterbody[i], dt, dt_forcing)
            river_v.q[i] = reservoir.variables.outflow[v]
            # average river discharge (here accumulated for model timestep Δt)
            river_v.q_av[i] += river_v.q[i] * dt
        end
    end

    (; lake) = model.boundary_conditions
    # inds_lake = domain.lake.network.river_indices
    if !isnothing(lake)  # should be declared before GPU kernel, otherwise it can't compile.
        AK.foreachindex(
            inds_lake;
            scheduler = :polyester,  # Use Polyester.jl on cpu
            min_elems = 1000,  # Same arg as `minbatch = 1000` in original Wflow.jl code
        ) do v
            i = inds_lake[v]
            # q_in_lake = get_inflow_waterbody(model, edges_at_node.src[i, :])
            q_in_lake = 0.0f0 # get_inflow_waterbody depends on !isnothing(floodplain). Not GPU compat.
            update!(lake, v, q_in_lake + inflow_waterbody[i], doy, dt, dt_forcing)
            river_v.q[i] = max(lake.variables.outflow[v], 0.0f0)
            # average river discharge (here accumulated for model timestep Δt)
            river_v.q_av[i] += river_v.q[i] * dt
        end
    end

    if update_h
        AK.foreachindex(river_p.active_n; scheduler = :polyester, min_elems = 1000) do i
            q_src = 0.0f0
            for j in axes(edges_at_node.src, 2)
                n = edges_at_node.src[i, j]
                if n ≠ Int(0)
                    q_src += river_v.q[n]
                end
            end
            m = edges_at_node.dst[i]
            if m ≠ Int(0)
                q_dst = river_v.q[m]
            else
                q_dst = 0.0f0
            end

            # internal abstraction (water demand) is limited by river storage and negative
            # external inflow as part of water allocation computations.
            river_v.storage[i] =
                river_v.storage[i] + (q_src - q_dst + inwater[i] - abstraction[i]) * dt

            if river_v.storage[i] < Float(0.0)
                river_v.error[i] = river_v.error[i] + abs(river_v.storage[i])
                river_v.storage[i] = Float(0.0) # set storage to zero
            end
            # limit negative external inflow
            if inflow[i] < Float(0.0)
                _inflow = max(-(Float(0.80) * river_v.storage[i] / dt), inflow[i])
            else
                _inflow = inflow[i]
            end

            river_v.storage[i] += _inflow * dt # add external inflow
            river_v.h[i] = river_v.storage[i] / (flow_length[i] * flow_width[i])

            ## TODO: make compilation for GPU work if floodplain is not defined.
            # if !isnothing(model.floodplain)
            #     floodplain_v = model.floodplain.variables
            #     floodplain_p = model.floodplain.parameters
            #     q_src = 0.0f0
            #     for j in axes(edges_at_node.src, 2)
            #         n = edges_at_node.src[i, j]
            #         if n ≠ Int(0)
            #             q_src += floodplain_v.q[n]
            #         end
            #     end
            #     m = edges_at_node.dst[i]
            #     if m ≠ Int(0)
            #         q_dst = floodplain_v.q[m]
            #     else
            #         q_dst = 0.0f0
            #     end

            #     floodplain_v.storage[i] = floodplain_v.storage[i] + (q_src - q_dst) * dt
            #     if floodplain_v.storage[i] < Float(0.0)
            #         floodplain_v.error[i] =
            #             floodplain_v.error[i] + abs(floodplain_v.storage[i])
            #         floodplain_v.storage[i] = Float(0.0)
            #     end
            #     storage_total = river_v.storage[i] + floodplain_v.storage[i]
            #     if storage_total > river_p.bankfull_storage[i]
            #         flood_storage = storage_total - river_p.bankfull_storage[i]
            #         h = flood_depth(floodplain_p.profile, flood_storage, flow_length[i], i)
            #         river_v.h[i] = river_p.bankfull_depth[i] + h
            #         river_v.storage[i] = river_v.h[i] * flow_width[i] * flow_length[i]
            #         floodplain_v.storage[i] =
            #             max(storage_total - river_v.storage[i], Float(0.0))
            #         floodplain_v.h[i] =
            #             floodplain_v.storage[i] > Float(0.0) ? h : Float(0.0)
            #     else
            #         river_v.h[i] = storage_total / (flow_length[i] * flow_width[i])
            #         river_v.storage[i] = storage_total
            #         floodplain_v.h[i] = Float(0.0)
            #         floodplain_v.storage[i] = Float(0.0)
            #     end
            #     # average variables (here accumulated for model timestep Δt)
            #     floodplain_v.storage_av[i] += floodplain_v.storage[i] * dt
            #     floodplain_v.h_av[i] += floodplain_v.h_av[i] * dt
            # end

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
    doy::Int,
    dt::Float;
    update_h = true,
)
    (; reservoir, lake) = model.boundary_conditions
    (; flow_length) = domain.river.parameters

    flow_length = adapt(BackendArray, domain.river.parameters.flow_length)
    flow_width = adapt(BackendArray, domain.river.parameters.flow_width)
    nodes_at_edge = adapt(BackendArray, domain.river.network.nodes_at_edge)
    edges_at_node = adapt(BackendArray, domain.river.network.edges_at_node)
    inds_reservoir = adapt(BackendArray, domain.reservoir.network.river_indices)
    inds_lake = adapt(BackendArray, domain.lake.network.river_indices)

    set_waterbody_vars!(reservoir)
    set_waterbody_vars!(lake)

    if !isnothing(model.floodplain)
        set_flow_vars!(model.floodplain.variables)
    end
    set_flow_vars!(model.variables)

    t = Float(0.0)
    steps = 100
    river = adapt(BackendArray, model)  # adapt to GPU
    while t < dt
        dt_s = stable_timestep(river, flow_length)
        dt_s *= Float(0.9)  # safety margin
        if t + steps * dt_s > dt
            dt_s = (dt - t) / steps
        end
        t = t + steps * dt_s
        for i in 1:steps
            @inbounds local_inertial_river_update!(
                river,
                # domain,
                Float(dt_s),
                dt,
                doy,
                update_h,
                nodes_at_edge,
                edges_at_node,
                flow_length,
                flow_width,
                inds_lake,
                inds_reservoir,
            )
        end
    end

    river_cpu = adapt(Array, river)
    copy!(model.variables.q_av, river_cpu.variables.q_av)
    copy!(model.variables.q, river_cpu.variables.q)
    copy!(model.variables.q0, river_cpu.variables.q0)
    copy!(model.variables.q_channel_av, river_cpu.variables.q_channel_av)
    copy!(model.variables.h, river_cpu.variables.h)
    copy!(model.variables.zs_max, river_cpu.variables.zs_max)
    copy!(model.variables.zs_src, river_cpu.variables.zs_src)
    copy!(model.variables.zs_dst, river_cpu.variables.zs_dst)
    copy!(model.variables.hf, river_cpu.variables.hf)
    copy!(model.variables.h_av, river_cpu.variables.h_av)
    copy!(model.variables.a, river_cpu.variables.a)
    copy!(model.variables.r, river_cpu.variables.r)
    copy!(model.variables.storage, river_cpu.variables.storage)
    copy!(model.variables.storage_av, river_cpu.variables.storage_av)
    copy!(model.variables.error, river_cpu.variables.error)

    # first copy! doesn't work for q_av, so we try again
    while !all(model.variables.q_av .=== river_cpu.variables.q_av)
        copy!(model.variables.q_av, river_cpu.variables.q_av)
    end

    # Set vars to nothing so gc can clean up. Wasn't happening automatically...
    river = nothing
    river_cpu = nothing
    flow_length = nothing
    flow_width = nothing
    nodes_at_edge = nothing
    edges_at_node = nothing
    inds_reservoir = nothing
    inds_lake = nothing
    GC.gc()  # force garbage collection

    average_flow_vars!(model.variables, dt)
    average_waterbody_vars!(reservoir, dt)
    average_waterbody_vars!(lake, dt)

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
    qy0::Vector{Float}              # flow in y direction at edge at previous time step [m³ s⁻¹]
    qx0::Vector{Float}              # flow in x direction at edge at previous time step [m³ s⁻¹]
    qx::Vector{Float}               # flow in x direction at egde [m³ s⁻¹]
    qy::Vector{Float}               # flow in y direction at edge [m³ s⁻¹]
    storage::Vector{Float}          # total storage of cell [m³] (including river storage for river cells) 
    storage_av::Vector{Float}       # average total storage of cell [m³] (including river storage for river cells) (model timestep Δt)
    error::Vector{Float}            # error storage [m³]
    h::Vector{Float}                # water depth of cell [m] (for river cells the reference is the river bed elevation `zb`)
    h_av::Vector{Float}             # average water depth [m] (for river cells the reference is the river bed elevation `zb`) (model timestep Δt)
end

"Initialize local inertial overland flow model variables"
function LocalInertialOverlandFlowVariables(n::Int)
    variables = LocalInertialOverlandFlowVariables(;
        qx0 = zeros(Float, n + 1),
        qy0 = zeros(Float, n + 1),
        qx = zeros(Float, n + 1),
        qy = zeros(Float, n + 1),
        storage = zeros(Float, n),
        storage_av = zeros(Float, n),
        error = zeros(Float, n),
        h = zeros(Float, n),
        h_av = zeros(Float, n),
    )
    return variables
end

"Struct to store local inertial overland flow model parameters"
@with_kw struct LocalInertialOverlandFlowParameters
    n::Int                            # number of cells [-]
    xwidth::Vector{Float}             # effective flow width x direction at edge (floodplain) [m]
    ywidth::Vector{Float}             # effective flow width y direction at edge (floodplain) [m]
    g::Float                          # acceleration due to gravity [m s⁻²]
    theta::Float                      # weighting factor (de Almeida et al., 2012) [-]
    h_thresh::Float                   # depth threshold for calculating flow [m]
    zx_max::Vector{Float}             # maximum cell elevation at edge [m] (x direction)
    zy_max::Vector{Float}             # maximum cell elevation at edge [m] (y direction)
    mannings_n_sq::Vector{Float}      # Manning's roughness squared at edge [(s m-1/3)2]
    z::Vector{Float}                  # elevation [m] of cell
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
    alpha =
        get(config.model, "land_local_inertial_flow__alpha_coefficient", Float(0.7))::Float # stability coefficient for model time step (0.2-0.7)
    theta =
        get(config.model, "land_local_inertial_flow__theta_coefficient", Float(0.8))::Float # weighting factor
    waterdepth_threshold =
        get(config.model, "land_surface_water_flow_threshold__depth", Float(1.0e-03))::Float # depth threshold for flow at edge

    (; edge_indices, indices) = domain.land.network
    (; x_length, y_length) = domain.land.parameters

    @info "Local inertial approach is used for overland flow." alpha theta waterdepth_threshold froude_limit

    lens = lens_input_parameter(config, "land_surface_water_flow__manning_n_parameter")
    mannings_n =
        ncread(dataset, config, lens; sel = indices, defaults = 0.072, type = Float)
    lens = lens_input_parameter(
        config,
        "land_surface_water_flow__ground_elevation";
        optional = false,
    )
    elevation_2d = ncread(dataset, config, lens; type = Float, fill = 0)
    elevation = elevation_2d[indices]
    n = length(domain.land.network.indices)

    zx_max = fill(Float(0), n)
    zy_max = fill(Float(0), n)
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
    runoff::Vector{Float}           # runoff from hydrological model [m³ s⁻¹]
    inflow_waterbody::Vector{Float} # inflow to water body from hydrological model [m³ s⁻¹]
end

"Struct to store shallow water overland flow model boundary conditions"
function LocalInertialOverlandFlowBC(n::Int)
    bc = LocalInertialOverlandFlowBC(;
        runoff = zeros(Float, n),
        inflow_waterbody = zeros(Float, n),
    )
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
    cfl =
        get(config.model, "land_local_inertial_flow__alpha_coefficient", Float(0.7))::Float # stability coefficient for model time step (0.2-0.7)
    timestepping = TimeStepping(; cfl)

    n = Int(length(domain.land.network.indices))
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
    stable_timestep(model::LocalInertialRiverFlow, flow_length::Vector{Float})
    stable_timestep(model::LocalInertialOverlandFlow, parameters::LandParameters)

Compute a stable timestep size for the local inertial approach, based on Bates et al. (2010).

dt = cfl * (Δx / sqrt(g max(h))
"""
function stable_timestep(
    model::LocalInertialRiverFlow,
    flow_length::T,
) where {T <: DenseArray{Float}}
    (; cfl) = model.timestepping
    (; n, g) = model.parameters
    (; h) = model.variables
    dt = array_from_host(Array{Float}(undef, n))
    AK.foreachindex(dt) do i
        @inbounds dt[i] = cfl * flow_length[i] / sqrt(g * h[i])
    end
    dt_min = AK.minimum(dt; init = Float(Inf))
    dt_min = (isinf(dt_min) || isnan(dt_min)) ? Float(60.0) : dt_min
    dt_min = return dt_min
end

function stable_timestep(model::LocalInertialOverlandFlow, parameters::LandParameters)
    dt_min = Inf
    (; cfl) = model.timestepping
    (; n, g) = model.parameters
    (; x_length, y_length, river_location) = parameters
    (; h) = model.variables
    # dt = array_from_host(Array{Float}(undef, n))
    dt = Array{Float}(undef, n)  # temporary to force on cpu
    @inbounds AK.foreachindex(dt; scheduler = :polyester, min_elems = 1000) do i
        if river_location[i] == 0
            dt[i] = cfl * min(x_length[i], y_length[i]) / sqrt(g * h[i])
        else
            dt[i] = Float(Inf)
        end
    end
    dt_min = AK.minimum(dt; init = Float(Inf))
    dt_min = (isinf(dt_min) || isnan(dt_min)) ? Float(60.0) : dt_min
    return dt_min
end

"""
Update boundary conditions `runoff` and inflow to a waterbody from land `inflow_waterbody` for
overland flow model `LocalInertialOverlandFlow` for a single timestep.
"""
function update_boundary_conditions!(
    model::LocalInertialOverlandFlow,
    external_models::NamedTuple,
    domain::Domain,
    dt::Float,
)
    (; river_flow, soil, subsurface_flow, runoff) = external_models
    (; inflow_waterbody) = model.boundary_conditions
    (; reservoir, lake) = river_flow.boundary_conditions
    (; net_runoff) = soil.variables
    (; net_runoff_river) = runoff.variables

    (; area) = domain.land.parameters
    (; network) = domain.river

    model.boundary_conditions.runoff .=
        net_runoff ./ 1000.0 .* area ./ dt .+ get_flux_to_river(subsurface_flow) .+
        net_runoff_river .* area .* 0.001 ./ dt

    if !isnothing(reservoir) || !isnothing(lake)
        inflow_subsurface = get_inflow_waterbody(river_flow, subsurface_flow)

        @. inflow_waterbody[network.land_indices] = inflow_subsurface[network.land_indices]
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
function average_flow_vars!(variables::LocalInertialOverlandFlowVariables, dt::Float)
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
    doy::Int,
    dt::Float;
    update_h = false,
)
    (; reservoir, lake) = river.boundary_conditions
    (; parameters) = domain.land
    # (; flow_length) = domain.river.parameters
    # flow_length = adapt(BackendArray, domain.river.parameters.flow_length)
    # flow_width = adapt(BackendArray, domain.river.parameters.flow_width)
    # nodes_at_edge = adapt(BackendArray, domain.river.network.nodes_at_edge)
    # edges_at_node = adapt(BackendArray, domain.river.network.edges_at_node)
    # inds_reservoir = adapt(BackendArray, domain.reservoir.network.river_indices)
    # inds_lake = adapt(BackendArray, domain.lake.network.river_indices)
    flow_length = domain.river.parameters.flow_length
    flow_width = domain.river.parameters.flow_width
    nodes_at_edge = domain.river.network.nodes_at_edge
    edges_at_node = domain.river.network.edges_at_node
    inds_reservoir = domain.reservoir.network.river_indices
    inds_lake = domain.lake.network.river_indices

    set_waterbody_vars!(reservoir)
    set_waterbody_vars!(lake)
    set_flow_vars!(river.variables)
    set_flow_vars!(land.variables)

    t = Float(0.0)
    steps = 100

    while t < dt
        # dt_s = stable_timestep(river, flow_length)
        dt_river = stable_timestep(river, flow_length)
        dt_land = stable_timestep(land, parameters)
        dt_s = min(dt_river, dt_land)

        dt_s *= Float(0.9)  # safety margin
        if t + steps * dt_s > dt
            dt_s = (dt - t) / steps
        end
        t = t + steps * dt_s

        for i in 1:steps
            @inbounds local_inertial_river_update!(
                river,
                # domain,
                Float(dt_s),
                dt,
                doy,
                update_h,
                nodes_at_edge,
                edges_at_node,
                flow_length,
                flow_width,
                inds_lake,
                inds_reservoir,
            )
            @inbounds local_inertial_update!(land, river, domain, Float(dt_s))
        end
    end

    average_flow_vars!(river.variables, dt)
    average_flow_vars!(land.variables, dt)

    average_waterbody_vars!(reservoir, dt)
    average_waterbody_vars!(lake, dt)

    return nothing
end

"""
Update combined river `LocalInertialRiverFlow`and overland flow `LocalInertialOverlandFlow`
models for a single timestep `dt`.
"""
function local_inertial_update!(
    land::LocalInertialOverlandFlow,
    river::LocalInertialRiverFlow,
    domain::Domain,
    dt::Float,
)
    indices = domain.land.network.edge_indices
    inds_river = domain.land.network.river_indices

    (; flow_width, flow_length) = domain.river.parameters
    (; x_length, y_length) = domain.land.parameters

    (; edges_at_node) = domain.river.network

    river_bc = river.boundary_conditions
    river_v = river.variables
    river_p = river.parameters
    land_bc = land.boundary_conditions
    land_v = land.variables
    land_p = land.parameters

    land_v.qx0 .= land_v.qx
    land_v.qy0 .= land_v.qy

    # update qx
    @inbounds AK.foreachindex(land_v.h; scheduler = :polyester, min_elems = 6000) do i
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

    # change in storage and water levels based on horizontal fluxes for river and land cells
    @inbounds AK.foreachindex(land_v.h; scheduler = :polyester, min_elems = 6000) do i
        yd = indices.yd[i]
        xd = indices.xd[i]
        j = inds_river[i]
        if domain.land.parameters.river_location[i]
            if domain.river.parameters.waterbody_outlet[j]
                # for reservoir or lake set inflow from land part, these are boundary points
                # and update of storage and h is not required
                river_bc.inflow_waterbody[j] =
                    land_bc.inflow_waterbody[i] +
                    land_bc.runoff[i] +
                    (land_v.qx[xd] - land_v.qx[i] + land_v.qy[yd] - land_v.qy[i])
            else
                # internal abstraction (water demand) is limited by river storage and negative
                # external inflow as part of water allocation computations.
                land_v.storage[i] +=
                    (
                        sum_at(river_v.q, edges_at_node.src[j, :]) -
                        sum_at(river_v.q, edges_at_node.dst[j, :]) + land_v.qx[xd] -
                        land_v.qx[i] + land_v.qy[yd] - land_v.qy[i] + land_bc.runoff[i] -
                        river_bc.abstraction[j]
                    ) * dt
                if land_v.storage[i] < 0.0
                    land_v.error[i] = land_v.error[i] + abs(land_v.storage[i])
                    land_v.storage[i] = 0.0 # set storage to zero
                end
                # limit negative external inflow
                if river_bc.inflow[j] < 0.0
                    available_volume = if land_v.storage[i] >= river_p.bankfull_storage[j]
                        river_p.bankfull_depth[j]
                    else
                        river_v.storage[j]
                    end
                    _inflow = max(-(0.80 * available_volume / dt), river_bc.inflow[j])
                else
                    _inflow = river_bc.inflow[j]
                end
                land_v.storage[i] += _inflow * dt
                if land_v.storage[i] >= river_p.bankfull_storage[j]
                    river_v.h[j] =
                        river_p.bankfull_depth[j] +
                        (land_v.storage[i] - river_p.bankfull_storage[j]) /
                        (x_length[i] * y_length[i])
                    land_v.h[i] = river_v.h[j] - river_p.bankfull_depth[j]
                    river_v.storage[j] = river_v.h[j] * flow_length[j] * flow_width[j]
                else
                    river_v.h[j] = land_v.storage[i] / (flow_length[j] * flow_width[j])
                    land_v.h[i] = 0.0
                    river_v.storage[j] = land_v.storage[i]
                end
                # average variables (here accumulated for model timestep Δt)
                river_v.h_av[j] += river_v.h[j] * dt
                river_v.storage_av[j] += river_v.storage[j] * dt
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
@with_kw struct FloodPlainProfile{N, T <: DenseArray{Float}}
    depth::T                      # Flood depth [m]
    storage::T                    # Flood storage (cumulative) [m³]
    width::T                      # Flood width [m]
    a::T                          # Flow area (cumulative) [m²]
    p::T                          # Wetted perimeter (cumulative) [m]
end

function Adapt.adapt_structure(to, from::FloodPlainProfile)
    return FloodPlainProfile(
        adapt(to, from.depth),
        adapt(to, from.storage),
        adapt(to, from.width),
        adapt(to, from.a),
        adapt(to, from.p),
    )
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
        ncread(dataset, config, lens; sel = indices, type = Float, dimname = :flood_depth)
    n = length(indices)

    # for convenience (interpolation) flood depth 0.0 m is added, with associated area (a),
    # storage, width (river width) and wetted perimeter (p).
    storage = vcat(fill(Float(0), n)', storage)
    start_storage = storage
    flood_depths = Float.(dataset["flood_depth"][:])
    pushfirst!(flood_depths, 0.0)
    n_depths = length(flood_depths)

    p = zeros(Float, n_depths, n)
    a = zeros(Float, n_depths, n)
    segment_storage = zeros(Float, n_depths, n)
    width = zeros(Float, n_depths, n)
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
@with_kw struct FloodPlainParameters{P, T <: DenseArray{Float}}
    profile::P                          # floodplain profile
    mannings_n::T                       # manning's roughness [s m-1/3]
    mannings_n_sq::T                    # manning's roughness squared at edge [(s m-1/3)2]
    zb_max::T                           # maximum bankfull elevation at edge [m]
end

function Adapt.adapt_structure(to, from::FloodPlainParameters)
    return FloodPlainParameters(
        adapt(to, from.profile),
        adapt(to, from.mannings_n),
        adapt(to, from.mannings_n_sq),
        adapt(to, from.zb_max),
    )
end

"Initialize floodplain flow model parameters"
function FloodPlainParameters(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    zb_floodplain::Vector{Float},
    index_pit::Vector{Int},
)
    (; indices, nodes_at_edge, graph) = domain.network
    (; flow_length) = domain.parameters
    n_edges = ne(graph)
    profile = FloodPlainProfile(dataset, config, domain, index_pit)

    lens = lens_input_parameter(config, "floodplain_water_flow__manning_n_parameter")
    mannings_n =
        ncread(dataset, config, lens; sel = indices, defaults = 0.072, type = Float)
    # manning roughness at edges
    append!(mannings_n, mannings_n[index_pit]) # copy to ghost nodes
    mannings_n_sq = fill(Float(0), n_edges)
    zb_max = fill(Float(0), n_edges)
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
@with_kw struct FloodPlainVariables{T <: DenseArray{Float}, I <: DenseArray{Int}}
    storage::T                      # storage [m³]
    storage_av::T                   # average storage for model timestep Δt [m³]
    h::T                            # water depth [m]
    h_av::T                         # average water depth [m] for model timestep Δt
    error::T                        # error storage [m³]
    a::T                            # flow area at egde [m²]
    r::T                            # hydraulic radius at edge [m]
    hf::T                           # water depth at edge [m]
    q0::T                           # discharge at edge at previous time step
    q::T                            # discharge at edge  [m³ s⁻¹]
    q_av::T                         # average river discharge at edge  [m³ s⁻¹] for model timestep Δt
    hf_index::I                     # edge index with `hf` [-] above depth threshold
end

function Adapt.adapt_structure(to, from::FloodPlainVariables)
    return FloodPlainVariables(
        adapt(to, from.storage),
        adapt(to, from.storage_av),
        adapt(to, from.h),
        adapt(to, from.h_av),
        adapt(to, from.error),
        adapt(to, from.a),
        adapt(to, from.r),
        adapt(to, from.hf),
        adapt(to, from.q0),
        adapt(to, from.q),
        adapt(to, from.q_av),
        adapt(to, from.hf_index),
    )
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
@with_kw struct FloodPlain{P, T <: DenseArray{Float}, I <: DenseArray{Int}}
    parameters::FloodPlainParameters{P, T}
    variables::FloodPlainVariables{T, I}
end

"Determine the initial floodplain storage"
function initialize_storage!(river, nriv::Int)
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
    p = p + (2.0f0 * dh) # p at i1
    return p
end

"Compute flood depth by interpolating flood storage `flood_storage` using flood depth intervals."
function flood_depth(
    profile::FloodPlainProfile,
    flood_storage::Float,
    flow_length::Float,
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
    zb_floodplain::Vector{Float},
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
