"Struct for storing variables for river flow model"
@with_kw struct RiverFlowVariables <: AbstractRiverFlowVariables
    n::Int
    q::Vector{Float64} = zeros(n)            # River discharge [m³ s⁻¹]
    qlat::Vector{Float64} = zeros(n)         # Lateral inflow per unit length [m² s⁻¹]
    qin::Vector{Float64} = zeros(n)          # River inflow from upstream cells [m³ s⁻¹]
    qin_av::Vector{Float64} = zeros(n)       # Average river inflow from upstream cells [m³ s⁻¹] for model timestep Δt
    q_av::Vector{Float64}                    # Average river channel (+ floodplain) discharge discharge [m³ s⁻¹] for model timestep Δt
    q_channel_av::Vector{Float64}            # Average river channel discharge [m³ s⁻¹] (for model timestep Δt)
    storage::Vector{Float64} = zeros(n)      # River kinematic wave storage [m³] (based on water depth h)
    h::Vector{Float64} = zeros(n)            # River water depth [m]
end

"Initialize Manning flow parameters"
function ManningFlowParameters(
    mannings_n::Vector{Float64},
    slope::Vector{Float64},
    wetted_perimeter::Vector{Float64},
)
    beta = Float64(0.6)
    hydraulic_radius_pow = Float64(2.0 / 3.0)
    alpha_term = @. pow(mannings_n / sqrt(slope), beta)
    alpha_pow = hydraulic_radius_pow * beta
    alpha = @. alpha_term * pow(wetted_perimeter, alpha_pow)

    parameters =
        ManningFlowParameters(; beta, slope, mannings_n, alpha_pow, alpha_term, alpha)
    return parameters
end

"Struct for storing river flow model parameters"
@with_kw struct RiverFlowParameters <: AbstractRiverFlowParameters
    flow::ManningFlowParameters
    bankfull_depth::Vector{Float64}     # Bankfull water level [m]
    bankfull_storage::Vector{Float64}   # Bankfull storage [m³]
    bankfull_flow::Vector{Float64}      # Bankfull discharge [m³ s⁻¹]
end

"Overload `getproperty` for river flow model parameters"
function Base.getproperty(v::RiverFlowParameters, s::Symbol)
    if s === :bankfull_depth
        getfield(v, s)
    elseif s === :bankfull_storage
        getfield(v, s)
    elseif s === :bankfull_flow
        getfield(v, s)
    elseif s === :flow
        getfield(v, :flow)
    else
        getfield(getfield(v, :flow), s)
    end
end

"Initialize river flow model parameters"
function RiverFlowParameters(dataset::NCDataset, config::Config, domain::DomainRiver)
    (; indices) = domain.network
    (; slope, flow_length, flow_width) = domain.parameters
    mannings_n = ncread(
        dataset,
        config,
        "river_water_flow__manning_n_parameter";
        sel = indices,
        defaults = 0.036,
        type = Float64,
    )
    bankfull_depth = ncread(
        dataset,
        config,
        "river_bank_water__depth";
        sel = indices,
        defaults = 1.0,
        type = Float64,
    )

    # use fixed wetted perimeter based on 0.5 * bankfull_depth
    wetted_perimeter = flow_width + bankfull_depth
    flow_params = ManningFlowParameters(mannings_n, slope, wetted_perimeter)
    bankfull_storage = bankfull_depth .* flow_length .* flow_width
    if config.model.floodplain_1d__flag
        wetted_perimeter = @. wetted_perimeter_channel(bankfull_depth, flow_width)
        alpha = @. flow_params.alpha_term * pow(wetted_perimeter, flow_params.alpha_pow)
        bankfull_flow = @. pow(bankfull_depth * flow_width / alpha, 1.0 / flow_params.beta)
    else
        bankfull_flow = Float64[]
    end
    parameters = RiverFlowParameters(;
        flow = flow_params,
        bankfull_depth,
        bankfull_storage,
        bankfull_flow,
    )
    return parameters
end

"Initialize river flow model boundary conditions"
function RiverFlowBC(
    dataset::NCDataset,
    config::Config,
    network::NetworkRiver,
    reservoir::Union{ReservoirModel, Nothing},
)
    (; indices) = network
    external_inflow = ncread(
        dataset,
        config,
        "river_water__external_inflow_volume_flow_rate";
        sel = indices,
        defaults = 0.0,
        type = Float64,
    )
    n = length(indices)
    bc = RiverFlowBC(; n, external_inflow, reservoir)
    return bc
end

"Initialize kinematic wave river flow model"
function init_kinematic_wave_river_flow(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    reservoir::Union{ReservoirModel, Nothing},
)
    (; indices) = domain.network
    n = length(indices)

    timestepping = init_kinematic_wave_timestepping(config, n; domain = "river")

    allocation =
        do_water_demand(config) ? AllocationRiverModel(; n) : NoAllocationRiverModel(n)

    q_av = zeros(n)
    variables = RiverFlowVariables(;
        n,
        q_av,
        q_channel_av = config.model.floodplain_1d__flag ? zeros(n) : q_av,
    )
    parameters = RiverFlowParameters(dataset, config, domain)
    boundary_conditions = RiverFlowBC(dataset, config, domain.network, reservoir)
    if config.model.floodplain_1d__flag
        floodplain = FloodPlainModel(dataset, config, domain)
    else
        floodplain = nothing
    end
    routing_method = KinematicWave()

    river_flow = RiverFlowModel(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain,
        allocation,
        routing_method,
    )

    return river_flow
end

"Struct for storing overland flow model variables"
@with_kw struct OverLandFlowVariables <: AbstractOverlandFlowVariables
    n::Int
    q::Vector{Float64} = zeros(n)            # Overland discharge [m³ s⁻¹]
    qlat::Vector{Float64} = zeros(n)         # Lateral inflow per unit length [m² s⁻¹]
    qin::Vector{Float64} = zeros(n)          # Overland inflow from upstream cells [m³ s⁻¹]
    qin_av::Vector{Float64} = zeros(n)       # Average overland inflow from upstream cells [m³ s⁻¹] for model timestep Δt
    q_av::Vector{Float64} = zeros(n)         # Average overland discharge [m³ s⁻¹] for model timestep Δt
    storage::Vector{Float64} = zeros(n)      # Overland kinematic wave storage [m³] (based on water depth h)
    h::Vector{Float64} = zeros(n)            # Overland water depth [m]
    to_river::Vector{Float64} = zeros(n)     # Part of overland flow [m³ s⁻¹] that flows to the river
end

"Struct for storing overland flow model boundary conditions"
@with_kw struct LandFlowBC <: AbstractOverlandFlowBC
    n::Int
    inwater::Vector{Float64} = zeros(n) # Lateral inflow [m³ s⁻¹]
end

"Initialize kinematic wave overland flow model"
function init_kinematic_wave_overland_flow(
    dataset::NCDataset,
    config::Config,
    domain::DomainLand,
)
    (; indices) = domain.network
    (; slope, surface_flow_width) = domain.parameters
    mannings_n = ncread(
        dataset,
        config,
        "land_surface_water_flow__manning_n_parameter";
        sel = indices,
        defaults = 0.072,
        type = Float64,
    )

    n = length(indices)
    timestepping = init_kinematic_wave_timestepping(config, n; domain = "land")

    variables = OverLandFlowVariables(; n)
    parameters = ManningFlowParameters(mannings_n, slope, surface_flow_width)
    boundary_conditions = LandFlowBC(; n)
    routing_method = KinematicWave()

    overland_flow = OverlandFlowModel(;
        timestepping,
        boundary_conditions,
        variables,
        parameters,
        routing_method,
    )

    return overland_flow
end

"""
Helper function to set reservoir variables to zero. This is done at the start of each
simulation timestep, during the timestep the total (weighted) sum is computed from values at
each internal timestep.
"""
function set_reservoir_vars!(reservoir_model::ReservoirModel)
    reservoir_model.boundary_conditions.inflow .= 0.0
    reservoir_model.boundary_conditions.actual_external_abstraction_av .= 0.0
    reservoir_model.variables.outflow_av .= 0.0
    reservoir_model.variables.actevap .= 0.0

    return nothing
end
set_reservoir_vars!(reservoir_model) = nothing

"""
Helper function to compute the average of reservoir variables. This is done at the end of
each simulation timestep.
"""
function average_reservoir_vars!(reservoir_model::ReservoirModel, dt::Float64)
    reservoir_model.variables.outflow_av ./= dt
    reservoir_model.boundary_conditions.inflow ./= dt
    reservoir_model.boundary_conditions.actual_external_abstraction_av ./= dt

    return nothing
end
average_reservoir_vars!(reservoir_model, dt) = nothing

"""
    set_flow_vars!(river_flow_model::AbstractRiverFlowModel)

Helper functions to set river flow routing variables discharge and actual abstraction (based
on external negative inflow) from river to zero. This is done at the start of each
simulation timestep, during the timestep the total (weighted) sum is computed from values at
each sub timestep.
"""
function set_flow_vars!(river_flow_model::AbstractRiverFlowModel)
    (; q_av) = river_flow_model.variables
    (; actual_external_abstraction_av) = river_flow_model.boundary_conditions
    q_av .= 0.0
    actual_external_abstraction_av .= 0.0
    return nothing
end

"""
    average_flow_vars!(river_flow::AbstractRiverFlowModel, dt::Float64)

Helper functions to compute average river flow routing variables. This is done at the end of
each simulation timestep.
"""
function average_flow_vars!(river_flow_model::AbstractRiverFlowModel, dt::Float64)
    (; q_av) = river_flow_model.variables
    (; actual_external_abstraction_av) = river_flow_model.boundary_conditions
    q_av ./= dt
    actual_external_abstraction_av ./= dt
    return nothing
end

"Update overland flow model `OverlandFlowModel{<:KinematicWave}` for a single timestep"
function kinwave_land_update!(
    overland_flow_model::OverlandFlowModel{<:KinematicWave},
    domain::DomainLand,
    dt::Float64,
)
    (; order_of_subdomains, order_subdomain, subdomain_indices, upstream_nodes) =
        domain.network

    (; beta, alpha) = overland_flow_model.parameters
    (; h, q, q_av, storage, qin, qin_av, qlat, to_river) = overland_flow_model.variables
    (; surface_flow_width, flow_length, flow_fraction_to_river) = domain.parameters

    ns = length(order_of_subdomains)
    qin .= 0.0
    for k in 1:ns
        threaded_foreach(eachindex(order_of_subdomains[k]); basesize = 1) do i
            m = order_of_subdomains[k][i]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                # for a river cell without a reservoir part of the upstream surface flow
                # goes to the river (flow_fraction_to_river) and part goes to the surface
                # flow reservoir (1.0 - flow_fraction_to_river), upstream nodes with a
                # reservoir are excluded
                to_river[v] +=
                    sum_at(i -> q[i] * flow_fraction_to_river[i], upstream_nodes[n]) * dt
                if surface_flow_width[v] > 0.0
                    qin[v] = sum_at(
                        i -> q[i] * (1.0 - flow_fraction_to_river[i]),
                        upstream_nodes[n],
                    )
                end

                q[v] = kinematic_wave(
                    qin[v],
                    q[v],
                    qlat[v],
                    alpha[v],
                    beta,
                    dt,
                    flow_length[v],
                )

                # update h, only if flow width > 0.0
                if surface_flow_width[v] > 0.0
                    crossarea = alpha[v] * pow(q[v], beta)
                    h[v] = crossarea / surface_flow_width[v]
                end
                storage[v] = flow_length[v] * surface_flow_width[v] * h[v]

                # average flow (here accumulated for model timestep Δt)
                q_av[v] += q[v] * dt
                qin_av[v] += qin[v] * dt
            end
        end
    end
end

"""
Update overland flow model `OverlandFlowModel{<:KinematicWave}` for a single timestep `dt`.
Timestepping within `dt` is either with a fixed timestep `dt_fixed` or adaptive.
"""
function update_overland_flow_model!(
    overland_flow_model::OverlandFlowModel{<:KinematicWave},
    domain::DomainLand,
    dt::Float64,
)
    (; inwater) = overland_flow_model.boundary_conditions
    (; flow_length) = domain.parameters
    (; q_av, qlat, qin_av, to_river) = overland_flow_model.variables
    (; adaptive) = overland_flow_model.timestepping

    @. qlat = inwater / flow_length

    q_av .= 0.0
    to_river .= 0.0
    qin_av .= 0.0

    t = 0.0
    while t < dt
        dt_s =
            adaptive ? stable_timestep(overland_flow_model, flow_length, 0.02) :
            overland_flow_model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        kinwave_land_update!(overland_flow_model, domain, dt_s)
        t += dt_s
    end
    q_av ./= dt
    to_river ./= dt
    qin_av ./= dt
    return nothing
end

function update_floodplain_model!(
    river_flow_model::RiverFlowModel{T, F},
    domain::DomainRiver,
    v::Int,
    dt::Float64,
) where {T <: KinematicWave, F <: FloodPlainModel}
    (; floodplain) = river_flow_model
    (; beta, alpha, bankfull_depth, bankfull_flow) = river_flow_model.parameters
    (; flow_length, flow_width) = domain.parameters
    (; q, h, storage) = river_flow_model.variables

    flow_area = alpha[v] * pow(q[v], beta)
    if q[v] > bankfull_flow[v]
        bankfull_area = flow_width[v] * bankfull_depth[v]
        flood_area = flow_area - bankfull_area
        flood_depth, flood_storage = compute_flood_depth_storage(
            floodplain.parameters.profile,
            flood_area,
            flow_length[v],
            v,
        )
        floodplain_flow = q[v] - bankfull_flow[v]
        # update total river and floodplain depth and storage
        h[v] = bankfull_depth[v] + flood_depth
        storage[v] = bankfull_area * flow_length[v] + flood_storage
    else
        flood_storage = 0.0
        flood_depth = 0.0
        floodplain_flow = 0.0
        # update river depth and storage
        h[v] = flow_area / flow_width[v]
        storage[v] = flow_length[v] * flow_width[v] * h[v]
    end
    floodplain.variables.storage[v] = flood_storage
    floodplain.variables.h[v] = flood_depth
    floodplain.variables.q[v] = floodplain_flow
    floodplain.variables.q_av[v] += floodplain_flow * dt
end

"Update river flow model `RiverFlowModel{<:KinematicWave}` for a single timestep"
function kinwave_river_update!(
    river_flow_model::RiverFlowModel{<:KinematicWave},
    domain::DomainRiver,
    dt::Float64,
    dt_forcing::Float64,
)
    (;
        graph,
        order_of_subdomains,
        order_subdomain,
        subdomain_indices,
        upstream_nodes,
        reservoir_indices,
    ) = domain.network

    (; reservoir, external_inflow, actual_external_abstraction_av, abstraction) =
        river_flow_model.boundary_conditions

    (; beta, alpha) = river_flow_model.parameters
    (; flow_width, flow_length) = domain.parameters
    (; h, q, q_av, storage, qin, qin_av, qlat) = river_flow_model.variables
    (; floodplain) = river_flow_model

    if !isnothing(reservoir)
        res_bc = reservoir.boundary_conditions
    end

    ns = length(order_of_subdomains)
    qin .= 0.0
    for k in 1:ns
        threaded_foreach(eachindex(order_of_subdomains[k]); basesize = 1) do i
            m = order_of_subdomains[k][i]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                # qin by outflow from upstream reservoir location is added
                qin[v] += sum_at(q, upstream_nodes[n])
                # Inflow supply/abstraction is added to qlat (divide by flow length)
                # If external_inflow < 0, abstraction is limited
                if external_inflow[v] < 0.0
                    _abstraction = min(-external_inflow[v], (storage[v] / dt) * 0.80)
                    actual_external_abstraction_av[v] += _abstraction * dt
                    _inflow = -_abstraction / flow_length[v]
                else
                    _inflow = external_inflow[v] / flow_length[v]
                end
                # internal abstraction (water demand) is limited by river storage and
                # negative external inflow as part of water allocation computations.
                _inflow -= abstraction[v] / flow_length[v]

                q[v] = kinematic_wave(
                    qin[v],
                    q[v],
                    qlat[v] + _inflow,
                    alpha[v],
                    beta,
                    dt,
                    flow_length[v],
                )

                if !isnothing(reservoir) && reservoir_indices[v] != 0
                    # run reservoir model and copy reservoir outflow to inflow (qin) of
                    # downstream river cell
                    i = reservoir_indices[v]
                    # If external_inflow < 0, abstraction is limited
                    if res_bc.external_inflow[i] < 0.0
                        _abstraction = min(
                            -res_bc.external_inflow[i],
                            (reservoir.variables.storage[i] / dt) * 0.98,
                        )
                        res_bc.actual_external_abstraction_av[i] += _abstraction * dt
                        _inflow = -_abstraction
                    else
                        _inflow = res_bc.external_inflow[i]
                    end
                    net_inflow =
                        q[v] +
                        res_bc.inflow_overland[i] +
                        res_bc.inflow_subsurface[i] +
                        _inflow
                    update_reservoir_model!(reservoir, i, net_inflow, dt, dt_forcing)

                    downstream_nodes = outneighbors(graph, v)
                    n_downstream = length(downstream_nodes)
                    if n_downstream == 1
                        j = only(downstream_nodes)
                        qin[j] = reservoir.variables.outflow[i]
                    elseif n_downstream == 0
                        error(
                            """A reservoir without a downstream river node is not supported.
                            Add a downstream river node or move the reservoir to an upstream node (model schematization).
                            """,
                        )
                    else
                        error("bifurcations not supported")
                    end
                end
                # update h and storage
                flow_area = alpha[v] * pow(q[v], beta)
                if isnothing(floodplain)
                    h[v] = flow_area / flow_width[v]
                    storage[v] = flow_length[v] * flow_width[v] * h[v]
                else
                    update_floodplain_model!(river_flow_model, domain, v, dt)
                end
                # average variables (here accumulated for model timestep Δt)
                q_av[v] += q[v] * dt
                qin_av[v] += qin[v] * dt
            end
        end
    end
end

function update_alpha_parameter!(
    river_flow_model::RiverFlowModel{T, F},
    parameters::RiverParameters,
) where {T <: KinematicWave, F <: FloodPlainModel}
    (; h) = river_flow_model.variables
    (; bankfull_depth, alpha, mannings_n, beta, alpha_term, alpha_pow) =
        river_flow_model.parameters
    (; flow_width, slope) = parameters

    floodplain_p = river_flow_model.floodplain.parameters
    floodplain_v = river_flow_model.floodplain.variables

    for i in eachindex(h)
        if h[i] > bankfull_depth[i]
            wp_channel = wetted_perimeter_channel(bankfull_depth[i], flow_width[i])
            i1, i2 = interpolation_indices(floodplain_v.h[i], floodplain_p.profile.depth)
            wp_floodplain =
                compute_wetted_perimeter(floodplain_p.profile, floodplain_v.h[i], i, i1)
            wp_combined = wp_channel + wp_floodplain
            mannings_n_combined = pow(
                wp_floodplain/wp_combined * pow(floodplain_p.mannings_n[i], 1.5) +
                wp_channel/wp_combined * pow(mannings_n[i], 1.5),
                2.0/3.0,
            )
            alpha[i] = alpha_term[i] * pow(wp_combined, alpha_pow)
        else
            wp_channel = wetted_perimeter_channel(h[i], flow_width[i])
            alpha[i] = alpha_term[i] * pow(wp_channel, alpha_pow)
        end
    end
end

update_alpha_parameter!(
    river_flow_model::RiverFlowModel{T, F},
    parameters::RiverParameters,
) where {T <: KinematicWave, F <: Nothing} = nothing

"""
Update river flow model `RiverFlowModel{<:KinematicWave}` for a single timestep `dt`.
Timestepping within `dt` is either with a fixed timestep `dt_fixed` or adaptive.
"""
function update_river_flow_model!(
    river_flow_model::RiverFlowModel{<:KinematicWave},
    domain::Domain,
    clock::Clock,
)
    (; floodplain) = river_flow_model
    (; reservoir, inwater) = river_flow_model.boundary_conditions
    (; flow_length) = domain.river.parameters
    (; qlat, qin_av) = river_flow_model.variables
    (; adaptive) = river_flow_model.timestepping

    @. qlat = inwater / flow_length
    set_flow_vars!(river_flow_model)
    qin_av .= 0.0
    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)
    if !isnothing(floodplain)
        floodplain.variables.q_av .= 0.0
    end

    dt = tosecond(clock.dt)
    t = 0.0
    while t < dt
        update_alpha_parameter!(river_flow_model, domain.river.parameters)
        dt_s =
            adaptive ? stable_timestep(river_flow_model, flow_length, 0.05) :
            river_flow_model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        kinwave_river_update!(river_flow_model, domain.river, dt_s, dt)
        t += dt_s
    end

    average_reservoir_vars!(reservoir, dt)
    average_flow_vars!(river_flow_model, dt)
    qin_av ./= dt
    if !isnothing(floodplain)
        floodplain.variables.q_av ./= dt
        river_flow_model.variables.q_channel_av .=
            river_flow_model.variables.q_av .- floodplain.variables.q_av
    end
    return nothing
end

"""
Compute a stable timestep size for the kinematice wave method for a river or overland flow
model using a nonlinear scheme (Chow et al., 1988).

A stable time step is computed for each vector element based on the Courant timestep size
criterion. A quantile of the vector is computed based on probability `p` to remove potential
very small timestep sizes. Li et al. (1975) found that the nonlinear scheme is
unconditionally stable and that a wide range of dt/dx values can be used without loss of
accuracy.
"""
function stable_timestep(
    flow_model::S,
    flow_length::Vector{Float64},
    p::Float64,
) where {S <: Union{RiverFlowModel{<:KinematicWave}, OverlandFlowModel{<:KinematicWave}}}
    (; q) = flow_model.variables
    (; alpha, beta) = flow_model.parameters
    (; stable_timesteps) = flow_model.timestepping

    n = length(q)
    stable_timesteps .= Inf
    k = 0
    for i in 1:n
        if q[i] > 0.0
            k += 1
            c = 1.0 / (alpha[i] * beta * pow(q[i], (beta - 1.0)))
            stable_timesteps[k] = (flow_length[i] / c)
        end
    end

    dt_min = if k == 1
        stable_timesteps[k]
    elseif k > 0
        quantile!(@view(stable_timesteps[1:k]), p)
    else
        600.0
    end

    return dt_min
end

"""
Update boundary condition lateral inflow `inwater` of a river flow model for a single
timestep.
"""
function update_lateral_inflow!(
    river_flow_model::AbstractRiverFlowModel,
    external_models::NamedTuple,
    domain::Domain,
    dt::Float64,
)
    (; allocation, runoff, overland_flow, subsurface_flow) = external_models
    (; inwater) = river_flow_model.boundary_conditions
    (; net_runoff_river) = runoff.variables

    (; land_indices) = domain.river.network
    (; cell_area) = domain.river.parameters
    (; area) = domain.land.parameters

    inwater .= (
        get_flux_to_river(subsurface_flow, land_indices) .+
        overland_flow.variables.to_river[land_indices] .+
        (net_runoff_river[land_indices] .* area[land_indices] .* 0.001) ./ dt .+
        (get_nonirrigation_returnflow(allocation) .* 0.001 .* cell_area) ./ dt
    )
    return nothing
end

"""
Update boundary condition lateral inflow `inwater` of a kinematic wave overland flow model
`OverlandFlowModel{<:KinematicWave}` for a single timestep.
"""
function update_lateral_inflow!(
    overland_flow_model::OverlandFlowModel{<:KinematicWave},
    external_models::NamedTuple,
    domain::Domain,
    config::Config,
    dt::Float64,
)
    (; soil, subsurface_flow, allocation) = external_models
    (; net_runoff) = soil.variables
    (; inwater) = overland_flow_model.boundary_conditions

    (; area) = domain.land.parameters
    (; land_indices) = domain.drain.network

    if config.model.drain__flag
        drain = subsurface_flow.boundary_conditions.drain
        drainflux = zeros(length(net_runoff))
        drainflux[land_indices] = -drain.variables.flux ./ tosecond(BASETIMESTEP)
    else
        drainflux = 0.0
    end
    inwater .=
        (net_runoff .+ get_nonirrigation_returnflow(allocation)) .* area * 0.001 ./ dt .+
        drainflux

    return nothing
end

"""
Update overland and subsurface flow contribution to inflow of a reservoir model for a river
flow model `AbstractRiverFlowModel` for a single timestep.
"""
function update_inflow!(
    reservoir_model::Union{ReservoirModel, Nothing},
    river_flow_model::AbstractRiverFlowModel,
    external_models::NamedTuple,
    network::NetworkReservoir,
)
    (; overland_flow, subsurface_flow) = external_models
    (; land_indices) = network
    if !isnothing(reservoir_model)
        (; inflow_overland, inflow_subsurface) = reservoir_model.boundary_conditions
        inflow_overland .=
            get_inflow_reservoir(river_flow_model, overland_flow, land_indices)
        inflow_subsurface .=
            get_inflow_reservoir(river_flow_model, subsurface_flow, land_indices)
    end
    return nothing
end

# For the river kinematic wave, the variable `to_river` can be excluded, because this part
# is added to the river kinematic wave.
get_inflow_reservoir(
    ::RiverFlowModel{<:KinematicWave},
    overland_flow_model::OverlandFlowModel{<:KinematicWave},
    inds::Vector{Int},
) = overland_flow_model.variables.q_av[inds]
get_inflow_reservoir(
    ::RiverFlowModel{<:KinematicWave},
    subsurface_flow_model::LateralSSFModel,
    inds::Vector{Int},
) = subsurface_flow_model.variables.q_av[inds] ./ tosecond(BASETIMESTEP)

# Exclude subsurface flow from `GroundwaterFlowModel`.
get_inflow_reservoir(::AbstractRiverFlowModel, ::GroundwaterFlowModel, inds::Vector{Int}) =
    zeros(length(inds))
