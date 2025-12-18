"Struct for storing (shared) variables for river and overland flow models"
@with_kw struct FlowVariables
    n::Int
    q::Vector{Float64} = zeros(n)            # Discharge [m³ s⁻¹]
    qlat::Vector{Float64} = zeros(n)         # Lateral inflow per unit length [m² s⁻¹]
    qin::Vector{Float64} = zeros(n)          # Inflow from upstream cells [m³ s⁻¹]
    qin_av::Vector{Float64} = zeros(n)       # Average inflow from upstream cells  [m³ s⁻¹] for model timestep Δt
    q_av::Vector{Float64} = zeros(n)         # Average discharge [m³ s⁻¹] for model timestep Δt
    storage::Vector{Float64} = zeros(n)      # Kinematic wave storage [m³] (based on water depth h)
    h::Vector{Float64} = zeros(n)            # Water depth [m]
end

"Initialize timestepping for kinematic wave (river and overland flow models)"
function init_kinematic_wave_timestepping(config::Config, n::Int; domain::String)
    adaptive = config.model.kinematic_wave__adaptive_time_step_flag
    @info "Kinematic wave approach is used for $domain flow, adaptive timestepping = $adaptive."

    if adaptive
        stable_timesteps = zeros(n)
        timestepping = TimeStepping(; stable_timesteps, adaptive)
    else
        dt_fixed = getfield(config.model, Symbol("$(domain)_kinematic_wave__time_step"))
        @info "Using a fixed internal timestep (seconds) $dt_fixed for kinematic wave $domain flow."
        timestepping = TimeStepping(; dt_fixed, adaptive)
    end
    return timestepping
end

"Struct for storing Manning flow parameters"
@with_kw struct ManningFlowParameters
    n::Int
    beta::Float64                 # constant in Manning's equation [-]
    slope::Vector{Float64}        # Slope [m m⁻¹]
    mannings_n::Vector{Float64}   # Manning's roughness [s m⁻⅓]
    alpha_pow::Float64            # Used in the power part of alpha [-]
    alpha_term::Vector{Float64} = fill(MISSING_VALUE, n)   # Term used in computation of alpha [-]
    alpha::Vector{Float64} = fill(MISSING_VALUE, n)        # Constant in momentum equation A = alpha*Q^beta, based on Manning's equation [s3/5 m1/5]
end

"Initialize Manning flow parameters"
function ManningFlowParameters(mannings_n::Vector{Float64}, slope::Vector{Float64})
    n = length(slope)
    parameters = ManningFlowParameters(;
        n,
        beta = Float64(0.6),
        slope,
        mannings_n,
        alpha_pow = Float64((2.0 / 3.0) * 0.6),
    )
    return parameters
end

"Struct for storing river flow model parameters"
@with_kw struct RiverFlowParameters
    flow::ManningFlowParameters
    bankfull_depth::Vector{Float64} # Bankfull water level [m]
end

"Overload `getproperty` for river flow model parameters"
function Base.getproperty(v::RiverFlowParameters, s::Symbol)
    if s === :bankfull_depth
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
    (; slope) = domain.parameters
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

    flow_params = ManningFlowParameters(mannings_n, slope)
    parameters = RiverFlowParameters(; flow = flow_params, bankfull_depth)
    return parameters
end

"Struct for storing river flow model boundary conditions"
@with_kw struct RiverFlowBC{R}
    n::Int
    inwater::Vector{Float64} = zeros(n)                         # Lateral inflow [m³ s⁻¹]
    external_inflow::Vector{Float64} = zeros(n)                 # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    actual_external_abstraction_av::Vector{Float64} = zeros(n)  # Actual abstraction from external negative inflow [m³ s⁻¹]
    abstraction::Vector{Float64} = zeros(n)                     # Abstraction (computed as part of water demand and allocation) [m³ s⁻¹]
    reservoir::R                                                # Reservoir model struct of arrays
end

"Initialize river flow model boundary conditions"
function RiverFlowBC(
    dataset::NCDataset,
    config::Config,
    network::NetworkRiver,
    reservoir::Union{Reservoir, Nothing},
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

"River flow model using the kinematic wave method and the Manning flow equation"
@with_kw struct KinWaveRiverFlow{R <: RiverFlowBC, A <: AbstractAllocationModel} <:
                AbstractRiverFlowModel
    timestepping::TimeStepping
    boundary_conditions::R
    parameters::RiverFlowParameters
    variables::FlowVariables
    allocation::A   # Water allocation
end

"Initialize river flow model `KinWaveRiverFlow`"
function KinWaveRiverFlow(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    reservoir::Union{Reservoir, Nothing},
)
    (; indices) = domain.network
    n = length(indices)

    timestepping = init_kinematic_wave_timestepping(config, n; domain = "river")

    allocation = do_water_demand(config) ? AllocationRiver(; n) : NoAllocationRiver(n)

    variables = FlowVariables(; n)
    parameters = RiverFlowParameters(dataset, config, domain)
    boundary_conditions = RiverFlowBC(dataset, config, domain.network, reservoir)

    river_flow = KinWaveRiverFlow(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        allocation,
    )

    return river_flow
end

"Struct for storing overland flow model variables"
@with_kw struct OverLandFlowVariables
    n::Int
    flow::FlowVariables = FlowVariables(; n)
    to_river::Vector{Float64} = zeros(n) # Part of overland flow [m³ s⁻¹] that flows to the river
end

"Overload `getproperty` for overland flow model variables"
function Base.getproperty(v::OverLandFlowVariables, s::Symbol)
    if s === :to_river
        getfield(v, s)
    elseif s === :flow
        getfield(v, :flow)
    else
        getfield(getfield(v, :flow), s)
    end
end

"Struct for storing overland flow model boundary conditions"
@with_kw struct LandFlowBC
    n::Int
    inwater::Vector{Float64} = zeros(n) # Lateral inflow [m³ s⁻¹]
end

"Overland flow model using the kinematic wave method and the Manning flow{ equation"
@with_kw struct KinWaveOverlandFlow <: AbstractOverlandFlowModel
    timestepping::TimeStepping
    boundary_conditions::LandFlowBC
    parameters::ManningFlowParameters
    variables::OverLandFlowVariables
end

"Initialize Overland flow model `KinWaveOverlandFlow`"
function KinWaveOverlandFlow(dataset::NCDataset, config::Config, domain::DomainLand)
    (; indices) = domain.network
    (; slope) = domain.parameters
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
    parameters = ManningFlowParameters(mannings_n, slope)
    boundary_conditions = LandFlowBC(; n)

    overland_flow =
        KinWaveOverlandFlow(; timestepping, boundary_conditions, variables, parameters)

    return overland_flow
end

"""
Helper function to set reservoir variables to zero. This is done at the start of each
simulation timestep, during the timestep the total (weighted) sum is computed from values at
each internal timestep.
"""
function set_reservoir_vars!(reservoir::Reservoir)
    reservoir.boundary_conditions.inflow .= 0.0
    reservoir.boundary_conditions.actual_external_abstraction_av .= 0.0
    reservoir.variables.outflow_av .= 0.0
    reservoir.variables.actevap .= 0.0

    return nothing
end
set_reservoir_vars!(reservoir) = nothing

"""
Helper function to compute the average of reservoir variables. This is done at the end of
each simulation timestep.
"""
function average_reservoir_vars!(reservoir::Reservoir, dt::Float64)
    reservoir.variables.outflow_av ./= dt
    reservoir.boundary_conditions.inflow ./= dt
    reservoir.boundary_conditions.actual_external_abstraction_av ./= dt

    return nothing
end
average_reservoir_vars!(reservoir, dt) = nothing

"""
    set_flow_vars!(model::AbstractRiverFlowModel)

Helper functions to set river flow routing variables discharge and actual abstraction (based
on external negative inflow) from river to zero. This is done at the start of each
simulation timestep, during the timestep the total (weighted) sum is computed from values at
each sub timestep.
"""
function set_flow_vars!(model::AbstractRiverFlowModel)
    (; q_av) = model.variables
    (; actual_external_abstraction_av) = model.boundary_conditions
    q_av .= 0.0
    actual_external_abstraction_av .= 0.0
    return nothing
end

"""
    average_flow_vars!(model::AbstractRiverFlowModel, dt::Float64)

Helper functions to compute average river flow routing variables. This is done at the end of
each simulation timestep.
"""
function average_flow_vars!(model::AbstractRiverFlowModel, dt::Float64)
    (; q_av) = model.variables
    (; actual_external_abstraction_av) = model.boundary_conditions
    q_av ./= dt
    actual_external_abstraction_av ./= dt
    return nothing
end

"Update overland flow model `KinWaveOverlandFlow` for a single timestep"
function kinwave_land_update!(model::KinWaveOverlandFlow, domain::DomainLand, dt::Float64)
    (; order_of_subdomains, order_subdomain, subdomain_indices, upstream_nodes) =
        domain.network

    (; beta, alpha) = model.parameters
    (; h, q, q_av, storage, qin, qin_av, qlat, to_river) = model.variables
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
Update overland flow model `KinWaveOverlandFlow` for a single timestep `dt`. Timestepping within
`dt` is either with a fixed timestep `dt_fixed` or adaptive.
"""
function update!(model::KinWaveOverlandFlow, domain::DomainLand, dt::Float64)
    (; inwater) = model.boundary_conditions
    (; alpha_term, mannings_n, beta, alpha_pow, alpha) = model.parameters
    (; surface_flow_width, flow_length, slope) = domain.parameters
    (; q_av, qlat, qin_av, to_river) = model.variables
    (; adaptive) = model.timestepping

    @. alpha_term = pow(mannings_n / sqrt(slope), beta)
    # use fixed alpha value based flow width
    @. alpha = alpha_term * pow(surface_flow_width, alpha_pow)
    @. qlat = inwater / flow_length

    q_av .= 0.0
    to_river .= 0.0
    qin_av .= 0.0

    t = 0.0
    while t < dt
        dt_s =
            adaptive ? stable_timestep(model, flow_length, 0.02) :
            model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        kinwave_land_update!(model, domain, dt_s)
        t += dt_s
    end
    q_av ./= dt
    to_river ./= dt
    qin_av ./= dt
    return nothing
end

"Update river flow model `KinWaveRiverFlow` for a single timestep"
function kinwave_river_update!(
    model::KinWaveRiverFlow,
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
        model.boundary_conditions

    (; beta, alpha) = model.parameters
    (; flow_width, flow_length) = domain.parameters
    (; h, q, q_av, storage, qin, qin_av, qlat) = model.variables

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
                    update!(reservoir, i, net_inflow, dt, dt_forcing)

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
                crossarea = alpha[v] * pow(q[v], beta)
                h[v] = crossarea / flow_width[v]
                storage[v] = flow_length[v] * flow_width[v] * h[v]

                # average variables (here accumulated for model timestep Δt)
                q_av[v] += q[v] * dt
                qin_av[v] += qin[v] * dt
            end
        end
    end
end

"""
Update river flow model `KinWaveRiverFlow` for a single timestep `dt`. Timestepping within
`dt` is either with a fixed timestep `dt_fixed` or adaptive.
"""
function update!(model::KinWaveRiverFlow, domain::Domain, clock::Clock)
    (; reservoir, inwater) = model.boundary_conditions
    (; alpha_term, mannings_n, beta, alpha_pow, alpha, bankfull_depth) = model.parameters
    (; slope, flow_width, flow_length) = domain.river.parameters
    (; qlat, qin_av) = model.variables
    (; adaptive) = model.timestepping

    @. alpha_term = pow(mannings_n / sqrt(slope), beta)
    # use fixed alpha value based on 0.5 * bankfull_depth
    @. alpha = alpha_term * pow(flow_width + bankfull_depth, alpha_pow)
    @. qlat = inwater / flow_length

    set_flow_vars!(model)
    qin_av .= 0.0
    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)

    dt = tosecond(clock.dt)
    t = 0.0
    while t < dt
        dt_s =
            adaptive ? stable_timestep(model, flow_length, 0.05) :
            model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        kinwave_river_update!(model, domain.river, dt_s, dt)
        t += dt_s
    end

    average_reservoir_vars!(reservoir, dt)
    average_flow_vars!(model, dt)
    qin_av ./= dt
    return nothing
end

"""
Compute a stable timestep size for the kinematice wave method for a river or overland flow
model using a nonlinear scheme (Chow et al., 1988).

A stable time step is computed for each vector element based on the Courant timestep size
criterion. A quantile of the vector is computed based on probability `p` to remove potential
very low timestep sizes. Li et al. (1975) found that the nonlinear scheme is unconditonally
stable and that a wide range of dt/dx values can be used without loss of accuracy.
"""
function stable_timestep(
    model::S,
    flow_length::Vector{Float64},
    p::Float64,
) where {S <: Union{KinWaveOverlandFlow, KinWaveRiverFlow}}
    (; q) = model.variables
    (; alpha, beta) = model.parameters
    (; stable_timesteps) = model.timestepping

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
    model::AbstractRiverFlowModel,
    external_models::NamedTuple,
    domain::Domain,
    dt::Float64,
)
    (; allocation, runoff, overland_flow, subsurface_flow) = external_models
    (; inwater) = model.boundary_conditions
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
`KinWaveOverlandFlow` for a single timestep.
"""
function update_lateral_inflow!(
    model::KinWaveOverlandFlow,
    external_models::NamedTuple,
    area::Vector{Float64},
    config::Config,
    dt::Float64,
)
    (; soil, subsurface_flow, allocation) = external_models
    (; net_runoff) = soil.variables
    (; inwater) = model.boundary_conditions

    if config.model.drain__flag
        drain = subsurface_flow.boundaries.drain
        drainflux = zeros(length(net_runoff))
        drainflux[drain.index] = -drain.variables.flux ./ tosecond(BASETIMESTEP)
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
    model::Union{Reservoir, Nothing},
    river_flow::AbstractRiverFlowModel,
    external_models::NamedTuple,
    network::NetworkReservoir,
)
    (; overland_flow, subsurface_flow) = external_models
    (; land_indices) = network
    if !isnothing(model)
        (; inflow_overland, inflow_subsurface) = model.boundary_conditions
        inflow_overland .= get_inflow_reservoir(river_flow, overland_flow, land_indices)
        inflow_subsurface .= get_inflow_reservoir(river_flow, subsurface_flow, land_indices)
    end
    return nothing
end

# For the river kinematic wave, the variable `to_river` can be excluded, because this part
# is added to the river kinematic wave.
get_inflow_reservoir(::KinWaveRiverFlow, model::KinWaveOverlandFlow, inds::Vector{Int}) =
    model.variables.q_av[inds]
get_inflow_reservoir(::KinWaveRiverFlow, model::LateralSSF, inds::Vector{Int}) =
    model.variables.ssf[inds] ./ tosecond(BASETIMESTEP)

# Exclude subsurface flow from `GroundwaterFlow`.
get_inflow_reservoir(::AbstractRiverFlowModel, model::GroundwaterFlow, inds::Vector{Int}) =
    zeros(length(inds))
