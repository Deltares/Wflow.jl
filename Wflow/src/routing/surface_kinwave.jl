"Struct for storing (shared) variables for river and overland flow models"
@with_kw struct FlowVariables
    n::Int
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = zeros(n)
    # Lateral inflow per unit length [m² s⁻¹]
    qlat::Vector{Float64} = zeros(n)
    # Inflow from upstream cells [m³ s⁻¹]
    qin::Vector{Float64} = zeros(n)
    # Average inflow from upstream cells  [m³ s⁻¹] for model time step dt
    qin_av::AverageVector = AverageVector(; n)
    # Average discharge [m³ s⁻¹] for model time step dt
    q_av::AverageVector = AverageVector(; n)
    # Kinematic wave storage [m³] (based on water depth h)
    storage::Vector{Float64} = zeros(n)
    # Water depth [m]
    h::Vector{Float64} = zeros(n)
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
    # constant in Manning's equation [-]
    beta::Float64 = 0.6
    # Slope [m m⁻¹]
    slope::Vector{Float64}
    # Manning's roughness [s m⁻⅓]
    mannings_n::Vector{Float64}
    # Used in the power part of alpha [-]
    alpha_pow::Float64 = (2 // 3) * 0.6
    # Term used in computation of alpha [s^3/5 m^-1/5]
    alpha_term::Vector{Float64} = fill(MISSING_VALUE, length(slope))
    # Constant in momentum equation A = alpha*Q^beta, based on Manning's equation [s^3/5 m^1/5]
    alpha::Vector{Float64} = fill(MISSING_VALUE, length(slope))
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
        "river_water_flow__manning_n_parameter",
        Routing;
        sel = indices,
        defaults = 0.036,
        type = Float64,
    )
    bankfull_depth = ncread(
        dataset,
        config,
        "river_bank_water__depth",
        Routing;
        sel = indices,
        defaults = 1.0,
        type = Float64,
    )

    flow_params = ManningFlowParameters(; mannings_n, slope)
    parameters = RiverFlowParameters(; flow = flow_params, bankfull_depth)
    return parameters
end

"Struct for storing river flow model boundary conditions"
@with_kw struct RiverFlowBC{R}
    # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    external_inflow::Vector{Float64}
    # Lateral inflow [m³ s⁻¹]
    inwater::Vector{Float64} = zeros(length(external_inflow))
    # Actual abstraction from external negative inflow [m³ s⁻¹]
    actual_external_abstraction_av::AverageVector =
        AverageVector(; n = length(external_inflow))
    # Abstraction (computed as part of water demand and allocation) [m³ s⁻¹]
    abstraction::Vector{Float64} = zeros(length(external_inflow))
    # Reservoir model struct of arrays
    reservoir::R
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
        "river_water__external_inflow_volume_flow_rate",
        Routing;
        sel = indices,
        defaults = 0.0,
        type = Float64,
    )
    bc = RiverFlowBC(; external_inflow, reservoir)
    return bc
end

"River flow model using the kinematic wave method and the Manning flow equation"
@with_kw struct KinWaveRiverFlow{R, A} <: AbstractRiverFlowModel
    timestepping::TimeStepping
    boundary_conditions::RiverFlowBC{R}
    parameters::RiverFlowParameters
    variables::FlowVariables
    # Water allocation
    allocation::A
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
    # Part of overland flow [m³ s⁻¹] that flows to the river
    to_river::AverageVector = AverageVector(; n)
end

"Overload `getproperty` for overland flow model variables"
function Base.getproperty(v::OverLandFlowVariables, s::Symbol)
    if s === :to_river
        getfield(v, s)
    elseif s === :flow
        getfield(v, s)
    else
        getfield(getfield(v, :flow), s)
    end
end

"Struct for storing overland flow model boundary conditions"
@with_kw struct LandFlowBC
    n::Int
    # Lateral inflow [m³ s⁻¹]
    inwater::Vector{Float64} = zeros(n)
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
        "land_surface_water_flow__manning_n_parameter",
        Routing;
        sel = indices,
        defaults = 0.072,
        type = Float64,
    )

    n = length(indices)
    timestepping = init_kinematic_wave_timestepping(config, n; domain = "land")

    variables = OverLandFlowVariables(; n)
    parameters = ManningFlowParameters(; mannings_n, slope)
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
    (; boundary_conditions, variables) = reservoir
    (; inflow, actual_external_abstraction_av) = boundary_conditions
    (; outflow_av, actevap) = variables

    zero!(inflow)
    zero!(actual_external_abstraction_av)
    zero!(outflow_av)
    zero!(actevap)
    return nothing
end
set_reservoir_vars!(reservoir) = nothing

"""
Helper function to compute the average of reservoir variables. This is done at the end of
each simulation timestep.
"""
function average_reservoir_vars!(reservoir::Reservoir, dt::Float64)
    (; variables, boundary_conditions) = reservoir
    (; outflow_av) = variables
    (; inflow, actual_external_abstraction_av) = boundary_conditions
    average!(outflow_av, dt)
    average!(inflow, dt)
    average!(actual_external_abstraction_av, dt)
    return nothing
end
average_reservoir_vars!(reservoir, dt) = nothing

"""
    set_flow_vars!(model::AbstractRiverFlowModel)

Helper functions to set cumulative river flow routing variables discharge and actual abstraction (based
on external negative inflow) from river to zero. This is done at the start of each
simulation timestep, during the timestep the total (weighted) sum is computed from values at
each sub timestep.
"""
function set_flow_vars!(model::AbstractRiverFlowModel)
    (; q_av) = model.variables
    (; actual_external_abstraction_av) = model.boundary_conditions
    zero!(q_av)
    zero!(actual_external_abstraction_av)
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
    average!(q_av, dt)
    average!(actual_external_abstraction_av, dt)
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
                # [m³ s⁻¹] = ∑ [m³ s⁻¹] * [-]
                add_to_cumulative!(
                    to_river,
                    v,
                    sum_at(i -> q[i] * flow_fraction_to_river[i], upstream_nodes[n]),
                )
                if surface_flow_width[v] > 0.0
                    # [m³ s⁻¹] = ∑ [m³ s⁻¹] * [-]
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
                    # [m²] = [s^3/5 m^1/5] * [m³ s⁻¹] ^ (3/5)
                    crossarea = alpha[v] * pow(q[v], beta)
                    # [m] = [m²] / [m]
                    h[v] = crossarea / surface_flow_width[v]
                end
                # [m³] = [m] * [m] * [m]
                storage[v] = flow_length[v] * surface_flow_width[v] * h[v]

                # average flow
                add_to_cumulative!(q_av, v, q[v])
                add_to_cumulative!(qin_av, v, qin[v])
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

    zero!(q_av)
    zero!(qin_av)
    zero!(to_river)

    t = 0.0
    while t < dt
        dt_s =
            adaptive ? stable_timestep(model, flow_length, 0.02) :
            model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        kinwave_land_update!(model, domain, dt_s)
        t += dt_s
    end
    average!(q_av, dt)
    average!(to_river, dt)
    average!(qin_av, dt)
    return nothing
end

"Update river flow model `KinWaveRiverFlow` for a single timestep"
function kinwave_river_update!(model::KinWaveRiverFlow, domain::DomainRiver, dt::Float64)
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
                    # [m³ s⁻¹] = min([m³ s⁻¹], [m³] / [s])
                    _abstraction = min(-external_inflow[v], (storage[v] / dt) * 0.80)
                    # [m³] += [m³ s⁻¹] * [s]
                    add_to_cumulative!(actual_external_abstraction_av, v, _abstraction * dt)
                    # [m² s⁻¹] = [m³ s⁻¹] / [m]
                    _inflow = -_abstraction / flow_length[v]
                else
                    # [m² s⁻¹] = [m³ s⁻¹] / [m]
                    _inflow = external_inflow[v] / flow_length[v]
                end
                # internal abstraction (water demand) is limited by river storage and
                # negative external inflow as part of water allocation computations.
                # [m² s⁻¹] = [m³ s⁻¹] / [m]
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
                        # [?] = min([m³ s⁻¹], [m³] / [s])
                        _abstraction = min(
                            -res_bc.external_inflow[i],
                            (reservoir.variables.storage[i] / dt) * 0.98,
                        )

                        res_bc.actual_external_abstraction_av[i] += _abstraction
                        _inflow = -_abstraction
                    else
                        _inflow = res_bc.external_inflow[i]
                    end
                    net_inflow =
                        q[v] +
                        res_bc.inflow_overland[i] +
                        res_bc.inflow_subsurface[i] +
                        _inflow
                    update!(reservoir, i, net_inflow, dt)

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
                # [m²] = [s^3/5 m^1/5] * [m³ s⁻¹] ^ (3/5)
                crossarea = alpha[v] * pow(q[v], beta)
                # [m] = [m²] / [m]
                h[v] = crossarea / flow_width[v]
                # [m³] = [m] * [m] * [m]
                storage[v] = flow_length[v] * flow_width[v] * h[v]

                # average variables
                add_to_cumulative!(q_av, v, q[v])
                add_to_cumulative!(qin_av, v, qin[v])
            end
        end
    end
end

"""
Update river flow model `KinWaveRiverFlow` for a single timestep `dt`. Timestepping within
`dt` is either with a fixed timestep `dt_fixed` or adaptive.
"""
function update!(model::KinWaveRiverFlow, domain::Domain, clock::Clock, dt::Number)
    (; reservoir, inwater) = model.boundary_conditions
    (; alpha_term, mannings_n, beta, alpha_pow, alpha, bankfull_depth) = model.parameters
    (; slope, flow_width, flow_length) = domain.river.parameters
    (; qlat, qin_av) = model.variables
    (; adaptive) = model.timestepping

    # [s³ᐟ⁵ m⁻¹ᐟ⁵] = ([s m⁻¹ᐟ³] / [-])³ᐟ⁵
    @. alpha_term = pow(mannings_n / sqrt(slope), beta)
    # use fixed alpha value based on 0.5 * bankfull_depth
    # [s³ᐟ⁵ m¹ᐟ⁵] =  [s³ᐟ⁵ m⁻¹ᐟ⁵] * ([m] + [m])²ᐟ⁵
    @. alpha = alpha_term * pow(flow_width + bankfull_depth, alpha_pow)
    # [m² s⁻¹] = [m³ s⁻¹] / [m]
    @. qlat = inwater / flow_length

    set_flow_vars!(model)
    zero!(qin_av)
    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)

    t = 0.0
    while t < dt
        dt_s =
            adaptive ? stable_timestep(model, flow_length, 0.05) :
            model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        kinwave_river_update!(model, domain.river, dt_s)
        t += dt_s
    end

    average_reservoir_vars!(reservoir, dt)
    average_flow_vars!(model, dt)
    average!(qin_av, dt)
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
            # [s^-1/5 m^-7/5] = ([s^3/5 m^1/5] * [-] * [m³ s⁻¹] ^ 2/5)⁻¹ ??
            c = inv(alpha[i] * beta * pow(q[i], (beta - 1.0)))
            stable_timesteps[k] = (flow_length[i] / c)
        end
    end

    dt_min = if isone(k)
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

    nonirrigation_returnflow = get_nonirrigation_returnflow(allocation)
    flux_to_river = get_flux_to_river(subsurface_flow, land_indices)
    # [m³ s⁻¹] = [m³ s⁻¹] + [m³ s⁻¹] + [m s⁻¹] * [m²] + [m s⁻¹] * [m²]
    @. inwater = (
        flux_to_river +
        overland_flow.variables.to_river.average[land_indices] +
        net_runoff_river[land_indices] * area[land_indices] +
        nonirrigation_returnflow * cell_area
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

    nonirrigation_returnflow = get_nonirrigation_returnflow(allocation)
    # [m³ s⁻¹] = ([m s⁻¹] + [m s⁻¹]) * [m²]
    @. inwater = (net_runoff + nonirrigation_returnflow) * area

    if config.model.drain__flag
        (; drain) = subsurface_flow.boundaries
        for (i, index) in drain.index
            inwater[index] -= drain_variables.flux[i]
        end
    end
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
    model.variables.q_av.average[inds]
get_inflow_reservoir(::KinWaveRiverFlow, model::LateralSSF, inds::Vector{Int}) =
    model.variables.ssf[inds]

# Exclude subsurface flow from `GroundwaterFlow`.
get_inflow_reservoir(::AbstractRiverFlowModel, model::GroundwaterFlow, inds::Vector{Int}) =
    zeros(length(inds))
