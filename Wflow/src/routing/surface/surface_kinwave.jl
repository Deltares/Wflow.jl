"Struct for storing (shared) variables for river and overland flow models"
@with_kw struct FlowVariables
    n::Int
    # Discharge [m³ s⁻¹]
    q::Vector{Float64} = zeros(n)
    # Lateral inflow per unit length [m² s⁻¹]
    qlat::Vector{Float64} = zeros(n)
    # Inflow from upstream cells [m³ s⁻¹]
    qin::Vector{Float64} = zeros(n)
    # Cumulative inflow from upstream cells [m³] for model time step dt
    qin_cumulative::Vector{Float64} = zeros(n)
    # Average inflow from upstream cells  [m³ s⁻¹] for model time step dt
    qin_average::Vector{Float64} = zeros(n)
    # Cumulative discharge [m³] for model timestep dt
    q_cumulative::Vector{Float64} = zeros(n)
    # Average discharge [m³ s⁻¹] for model time step dt
    q_average::Vector{Float64} = zeros(n)
    # Kinematic wave storage [m³] (based on water depth h)
    storage::Vector{Float64} = zeros(n)
    # Water depth [m]
    h::Vector{Float64} = zeros(n)
end

"Struct for storing Manning flow parameters"
@with_kw struct ManningFlowParameters
    # Slope [m m⁻¹]
    slope::Vector{Float64}
    # Manning's roughness [s m⁻⅓]
    mannings_n::Vector{Float64}
    # Used in the power part of alpha [-]
    alpha_pow::Float64 = (2 // 3) * 0.6
    # Term used in computation of alpha [s^3/5 m^-1/5]
    alpha_term::Vector{Float64} = fill(MISSING_VALUE, length(slope))
    # Constant in momentum equation A = alpha*Q^beta, based on Manning's equation [s³ᐟ⁵ m¹ᐟ⁵]
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
    )
    bankfull_depth =
        ncread(dataset, config, "river_bank_water__depth", Routing; sel = indices)

    flow_params = ManningFlowParameters(; mannings_n, slope)
    parameters = RiverFlowParameters(; flow = flow_params, bankfull_depth)
    return parameters
end

"Struct for storing river flow model boundary conditions"
@with_kw struct RiverFlowBC{R <: Union{ReservoirModel, Nothing}}
    n::Int
    # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    external_inflow::Vector{Float64}
    # Lateral inflow [m³ s⁻¹]
    inwater::Vector{Float64} = zeros(n)
    # Actual cumulative abstraction from external negative inflow [m³]
    actual_external_abstraction_cumulative::Vector{Float64} = zeros(n)
    # Actual average abstraction from external negative inflow [m³ s⁻¹]
    actual_external_abstraction_average::Vector{Float64} = zeros(n)
    # Abstraction (computed as part of water demand and allocation) [m³ s⁻¹]
    abstraction::Vector{Float64} = zeros(n)
    # Reservoir model struct of arrays
    reservoir::R
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
        "river_water__external_inflow_volume_flow_rate",
        Routing;
        sel = indices,
    )
    n = length(indices)
    bc = RiverFlowBC(; n, external_inflow, reservoir)
    return bc
end

"River flow model using the kinematic wave method and the Manning flow equation"
@with_kw struct KinWaveRiverFlowModel{R <: RiverFlowBC, A <: AbstractAllocationModel} <:
                AbstractRiverFlowModel
    timestepping::TimeStepping
    boundary_conditions::R
    parameters::RiverFlowParameters
    variables::FlowVariables
    # Water allocation
    allocation::A
end

"Initialize river flow model `KinWaveRiverFlowModel`"
function KinWaveRiverFlowModel(
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

    variables = FlowVariables(; n)
    parameters = RiverFlowParameters(dataset, config, domain)
    boundary_conditions = RiverFlowBC(dataset, config, domain.network, reservoir)

    river_flow = KinWaveRiverFlowModel(;
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
    # Part of cumulative overland flow [m³ s⁻¹] that flows to the river
    to_river_cumulative::Vector{Float64} = zeros(n)
    # Part of average overland flow [m³ s⁻¹] that flows to the river
    to_river_average::Vector{Float64} = zeros(n)
end

"Overload `getproperty` for overland flow model variables"
function Base.getproperty(v::OverLandFlowVariables, s::Symbol)
    if s ∈ (:to_river_average, :to_river_cumulative, :flow)
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
@with_kw struct KinWaveOverlandFlowModel <: AbstractOverlandFlowModel
    n::Int
    timestepping::TimeStepping
    boundary_conditions::LandFlowBC = LandFlowBC(; n)
    parameters::ManningFlowParameters
    variables::OverLandFlowVariables = OverLandFlowVariables(; n)
end

"Initialize Overland flow model `KinWaveOverlandFlowModel`"
function KinWaveOverlandFlowModel(dataset::NCDataset, config::Config, domain::DomainLand)
    (; indices) = domain.network
    (; slope) = domain.parameters
    mannings_n = ncread(
        dataset,
        config,
        "land_surface_water_flow__manning_n_parameter",
        Routing;
        sel = indices,
    )

    n = length(indices)
    timestepping = init_kinematic_wave_timestepping(config, n; domain = "land")
    parameters = ManningFlowParameters(; mannings_n, slope)
    overland_flow = KinWaveOverlandFlowModel(; n, timestepping, parameters)

    return overland_flow
end

"""
Helper function to set reservoir variables to zero. This is done at the start of each
simulation timestep, during the timestep the total (weighted) sum is computed from values at
each internal timestep.
"""
function set_reservoir_vars!(reservoir_model::ReservoirModel)
    (; boundary_conditions, variables) = reservoir_model
    (; inflow_cumulative, actual_external_abstraction_cumulative) = boundary_conditions
    (; outflow_cumulative, actevap_cumulative) = variables

    inflow_cumulative .= 0.0
    actual_external_abstraction_cumulative .= 0.0
    outflow_cumulative .= 0.0
    actevap_cumulative .= 0.0
    return nothing
end
set_reservoir_vars!(reservoir_model::Any) = nothing

"""
Helper function to compute the average of reservoir variables. This is done at the end of
each simulation timestep.
"""
function average_reservoir_vars!(reservoir_model::ReservoirModel, dt::Float64)
    (; variables, boundary_conditions) = reservoir_model
    (; outflow_average, outflow_cumulative) = variables
    (;
        inflow_average,
        inflow_cumulative,
        actual_external_abstraction_average,
        actual_external_abstraction_cumulative,
    ) = boundary_conditions

    @. outflow_average = outflow_cumulative / dt
    @. inflow_average = inflow_cumulative / dt
    @. actual_external_abstraction_average = actual_external_abstraction_cumulative / dt
    return nothing
end
average_reservoir_vars!(reservoir::Any, dt::Float64) = nothing

"""
    set_flow_vars!(river_flow_model::AbstractRiverFlowModel)

Helper functions to set cumulative river flow routing variables discharge and actual abstraction (based
on external negative inflow) from river to zero. This is done at the start of each
simulation timestep, during the timestep the total (weighted) sum is computed from values at
each sub timestep.
"""
function set_flow_vars!(river_flow_model::AbstractRiverFlowModel)
    (; q_cumulative) = river_flow_model.variables
    (; actual_external_abstraction_cumulative) = river_flow_model.boundary_conditions
    q_cumulative .= 0.0
    actual_external_abstraction_cumulative .= 0.0
    return nothing
end

"""
    average_flow_vars!(river_flow_model::AbstractRiverFlowModel, dt::Float64)

Helper functions to compute average river flow routing variables. This is done at the end of
each simulation timestep.
"""
function average_flow_vars!(river_flow_model::AbstractRiverFlowModel, dt::Float64)
    (; q_average, q_cumulative) = river_flow_model.variables
    (; actual_external_abstraction_average, actual_external_abstraction_cumulative) =
        river_flow_model.boundary_conditions
    @. q_average = q_cumulative / dt
    @. actual_external_abstraction_average = actual_external_abstraction_cumulative / dt
    return nothing
end

"Update overland flow model `KinWaveOverlandFlowModel` for a single timestep"
function kinwave_land_update!(
    overland_flow_model::KinWaveOverlandFlowModel,
    domain::DomainLand,
    dt::Float64,
)
    (; order_of_subdomains, order_subdomain, subdomain_indices, upstream_nodes) =
        domain.network

    (; alpha) = overland_flow_model.parameters
    (; h, q, q_cumulative, storage, qin, qin_cumulative, qlat, to_river_cumulative) =
        overland_flow_model.variables
    (; surface_flow_width, flow_length, flow_fraction_to_river) = domain.parameters

    n_subdomain_levels = length(order_of_subdomains)
    qin .= 0.0
    for level_idx in 1:n_subdomain_levels
        threaded_foreach(
            eachindex(order_of_subdomains[level_idx]);
            basesize = 1,
        ) do subdomain_idx
            m = order_of_subdomains[level_idx][subdomain_idx]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                # for a river cell without a reservoir part of the upstream surface flow
                # goes to the river (flow_fraction_to_river) and part goes to the surface
                # flow reservoir (1.0 - flow_fraction_to_river), upstream nodes with a
                # reservoir are excluded
                to_river_cumulative[v] +=
                    sum_at(
                        upstream_idx ->
                            q[upstream_idx] * flow_fraction_to_river[upstream_idx],
                        upstream_nodes[n],
                    ) * dt

                if surface_flow_width[v] > 0.0
                    qin[v] = sum_at(
                        upstream_idx ->
                            q[upstream_idx] * (1.0 - flow_fraction_to_river[upstream_idx]),
                        upstream_nodes[n],
                    )
                end

                q[v], crossarea =
                    kinematic_wave(qin[v], q[v], qlat[v], alpha[v], dt, flow_length[v])

                # update h, only if flow width > 0.0
                if surface_flow_width[v] > 0.0
                    h[v] = crossarea / surface_flow_width[v]
                end
                storage[v] = flow_length[v] * surface_flow_width[v] * h[v]

                # average flow
                q_cumulative[v] += q[v] * dt
                qin_cumulative[v] += qin[v] * dt
            end
        end
    end
end

"""
Update overland flow model `KinWaveOverlandFlowModel` for a single timestep `dt`. Timestepping within
`dt` is either with a fixed timestep `dt_fixed` or adaptive.
"""
function update_overland_flow_model!(
    overland_flow_model::KinWaveOverlandFlowModel,
    domain::DomainLand,
    dt::Float64,
)
    (; inwater) = overland_flow_model.boundary_conditions
    (; alpha_term, mannings_n, alpha_pow, alpha) = overland_flow_model.parameters
    (; surface_flow_width, flow_length, slope) = domain.parameters
    (;
        q_average,
        q_cumulative,
        qlat,
        qin_average,
        qin_cumulative,
        to_river_average,
        to_river_cumulative,
    ) = overland_flow_model.variables
    (; adaptive) = overland_flow_model.timestepping

    @. alpha_term = pow(mannings_n / sqrt(slope), BETA_KINWAVE)
    # use fixed alpha value based flow width
    @. alpha = alpha_term * pow(surface_flow_width, alpha_pow)
    @. qlat = inwater / flow_length

    q_cumulative .= 0.0
    qin_cumulative .= 0.0
    to_river_cumulative .= 0.0

    t = 0.0
    while t < dt
        dt_s =
            adaptive ? stable_timestep(overland_flow_model, flow_length, 0.02) :
            overland_flow_model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        kinwave_land_update!(overland_flow_model, domain, dt_s)
        t += dt_s
    end

    @. q_average = q_cumulative / dt
    @. to_river_average = to_river_cumulative / dt
    @. qin_average = qin_cumulative / dt
    return nothing
end

"""
run reservoir model and copy reservoir outflow to inflow (qin) of
downstream river cell
"""
function update_reservoir_model!(
    reservoir_model::ReservoirModel,
    river_flow_vars::FlowVariables,
    network::NetworkRiver,
    v::Int,
    dt::Float64,
)
    (; boundary_conditions, variables) = reservoir_model
    (; storage, outflow) = variables
    (;
        external_inflow,
        actual_external_abstraction_cumulative,
        inflow_overland,
        inflow_subsurface,
    ) = boundary_conditions
    (; q, qin) = river_flow_vars
    (; reservoir_indices, graph) = network

    i = reservoir_indices[v]
    iszero(i) && return nothing

    # If the external inflow is negative, the abstraction is limited
    inflow_ext = external_inflow[i]
    if inflow_ext < 0.0
        abstraction = min(-inflow_ext, (storage[i] / dt) * 0.98)
        actual_external_abstraction_cumulative[i] += abstraction * dt
        inflow = -abstraction
    else
        inflow = inflow_ext
    end
    net_inflow = q[v] + inflow_overland[i] + inflow_subsurface[i] + inflow
    update_reservoir_model!(reservoir_model, i, net_inflow, dt)

    downstream_nodes = outneighbors(graph, v)
    n_downstream = length(downstream_nodes)
    if n_downstream == 1
        j = only(downstream_nodes)
        qin[j] = outflow[i]
    elseif n_downstream == 0
        error(
            """A reservoir without a downstream river node is not supported.
            Add a downstream river node or move the reservoir to an upstream node (model schematization).
            """,
        )
    else
        error("bifurcations not supported")
    end
    return nothing
end

"Update river flow model `KinWaveRiverFlowModel` for a single timestep"
function kinwave_river_update!(
    river_flow_model::KinWaveRiverFlowModel,
    domain::DomainRiver,
    dt::Float64,
)
    (; order_of_subdomains, order_subdomain, subdomain_indices, upstream_nodes) =
        domain.network

    (; reservoir, external_inflow, actual_external_abstraction_cumulative, abstraction) =
        river_flow_model.boundary_conditions

    (; alpha) = river_flow_model.parameters
    (; flow_width, flow_length) = domain.parameters
    (; h, q, q_cumulative, storage, qin, qin_cumulative, qlat) = river_flow_model.variables

    n_subdomain_levels = length(order_of_subdomains)
    qin .= 0.0
    for level_idx in 1:n_subdomain_levels
        threaded_foreach(
            eachindex(order_of_subdomains[level_idx]);
            basesize = 1,
        ) do subdomain_idx
            m = order_of_subdomains[level_idx][subdomain_idx]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                # qin by outflow from upstream reservoir location is added
                qin[v] += sum_at(q, upstream_nodes[n])
                # Inflow supply/abstraction is added to qlat (divide by flow length)
                # If external_inflow < 0, abstraction is limited
                if external_inflow[v] < 0.0
                    _abstraction = min(-external_inflow[v], (storage[v] / dt) * 0.80)
                    actual_external_abstraction_cumulative[v] += _abstraction * dt
                    _inflow = -_abstraction / flow_length[v]
                else
                    _inflow = external_inflow[v] / flow_length[v]
                end
                # internal abstraction (water demand) is limited by river storage and
                # negative external inflow as part of water allocation computations.
                _inflow -= abstraction[v] / flow_length[v]

                q[v], crossarea = kinematic_wave(
                    qin[v],
                    q[v],
                    qlat[v] + _inflow,
                    alpha[v],
                    dt,
                    flow_length[v],
                )

                if !isnothing(reservoir)
                    update_reservoir_model!(
                        reservoir,
                        river_flow_model.variables,
                        domain.network,
                        v,
                        dt,
                    )
                end
                # update h and storage
                h[v] = crossarea / flow_width[v]
                storage[v] = flow_length[v] * flow_width[v] * h[v]

                # average variables
                q_cumulative[v] += q[v] * dt
                qin_cumulative[v] += qin[v] * dt
            end
        end
    end
end

"""
Update river flow model `KinWaveRiverFlowModel` for a single timestep `dt`. Timestepping within
`dt` is either with a fixed timestep `dt_fixed` or adaptive.
"""
function update_river_flow_model!(
    river_flow_model::KinWaveRiverFlowModel,
    domain::Domain,
    clock::Clock,
    dt::Float64,
)
    (; reservoir, inwater) = river_flow_model.boundary_conditions
    (; alpha_term, mannings_n, alpha_pow, alpha, bankfull_depth) =
        river_flow_model.parameters
    (; slope, flow_width, flow_length) = domain.river.parameters
    (; qlat, qin_average, qin_cumulative) = river_flow_model.variables
    (; adaptive) = river_flow_model.timestepping

    @. alpha_term = pow(mannings_n / sqrt(slope), BETA_KINWAVE)
    # use fixed alpha value based on 0.5 * bankfull_depth
    @. alpha = alpha_term * pow(flow_width + bankfull_depth, alpha_pow)
    @. qlat = inwater / flow_length

    set_flow_vars!(river_flow_model)
    qin_cumulative .= 0.0
    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)

    t = 0.0
    while t < dt
        dt_s =
            adaptive ? stable_timestep(river_flow_model, flow_length, 0.05) :
            river_flow_model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        kinwave_river_update!(river_flow_model, domain.river, dt_s)
        t += dt_s
    end

    average_reservoir_vars!(reservoir, dt)
    average_flow_vars!(river_flow_model, dt)
    @. qin_average = qin_cumulative / dt

    return nothing
end

# constant in Manning's equation [-]
const BETA_KINWAVE = 0.6

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
) where {S <: Union{KinWaveOverlandFlowModel, KinWaveRiverFlowModel}}
    (; q) = flow_model.variables
    (; alpha) = flow_model.parameters
    (; stable_timesteps) = flow_model.timestepping

    n_cells = length(q)
    stable_timesteps .= Inf
    k = 0
    for cell_idx in 1:n_cells
        if q[cell_idx] > KIN_WAVE_MIN_FLOW
            k += 1
            c = inv(alpha[cell_idx] * BETA_KINWAVE * pow(q[cell_idx], (BETA_KINWAVE - 1.0)))
            stable_timesteps[k] = (flow_length[cell_idx] / c)
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

    nonirrigation_returnflow = get_nonirrigation_returnflow(allocation)
    flux_subsurface_to_river = get_flux_to_river(subsurface_flow, land_indices)
    flux_overland_to_river = overland_flow.variables.to_river_average
    @. inwater = (
        flux_subsurface_to_river +
        flux_overland_to_river[land_indices] +
        net_runoff_river[land_indices] * area[land_indices] +
        nonirrigation_returnflow * cell_area
    )
    return nothing
end

"""
Update boundary condition lateral inflow `inwater` of a kinematic wave overland flow model
`KinWaveOverlandFlowModel` for a single timestep.
"""
function update_lateral_inflow!(
    overland_flow_model::KinWaveOverlandFlowModel,
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
        drainflux[land_indices] = -drain.variables.flux
    else
        drainflux = 0.0
    end

    nonirrigation_returnflow = get_nonirrigation_returnflow(allocation)
    @. inwater = (net_runoff + nonirrigation_returnflow) * area + drainflux

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
    river_flow_model::KinWaveRiverFlowModel,
    overland_flow_model::KinWaveOverlandFlowModel,
    inds::Vector{Int},
) = overland_flow_model.variables.q_average[inds]
get_inflow_reservoir(
    river_flow_model::KinWaveRiverFlowModel,
    subsurface_flow_model::LateralSSFModel,
    inds::Vector{Int},
) = subsurface_flow_model.variables.q_average[inds]

# Exclude subsurface flow from `GroundwaterFlowModel`.
get_inflow_reservoir(::AbstractRiverFlowModel, ::GroundwaterFlowModel, inds::Vector{Int}) =
    zeros(length(inds))
