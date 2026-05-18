"Struct for storing (shared) variables for river and overland flow models"
@with_kw struct FlowVariables
    n_cells::Int
    q::Vector{Float64} = zeros(n_cells)            # Discharge [m³ s⁻¹]
    qlat::Vector{Float64} = zeros(n_cells)         # Lateral inflow per unit length [m² s⁻¹]
    qin::Vector{Float64} = zeros(n_cells)          # Inflow from upstream cells [m³ s⁻¹]
    qin_av::Vector{Float64} = zeros(n_cells)       # Average inflow from upstream cells  [m³ s⁻¹] for model timestep Δt
    q_av::Vector{Float64} = zeros(n_cells)         # Average discharge [m³ s⁻¹] for model timestep Δt
    storage::Vector{Float64} = zeros(n_cells)      # Kinematic wave storage [m³] (based on water depth h)
    h::Vector{Float64} = zeros(n_cells)            # Water depth [m]
end

"Struct for storing Manning flow parameters"
@with_kw struct ManningFlowParameters
    n_river_cells::Int
    slope::Vector{Float64}        # Slope [m m⁻¹]
    mannings_n::Vector{Float64}   # Manning's roughness [s m⁻⅓]
    alpha_pow::Float64            # Used in the power part of alpha [-]
    alpha_term::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)   # Term used in computation of alpha [-]
    alpha::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)        # Constant in momentum equation A = alpha*Q^beta, based on Manning's equation [s3/5 m1/5]
end

"Initialize Manning flow parameters"
function ManningFlowParameters(mannings_n::Vector{Float64}, slope::Vector{Float64})
    n_river_cells = length(slope)
    parameters = ManningFlowParameters(;
        n_river_cells,
        slope,
        mannings_n,
        alpha_pow=Float64((2.0 / 3.0) * 0.6),
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
    (; river_indices_2d) = domain.network
    (; slope) = domain.parameters
    mannings_n = ncread(
        dataset,
        config,
        "river_water_flow__manning_n_parameter",
        Routing;
        sel=river_indices_2d,
    )
    bankfull_depth =
        ncread(dataset, config, "river_bank_water__depth", Routing; sel=river_indices_2d)

    flow_params = ManningFlowParameters(mannings_n, slope)
    parameters = RiverFlowParameters(; flow=flow_params, bankfull_depth)
    return parameters
end

"Struct for storing river flow model boundary conditions"
@with_kw struct RiverFlowBC{R<:Union{ReservoirModel,Nothing}}
    n_river_cells::Int
    inwater::Vector{Float64} = zeros(n_river_cells)                         # Lateral inflow [m³ s⁻¹]
    external_inflow::Vector{Float64} = zeros(n_river_cells)                 # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    actual_external_abstraction_av::Vector{Float64} = zeros(n_river_cells)  # Actual abstraction from external negative inflow [m³ s⁻¹]
    abstraction::Vector{Float64} = zeros(n_river_cells)                     # Abstraction (computed as part of water demand and allocation) [m³ s⁻¹]
    reservoir::R                                                # ReservoirModel model struct of arrays
end

"Initialize river flow model boundary conditions"
function RiverFlowBC(
    dataset::NCDataset,
    config::Config,
    network::NetworkRiver,
    reservoir::Union{ReservoirModel,Nothing},
)
    (; river_indices_2d) = network
    external_inflow = ncread(
        dataset,
        config,
        "river_water__external_inflow_volume_flow_rate",
        Routing;
        sel=river_indices_2d,
    )
    n_river_cells = length(river_indices_2d)
    bc = RiverFlowBC(; n_river_cells, external_inflow, reservoir)
    return bc
end

"River flow model using the kinematic wave method and the Manning flow equation"
@with_kw struct KinWaveRiverFlowModel{R<:RiverFlowBC,A<:AbstractAllocationModel} <:
                AbstractRiverFlowModel
    timestepping::TimeStepping
    boundary_conditions::R
    parameters::RiverFlowParameters
    variables::FlowVariables
    allocation::A   # Water allocation
end

"Initialize river flow model `KinWaveRiverFlowModel`"
function KinWaveRiverFlowModel(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    reservoir::Union{ReservoirModel,Nothing},
)
    (; river_indices_2d) = domain.network
    n_river_cells = length(river_indices_2d)

    timestepping = init_kinematic_wave_timestepping(config, n_river_cells; domain="river")

    allocation =
        do_water_demand(config) ? AllocationRiverModel(; n_river_cells) :
        NoAllocationRiverModel(n_river_cells)

    variables = FlowVariables(; n_cells=n_river_cells)
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
    n_land_cells::Int
    flow::FlowVariables = FlowVariables(; n_cells=n_land_cells)
    to_river::Vector{Float64} = zeros(n_land_cells) # Part of overland flow [m³ s⁻¹] that flows to the river
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
    n_land_cells::Int
    inwater::Vector{Float64} = zeros(n_land_cells) # Lateral inflow [m³ s⁻¹]
end

"Overland flow model using the kinematic wave method and the Manning flow{ equation"
@with_kw struct KinWaveOverlandFlowModel <: AbstractOverlandFlowModel
    timestepping::TimeStepping
    boundary_conditions::LandFlowBC
    parameters::ManningFlowParameters
    variables::OverLandFlowVariables
end

"Initialize Overland flow model `KinWaveOverlandFlowModel`"
function KinWaveOverlandFlowModel(dataset::NCDataset, config::Config, domain::DomainLand)
    (; land_indices_2d) = domain.network
    (; slope) = domain.parameters
    mannings_n = ncread(
        dataset,
        config,
        "land_surface_water_flow__manning_n_parameter",
        Routing;
        sel=land_indices_2d,
    )

    n_land_cells = length(land_indices_2d)
    timestepping = init_kinematic_wave_timestepping(config, n_land_cells; domain="land")

    variables = OverLandFlowVariables(; n_land_cells)
    parameters = ManningFlowParameters(mannings_n, slope)
    boundary_conditions = LandFlowBC(; n_land_cells)

    overland_flow =
        KinWaveOverlandFlowModel(; timestepping, boundary_conditions, variables, parameters)

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

"Update overland flow model `KinWaveOverlandFlowModel` for a single timestep"
function kinwave_land_update!(
    overland_flow_model::KinWaveOverlandFlowModel,
    domain::DomainLand,
    dt::Float64,
)
    (; order_of_subdomains, subdomain_global_order, order_subdomain, upstream_nodes) =
        domain.network

    (; alpha) = overland_flow_model.parameters
    (; h, q, q_av, storage, qin, qin_av, qlat, to_river) = overland_flow_model.variables
    (; surface_flow_width, flow_length, flow_fraction_to_river) = domain.parameters

    n_subdomain_sets = length(order_of_subdomains)
    qin .= 0.0
    for subdomain_set_idx in 1:n_subdomain_sets
        threaded_foreach(
            eachindex(order_of_subdomains[subdomain_set_idx]);
            basesize=1,
        ) do in_subdomain_set_idx
            subdomain_idx = order_of_subdomains[subdomain_set_idx][in_subdomain_set_idx]
            for (river_global_traversion_idx, land_cell_idx) in
                zip(subdomain_global_order[subdomain_idx], order_subdomain[subdomain_idx])
                # for a river cell without a reservoir part of the upstream surface flow
                # goes to the river (flow_fraction_to_river) and part goes to the surface
                # flow reservoir (1.0 - flow_fraction_to_river), upstream nodes with a
                # reservoir are excluded
                to_river[land_cell_idx] +=
                    sum_at(
                        land_cell_idx_other ->
                            q[land_cell_idx_other] *
                            flow_fraction_to_river[land_cell_idx_other],
                        upstream_nodes[river_global_traversion_idx],
                    ) * dt
                if surface_flow_width[land_cell_idx] > 0.0
                    qin[land_cell_idx] = sum_at(
                        land_cell_idx_other ->
                            q[land_cell_idx_other] *
                            (1.0 - flow_fraction_to_river[land_cell_idx_other]),
                        upstream_nodes[river_global_traversion_idx],
                    )
                end

                q[land_cell_idx], crossarea = kinematic_wave(
                    qin[land_cell_idx],
                    q[land_cell_idx],
                    qlat[land_cell_idx],
                    alpha[land_cell_idx],
                    dt,
                    flow_length[land_cell_idx],
                )

                # update h, only if flow width > 0.0
                if surface_flow_width[land_cell_idx] > 0.0
                    h[land_cell_idx] = crossarea / surface_flow_width[land_cell_idx]
                end
                storage[land_cell_idx] =
                    flow_length[land_cell_idx] *
                    surface_flow_width[land_cell_idx] *
                    h[land_cell_idx]

                # average flow (here accumulated for model timestep Δt)
                q_av[land_cell_idx] += q[land_cell_idx] * dt
                qin_av[land_cell_idx] += qin[land_cell_idx] * dt
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
    (; q_av, qlat, qin_av, to_river) = overland_flow_model.variables
    (; adaptive) = overland_flow_model.timestepping

    @. alpha_term = pow(mannings_n / sqrt(slope), BETA_KINWAVE)
    # use fixed alpha value based flow width
    @. alpha = alpha_term * pow(surface_flow_width, alpha_pow)
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

"Run reservoir model and copy reservoir outflow to inflow (qin) of downstream river cell"
function update_reservoir_model!(
    reservoir_model::ReservoirModel,
    river_flow_vars::FlowVariables,
    network::NetworkRiver,
    river_cell_idx::Int,
    dt::Float64,
    dt_forcing::Float64,
)
    (; boundary_conditions, variables) = reservoir_model
    (; storage, outflow) = variables
    (;
        external_inflow,
        actual_external_abstraction_av,
        inflow_overland,
        inflow_subsurface,
    ) = boundary_conditions
    (; q, qin) = river_flow_vars
    (; reservoir_indices, graph) = network

    reservoir_idx = reservoir_indices[river_cell_idx]
    iszero(reservoir_idx) && return nothing

    # If the external inflow is negative, the abstraction is limited
    inflow_ext = external_inflow[reservoir_idx]
    if inflow_ext < 0.0
        abstraction = min(-inflow_ext, (storage[reservoir_idx] / dt) * 0.98)
        actual_external_abstraction_av[reservoir_idx] += abstraction * dt
        inflow = -abstraction
    else
        inflow = inflow_ext
    end

    net_inflow =
        q[river_cell_idx] +
        inflow_overland[reservoir_idx] +
        inflow_subsurface[reservoir_idx] +
        inflow
    update_reservoir_model!(reservoir_model, reservoir_idx, net_inflow, dt, dt_forcing)

    downstream_nodes = outneighbors(graph, river_cell_idx)
    n_downstream = length(downstream_nodes)
    if n_downstream == 1
        downstream_cell_idx = only(downstream_nodes)
        qin[downstream_cell_idx] = outflow[reservoir_idx]
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
    dt_forcing::Float64,
)
    (;
        graph,
        order_of_subdomains,
        order_subdomain,
        subdomain_global_order,
        upstream_nodes,
        reservoir_indices,
    ) = domain.network

    (; reservoir, external_inflow, actual_external_abstraction_av, abstraction) =
        river_flow_model.boundary_conditions

    (; alpha) = river_flow_model.parameters
    (; flow_width, flow_length) = domain.parameters
    (; h, q, q_av, storage, qin, qin_av, qlat) = river_flow_model.variables

    n_subdomain_sets = length(order_of_subdomains)
    qin .= 0.0
    for subdomain_set_idx in 1:n_subdomain_sets
        threaded_foreach(
            eachindex(order_of_subdomains[subdomain_set_idx]);
            basesize=1,
        ) do in_subdomain_set_idx
            subdomain_idx = order_of_subdomains[subdomain_set_idx][in_subdomain_set_idx]
            for (river_global_traversion_idx, river_cell_idx) in
                zip(subdomain_global_order[subdomain_idx], order_subdomain[subdomain_idx])
                # qin by outflow from upstream reservoir location is added
                qin[river_cell_idx] +=
                    sum_at(q, upstream_nodes[river_global_traversion_idx])
                # Inflow supply/abstraction is added to qlat (divide by flow length)
                # If external_inflow < 0, abstraction is limited
                if external_inflow[river_cell_idx] < 0.0
                    _abstraction = min(
                        -external_inflow[river_cell_idx],
                        (storage[river_cell_idx] / dt) * 0.80,
                    )
                    actual_external_abstraction_av[river_cell_idx] += _abstraction * dt
                    _inflow = -_abstraction / flow_length[river_cell_idx]
                else
                    _inflow = external_inflow[river_cell_idx] / flow_length[river_cell_idx]
                end
                # internal abstraction (water demand) is limited by river storage and
                # negative external inflow as part of water allocation computations.
                _inflow -= abstraction[river_cell_idx] / flow_length[river_cell_idx]

                q[river_cell_idx], crossarea = kinematic_wave(
                    qin[river_cell_idx],
                    q[river_cell_idx],
                    qlat[river_cell_idx] + _inflow,
                    alpha[river_cell_idx],
                    dt,
                    flow_length[river_cell_idx],
                )

                if !isnothing(reservoir)
                    update_reservoir_model!(
                        reservoir,
                        river_flow_model.variables,
                        domain.network,
                        river_cell_idx,
                        dt,
                        dt_forcing,
                    )
                end
                # update h and storage
                h[river_cell_idx] = crossarea / flow_width[river_cell_idx]
                storage[river_cell_idx] =
                    flow_length[river_cell_idx] *
                    flow_width[river_cell_idx] *
                    h[river_cell_idx]

                # average variables (here accumulated for model timestep Δt)
                q_av[river_cell_idx] += q[river_cell_idx] * dt
                qin_av[river_cell_idx] += qin[river_cell_idx] * dt
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
)
    (; reservoir, inwater) = river_flow_model.boundary_conditions
    (; alpha_term, mannings_n, alpha_pow, alpha, bankfull_depth) =
        river_flow_model.parameters
    (; slope, flow_width, flow_length) = domain.river.parameters
    (; qlat, qin_av) = river_flow_model.variables
    (; adaptive) = river_flow_model.timestepping

    @. alpha_term = pow(mannings_n / sqrt(slope), BETA_KINWAVE)
    # use fixed alpha value based on 0.5 * bankfull_depth
    @. alpha = alpha_term * pow(flow_width + bankfull_depth, alpha_pow)
    @. qlat = inwater / flow_length

    set_flow_vars!(river_flow_model)
    qin_av .= 0.0
    set_reservoir_vars!(reservoir)
    update_index_hq!(reservoir, clock)

    dt = tosecond(clock.dt)
    t = 0.0
    while t < dt
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
) where {S<:Union{KinWaveOverlandFlowModel,KinWaveRiverFlowModel}}
    (; q) = flow_model.variables
    (; alpha) = flow_model.parameters
    (; stable_timesteps) = flow_model.timestepping

    n_cells = length(q)
    stable_timesteps .= Inf
    stable_timestep_idx = 0
    for cell_idx in 1:n_cells
        if q[cell_idx] > 0.0
            stable_timestep_idx += 1
            c = inv(alpha[cell_idx] * BETA_KINWAVE * pow(q[cell_idx], (BETA_KINWAVE - 1.0)))
            stable_timesteps[stable_timestep_idx] = (flow_length[cell_idx] / c)
        end
    end

    dt_min = if stable_timestep_idx == 1
        stable_timesteps[stable_timestep_idx]
    elseif stable_timestep_idx > 0
        quantile!(@view(stable_timesteps[1:stable_timestep_idx]), p)
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

    (; land_cell_indices_containing_river) = domain.river.network
    (; cell_area) = domain.river.parameters
    (; area) = domain.land.parameters

    inwater .= (
        get_flux_to_river(subsurface_flow, land_cell_indices_containing_river) .+
        overland_flow.variables.to_river[land_cell_indices_containing_river] .+
        (
            net_runoff_river[land_cell_indices_containing_river] .*
            area[land_cell_indices_containing_river] .* 0.001
        ) ./ dt .+
        (get_nonirrigation_returnflow(allocation) .* 0.001 .* cell_area) ./ dt
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
    (; land_cell_indices_containing_drainage) = domain.drain.network

    if config.model.drain__flag
        drain = subsurface_flow.boundary_conditions.drain
        drainflux = zeros(length(net_runoff))
        drainflux[land_cell_indices_containing_drainage] =
            -drain.variables.flux ./ tosecond(BASETIMESTEP)
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
    reservoir_model::Union{ReservoirModel,Nothing},
    river_flow_model::AbstractRiverFlowModel,
    external_models::NamedTuple,
    network::NetworkReservoir,
)
    (; overland_flow, subsurface_flow) = external_models
    (; land_cell_indices_containing_reservoir) = network
    if !isnothing(reservoir_model)
        (; inflow_overland, inflow_subsurface) = reservoir_model.boundary_conditions
        inflow_overland .= get_inflow_reservoir(
            river_flow_model,
            overland_flow,
            land_cell_indices_containing_reservoir,
        )
        inflow_subsurface .= get_inflow_reservoir(
            river_flow_model,
            subsurface_flow,
            land_cell_indices_containing_reservoir,
        )
    end
    return nothing
end

# For the river kinematic wave, the variable `to_river` can be excluded, because this part
# is added to the river kinematic wave.
get_inflow_reservoir(
    ::KinWaveRiverFlowModel,
    overland_flow_model::KinWaveOverlandFlowModel,
    inds::Vector{Int},
) = overland_flow_model.variables.q_av[inds]
get_inflow_reservoir(
    ::KinWaveRiverFlowModel,
    subsurface_flow_model::LateralSSFModel,
    inds::Vector{Int},
) = subsurface_flow_model.variables.q_av[inds] ./ tosecond(BASETIMESTEP)

# Exclude subsurface flow from `GroundwaterFlowModel`.
get_inflow_reservoir(::AbstractRiverFlowModel, ::GroundwaterFlowModel, inds::Vector{Int}) =
    zeros(length(inds))
