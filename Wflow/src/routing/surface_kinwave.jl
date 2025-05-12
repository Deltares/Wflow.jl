"Struct for storing (shared) variables for river and overland flow models"
@with_kw struct FlowVariables
    q::Vector{Float}            # Discharge [m³ s⁻¹]
    qlat::Vector{Float}         # Lateral inflow per unit length [m² s⁻¹]
    qin::Vector{Float}          # Inflow from upstream cells [m³ s⁻¹]
    q_av::Vector{Float}         # Average discharge [m³ s⁻¹] for model timestep Δt
    storage::Vector{Float}      # Kinematic wave storage [m³] (based on water depth h)
    storage_av::Vector{Float}   # Average kinematic wave storage [m³] for model timestep Δt
    h::Vector{Float}            # Water depth [m]
    h_av::Vector{Float}         # Average water depth [m] for model timestep Δt
end

"Initialize timestepping for kinematic wave (river and overland flow models)"
function init_kinematic_wave_timestepping(
    config::Config,
    n::Int;
    domain::String,
    dt_fixed::Float,
)
    adaptive = get(config.model, "kinematic_wave__adaptive_time_step_flag", false)::Bool
    @info "Kinematic wave approach is used for $domain flow, adaptive timestepping = $adaptive."

    if adaptive
        stable_timesteps = zeros(Float, n)
        timestepping = TimeStepping(; stable_timesteps, adaptive)
    else
        dt_fixed = get(config.model, "$(domain)_kinematic_wave__time_step", dt_fixed)
        @info "Using a fixed internal timestep (seconds) $dt_fixed for kinematic wave $domain flow."
        timestepping = TimeStepping(; dt_fixed, adaptive)
    end
    return timestepping
end

"Initialize variables for river or overland flow models"
function FlowVariables(n::Int)
    variables = FlowVariables(;
        q = zeros(Float, n),
        qlat = zeros(Float, n),
        qin = zeros(Float, n),
        q_av = zeros(Float, n),
        storage = zeros(Float, n),
        storage_av = zeros(Float, n),
        h = zeros(Float, n),
        h_av = zeros(Float, n),
    )
    return variables
end

"Struct for storing Manning flow parameters"
@with_kw struct ManningFlowParameters
    beta::Float                 # constant in Manning's equation [-]
    slope::Vector{Float}        # Slope [m m⁻¹]
    mannings_n::Vector{Float}   # Manning's roughness [s m⁻⅓]
    alpha_pow::Float            # Used in the power part of alpha [-]
    alpha_term::Vector{Float}   # Term used in computation of alpha [-]
    alpha::Vector{Float}        # Constant in momentum equation A = alpha*Q^beta, based on Manning's equation [s3/5 m1/5]
end

"Initialize Manning flow parameters"
function ManningFlowParameters(mannings_n::Vector{Float}, slope::Vector{Float})
    n = length(slope)
    parameters = ManningFlowParameters(;
        beta = Float(0.6),
        slope,
        mannings_n,
        alpha_pow = Float((2.0 / 3.0) * 0.6),
        alpha_term = fill(MISSING_VALUE, n),
        alpha = fill(MISSING_VALUE, n),
    )
    return parameters
end

"Struct for storing river flow model parameters"
@with_kw struct RiverFlowParameters
    flow::ManningFlowParameters
    bankfull_depth::Vector{Float} # Bankfull water level [m]
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

    lens = lens_input_parameter(config, "river_water_flow__manning_n_parameter")
    mannings_n =
        ncread(dataset, config, lens; sel = indices, defaults = 0.036, type = Float)
    lens = lens_input_parameter(config, "river_bank_water__depth")
    bankfull_depth =
        ncread(dataset, config, lens; sel = indices, defaults = 1.0, type = Float)

    flow_params = ManningFlowParameters(mannings_n, slope)
    parameters = RiverFlowParameters(; flow = flow_params, bankfull_depth)
    return parameters
end

"Struct for storing river flow model boundary conditions"
@with_kw struct RiverFlowBC{T <: DenseArray{Float}, R, L}
    inwater::T                          # Lateral inflow [m³ s⁻¹]
    inflow::T                           # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    inflow_waterbody::T                 # inflow waterbody (lake or reservoir model) from land part [m³ s⁻¹]
    abstraction::T                      # Abstraction (computed as part of water demand and allocation) [m³ s⁻¹]
    reservoir::R                        # Reservoir model struct of arrays
    lake::L                             # Lake model struct of arrays
end

function Adapt.adapt_structure(to, from::RiverFlowBC)
    return RiverFlowBC(
        adapt(to, from.inwater),
        adapt(to, from.inflow),
        adapt(to, from.inflow_waterbody),
        adapt(to, from.abstraction),
        adapt(to, from.reservoir),
        adapt(to, from.lake),
    )
end

"Initialize river flow model boundary conditions"
function RiverFlowBC(
    n::Int,
    reservoir::Union{SimpleReservoir, Nothing},
    lake::Union{Lake, Nothing},
)
    bc = RiverFlowBC(;
        inwater = zeros(Float, n),
        inflow = zeros(Float, n),
        inflow_waterbody = zeros(Float, n),
        abstraction = zeros(Float, n),
        reservoir = reservoir,
        lake = lake,
    )
    return bc
end

"River flow model using the kinematic wave method and the Manning flow equation"
@with_kw struct KinWaveRiverFlow{R, L, A} <: AbstractRiverFlowModel
    timestepping::TimeStepping
    boundary_conditions::RiverFlowBC{R, L}
    parameters::RiverFlowParameters
    variables::FlowVariables
    allocation::A   # Water allocation
end

"Initialize river flow model `KinWaveRiverFlow`"
function KinWaveRiverFlow(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    reservoir::Union{SimpleReservoir, Nothing},
    lake::Union{Lake, Nothing},
)
    (; indices) = domain.network
    n = length(indices)

    timestepping = init_kinematic_wave_timestepping(
        config,
        n;
        domain = "river",
        dt_fixed = Float(900.0),
    )

    do_water_demand = haskey(config.model, "water_demand")
    allocation = do_water_demand ? AllocationRiver(n) : NoAllocationRiver()

    variables = FlowVariables(n)
    parameters = RiverFlowParameters(dataset, config, domain)
    boundary_conditions = RiverFlowBC(n, reservoir, lake)

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
    flow::FlowVariables
    to_river::Vector{Float} # Part of overland flow [m³ s⁻¹] that flows to the river
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
    inwater::Vector{Float} # Lateral inflow [m³ s⁻¹]
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
    lens = lens_input_parameter(config, "land_surface_water_flow__manning_n_parameter")
    mannings_n =
        ncread(dataset, config, lens; sel = indices, defaults = 0.072, type = Float)

    n = length(indices)
    timestepping = init_kinematic_wave_timestepping(
        config,
        n;
        domain = "land",
        dt_fixed = Float(3600.0),
    )

    variables = OverLandFlowVariables(; flow = FlowVariables(n), to_river = zeros(Float, n))
    parameters = ManningFlowParameters(mannings_n, slope)
    boundary_conditions = LandFlowBC(; inwater = zeros(Float, n))

    overland_flow =
        KinWaveOverlandFlow(; timestepping, boundary_conditions, variables, parameters)

    return overland_flow
end

"""
Helper function to set waterbody variables `inflow`,`outflow_av` and `actevap` to zero. This
is done at the start of each simulation timestep, during the timestep the total (weighted)
sum is computed from values at each sub timestep.
"""
function set_waterbody_vars!(waterbody::W) where {W <: Union{SimpleReservoir, Lake}}
    waterbody.boundary_conditions.inflow .= 0.0
    waterbody.variables.outflow_av .= 0.0
    waterbody.variables.storage_av .= 0.0
    waterbody.variables.actevap .= 0.0
    if isa(waterbody, Lake)
        waterbody.variables.waterlevel_av .= 0.0
    end
    return nothing
end
set_waterbody_vars!(waterbody) = nothing

"""
Helper function to compute the average of waterbody variables inflow, storage, outflow and
water level. This is done at the end of each simulation timestep.
"""
function average_waterbody_vars!(
    waterbody::W,
    dt::Float,
) where {W <: Union{SimpleReservoir, Lake}}
    waterbody.variables.outflow_av ./= dt
    waterbody.variables.storage_av ./= dt
    waterbody.boundary_conditions.inflow ./= dt
    if isa(waterbody, Lake)
        waterbody.variables.waterlevel_av ./= dt
    end
    return nothing
end
average_waterbody_vars!(waterbody, dt) = nothing

"""
Helper function to set flow routing variables discharge, water depth and storage to zero.
This is done at the start of each simulation timestep, during the timestep the total
(weighted) sum is computed from values at each sub timestep.
"""
function set_flow_vars!(variables)
    (; q_av, h_av, storage_av) = variables
    q_av .= 0.0
    h_av .= 0.0
    storage_av .= 0.0
    return nothing
end

"""
Helper function to compute average flow routing variables. This is done at the end of each
simulation timestep.
"""
function average_flow_vars!(variables, dt::Float)
    (; q_av, h_av, storage_av) = variables
    q_av ./= dt
    h_av ./= dt
    storage_av ./= dt
    return nothing
end

"Update overland flow model `KinWaveOverlandFlow` for a single timestep"
function kinwave_land_update!(model::KinWaveOverlandFlow, domain::DomainLand, dt::Float)
    (; order_of_subdomains, order_subdomain, subdomain_indices, upstream_nodes) =
        domain.network

    (; beta, alpha) = model.parameters
    (; h, h_av, q, q_av, storage, storage_av, qin, qlat, to_river) = model.variables
    (; surface_flow_width, flow_length, flow_fraction_to_river) = domain.parameters

    ns = length(order_of_subdomains)
    qin .= 0.0
    for k in 1:ns
        threaded_foreach(eachindex(order_of_subdomains[k]); basesize = 1) do i
            m = order_of_subdomains[k][i]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                # for a river cell without a reservoir or lake part of the upstream surface
                # flow goes to the river (flow_fraction_to_river) and part goes to the
                # surface flow reservoir (1.0 - flow_fraction_to_river), upstream nodes with
                # a reservoir or lake are excluded
                to_river[v] +=
                    sum_at(
                        i -> q[i] * flow_fraction_to_river[i],
                        upstream_nodes[n],
                        eltype(to_river),
                    ) * dt
                if surface_flow_width[v] > 0.0
                    qin[v] = sum_at(
                        i -> q[i] * (1.0 - flow_fraction_to_river[i]),
                        upstream_nodes[n],
                        eltype(q),
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

                # average variables (here accumulated for model timestep Δt)
                storage_av[v] += storage[v] * dt
                h_av[v] += h[v] * dt
                q_av[v] += q[v] * dt
            end
        end
    end
end

"""
Update overland flow model `KinWaveOverlandFlow` for a single timestep `dt`. Timestepping within
`dt` is either with a fixed timestep `dt_fixed` or adaptive.
"""
function update!(model::KinWaveOverlandFlow, domain::DomainLand, dt::Float)
    (; inwater) = model.boundary_conditions
    (; alpha_term, mannings_n, beta, alpha_pow, alpha) = model.parameters
    (; surface_flow_width, flow_length, slope) = domain.parameters
    (; qlat, to_river) = model.variables
    (; adaptive) = model.timestepping

    @. alpha_term = pow(mannings_n / sqrt(slope), beta)
    # use fixed alpha value based flow width
    @. alpha = alpha_term * pow(surface_flow_width, alpha_pow)
    @. qlat = inwater / flow_length

    set_flow_vars!(model.variables)
    to_river .= 0.0

    t = 0.0
    while t < dt
        dt_s =
            adaptive ? stable_timestep(model, flow_length, Float(0.02)) :
            model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        kinwave_land_update!(model, domain, dt_s)
        t = t + dt_s
    end
    average_flow_vars!(model.variables, dt)
    to_river ./= dt
    return nothing
end

"Update river flow model `KinWaveRiverFlow` for a single timestep"
function kinwave_river_update!(
    model::KinWaveRiverFlow,
    domain::DomainRiver,
    doy::Int,
    dt::Float,
    dt_forcing::Float,
)
    (;
        graph,
        order_of_subdomains,
        order_subdomain,
        subdomain_indices,
        upstream_nodes,
        reservoir_indices,
        lake_indices,
    ) = domain.network

    (; reservoir, lake, inwater, inflow, abstraction, inflow_waterbody) =
        model.boundary_conditions

    (; beta, alpha) = model.parameters
    (; flow_width, flow_length) = domain.parameters
    (; h, h_av, q, q_av, storage, storage_av, qin, qlat) = model.variables

    ns = length(order_of_subdomains)
    qin .= 0.0
    for k in 1:ns
        threaded_foreach(eachindex(order_of_subdomains[k]); basesize = 1) do i
            m = order_of_subdomains[k][i]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                # qin by outflow from upstream reservoir or lake location is added
                qin[v] += sum_at(q, upstream_nodes[n])
                # Inflow supply/abstraction is added to qlat (divide by flow length)
                # If inflow < 0, abstraction is limited
                if inflow[v] < 0.0
                    _inflow =
                        max(-((inwater[v] + qin[v] + storage[v] / dt) * 0.80), inflow[v])
                    _inflow = _inflow / flow_length[v]
                else
                    _inflow = inflow[v] / flow_length[v]
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
                    update!(reservoir, i, q[v] + inflow_waterbody[v], dt, dt_forcing)

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

                elseif !isnothing(lake) && lake_indices[v] != 0
                    # run lake model and copy lake outflow to inflow (qin) of downstream river
                    # cell
                    i = lake_indices[v]
                    update!(lake, i, q[v] + inflow_waterbody[v], doy, dt, dt_forcing)

                    downstream_nodes = outneighbors(graph, v)
                    n_downstream = length(downstream_nodes)
                    if n_downstream == 1
                        j = only(downstream_nodes)
                        qin[j] = max(lake.variables.outflow[i], 0.0)
                    elseif n_downstream == 0
                        error(
                            """A lake without a downstream river node is not supported.
                            Add a downstream river node or move the lake to an upstream node (model schematization).
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
                storage_av[v] += storage[v] * dt
                h_av[v] += h[v] * dt
                q_av[v] += q[v] * dt
            end
        end
    end
end

"""
Update river flow model `KinWaveRiverFlow` for a single timestep `dt`. Timestepping within
`dt` is either with a fixed timestep `dt_fixed` or adaptive.
"""
function update!(model::KinWaveRiverFlow, domain::Domain, doy::Int, dt::Float)
    (; reservoir, lake, inwater) = model.boundary_conditions
    (; alpha_term, mannings_n, beta, alpha_pow, alpha, bankfull_depth) = model.parameters
    (; slope, flow_width, flow_length) = domain.river.parameters
    (; qlat) = model.variables
    (; adaptive) = model.timestepping

    @. alpha_term = pow(mannings_n / sqrt(slope), beta)
    # use fixed alpha value based on 0.5 * bankfull_depth
    @. alpha = alpha_term * pow(flow_width + bankfull_depth, alpha_pow)
    @. qlat = inwater / flow_length

    set_flow_vars!(model.variables)
    set_waterbody_vars!(reservoir)
    set_waterbody_vars!(lake)

    t = 0.0
    while t < dt
        dt_s =
            adaptive ? stable_timestep(model, flow_length, Float(0.05)) :
            model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        kinwave_river_update!(model, domain.river, doy, dt_s, dt)
        t = t + dt_s
    end

    average_waterbody_vars!(reservoir, dt)
    average_waterbody_vars!(lake, dt)
    average_flow_vars!(model.variables, dt)
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
    flow_length::Vector{Float},
    p::Float,
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
    dt::Float,
)
    (; allocation, runoff, overland_flow, subsurface_flow) = external_models
    (; inwater) = model.boundary_conditions
    (; net_runoff_river) = runoff.variables

    (; land_indices) = domain.river.network
    (; cell_area) = domain.river.parameters
    (; area) = domain.land.parameters

    inwater .= (
        get_flux_to_river(subsurface_flow)[land_indices] .+
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
    area::Vector{Float},
    config::Config,
    dt::Float,
)
    (; soil, subsurface_flow, allocation) = external_models
    (; net_runoff) = soil.variables
    (; inwater) = model.boundary_conditions

    do_drains = get(config.model, "drain__flag", false)::Bool
    if do_drains
        drain = subsurface_flow.boundaries.drain
        drainflux = zeros(Float, length(net_runoff))
        drainflux[drain.index] = -drain.variables.flux ./ Float(tosecond(BASETIMESTEP))
    else
        drainflux = 0.0
    end
    inwater .=
        (net_runoff .+ get_nonirrigation_returnflow(allocation)) .* area * 0.001 ./ dt .+
        drainflux

    return nothing
end

"""
Update boundary condition inflow to a waterbody from land `inflow_waterbody` of a model
`AbstractRiverFlowModel` for a single timestep.
"""
function update_inflow_waterbody!(
    model::AbstractRiverFlowModel,
    external_models::NamedTuple,
    river_indices::Vector{Int},
)
    (; overland_flow, subsurface_flow) = external_models
    (; reservoir, lake, inflow_waterbody) = model.boundary_conditions

    if !isnothing(reservoir) || !isnothing(lake)
        inflow_land = get_inflow_waterbody(model, overland_flow)
        inflow_subsurface = get_inflow_waterbody(model, subsurface_flow)

        @. inflow_waterbody = inflow_land[river_indices] + inflow_subsurface[river_indices]
    end
    return nothing
end

# For the river kinematic wave, the variable `to_river` can be excluded, because this part
# is added to the river kinematic wave.
get_inflow_waterbody(::KinWaveRiverFlow, model::KinWaveOverlandFlow) = model.variables.q_av
get_inflow_waterbody(::KinWaveRiverFlow, model::LateralSSF) =
    model.variables.ssf ./ Float(tosecond(BASETIMESTEP))

# Exclude subsurface flow from `GroundwaterFlow`.
get_inflow_waterbody(::AbstractRiverFlowModel, model::GroundwaterFlow) =
    zeros(model.connectivity.ncell)
