abstract type AbstractRiverFlowModel end

@get_units @grid_loc @with_kw struct FlowVariables{T}
    q::Vector{T} | "m3 s-1"                 # Discharge [m³ s⁻¹]
    qlat::Vector{T} | "m2 s-1"              # Lateral inflow per unit length [m² s⁻¹]
    qin::Vector{T} | "m3 s-1"               # Inflow from upstream cells [m³ s⁻¹]
    q_av::Vector{T} | "m3 s-1"              # Average discharge [m³ s⁻¹]
    volume::Vector{T} | "m3"                # Kinematic wave volume [m³] (based on water depth h)
    h::Vector{T} | "m"                      # Water depth [m]
    h_av::Vector{T} | "m"                   # Average water depth [m]
end

function init_kinematic_wave_timestepping(config, n; domain, dt_fixed)
    adaptive = get(config.model, "kin_wave_iteration", false)::Bool
    @info "Kinematic wave approach is used for $domain flow, adaptive timestepping = $adaptive."

    if adaptive
        stable_timesteps = zeros(n)
        timestepping = TimeStepping(; stable_timesteps, adaptive)
    else
        dt_fixed = get(config.model, "kw_$(domain)_tstep", dt_fixed)
        @info "Using a fixed sub-timestep (seconds) $dt_fixed for kinematic wave $domain flow."
        timestepping = TimeStepping(; dt_fixed, adaptive)
    end
    return timestepping
end

function FlowVariables(n)
    variables = FlowVariables(;
        q = zeros(Float, n),
        qlat = zeros(Float, n),
        qin = zeros(Float, n),
        q_av = zeros(Float, n),
        volume = zeros(Float, n),
        h = zeros(Float, n),
        h_av = zeros(Float, n),
    )
    return variables
end

@get_units @grid_loc @with_kw struct ManningFlowParameters{T}
    beta::T                                 # constant in Manning's equation [-]
    slope::Vector{T} | "m m-1"              # Slope [m m⁻¹]
    mannings_n::Vector{T} | "s m-1/3"       # Manning's roughness [s m⁻⅓]
    flow_length::Vector{T} | "m"            # Flow length [m]
    flow_width::Vector{T} | "m"             # Flow width [m]
    alpha_pow::T                            # Used in the power part of alpha [-]
    alpha_term::Vector{T} | "-"             # Term used in computation of alpha [-]
    alpha::Vector{T} | "s3/5 m1/5"          # Constant in momentum equation A = alpha*Q^beta, based on Manning's equation
end

function ManningFlowParameters(slope, mannings_n, flow_length, flow_width)
    n = length(slope)
    parameters = ManningFlowParameters(;
        beta = Float(0.6),
        slope,
        mannings_n,
        flow_length,
        flow_width,
        alpha_pow = Float((2.0 / 3.0) * 0.6),
        alpha_term = fill(mv, n),
        alpha = fill(mv, n),
    )
    return parameters
end

@get_units @grid_loc @with_kw struct RiverFlowParameters{T}
    flow::ManningFlowParameters{T}
    bankfull_depth::Vector{T} | "m"         # Bankfull water level [m]
end

function Base.getproperty(v::RiverFlowParameters, s::Symbol)
    if s === :bankfull_depth
        getfield(v, s)
    elseif s === :flow
        getfield(v, :flow)
    else
        getfield(getfield(v, :flow), s)
    end
end

function RiverFlowParameters(dataset, config, indices, river_length, river_width)
    mannings_n = ncread(
        dataset,
        config,
        "lateral.river.mannings_n";
        sel = indices,
        defaults = 0.036,
        type = Float,
    )
    bankfull_depth = ncread(
        dataset,
        config,
        "lateral.river.bankfull_depth";
        alias = "lateral.river.h_bankfull",
        sel = indices,
        defaults = 1.0,
        type = Float,
    )
    if haskey(config.input.lateral.river, "h_bankfull")
        @warn string(
            "The `h_bankfull` key in `[input.lateral.river]` is now called ",
            "`bankfull_depth`. Please update your TOML file.",
        )
    end
    slope = ncread(
        dataset,
        config,
        "lateral.river.slope";
        optional = false,
        sel = indices,
        type = Float,
    )
    clamp!(slope, 0.00001, Inf)

    flow_parameter_set = ManningFlowParameters(slope, mannings_n, river_length, river_width)
    parameters =
        RiverFlowParameters(; flow = flow_parameter_set, bankfull_depth = bankfull_depth)
    return parameters
end

@get_units @grid_loc @with_kw struct RiverFlowBC{T, R, L}
    inwater::Vector{T} | "m3 s-1"           # Lateral inflow [m³ s⁻¹]
    inflow::Vector{T} | "m3 s-1"            # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    inflow_waterbody::Vector{T} | "m3 s-1"  # inflow waterbody (lake or reservoir model) from land part [m³ s⁻¹]
    abstraction::Vector{T} | "m3 s-1"       # Abstraction (computed as part of water demand and allocation) [m³ s⁻¹]
    reservoir::R                            # Reservoir model struct of arrays
    lake::L                                 # Lake model struct of arrays
end

function RiverFlowBC(n, reservoir, lake)
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

@with_kw struct SurfaceFlowRiver{T, R, L, A} <: AbstractRiverFlowModel
    timestepping::TimeStepping{T}
    boundary_conditions::RiverFlowBC{T, R, L}
    parameters::RiverFlowParameters{T}
    variables::FlowVariables{T}
    allocation::A   # Water allocation
end

function SurfaceFlowRiver(
    dataset,
    config,
    indices;
    river_length,
    river_width,
    reservoir,
    lake,
)
    n = length(indices)

    timestepping =
        init_kinematic_wave_timestepping(config, n; domain = "river", dt_fixed = 900.0)

    do_water_demand = haskey(config.model, "water_demand")
    allocation = do_water_demand ? AllocationRiver(n) : NoAllocationRiver{Float}()

    variables = FlowVariables(n)
    parameters = RiverFlowParameters(dataset, config, indices, river_length, river_width)
    boundary_conditions = RiverFlowBC(n, reservoir, lake)

    sf_river = SurfaceFlowRiver(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        allocation,
    )

    return sf_river
end

@get_units @grid_loc @with_kw struct LandFlowVariables{T}
    flow::FlowVariables{T}
    to_river::Vector{T} | "m3 s-1"      # Part of overland flow [m³ s⁻¹] that flows to the river
end

function Base.getproperty(v::LandFlowVariables, s::Symbol)
    if s === :to_river
        getfield(v, s)
    elseif s === :flow
        getfield(v, :flow)
    else
        getfield(getfield(v, :flow), s)
    end
end

@get_units @grid_loc @with_kw struct LandFlowBC{T}
    inwater::Vector{T} | "m3 s-1"       # Lateral inflow [m³ s⁻¹]
end

@with_kw struct SurfaceFlowLand{T}
    timestepping::TimeStepping{T}
    boundary_conditions::LandFlowBC{T}
    parameters::ManningFlowParameters{T}
    variables::LandFlowVariables{T}
end

function SurfaceFlowLand(dataset, config, indices; slope, flow_length, flow_width)
    mannings_n = ncread(
        dataset,
        config,
        "lateral.land.mannings_n";
        sel = indices,
        defaults = 0.072,
        type = Float,
    )

    n = length(indices)
    timestepping =
        init_kinematic_wave_timestepping(config, n; domain = "land", dt_fixed = 3600.0)

    variables = LandFlowVariables(; flow = FlowVariables(n), to_river = zeros(Float, n))
    parameters = ManningFlowParameters(slope, mannings_n, flow_length, flow_width)
    boundary_conditions = LandFlowBC(; inwater = zeros(Float, n))
    sf_land = SurfaceFlowLand(; timestepping, boundary_conditions, variables, parameters)

    return sf_land
end

function surfaceflow_land_update!(model::SurfaceFlowLand, network, dt)
    (;
        order_of_subdomains,
        order_subdomain,
        subdomain_indices,
        upstream_nodes,
        frac_to_river,
    ) = network

    (; beta, alpha, flow_width, flow_length) = model.parameters
    (; h, h_av, q, q_av, qin, qlat, to_river) = model.variables

    ns = length(order_of_subdomains)
    qin .= 0.0
    for k in 1:ns
        threaded_foreach(eachindex(order_of_subdomains[k]); basesize = 1) do i
            m = order_of_subdomains[k][i]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                # for a river cell without a reservoir or lake part of the upstream
                # surface flow goes to the river (frac_to_river) and part goes to the
                # surface flow reservoir (1.0 - frac_to_river), upstream nodes with a
                # reservoir or lake are excluded
                to_river[v] +=
                    sum_at(
                        i -> q[i] * frac_to_river[i],
                        upstream_nodes[n],
                        eltype(to_river),
                    ) * dt
                if flow_width[v] > 0.0
                    qin[v] = sum_at(
                        i -> q[i] * (1.0 - frac_to_river[i]),
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
                if flow_width[v] > 0.0
                    crossarea = alpha[v] * pow(q[v], beta)
                    h[v] = crossarea / flow_width[v]
                end
                q_av[v] += q[v] * dt
                h_av[v] += h[v] * dt
            end
        end
    end
end

function update!(model::SurfaceFlowLand, network, dt)
    (; inwater) = model.boundary_conditions
    (; alpha_term, mannings_n, slope, beta, alpha_pow, alpha, flow_width, flow_length) =
        model.parameters
    (; h, h_av, q_av, qlat, volume, to_river) = model.variables
    (; adaptive) = model.timestepping

    @. alpha_term = pow(mannings_n / sqrt(slope), beta)
    # use fixed alpha value based flow width
    @. alpha = alpha_term * pow(flow_width, alpha_pow)
    @. qlat = inwater / flow_length

    q_av .= 0.0
    h_av .= 0.0
    to_river .= 0.0

    t = 0.0
    while t < dt
        dt_s = adaptive ? stable_timestep(model, 0.02) : model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        surfaceflow_land_update!(model, network, dt_s)
        t = t + dt_s
    end
    q_av ./= dt
    h_av ./= dt
    to_river ./= dt
    volume .= flow_length .* flow_width .* h
    return nothing
end

function surfaceflow_river_update!(model::SurfaceFlowRiver, network, doy, dt)
    (;
        graph,
        order_of_subdomains,
        order_subdomain,
        subdomain_indices,
        upstream_nodes,
        reservoir_indices,
        lake_indices,
    ) = network.river

    (; reservoir, lake, inwater, inflow, abstraction, inflow_waterbody) =
        model.boundary_conditions

    (; beta, alpha, flow_width, flow_length) = model.parameters
    (; h, h_av, q, q_av, qin, qlat, volume) = model.variables

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
                    max_abstract =
                        min((inwater[v] + qin[v] + volume[v] / dt) * 0.80, -inflow[v])
                    _inflow = -max_abstract / flow_length[v]
                else
                    _inflow = inflow[v] / flow_length[v]
                end
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
                    update!(reservoir, i, q[v] + inflow_waterbody[v], dt)

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
                    update!(lake, i, q[v] + inflow_waterbody[v], doy, dt)

                    downstream_nodes = outneighbors(graph, v)
                    n_downstream = length(downstream_nodes)
                    if n_downstream == 1
                        j = only(downstream_nodes)
                        qin[j] = lake.variables.outflow[i]
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
                # update h
                crossarea = alpha[v] * pow(q[v], beta)
                h[v] = crossarea / flow_width[v]
                volume[v] = flow_length[v] * flow_width[v] * h[v]
                q_av[v] += q[v] * dt
                h_av[v] += h[v] * dt
            end
        end
    end
end

function update!(model::SurfaceFlowRiver, network, doy, dt)
    (; reservoir, lake, inwater) = model.boundary_conditions

    (;
        alpha_term,
        mannings_n,
        slope,
        beta,
        alpha_pow,
        alpha,
        flow_width,
        flow_length,
        bankfull_depth,
    ) = model.parameters
    (; h, h_av, q_av, qlat, volume) = model.variables
    (; adaptive) = model.timestepping

    @. alpha_term = pow(mannings_n / sqrt(slope), beta)
    # use fixed alpha value based on 0.5 * bankfull_depth
    @. alpha = alpha_term * pow(flow_width + bankfull_depth, alpha_pow)
    @. qlat = inwater / flow_length

    q_av .= 0.0
    h_av .= 0.0
    # because of possible iterations set reservoir and lake inflow and total outflow at
    # start to zero, the total sum of inflow and outflow at each sub time step is calculated
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

    t = 0.0
    while t < dt
        dt_s = adaptive ? stable_timestep(model, 0.05) : model.timestepping.dt_fixed
        dt_s = check_timestepsize(dt_s, t, dt)
        surfaceflow_river_update!(model, network, doy, dt_s)
        t = t + dt_s
    end
    q_av ./= dt
    h_av ./= dt
    volume .= flow_length .* flow_width .* h
    return nothing
end

function stable_timestep(model::S, p) where {S <: Union{SurfaceFlowLand, SurfaceFlowRiver}}
    (; q) = model.variables
    (; alpha, beta, flow_length) = model.parameters
    (; stable_timesteps) = model.timestepping

    n = length(q)
    stable_timesteps .= Inf
    for i in 1:n
        if q[i] > 0.0
            c = 1.0 / (alpha[i] * beta * pow(q[i], (beta - 1.0)))
            stable_timesteps[i] = (flow_length[i] / c)
        end
    end
    _stable_timesteps = filter(x -> !isinf(x), stable_timesteps)

    if !isempty(_stable_timesteps)
        dt_s = quantile!(_stable_timesteps, p)
    else
        dt_s = 600.0
    end

    return dt_s
end

"""
Update boundary condition lateral inflow `inwater` of a `river` flow model for a single
timestep.
"""
function update_lateral_inflow!(
    model::AbstractRiverFlowModel,
    external_models::NamedTuple,
    river_cell_area,
    land_area,
    river_indices,
    dt,
)
    (; allocation, runoff, land, subsurface) = external_models
    (; inwater) = model.boundary_conditions
    (; net_runoff_river) = runoff.variables

    inwater .= (
        get_flux_to_river(subsurface)[river_indices] .+
        land.variables.to_river[river_indices] .+
        (net_runoff_river[river_indices] .* land_area[river_indices] .* 0.001) ./ dt .+
        (get_nonirrigation_returnflow(allocation) .* 0.001 .* river_cell_area) ./ dt
    )
    return nothing
end

"""
Update boundary condition lateral inflow `inwater` of kinematic wave overland flow model
`SurfaceFlowLand` for a single timestep.
"""
function update_lateral_inflow!(
    model::SurfaceFlowLand,
    external_models::NamedTuple,
    area,
    config,
    dt,
)
    (; soil, subsurface, allocation) = external_models
    (; net_runoff) = soil.variables
    (; inwater) = model.boundary_conditions

    do_drains = get(config.model, "drains", false)::Bool
    if do_drains
        drainflux = zeros(length(net_runoff))
        drainflux[subsurface.drain.index] =
            -subsurface.drain.variables.flux ./ tosecond(basetimestep)
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
    river_indices,
)
    (; land, subsurface) = external_models
    (; reservoir, lake, inflow_waterbody) = model.boundary_conditions

    if !isnothing(reservoir) || !isnothing(lake)
        inflow_land = get_inflow_waterbody(model, land)
        inflow_subsurface = get_inflow_waterbody(model, subsurface)

        @. inflow_waterbody = inflow_land[river_indices] + inflow_subsurface[river_indices]
    end
    return nothing
end

# For the river kinematic wave, the variable `to_river` can be excluded, because this part
# is added to the river kinematic wave.
get_inflow_waterbody(::SurfaceFlowRiver, model::SurfaceFlowLand) = model.variables.q_av
get_inflow_waterbody(::SurfaceFlowRiver, model::LateralSSF) =
    model.variables.ssf ./ tosecond(basetimestep)

# Exclude subsurface flow for other groundwater components than `LateralSSF`.
get_inflow_waterbody(::AbstractRiverFlowModel, model::GroundwaterFlow) =
    model.flow.connectivity.ncell .* 0.0
get_inflow_waterbody(::SurfaceFlowRiver, model) = model.variables.to_river .* 0.0