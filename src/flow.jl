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

@with_kw struct TimeStepping{T}
    stable_timesteps::Vector{T} = Float[]
    dt_fixed::T = 0.0
    adaptive::Bool = true
    cfl::T = 0.70
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

function check_timestepsize(timestepsize, currenttime, endtime)
    if currenttime + timestepsize > endtime
        timestepsize = endtime - currenttime
    end
    return timestepsize
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

@get_units @grid_loc struct KhExponential{T}
    # Horizontal hydraulic conductivity at soil surface [m d⁻¹]
    kh_0::Vector{T} | "m d-1"
    # A scaling parameter [m⁻¹] (controls exponential decline of kh_0)
    f::Vector{T} | "m-1"
end

@get_units @grid_loc struct KhExponentialConstant{T}
    # Exponential horizontal hydraulic conductivity profile type
    exponential::KhExponential
    # Depth [m] from soil surface for which exponential decline of kv_0 is valid
    z_exp::Vector{T} | "m"
end

@get_units @grid_loc struct KhLayered{T}
    # Horizontal hydraulic conductivity [m d⁻¹]
    kh::Vector{T} | "m d-1"
end

abstract type SubsurfaceFlow end

@get_units @grid_loc @with_kw struct LateralSsfParameters{T, Kh}
    kh_profile::Kh                         # Horizontal hydraulic conductivity profile type [-]  
    khfrac::Vector{T} | "-"                # A muliplication factor applied to vertical hydraulic conductivity `kv` [-]
    soilthickness::Vector{T} | "m"         # Soil thickness [m]
    theta_s::Vector{T} | "-"               # Saturated water content (porosity) [-]
    theta_r::Vector{T} | "-"               # Residual water content [-]   
    dt::T                                  # model time step [d]
    slope::Vector{T} | "m m-1"             # Slope [m m⁻¹]
    flow_length::Vector{T} | "m"           # Flow length [m]
    flow_width::Vector{T} | "m"            # Flow width [m]
end

function LateralSsfParameters(
    dataset,
    config,
    indices;
    soil,
    slope,
    flow_length,
    flow_width,
    dt,
)
    khfrac = ncread(
        dataset,
        config,
        "lateral.subsurface.ksathorfrac";
        sel = indices,
        defaults = 1.0,
        type = Float,
    )
    n_cells = length(khfrac)

    (; theta_s, theta_r, soilthickness) = soil
    soilthickness = soilthickness .* 0.001

    kh_profile_type = get(config.input.vertical, "ksat_profile", "exponential")::String
    if kh_profile_type == "exponential"
        (; kv_0, f) = soil.kv_profile
        kh_0 = khfrac .* kv_0 .* 0.001 .* dt
        kh_profile = KhExponential(kh_0, f .* 1000.0)
    elseif kh_profile_type == "exponential_constant"
        (; z_exp) = soil.kv_profile
        (; kv_0, f) = soil.kv_profile.exponential
        kh_0 = khfrac .* kv_0 .* 0.001 .* dt
        exp_profile = KhExponential(kh_0, f .* 1000.0)
        kh_profile = KhExponentialConstant(exp_profile, z_exp .* 0.001)
    elseif kh_profile_type == "layered" || kh_profile_type == "layered_exponential"
        kh_profile = KhLayered(fill(mv, n_cells))
    end
    parameters = LateralSsfParameters(
        kh_profile,
        khfrac,
        soilthickness,
        theta_s,
        theta_r,
        dt,
        slope,
        flow_length,
        flow_width,
    )
    return parameters
end

@get_units @grid_loc @with_kw struct LateralSsfVariables{T}
    zi::Vector{T} | "m"                    # Pseudo-water table depth [m] (top of the saturated zone)
    exfiltwater::Vector{T} | "m dt-1"      # Exfiltration [m Δt⁻¹] (groundwater above surface level, saturated excess conditions)
    recharge::Vector{T} | "m2 dt-1"        # Net recharge to saturated store [m² Δt⁻¹]
    ssf::Vector{T} | "m3 d-1"              # Subsurface flow [m³ d⁻¹]
    ssfin::Vector{T} | "m3 d-1"            # Inflow from upstream cells [m³ d⁻¹]
    ssfmax::Vector{T} | "m2 d-1"           # Maximum subsurface flow [m² d⁻¹]
    to_river::Vector{T} | "m3 d-1"         # Part of subsurface flow [m³ d⁻¹] that flows to the river
    volume::Vector{T} | "m3"               # Subsurface volume [m³]
end

function LateralSsfVariables(ssf, zi, xl, yl)
    n = length(zi)
    volume = @. (ssf.theta_s - ssf.theta_r) * (ssf.soilthickness - zi) * (xl * yl)
    variables = LateralSsfVariables(;
        zi,
        exfiltwater = fill(mv, n),
        recharge = fill(mv, n),
        ssf = fill(mv, n),
        ssfin = fill(mv, n),
        ssfmax = fill(mv, n),
        to_river = zeros(n),
        volume,
    )
    return variables
end

@get_units @grid_loc @with_kw struct LateralSsfBC{T}
    recharge::Vector{T} | "m2 dt-1"        # Net recharge to saturated store [m² Δt⁻¹]
end

@with_kw struct LateralSSF{T, Kh} <: SubsurfaceFlow
    boundary_conditions::LateralSsfBC{T}
    parameters::LateralSsfParameters{T, Kh}
    variables::LateralSsfVariables{T}
end

function LateralSSF(
    dataset,
    config,
    indices;
    soil,
    slope,
    flow_length,
    flow_width,
    x_length,
    y_length,
    dt,
)
    parameters = LateralSsfParameters(
        dataset,
        config,
        indices;
        soil = soil.parameters,
        slope,
        flow_length,
        flow_width,
        dt,
    )
    zi = 0.001 * soil.variables.zi
    variables = LateralSsfVariables(parameters, zi, x_length, y_length)
    boundary_conditions = LateralSsfBC(; recharge = fill(mv, length(zi)))
    ssf = LateralSSF(; boundary_conditions, parameters, variables)
    return ssf
end

function update!(model::LateralSSF, network)
    (;
        order_of_subdomains,
        order_subdomain,
        subdomain_indices,
        upstream_nodes,
        area,
        frac_to_river,
    ) = network

    (; recharge) = model.boundary_conditions
    (; ssfin, ssf, ssfmax, to_river, zi, exfiltwater, volume) = model.variables
    (; slope, theta_s, theta_r, soilthickness, flow_length, flow_width, dt, kh_profile) =
        model.parameters

    ns = length(order_of_subdomains)
    for k in 1:ns
        threaded_foreach(eachindex(order_of_subdomains[k]); basesize = 1) do i
            m = order_of_subdomains[k][i]
            for (n, v) in zip(subdomain_indices[m], order_subdomain[m])
                # for a river cell without a reservoir or lake part of the upstream
                # subsurface flow goes to the river (frac_to_river) and part goes to the
                # subsurface flow reservoir (1.0 - frac_to_river) upstream nodes with a
                # reservoir or lake are excluded
                ssfin[v] = sum_at(
                    i -> ssf[i] * (1.0 - frac_to_river[i]),
                    upstream_nodes[n],
                    eltype(ssfin),
                )
                to_river[v] = sum_at(
                    i -> ssf[i] * frac_to_river[i],
                    upstream_nodes[n],
                    eltype(to_river),
                )
                ssf[v], zi[v], exfiltwater[v] = kinematic_wave_ssf(
                    ssfin[v],
                    ssf[v],
                    zi[v],
                    recharge[v],
                    slope[v],
                    theta_s[v] - theta_r[v],
                    soilthickness[v],
                    dt,
                    flow_length[v],
                    flow_width[v],
                    ssfmax[v],
                    kh_profile,
                    v,
                )
                volume[v] = (theta_s[v] - theta_r[v]) * (soilthickness[v] - zi[v]) * area[v]
            end
        end
    end
    return nothing
end

@get_units@grid_loc @with_kw struct GroundwaterExchangeVariables{T}
    exfiltwater::Vector{T} | "m dt-1"   # Exfiltration [m Δt⁻¹]  (groundwater above surface level, saturated excess conditions)
    zi::Vector{T} | "m"                 # Pseudo-water table depth [m] (top of the saturated zone)
    to_river::Vector{T} | "m3 d-1"      # Part of subsurface flow [m³ d⁻¹] that flows to the river
    ssf::Vector{T} | "m3 d-1"           # Subsurface flow [m³ d⁻¹]
end

function GroundwaterExchangeVariables(n)
    variables = GroundwaterExchangeVariables{Float}(;
        exfiltwater = fill(mv, n),
        zi = fill(mv, n),
        to_river = fill(mv, n),
        ssf = zeros(n),
    )
    return variables
end

@with_kw struct GroundwaterExchangeParameters{T}
    dt::T       # model time step [d]
end

@with_kw struct GroundwaterExchange{T} <: SubsurfaceFlow
    parameters::GroundwaterExchangeParameters{T}
    variables::GroundwaterExchangeVariables{T}
end

function GroundwaterExchange(n, dt)
    parameters = GroundwaterExchangeParameters{Float}(; dt = dt / basetimestep)
    variables = GroundwaterExchangeVariables(n)
    ssf = GroundwaterExchange{Float}(; parameters, variables)
    return ssf
end

get_water_depth(subsurface::SubsurfaceFlow) = subsurface.variables.zi
get_exfiltwater(subsurface::SubsurfaceFlow) = subsurface.variables.exfiltwater

get_flux_to_river(subsurface::SubsurfaceFlow) =
    subsurface.variables.to_river ./ tosecond(basetimestep) # [m³ s⁻¹]

@get_units @grid_loc @with_kw struct ShallowWaterRiverParameters{T}
    n::Int                                              # number of cells [-]
    ne::Int                                             # number of edges/links [-]
    active_n::Vector{Int} | "-"                         # active nodes [-]
    active_e::Vector{Int} | "-" | "edge"                # active edges/links [-]
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
    flow_length_at_link::Vector{T} | "m" | "edge"       # flow (river) length at edge
    flow_width::Vector{T} | "m"                         # flow (river) width
    flow_width_at_link::Vector{T} | "m" | "edge"        # flow (river) width at edge
    waterbody::Vector{Bool} | "-"                       # water body cells (reservoir or lake)
end

function ShallowWaterRiverParameters(
    dataset,
    config,
    indices;
    river_length,
    river_width,
    waterbody,
    n_edges,
    nodes_at_link,
    index_pit,
    inds_pit,
)
    cfl = get(config.model, "inertial_flow_alpha", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    h_thresh = get(config.model, "h_thresh", 1.0e-03)::Float64 # depth threshold for flow at link
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
        "lateral.river.riverlength_bc";
        sel = inds_pit,
        defaults = 1.0e04,
        type = Float,
    )
    bankfull_elevation_2d = ncread(
        dataset,
        config,
        "lateral.river.bankfull_elevation";
        optional = false,
        type = Float,
        fill = 0,
    )
    bankfull_depth_2d = ncread(
        dataset,
        config,
        "lateral.river.bankfull_depth";
        optional = false,
        type = Float,
        fill = 0,
    )
    bankfull_depth = bankfull_depth_2d[indices]
    zb = bankfull_elevation_2d[indices] - bankfull_depth # river bed elevation

    bankfull_volume = bankfull_depth .* river_width .* river_length
    mannings_n = ncread(
        dataset,
        config,
        "lateral.river.mannings_n";
        sel = indices,
        defaults = 0.036,
        type = Float,
    )

    n = length(indices)

    # set ghost points for boundary condition (downstream river outlet): river width, bed
    # elevation, manning n is copied from the upstream cell.
    append!(river_length, riverlength_bc)
    append!(zb, zb[index_pit])
    append!(river_width, river_width[index_pit])
    append!(mannings_n, mannings_n[index_pit])
    append!(bankfull_depth, bankfull_depth[index_pit])

    # determine z, width, length and manning's n at links
    zb_max = fill(Float(0), n_edges)
    width_at_link = fill(Float(0), n_edges)
    length_at_link = fill(Float(0), n_edges)
    mannings_n_sq = fill(Float(0), n_edges)
    for i in 1:n_edges
        src_node = nodes_at_link.src[i]
        dst_node = nodes_at_link.dst[i]
        zb_max[i] = max(zb[src_node], zb[dst_node])
        width_at_link[i] = min(river_width[src_node], river_width[dst_node])
        length_at_link[i] = 0.5 * (river_length[dst_node] + river_length[src_node])
        mannings_n_i =
            (
                mannings_n[dst_node] * river_length[dst_node] +
                mannings_n[src_node] * river_length[src_node]
            ) / (river_length[dst_node] + river_length[src_node])
        mannings_n_sq[i] = mannings_n_i * mannings_n_i
    end
    active_index = findall(x -> x == 0, waterbody)

    parameters = ShallowWaterRiverParameters(;
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
        flow_length_at_link = length_at_link,
        flow_width = river_width,
        flow_width_at_link = width_at_link,
        waterbody,
    )
    return parameters
end

@get_units @grid_loc @with_kw struct ShallowWaterRiverVariables{T}
    q::Vector{T} | "m3 s-1" | "edge"                    # river discharge (subgrid channel)
    q0::Vector{T} | "m3 s-1" | "edge"                   # river discharge (subgrid channel) at previous time step
    q_av::Vector{T} | "m3 s-1" | "edge"                 # average river channel (+ floodplain) discharge [m³ s⁻¹]
    q_channel_av::Vector{T} | "m3 s-1"                  # average river channel discharge [m³ s⁻¹]
    h::Vector{T} | "m"                                  # water depth
    zs_max::Vector{T} | "m" | "edge"                    # maximum water elevation at edge
    zs_src::Vector{T} | "m"                             # water elevation of source node of edge
    zs_dst::Vector{T} | "m"                             # water elevation of downstream node of edge
    hf::Vector{T} | "m" | "edge"                        # water depth at edge/link
    h_av::Vector{T} | "m"                               # average water depth
    a::Vector{T} | "m2" | "edge"                        # flow area at edge/link
    r::Vector{T} | "m" | "edge"                         # wetted perimeter at edge/link
    volume::Vector{T} | "m3"                            # river volume
    error::Vector{T} | "m3"                             # error volume    
end

function ShallowWaterRiverVariables(dataset, config, indices, n_edges, inds_pit)
    floodplain_1d = get(config.model, "floodplain_1d", false)::Bool
    riverdepth_bc = ncread(
        dataset,
        config,
        "lateral.river.riverdepth_bc";
        sel = inds_pit,
        defaults = 0.0,
        type = Float,
    )

    n = length(indices)
    # set river depth h to zero (including reservoir and lake locations)
    h = zeros(n)
    q_av = zeros(n_edges)
    # set ghost points for boundary condition (downstream river outlet): river depth `h`
    append!(h, riverdepth_bc)
    variables = ShallowWaterRiverVariables(;
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

@with_kw struct ShallowWaterRiver{T, R, L, F, A} <: AbstractRiverFlowModel
    timestepping::TimeStepping{T}
    boundary_conditions::RiverFlowBC{T, R, L}
    parameters::ShallowWaterRiverParameters{T}
    variables::ShallowWaterRiverVariables{T}
    floodplain::F                                       # Floodplain (1D) schematization
    allocation::A                                       # Water allocation
end

function ShallowWaterRiver(
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
    # with nodes and links. This information is extracted from the directed graph of the
    # river. Discharge q is calculated at links between nodes and mapped to the source
    # nodes for gridded output (index of link is equal to source node index, e.g.:
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

    #waterbody = !=(0).(reservoir_indices .+ lake_indices)
    add_vertex_edge_graph!(graph_river, index_pit)
    nodes_at_link = adjacent_nodes_at_link(graph_river)
    n_edges = ne(graph_river)

    parameters = ShallowWaterRiverParameters(
        dataset,
        config,
        indices;
        river_length,
        river_width,
        waterbody,
        n_edges,
        nodes_at_link,
        index_pit,
        inds_pit,
    )
    variables = ShallowWaterRiverVariables(dataset, config, indices, n_edges, inds_pit)

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
            nodes_at_link,
        )
    else
        floodplain = nothing
    end

    do_water_demand = haskey(config.model, "water_demand")
    sw_river = ShallowWaterRiver(;
        timestepping,
        boundary_conditions,
        parameters,
        variables,
        floodplain,
        allocation = do_water_demand ? AllocationRiver(n) : NoAllocationRiver{Float}(),
    )
    return sw_river, nodes_at_link
end

"Return the upstream inflow for a waterbody in `ShallowWaterRiver`"
function get_inflow_waterbody(sw::ShallowWaterRiver, src_edge)
    q_in = sum_at(sw.variables.q, src_edge)
    if !isnothing(sw.floodplain)
        q_in = q_in + sum_at(sw.floodplain.variables.q, src_edge)
    end
    return q_in
end

function shallowwater_river_update!(model::ShallowWaterRiver, network, dt, doy, update_h)
    (; nodes_at_link, links_at_node) = network.river
    (; inwater, abstraction, inflow) = model.boundary_conditions
    river_v = model.variables
    river_p = model.parameters

    river_v.q0 .= river_v.q
    if !isnothing(model.floodplain)
        model.floodplain.variables.q0 .= model.floodplain.variables.q
    end
    @tturbo for j in eachindex(river_p.active_e)
        i = river_p.active_e[j]
        i_src = nodes_at_link.src[i]
        i_dst = nodes_at_link.dst[i]
        river_v.zs_src[i] = river_p.zb[i_src] + river_v.h[i_src]
        river_v.zs_dst[i] = river_p.zb[i_dst] + river_v.h[i_dst]

        river_v.zs_max[i] = max(river_v.zs_src[i], river_v.zs_dst[i])
        river_v.hf[i] = (river_v.zs_max[i] - river_p.zb_max[i])

        river_v.a[i] = river_p.flow_width_at_link[i] * river_v.hf[i] # flow area (rectangular channel)
        river_v.r[i] = river_v.a[i] / (river_p.flow_width_at_link[i] + 2.0 * river_v.hf[i]) # hydraulic radius (rectangular channel)

        river_v.q[i] = IfElse.ifelse(
            river_v.hf[i] > river_p.h_thresh,
            local_inertial_flow(
                river_v.q0[i],
                river_v.zs_src[i],
                river_v.zs_dst[i],
                river_v.hf[i],
                river_v.a[i],
                river_v.r[i],
                river_p.flow_length_at_link[i],
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
            i_src = nodes_at_link.src[i]
            i_dst = nodes_at_link.dst[i]

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
                    river_p.flow_length_at_link[i],
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

        q_in = get_inflow_waterbody(model, links_at_node.src[i])
        update!(reservoir, v, q_in + inflow_waterbody[i], dt)
        river_v.q[i] = reservoir.variables.outflow[v]
        river_v.q_av[i] += river_v.q[i] * dt
    end
    (; lake, inflow_waterbody) = model.boundary_conditions
    inds_lake = network.lake.river_indices
    for v in eachindex(inds_lake)
        i = inds_lake[v]

        q_in = get_inflow_waterbody(model, links_at_node.src[i])
        update!(lake, v, q_in + inflow_waterbody[i], doy, dt)
        river_v.q[i] = lake.variables.outflow[v]
        river_v.q_av[i] += river_v.q[i] * dt
    end
    if update_h
        @batch per = thread minbatch = 2000 for i in river_p.active_n
            q_src = sum_at(river_v.q, links_at_node.src[i])
            q_dst = sum_at(river_v.q, links_at_node.dst[i])
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
                q_src = sum_at(floodplain_v.q, links_at_node.src[i])
                q_dst = sum_at(floodplain_v.q, links_at_node.dst[i])
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

function update!(model::ShallowWaterRiver{T}, network, doy, dt; update_h = true) where {T}
    (; reservoir, lake) = model.boundary_conditions
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
        shallowwater_river_update!(model, network, dt_s, doy, update_h)
        t = t + dt_s
    end
    model.variables.q_av ./= dt
    model.variables.h_av ./= dt

    if !isnothing(model.floodplain)
        model.floodplain.variables.q_av ./= dt
        model.floodplain.variables.h_av ./= dt
        model.variables.q_channel_av .= model.variables.q_av
        model.variables.q_av .=
            model.variables.q_channel_av .+ model.floodplain.variables.q_av
    end

    return nothing
end

# Stores links in x and y direction between cells of a Vector with CartesianIndex(x, y), for
# staggered grid calculations.
@with_kw struct Indices
    xu::Vector{Int}     # index of neighbor cell in the (+1, 0) direction
    xd::Vector{Int}     # index of neighbor cell in the (-1, 0) direction
    yu::Vector{Int}     # index of neighbor cell in the (0, +1) direction
    yd::Vector{Int}     # index of neighbor cell in the (0, -1) direction
end

# maps the fields of struct Indices to the defined Wflow cartesian indices of const
# neigbors.
const dirs = (:yd, :xd, :xu, :yu)

@get_units @grid_loc @with_kw struct ShallowWaterLandVariables{T}
    qy0::Vector{T} | "m3 s-1" | "edge"      # flow in y direction at previous time step
    qx0::Vector{T} | "m3 s-1" | "edge"      # flow in x direction at previous time step
    qx::Vector{T} | "m3 s-1" | "edge"       # flow in x direction
    qy::Vector{T} | "m3 s-1" | "edge"       # flow in y direction
    volume::Vector{T} | "m3"                # total volume of cell (including river volume for river cells)
    error::Vector{T} | "m3"                 # error volume
    h::Vector{T} | "m"                      # water depth of cell (for river cells the reference is the river bed elevation `zb`)
    h_av::Vector{T} | "m"                   # average water depth (for river cells the reference is the river bed elevation `zb`)
end

function ShallowWaterLandVariables(n)
    variables = ShallowWaterLandVariables(;
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

@get_units @grid_loc @with_kw struct ShallowWaterLandParameters{T}
    n::Int                                              # number of cells [-]
    xl::Vector{T} | "m"                                 # cell length x direction [m]
    yl::Vector{T} | "m"                                 # cell length y direction [m]
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

function ShallowWaterLandParameters(
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
    h_thresh = get(config.model, "h_thresh", 1.0e-03)::Float64 # depth threshold for flow at link

    @info "Local inertial approach is used for overlandflow." cfl theta h_thresh froude_limit

    mannings_n = ncread(
        dataset,
        config,
        "lateral.land.mannings_n";
        sel = indices,
        defaults = 0.072,
        type = Float,
    )
    elevation_2d = ncread(
        dataset,
        config,
        "lateral.land.elevation";
        optional = false,
        type = Float,
        fill = 0,
    )
    elevation = elevation_2d[indices]
    n = length(indices)

    # initialize edges between cells in x and y direction.
    staggered_indices =
        Indices(; xu = zeros(n), xd = zeros(n), yu = zeros(n), yd = zeros(n))

    # edges without neigbors are handled by an extra index (at n + 1, with n links), which
    # is set to a value of 0.0 m³ s⁻¹ for qx and qy fields at initialization.
    # edges are defined as follows for the x and y direction, respectively:
    # node i => node xu (node i + CartesianIndex(1, 0))
    # node i => node yu (node i + CartesianIndex(0, 1))
    # where i is the index of indices
    nrow, ncol = modelsize_2d
    for (v, i) in enumerate(indices)
        for (m, neighbor) in enumerate(neighbors)
            j = i + neighbor
            dir = dirs[m]
            if (1 <= j[1] <= nrow) && (1 <= j[2] <= ncol) && (reverse_indices[j] != 0)
                getfield(staggered_indices, dir)[v] = reverse_indices[j]
            else
                getfield(staggered_indices, dir)[v] = n + 1
            end
        end
    end

    # determine z at links in x and y direction
    zx_max = fill(Float(0), n)
    zy_max = fill(Float(0), n)
    for i in 1:n
        xu = staggered_indices.xu[i]
        if xu <= n
            zx_max[i] = max(elevation[i], elevation[xu])
        end
        yu = staggered_indices.yu[i]
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
        staggered_indices,
        graph_river,
        river_width,
        ldd_river,
        waterbody,
        reverse_indices[inds_river],
    )
    parameters = ShallowWaterLandParameters(;
        n = n,
        xl = x_length,
        yl = y_length,
        xwidth = we_x,
        ywidth = we_y,
        g = 9.80665,
        theta = theta,
        h_thresh = h_thresh,
        zx_max = zx_max,
        zy_max = zy_max,
        mannings_n_sq = mannings_n .* mannings_n,
        z = elevation,
        froude_limit = froude_limit,
        rivercells = river_location,
    )
    return parameters, staggered_indices
end

@get_units @grid_loc @with_kw struct ShallowWaterLandBC{T}
    runoff::Vector{T} | "m3 s-1"               # runoff from hydrological model
    inflow_waterbody::Vector{T} | "m3 s-1"     # inflow to water body from hydrological model
end

function ShallowWaterLandBC(n)
    bc = ShallowWaterLandBC(; runoff = zeros(n), inflow_waterbody = zeros(n))
    return bc
end

@with_kw struct ShallowWaterLand{T}
    timestepping::TimeStepping{T}
    boundary_conditions::ShallowWaterLandBC{T}
    parameters::ShallowWaterLandParameters{T}
    variables::ShallowWaterLandVariables{T}
end

function ShallowWaterLand(
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
    boundary_conditions = ShallowWaterLandBC(n)
    parameters, staggered_indices = ShallowWaterLandParameters(
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
    variables = ShallowWaterLandVariables(n)

    sw_land =
        ShallowWaterLand{Float}(; timestepping, boundary_conditions, parameters, variables)

    return sw_land, staggered_indices
end

"""
    stable_timestep(sw::ShallowWaterRiver)
    stable_timestep(sw::ShallowWaterLand)

Compute a stable timestep size for the local inertial approach, based on Bates et al. (2010).

dt = cfl * (Δx / sqrt(g max(h))
"""
function stable_timestep(sw::ShallowWaterRiver{T})::T where {T}
    dt_min = T(Inf)
    (; cfl) = sw.timestepping
    (; n, flow_length, g) = sw.parameters
    (; h) = sw.variables
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(n)
        @fastmath @inbounds dt = cfl * flow_length[i] / sqrt(g * h[i])
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? T(10.0) : dt_min
    return dt_min
end

function stable_timestep(sw::ShallowWaterLand{T})::T where {T}
    dt_min = T(Inf)
    (; cfl) = sw.timestepping
    (; n, g, xl, yl, rivercells) = sw.parameters
    (; h) = sw.variables
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(n)
        @fastmath @inbounds dt = if rivercells[i] == 0
            cfl * min(xl[i], yl[i]) / sqrt(g * h[i])
        else
            T(Inf)
        end
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? T(10.0) : dt_min
    return dt_min
end

function update!(
    land::ShallowWaterLand{T},
    river::ShallowWaterRiver{T},
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
        shallowwater_river_update!(river, network, dt_s, doy, update_h)
        shallowwater_update!(land, river, network, dt_s)
        t = t + dt_s
    end
    river.variables.q_av ./= dt
    river.variables.h_av ./= dt
    land.variables.h_av ./= dt

    return nothing
end

function shallowwater_update!(
    land::ShallowWaterLand{T},
    river::ShallowWaterRiver{T},
    network,
    dt,
) where {T}
    indices = network.land.staggered_indices
    inds_river = network.land.river_indices

    (; links_at_node) = network.river

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
                length = T(0.5) * (land_p.xl[i] + land_p.xl[xu]) # can be precalculated
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
                length = T(0.5) * (land_p.yl[i] + land_p.yl[yu]) # can be precalculated
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
                        sum_at(river_v.q, links_at_node.src[inds_river[i]]) -
                        sum_at(river_v.q, links_at_node.dst[inds_river[i]]) +
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
                        (land_p.xl[i] * land_p.yl[i])
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
            land_v.h[i] = land_v.volume[i] / (land_p.xl[i] * land_p.yl[i])
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

function FloodPlainProfile(dataset, config, indices; river_width, river_length, index_pit)
    volume = ncread(
        dataset,
        config,
        "lateral.river.floodplain.volume";
        sel = indices,
        type = Float,
        dimname = :flood_depth,
    )
    n = length(indices)

    # for convenience (interpolation) flood depth 0.0 m is added, with associated area (a),
    # volume, width (river width) and wetted perimeter (p).
    volume = vcat(fill(Float(0), n)', volume)
    start_volume = volume
    flood_depths = Float.(dataset["flood_depth"][:])
    pushfirst!(flood_depths, 0.0)
    n_depths = length(flood_depths)

    p = zeros(Float, n_depths, n)
    a = zeros(Float, n_depths, n)
    segment_volume = zeros(Float, n_depths, n)
    width = zeros(Float, n_depths, n)
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
        FloodPlainProfile{Float, n_depths}(; volume, width, depth = flood_depths, a, p)
    return profile
end

@get_units @grid_loc @with_kw struct FloodPlainParameters{T, P}
    profile::P                                          # floodplain profile
    mannings_n::Vector{T} | "s m-1/3"                   # manning's roughness
    mannings_n_sq::Vector{T} | "(s m-1/3)2" | "edge"    # manning's roughness squared
    zb_max::Vector{T} | "m" | "edge"                    # maximum bankfull elevation (edge/link)
end

function FloodPlainParameters(
    dataset,
    config,
    indices;
    river_width,
    river_length,
    zb_floodplain,
    nodes_at_link,
    n_edges,
    index_pit,
)
    profile =
        FloodPlainProfile(dataset, config, indices; river_width, river_length, index_pit)

    mannings_n = ncread(
        dataset,
        config,
        "lateral.river.floodplain.mannings_n";
        sel = indices,
        defaults = 0.072,
        type = Float,
    )
    # manning roughness at edges
    append!(mannings_n, mannings_n[index_pit]) # copy to ghost nodes
    mannings_n_sq = fill(Float(0), n_edges)
    zb_max = fill(Float(0), n_edges)
    for i in 1:n_edges
        src_node = nodes_at_link.src[i]
        dst_node = nodes_at_link.dst[i]
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

@get_units @grid_loc @with_kw struct FloodPlainVariables{T}
    volume::Vector{T} | "m3"                            # volume
    h::Vector{T} | "m"                                  # water depth
    h_av::Vector{T} | "m"                               # average water depth
    error::Vector{T} | "m3"                             # error volume
    a::Vector{T} | "m2" | "edge"                        # flow area
    r::Vector{T} | "m" | "edge"                         # hydraulic radius
    hf::Vector{T} | "m" | "edge"                        # water depth at edge/link
    q0::Vector{T} | "m3 s-1" | "edge"                   # discharge at previous time step
    q::Vector{T} | "m3 s-1" | "edge"                    # discharge
    q_av::Vector{T} | "m" | "edge"                      # average river discharge
    hf_index::Vector{Int} | "-" | "edge"                # index with `hf` above depth threshold
end

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
    nodes_at_link,
)
    n = length(indices)
    parameters = FloodPlainParameters(
        dataset,
        config,
        indices;
        river_width,
        river_length,
        zb_floodplain,
        nodes_at_link,
        n_edges,
        index_pit,
    )
    variables = FloodPlainVariables(n, n_edges, index_pit)

    floodplain = FloodPlain(; parameters, variables)
    return floodplain
end

"""
Update boundary condition lateral inflow `inwater` of a `river` flow model for a single
timestep.
"""
function update_lateral_inflow!(
    model::AbstractRiverFlowModel,
    external_models::NamedTuple,
    river_area,
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
        (get_nonirrigation_returnflow(allocation) .* 0.001 .* river_area) ./ dt
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

"""
Update boundary conditions `runoff` and inflow to a waterbody from land `inflow_waterbody` for
overland flow model `ShallowWaterLand` for a single timestep.
"""
function update_boundary_conditions!(
    model::ShallowWaterLand,
    external_models::NamedTuple,
    network,
    dt,
)
    (; river, soil, subsurface, runoff) = external_models
    (; inflow_waterbody) = model.boundary_conditions
    (; reservoir, lake) = river.boundary_conditions
    (; net_runoff) = soil.variables
    (; net_runoff_river) = runoff.variables

    model.boundary_conditions.runoff .=
        net_runoff ./ 1000.0 .* network.land.area ./ dt .+ get_flux_to_river(subsurface) .+
        net_runoff_river .* network.land.area .* 0.001 ./ dt

    if !isnothing(reservoir) || !isnothing(lake)
        inflow_land = get_inflow_waterbody(river, model)
        inflow_subsurface = get_inflow_waterbody(river, subsurface)

        @. inflow_waterbody[network.river_indices] =
            inflow_land[network.river_indices] + inflow_subsurface[network.river_indices]
    end
    return nothing
end

# For the river kinematic wave, the variable `to_river` can be excluded, because this part
# is added to the river kinematic wave.
get_inflow_waterbody(::SurfaceFlowRiver, model::SurfaceFlowLand) = model.variables.q_av
get_inflow_waterbody(::SurfaceFlowRiver, model::LateralSSF) =
    model.variables.ssf ./ tosecond(basetimestep)

# Exclude subsurface flow for other groundwater components than `LateralSSF`.
get_inflow_waterbody(::Union{SurfaceFlowRiver, ShallowWaterRiver}, model::GroundwaterFlow) =
    model.flow.connectivity.ncell .* 0.0
get_inflow_waterbody(::SurfaceFlowRiver, model) = model.variables.to_river .* 0.0

# For local inertial river routing, `to_river` is included, as water body cells are excluded
# (boundary condition).
get_inflow_waterbody(::ShallowWaterRiver, model::SurfaceFlowLand) =
    model.variables.q_av .+ model.variables.to_river
get_inflow_waterbody(::ShallowWaterRiver, model::LateralSSF) =
    (model.variables.ssf .+ model.variables.to_river) ./ tosecond(basetimestep)

"""
    surface_routing!(model)

Run surface routing (land and river) for a single timestep. Kinematic wave for overland flow
and kinematic wave or local inertial model for river flow.
"""
function surface_routing!(model)
    (; vertical, lateral, network, config, clock) = model
    (; soil, runoff, allocation) = vertical
    (; land, river, subsurface) = lateral

    dt = tosecond(clock.dt)
    # update lateral inflow for kinematic wave overland flow
    update_lateral_inflow!(
        land,
        (; soil, allocation, subsurface),
        network.land.area,
        config,
        dt,
    )
    # run kinematic wave overland flow
    update!(land, network.land, dt)

    # update lateral inflow river flow
    update_lateral_inflow!(
        river,
        (; allocation = river.allocation, runoff, land, subsurface),
        network.river.area,
        network.land.area,
        network.river.land_indices,
        dt,
    )
    update_inflow_waterbody!(river, (; land, subsurface), network.river.land_indices)
    update!(river, network, julian_day(clock.time - clock.dt), dt)
    return nothing
end

"""
    surface_routing!(
        model::Model{N,L,V,R,W,T}
    ) where {N,L<:NamedTuple{<:Any,<:Tuple{Any,ShallowWaterLand,ShallowWaterRiver}},V,R,W,T}

Run surface routing (land and river) for a model type that contains the lateral components
`ShallowWaterLand` and `ShallowWaterRiver` for a single timestep.
"""
function surface_routing!(
    model::Model{N, L, V, R, W, T},
) where {
    N,
    L <: NamedTuple{<:Any, <:Tuple{Any, ShallowWaterLand, ShallowWaterRiver}},
    V,
    R,
    W,
    T,
}
    (; lateral, vertical, network, clock) = model
    (; land, river, subsurface) = lateral
    (; soil, runoff) = vertical

    dt = tosecond(clock.dt)
    update_boundary_conditions!(land, (; river, subsurface, soil, runoff), network, dt)

    update!(land, river, network, julian_day(clock.time - clock.dt), dt)

    return nothing
end
