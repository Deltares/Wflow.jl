@get_units @grid_loc @with_kw struct FlowVariables{T}
    q::Vector{T} | "m3 s-1"                 # Discharge [m³ s⁻¹]
    qin::Vector{T} | "m3 s-1"               # Inflow from upstream cells [m³ s⁻¹]
    q_av::Vector{T} | "m3 s-1"              # Average discharge [m³ s⁻¹]
    volume::Vector{T} | "m3"                # Kinematic wave volume [m³] (based on water depth h)
    h::Vector{T} | "m"                      # Water depth [m]
    h_av::Vector{T} | "m"                   # Average water depth [m]
    celerity::Vector{T} | "m s-1"           # Celerity of the kinematic wave
end

function FlowVariables(n)
    variables = FlowVariables(;
        q = zeros(Float, n),
        qin = zeros(Float, n),
        q_av = zeros(Float, n),
        volume = zeros(Float, n),
        h = zeros(Float, n),
        h_av = zeros(Float, n),
        celerity = zeros(Float, n),
    )
    return variables
end

@get_units @grid_loc @with_kw struct ManningFlowParameters{T}
    beta::T                                 # constant in Manning's equation [-]
    slope::Vector{T} | "m m-1"              # Slope [m m⁻¹]
    mannings_n::Vector{T} | "s m-1/3"       # Manning's roughness [s m⁻⅓]
    dl::Vector{T} | "m"                     # Drain length [m]
    width::Vector{T} | "m"                  # Flow width [m]
    alpha_pow::T                            # Used in the power part of alpha [-]
    alpha_term::Vector{T} | "-"             # Term used in computation of alpha [-]
    alpha::Vector{T} | "s3/5 m1/5"          # Constant in momentum equation A = alpha*Q^beta, based on Manning's equation
    dt::T                                   # Model time step [s]
    its::Int                                # Number of fixed iterations [-]
    kinwave_it::Bool                        # Boolean for iterations kinematic wave
end

function ManningFlowParameters(slope, mannings_n, dl, width, iterate, dt, tstep)
    n = length(slope)
    parameters = ManningFlowParameters(;
        beta = Float(0.6),
        slope,
        mannings_n,
        dl = dl,
        dt = Float(tosecond(dt)),
        its = tstep > 0 ? Int(cld(tosecond(dt), tstep)) : tstep,
        width,
        alpha_pow = Float((2.0 / 3.0) * 0.6),
        alpha_term = fill(mv, n),
        alpha = fill(mv, n),
        kinwave_it = iterate,
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
    else
        getfield(getfield(v, :flow), s)
    end
end

function RiverFlowParameters(nc, config, inds, dl, width, iterate, dt, tstep)
    n_river = ncread(
        nc,
        config,
        "lateral.river.mannings_n";
        sel = inds,
        defaults = 0.036,
        type = Float,
    )
    bankfull_depth = ncread(
        nc,
        config,
        "lateral.river.bankfull_depth";
        alias = "lateral.river.h_bankfull",
        sel = inds,
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
        nc,
        config,
        "lateral.river.slope";
        optional = false,
        sel = inds,
        type = Float,
    )
    clamp!(slope, 0.00001, Inf)

    flow_parameter_set =
        ManningFlowParameters(slope, n_river, dl, width, iterate, dt, tstep)
    parameters =
        RiverFlowParameters(; flow = flow_parameter_set, bankfull_depth = bankfull_depth)
    return parameters
end

@get_units @grid_loc @with_kw struct RiverFlowBC{T, R, L}
    qlat::Vector{T} | "m2 s-1"              # Lateral inflow per unit length [m² s⁻¹]
    inwater::Vector{T} | "m3 s-1"           # Lateral inflow [m³ s⁻¹]
    inflow::Vector{T} | "m3 s-1"            # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    inflow_wb::Vector{T} | "m3 s-1"         # inflow waterbody (lake or reservoir model) from land part [m³ s⁻¹]
    abstraction::Vector{T} | "m3 s-1"       # Abstraction (computed as part of water demand and allocation) [m³ s⁻¹]
    reservoir::R                            # Reservoir model struct of arrays
    lake::L                                 # Lake model struct of arrays
    reservoir_index::Vector{Int} | "-"      # map cell to 0 (no reservoir) or i (pick reservoir i in reservoir field)
    lake_index::Vector{Int} | "-"           # map cell to 0 (no lake) or i (pick lake i in lake field)
end

function RiverFlowBC(n, reservoir, reservoir_index, lake, lake_index)
    bc = RiverFlowBC(;
        qlat = zeros(Float, n),
        inwater = zeros(Float, n),
        inflow = zeros(Float, n),
        inflow_wb = zeros(Float, n),
        abstraction = zeros(Float, n),
        reservoir = reservoir,
        lake = lake,
        reservoir_index = reservoir_index,
        lake_index = lake_index,
    )
    return bc
end

@with_kw struct SurfaceFlowRiver{T, R, L, A}
    boundary_conditions::RiverFlowBC{T, R, L}
    parameters::RiverFlowParameters{T}
    variables::FlowVariables{T}
    allocation::A   # Water allocation
end

function SurfaceFlowRiver(
    nc,
    config,
    inds;
    dl,
    width,
    reservoir_index,
    reservoir,
    lake_index,
    lake,
    iterate,
    tstep,
    dt,
)
    @info "Kinematic wave approach is used for river flow." iterate
    if tstep > 0
        @info "Using a fixed sub-timestep (seconds) $tstep for kinematic wave river flow."
    end

    do_water_demand = haskey(config.model, "water_demand")
    allocation = do_water_demand ? AllocationRiver(n) : nothing

    n = length(inds)

    variables = FlowVariables(n)
    parameters = RiverFlowParameters(nc, config, inds, dl, width, iterate, dt, tstep)
    bc = RiverFlowBC(n, reservoir, reservoir_index, lake, lake_index)

    sf_river = SurfaceFlowRiver(;
        boundary_conditions = bc,
        parameters = parameters,
        variables = variables,
        allocation = allocation,
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
    else
        getfield(getfield(v, :flow), s)
    end
end

@get_units @grid_loc @with_kw struct LandFlowBC{T}
    qlat::Vector{T} | "m2 s-1"          # Lateral inflow per unit length [m² s⁻¹]
    inwater::Vector{T} | "m3 s-1"       # Lateral inflow [m³ s⁻¹]
end

@get_units @grid_loc @with_kw struct SurfaceFlowLand{T}
    boundary_conditions::LandFlowBC{T}
    parameters::ManningFlowParameters{T}
    variables::LandFlowVariables{T}
end

function SurfaceFlowLand(nc, config, inds; sl, dl, width, iterate, tstep, dt)
    @info "Kinematic wave approach is used for overland flow." iterate
    if tstep > 0
        @info "Using a fixed sub-timestep (seconds) $tstep for kinematic wave overland flow."
    end

    n_land = ncread(
        nc,
        config,
        "lateral.land.mannings_n";
        sel = inds,
        defaults = 0.072,
        type = Float,
    )
    n = length(inds)

    variables = LandFlowVariables(; flow = FlowVariables(n), to_river = zeros(Float, n))
    parameters = ManningFlowParameters(sl, n_land, dl, width, iterate, dt, tstep)
    bc = LandFlowBC(; qlat = zeros(Float, n), inwater = zeros(Float, n))
    sf_land = SurfaceFlowLand(;
        boundary_conditions = bc,
        variables = variables,
        parameters = parameters,
    )

    return sf_land
end

function update!(sf::SurfaceFlowLand, network, frac_toriver)
    (; subdomain_order, topo_subdomain, indices_subdomain, upstream_nodes) = network

    (; inwater, qlat) = sf.boundary_conditions
    (; alpha_term, mannings_n, slope, beta, alpha_pow, alpha, width, dl) = sf.parameters
    (; h, h_av, q, q_av, qin, volume, to_river) = sf.variables

    ns = length(subdomain_order)

    @. alpha_term = pow(mannings_n / sqrt(slope), beta)
    # use fixed alpha value based flow width
    @. alpha = alpha_term * pow(width, alpha_pow)
    @. qlat = inwater / dl

    q_av .= 0.0
    h_av .= 0.0
    to_river .= 0.0

    dt, its = stable_timestep(sf)
    for _ in 1:its
        qin .= 0.0
        for k in 1:ns
            threaded_foreach(eachindex(subdomain_order[k]); basesize = 1) do i
                m = subdomain_order[k][i]
                for (n, v) in zip(indices_subdomain[m], topo_subdomain[m])
                    # for a river cell without a reservoir or lake part of the upstream
                    # surface flow goes to the river (frac_toriver) and part goes to the
                    # surface flow reservoir (1.0 - frac_toriver), upstream nodes with a
                    # reservoir or lake are excluded
                    to_river[v] += sum_at(
                        i -> q[i] * frac_toriver[i],
                        upstream_nodes[n],
                        eltype(to_river),
                    )
                    if width[v] > 0.0
                        qin[v] = sum_at(
                            i -> q[i] * (1.0 - frac_toriver[i]),
                            upstream_nodes[n],
                            eltype(q),
                        )
                    end

                    q[v] = kinematic_wave(qin[v], q[v], qlat[v], alpha[v], beta, dt, dl[v])

                    # update h, only if surface width > 0.0
                    if width[v] > 0.0
                        crossarea = alpha[v] * pow(q[v], beta)
                        h[v] = crossarea / width[v]
                    end
                    q_av[v] += q[v]
                    h_av[v] += h[v]
                end
            end
        end
    end
    q_av ./= its
    h_av ./= its
    to_river ./= its
    volume .= dl .* width .* h
    return nothing
end

function update!(sf::SurfaceFlowRiver, network, doy)
    (; graph, subdomain_order, topo_subdomain, indices_subdomain, upstream_nodes) = network

    (;
        reservoir,
        reservoir_index,
        lake,
        lake_index,
        inwater,
        qlat,
        inflow,
        abstraction,
        inflow_wb,
    ) = sf.boundary_conditions

    (; alpha_term, mannings_n, slope, beta, alpha_pow, alpha, width, dl, bankfull_depth) =
        sf.parameters
    (; h, h_av, q, q_av, qin, volume) = sf.variables

    ns = length(subdomain_order)

    @. alpha_term = pow(mannings_n / sqrt(slope), beta)
    # use fixed alpha value based on 0.5 * bankfull_depth
    @. alpha = alpha_term * pow(width + bankfull_depth, alpha_pow)

    # move qlat to boundary conditions function?
    @. qlat = inwater / dl

    q_av .= 0.0
    h_av .= 0.0
    # because of possible iterations set reservoir and lake inflow and total outflow at
    # start to zero, the total sum of inflow and outflow at each sub time step is calculated
    if !isnothing(reservoir)
        reservoir.inflow .= 0.0
        reservoir.totaloutflow .= 0.0
        reservoir.actevap .= 0.0
    end
    if !isnothing(lake)
        lake.inflow .= 0.0
        lake.totaloutflow .= 0.0
        lake.actevap .= 0.0
    end

    dt, its = stable_timestep(sf)
    for _ in 1:its
        qin .= 0.0
        for k in 1:ns
            threaded_foreach(eachindex(subdomain_order[k]); basesize = 1) do i
                m = subdomain_order[k][i]
                for (n, v) in zip(indices_subdomain[m], topo_subdomain[m])
                    # sf.qin by outflow from upstream reservoir or lake location is added
                    qin[v] += sum_at(q, upstream_nodes[n])
                    # Inflow supply/abstraction is added to qlat (divide by flow length)
                    # If inflow < 0, abstraction is limited
                    if inflow[v] < 0.0
                        max_abstract =
                            min((inwater[v] + qin[v] + volume[v] / dt) * 0.80, -inflow[v])
                        _inflow = -max_abstract / sf.dl[v]
                    else
                        _inflow = inflow[v] / dl[v]
                    end
                    _inflow -= abstraction[v] / dl[v]

                    q[v] = kinematic_wave(
                        qin[v],
                        q[v],
                        qlat[v] + _inflow,
                        alpha[v],
                        beta,
                        dt,
                        dl[v],
                    )

                    if !isnothing(reservoir) && reservoir_index[v] != 0
                        # run reservoir model and copy reservoir outflow to inflow (qin) of
                        # downstream river cell
                        i = reservoir_index[v]
                        update!(reservoir, i, q[v] + inflow_wb[v], dt)

                        downstream_nodes = outneighbors(graph, v)
                        n_downstream = length(downstream_nodes)
                        if n_downstream == 1
                            j = only(downstream_nodes)
                            qin[j] = reservoir.outflow[i]
                        elseif n_downstream == 0
                            error(
                                """A reservoir without a downstream river node is not supported.
                                Add a downstream river node or move the reservoir to an upstream node (model schematization).
                                """,
                            )
                        else
                            error("bifurcations not supported")
                        end

                    elseif !isnothing(lake) && lake_index[v] != 0
                        # run lake model and copy lake outflow to inflow (qin) of downstream river
                        # cell
                        i = lake_index[v]
                        update!(lake, i, q[v] + inflow_wb[v], doy, dt)

                        downstream_nodes = outneighbors(graph, v)
                        n_downstream = length(downstream_nodes)
                        if n_downstream == 1
                            j = only(downstream_nodes)
                            qin[j] = lake.outflow[i]
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
                    h[v] = crossarea / width[v]
                    volume[v] = dl[v] * width[v] * h[v]
                    q_av[v] += q[v]
                    h_av[v] += h[v]
                end
            end
        end
    end
    q_av ./= its
    h_av ./= its
    volume .= dl .* width .* h
    return nothing
end

function stable_timestep(sf::S) where {S <: Union{SurfaceFlowLand, SurfaceFlowRiver}}
    (; q, celerity) = sf.variables
    (; alpha, beta, dl, dt, its, kinwave_it) = sf.parameters

    n = length(q)
    # two options for iteration, fixed or based on courant number.
    if kinwave_it
        if its > 0
            n_iterations = its
        else
            # calculate celerity
            courant = zeros(n)
            for v in 1:n
                if q[v] > 0.0
                    celerity[v] = 1.0 / (alpha[v] * beta * pow(q[v], (beta - 1.0)))
                    courant[v] = (dt / dl[v]) * celerity[v]
                end
            end
            filter!(x -> x ≠ 0.0, courant)
            n_iterations =
                isempty(courant) ? 1 : ceil(Int, (1.25 * quantile!(courant, 0.95)))
        end
    else
        n_iterations = 1
    end

    # sub time step
    dt_sub = dt / n_iterations
    return dt_sub, n_iterations
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

@get_units @grid_loc @with_kw struct LateralSSF{T, Kh} <: SubsurfaceFlow
    kh_profile::Kh                         # Horizontal hydraulic conductivity profile type [-]  
    khfrac::Vector{T} | "-"                # A muliplication factor applied to vertical hydraulic conductivity `kv` [-]
    soilthickness::Vector{T} | "m"         # Soil thickness [m]
    theta_s::Vector{T} | "-"               # Saturated water content (porosity) [-]
    theta_r::Vector{T} | "-"               # Residual water content [-]
    dt::T                                  # model time step [d]
    slope::Vector{T} | "m m-1"             # Slope [m m⁻¹]
    dl::Vector{T} | "m"                    # Drain length [m]
    dw::Vector{T} | "m"                    # Flow width [m]
    zi::Vector{T} | "m"                    # Pseudo-water table depth [m] (top of the saturated zone)
    exfiltwater::Vector{T} | "m dt-1"      # Exfiltration [m Δt⁻¹] (groundwater above surface level, saturated excess conditions)
    recharge::Vector{T} | "m2 dt-1"        # Net recharge to saturated store [m² Δt⁻¹]
    ssf::Vector{T} | "m3 d-1"              # Subsurface flow [m³ d⁻¹]
    ssfin::Vector{T} | "m3 d-1"            # Inflow from upstream cells [m³ d⁻¹]
    ssfmax::Vector{T} | "m2 d-1"           # Maximum subsurface flow [m² d⁻¹]
    to_river::Vector{T} | "m3 d-1"         # Part of subsurface flow [m³ d⁻¹] that flows to the river
    volume::Vector{T} | "m3"               # Subsurface volume [m³]

    function LateralSSF{T, Kh}(args...) where {T, Kh}
        equal_size_vectors(args)
        return new(args...)
    end
end

function update!(ssf::LateralSSF, network, frac_toriver)
    (; subdomain_order, topo_subdomain, indices_subdomain, upstream_nodes, area) = network

    ns = length(subdomain_order)
    for k in 1:ns
        threaded_foreach(eachindex(subdomain_order[k]); basesize = 1) do i
            m = subdomain_order[k][i]
            for (n, v) in zip(indices_subdomain[m], topo_subdomain[m])
                # for a river cell without a reservoir or lake part of the upstream
                # subsurface flow goes to the river (frac_toriver) and part goes to the
                # subsurface flow reservoir (1.0 - frac_toriver) upstream nodes with a
                # reservoir or lake are excluded
                ssf.ssfin[v] = sum_at(
                    i -> ssf.ssf[i] * (1.0 - frac_toriver[i]),
                    upstream_nodes[n],
                    eltype(ssf.ssfin),
                )
                ssf.to_river[v] = sum_at(
                    i -> ssf.ssf[i] * frac_toriver[i],
                    upstream_nodes[n],
                    eltype(ssf.to_river),
                )
                ssf.ssf[v], ssf.zi[v], ssf.exfiltwater[v] = kinematic_wave_ssf(
                    ssf.ssfin[v],
                    ssf.ssf[v],
                    ssf.zi[v],
                    ssf.recharge[v],
                    ssf.slope[v],
                    ssf.theta_s[v] - ssf.theta_r[v],
                    ssf.soilthickness[v],
                    ssf.dt,
                    ssf.dl[v],
                    ssf.dw[v],
                    ssf.ssfmax[v],
                    ssf.kh_profile,
                    v,
                )
                ssf.volume[v] =
                    (ssf.theta_s[v] - ssf.theta_r[v]) *
                    (ssf.soilthickness[v] - ssf.zi[v]) *
                    area[v]
            end
        end
    end
    return nothing
end

@get_units@grid_loc @with_kw struct GroundwaterExchange{T} <: SubsurfaceFlow
    dt::T                               # model time step [d]
    exfiltwater::Vector{T} | "m dt-1"   # Exfiltration [m Δt⁻¹]  (groundwater above surface level, saturated excess conditions)
    zi::Vector{T} | "m"                 # Pseudo-water table depth [m] (top of the saturated zone)
    to_river::Vector{T} | "m3 d-1"      # Part of subsurface flow [m³ d⁻¹] that flows to the river
    ssf::Vector{T} | "m3 d-1"           # Subsurface flow [m³ d⁻¹]
end

get_water_depth(subsurface::SubsurfaceFlow) = subsurface.zi
get_exfiltwater(subsurface::SubsurfaceFlow) = subsurface.exfiltwater

@get_units @grid_loc @with_kw struct ShallowWaterRiver{T, R, L, F, A}
    n::Int                                              # number of cells [-]
    ne::Int                                             # number of edges/links [-]
    active_n::Vector{Int} | "-"                         # active nodes [-]
    active_e::Vector{Int} | "-" | "edge"                # active edges/links [-]
    g::T                                                # acceleration due to gravity [m s⁻²]
    alpha::T                                            # stability coefficient (Bates et al., 2010) [-]
    h_thresh::T                                         # depth threshold for calculating flow [m]
    dt::T                                               # model time step [s]
    q::Vector{T} | "m3 s-1" | "edge"                    # river discharge (subgrid channel)
    q0::Vector{T} | "m3 s-1" | "edge"                   # river discharge (subgrid channel) at previous time step
    q_av::Vector{T} | "m3 s-1" | "edge"                 # average river channel (+ floodplain) discharge [m³ s⁻¹]
    q_channel_av::Vector{T} | "m3 s-1"                  # average river channel discharge [m³ s⁻¹]
    zb_max::Vector{T} | "m"                             # maximum channel bed elevation
    mannings_n_sq::Vector{T} | "(s m-1/3)2" | "edge"    # Manning's roughness squared at edge/link
    mannings_n::Vector{T} | "s m-1/3"                   # Manning's roughness at node
    h::Vector{T} | "m"                                  # water depth
    zs_max::Vector{T} | "m" | "edge"                    # maximum water elevation at edge
    zs_src::Vector{T} | "m"                             # water elevation of source node of edge
    zs_dst::Vector{T} | "m"                             # water elevation of downstream node of edge
    hf::Vector{T} | "m" | "edge"                        # water depth at edge/link
    h_av::Vector{T} | "m"                               # average water depth
    dl::Vector{T} | "m"                                 # river length
    dl_at_link::Vector{T} | "m" | "edge"                # river length at edge/link
    width::Vector{T} | "m"                              # river width
    width_at_link::Vector{T} | "m" | "edge"             # river width at edge/link
    a::Vector{T} | "m2" | "edge"                        # flow area at edge/link
    r::Vector{T} | "m" | "edge"                         # wetted perimeter at edge/link
    volume::Vector{T} | "m3"                            # river volume
    error::Vector{T} | "m3"                             # error volume
    inwater::Vector{T} | "m3 s-1"                       # lateral inflow [m³ s⁻¹]
    inflow::Vector{T} | "m3 s-1"                        # external inflow (abstraction/supply/demand) [m³ s⁻¹]
    abstraction::Vector{T} | "m3 s-1"                   # abstraction (computed as part of water demand and allocation) [m³ s⁻¹]
    inflow_wb::Vector{T} | "m3 s-1"                     # inflow waterbody (lake or reservoir model) from land part [m³ s⁻¹]
    bankfull_volume::Vector{T} | "m3"                   # bankfull volume
    bankfull_depth::Vector{T} | "m"                     # bankfull depth
    zb::Vector{T} | "m"                                 # river bed elevation
    froude_limit::Bool                                  # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
    reservoir_index::Vector{Int} | "-"                  # river cell index with a reservoir (each index of reservoir_index maps to reservoir i in reservoir field)
    lake_index::Vector{Int} | "-"                       # river cell index with a lake (each index of lake_index maps to lake i in lake field)
    waterbody::Vector{Bool} | "-"                       # water body cells (reservoir or lake)
    reservoir::R                                        # Reservoir model struct of arrays
    lake::L                                             # Lake model struct of arrays
    floodplain::F                                       # Floodplain (1D) schematization
    allocation::A                                       # Water allocation
end

function initialize_shallowwater_river(
    nc,
    config,
    inds;
    graph,
    ldd,
    dl,
    width,
    reservoir_index,
    reservoir,
    lake_index,
    lake,
    dt,
    floodplain,
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
    alpha = get(config.model, "inertial_flow_alpha", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    h_thresh = get(config.model, "h_thresh", 1.0e-03)::Float64 # depth threshold for flow at link
    froude_limit = get(config.model, "froude_limit", true)::Bool # limit flow to subcritical according to Froude number
    floodplain_1d = floodplain

    @info "Local inertial approach is used for river flow." alpha h_thresh froude_limit floodplain_1d
    @warn string(
        "Providing the boundary condition `riverlength_bc` as part of the `[model]` setting ",
        "in the TOML file has been deprecated as of Wflow v0.8.0.\n The boundary condition should ",
        "be provided as part of the file `$(config.input.path_static)`.",
    )
    # The following boundary conditions can be set at ghost nodes, downstream of river
    # outlets (pits): river length and river depth
    index_pit = findall(x -> x == 5, ldd)
    inds_pit = inds[index_pit]
    riverlength_bc = ncread(
        nc,
        config,
        "lateral.river.riverlength_bc";
        sel = inds_pit,
        defaults = 1.0e04,
        type = Float,
    )
    riverdepth_bc = ncread(
        nc,
        config,
        "lateral.river.riverdepth_bc";
        sel = inds_pit,
        defaults = 0.0,
        type = Float,
    )
    bankfull_elevation_2d = ncread(
        nc,
        config,
        "lateral.river.bankfull_elevation";
        optional = false,
        type = Float,
        fill = 0,
    )
    bankfull_depth_2d = ncread(
        nc,
        config,
        "lateral.river.bankfull_depth";
        optional = false,
        type = Float,
        fill = 0,
    )
    bankfull_depth = bankfull_depth_2d[inds]
    zb = bankfull_elevation_2d[inds] - bankfull_depth # river bed elevation

    bankfull_volume = bankfull_depth .* width .* dl

    n_river =
        ncread(nc, config, "lateral.river.n"; sel = inds, defaults = 0.036, type = Float)

    n = length(inds)

    # set river depth h to zero (including reservoir and lake locations)
    h = fill(0.0, n)

    # set ghost points for boundary condition (downstream river outlet): river width, bed
    # elevation, manning n is copied from the upstream cell.
    add_vertex_edge_graph!(graph, index_pit)
    append!(dl, riverlength_bc)
    append!(h, riverdepth_bc)
    append!(zb, zb[index_pit])
    append!(width, width[index_pit])
    append!(n_river, n_river[index_pit])
    append!(bankfull_depth, bankfull_depth[index_pit])

    # for each link the src and dst node is required
    nodes_at_link = adjacent_nodes_at_link(graph)
    _ne = ne(graph)

    if floodplain
        zb_floodplain = zb .+ bankfull_depth
        floodplain = initialize_floodplain_1d(
            nc,
            config,
            inds,
            width,
            dl,
            zb_floodplain,
            index_pit,
            _ne,
            nodes_at_link,
        )
    else
        floodplain = nothing
    end

    # determine z, width, length and manning's n at links
    zb_max = fill(Float(0), _ne)
    width_at_link = fill(Float(0), _ne)
    length_at_link = fill(Float(0), _ne)
    mannings_n_sq = fill(Float(0), _ne)
    for i in 1:_ne
        src_node = nodes_at_link.src[i]
        dst_node = nodes_at_link.dst[i]
        zb_max[i] = max(zb[src_node], zb[dst_node])
        width_at_link[i] = min(width[src_node], width[dst_node])
        length_at_link[i] = 0.5 * (dl[dst_node] + dl[src_node])
        mannings_n =
            (n_river[dst_node] * dl[dst_node] + n_river[src_node] * dl[src_node]) /
            (dl[dst_node] + dl[src_node])
        mannings_n_sq[i] = mannings_n * mannings_n
    end

    q_av = zeros(_ne)
    waterbody = !=(0).(reservoir_index .+ lake_index)
    active_index = findall(x -> x == 0, waterbody)

    do_water_demand = haskey(config.model, "water_demand")
    sw_river = ShallowWaterRiver(;
        n = n,
        ne = _ne,
        active_n = active_index,
        active_e = active_index,
        g = 9.80665,
        alpha = alpha,
        h_thresh = h_thresh,
        dt = tosecond(dt),
        q = zeros(_ne),
        q0 = zeros(_ne),
        q_av = q_av,
        q_channel_av = isnothing(floodplain) ? q_av : zeros(_ne),
        zb_max = zb_max,
        mannings_n_sq = mannings_n_sq,
        mannings_n = n_river,
        h = h,
        zs_max = zeros(_ne),
        zs_src = zeros(_ne),
        zs_dst = zeros(_ne),
        hf = zeros(_ne),
        h_av = zeros(n),
        width = width,
        width_at_link = width_at_link,
        a = zeros(_ne),
        r = zeros(_ne),
        volume = fill(0.0, n),
        error = zeros(n),
        inflow = zeros(n),
        abstraction = zeros(n),
        inflow_wb = zeros(n),
        inwater = zeros(n),
        dl = dl,
        dl_at_link = length_at_link,
        bankfull_volume = bankfull_volume,
        bankfull_depth = bankfull_depth,
        zb = zb,
        froude_limit = froude_limit,
        reservoir_index = findall(x -> x > 0, reservoir_index),
        lake_index = findall(x -> x > 0, lake_index),
        waterbody = waterbody,
        reservoir = reservoir,
        lake = lake,
        floodplain = floodplain,
        allocation = do_water_demand ? initialize_allocation_river(n) : nothing,
    )
    return sw_river, nodes_at_link
end

"Return the upstream inflow for a waterbody in `ShallowWaterRiver`"
function get_inflow_waterbody(sw::ShallowWaterRiver, src_edge)
    q_in = sum_at(sw.q, src_edge)
    if !isnothing(sw.floodplain)
        q_in = q_in + sum_at(sw.floodplain.q, src_edge)
    end
    return q_in
end

function shallowwater_river_update!(sw::ShallowWaterRiver, network, dt, doy, update_h)
    (; nodes_at_link, links_at_node) = network

    sw.q0 .= sw.q
    if !isnothing(sw.floodplain)
        sw.floodplain.q0 .= sw.floodplain.q
    end
    @tturbo for j in eachindex(sw.active_e)
        i = sw.active_e[j]
        i_src = nodes_at_link.src[i]
        i_dst = nodes_at_link.dst[i]
        sw.zs_src[i] = sw.zb[i_src] + sw.h[i_src]
        sw.zs_dst[i] = sw.zb[i_dst] + sw.h[i_dst]

        sw.zs_max[i] = max(sw.zs_src[i], sw.zs_dst[i])
        sw.hf[i] = (sw.zs_max[i] - sw.zb_max[i])

        sw.a[i] = sw.width_at_link[i] * sw.hf[i] # flow area (rectangular channel)
        sw.r[i] = sw.a[i] / (sw.width_at_link[i] + 2.0 * sw.hf[i]) # hydraulic radius (rectangular channel)

        sw.q[i] = IfElse.ifelse(
            sw.hf[i] > sw.h_thresh,
            local_inertial_flow(
                sw.q0[i],
                sw.zs_src[i],
                sw.zs_dst[i],
                sw.hf[i],
                sw.a[i],
                sw.r[i],
                sw.dl_at_link[i],
                sw.mannings_n_sq[i],
                sw.g,
                sw.froude_limit,
                dt,
            ),
            0.0,
        )

        # limit q in case water is not available
        sw.q[i] = IfElse.ifelse(sw.h[i_src] <= 0.0, min(sw.q[i], 0.0), sw.q[i])
        sw.q[i] = IfElse.ifelse(sw.h[i_dst] <= 0.0, max(sw.q[i], 0.0), sw.q[i])

        sw.q_av[i] += sw.q[i] * dt
    end
    if !isnothing(sw.floodplain)
        @tturbo @. sw.floodplain.hf = max(sw.zs_max - sw.floodplain.zb_max, 0.0)

        n = 0
        @inbounds for i in sw.active_e
            @inbounds if sw.floodplain.hf[i] > sw.h_thresh
                n += 1
                sw.floodplain.hf_index[n] = i
            else
                sw.floodplain.q[i] = 0.0
            end
        end

        @tturbo for j in 1:n
            i = sw.floodplain.hf_index[j]
            i_src = nodes_at_link.src[i]
            i_dst = nodes_at_link.dst[i]

            i0 = 0
            for k in eachindex(sw.floodplain.profile.depth)
                i0 += 1 * (sw.floodplain.profile.depth[k] <= sw.floodplain.hf[i])
            end
            i1 = max(i0, 1)
            i2 = IfElse.ifelse(i1 == length(sw.floodplain.profile.depth), i1, i1 + 1)

            a_src = flow_area(
                sw.floodplain.profile.width[i2, i_src],
                sw.floodplain.profile.a[i1, i_src],
                sw.floodplain.profile.depth[i1],
                sw.floodplain.hf[i],
            )
            a_src = max(a_src - (sw.floodplain.hf[i] * sw.width[i_src]), 0.0)

            a_dst = flow_area(
                sw.floodplain.profile.width[i2, i_dst],
                sw.floodplain.profile.a[i1, i_dst],
                sw.floodplain.profile.depth[i1],
                sw.floodplain.hf[i],
            )
            a_dst = max(a_dst - (sw.floodplain.hf[i] * sw.width[i_dst]), 0.0)

            sw.floodplain.a[i] = min(a_src, a_dst)

            sw.floodplain.r[i] = IfElse.ifelse(
                a_src < a_dst,
                a_src / wetted_perimeter(
                    sw.floodplain.profile.p[i1, i_src],
                    sw.floodplain.profile.depth[i1],
                    sw.floodplain.hf[i],
                ),
                a_dst / wetted_perimeter(
                    sw.floodplain.profile.p[i1, i_dst],
                    sw.floodplain.profile.depth[i1],
                    sw.floodplain.hf[i],
                ),
            )

            sw.floodplain.q[i] = IfElse.ifelse(
                sw.floodplain.a[i] > 1.0e-05,
                local_inertial_flow(
                    sw.floodplain.q0[i],
                    sw.zs_src[i],
                    sw.zs_dst[i],
                    sw.floodplain.hf[i],
                    sw.floodplain.a[i],
                    sw.floodplain.r[i],
                    sw.dl_at_link[i],
                    sw.floodplain.mannings_n_sq[i],
                    sw.g,
                    sw.froude_limit,
                    dt,
                ),
                0.0,
            )

            # limit floodplain q in case water is not available
            sw.floodplain.q[i] = IfElse.ifelse(
                sw.floodplain.h[i_src] <= 0.0,
                min(sw.floodplain.q[i], 0.0),
                sw.floodplain.q[i],
            )
            sw.floodplain.q[i] = IfElse.ifelse(
                sw.floodplain.h[i_dst] <= 0.0,
                max(sw.floodplain.q[i], 0.0),
                sw.floodplain.q[i],
            )

            sw.floodplain.q[i] =
                IfElse.ifelse(sw.floodplain.q[i] * sw.q[i] < 0.0, 0.0, sw.floodplain.q[i])
            sw.floodplain.q_av[i] += sw.floodplain.q[i] * dt
        end
    end
    # For reservoir and lake locations the local inertial solution is replaced by the
    # reservoir or lake model. These locations are handled as boundary conditions in the
    # local inertial model (fixed h).
    for v in eachindex(sw.reservoir_index)
        i = sw.reservoir_index[v]

        q_in = get_inflow_waterbody(sw, links_at_node.src[i])
        update!(sw.reservoir, v, q_in + sw.inflow_wb[i], dt)
        sw.q[i] = sw.reservoir.outflow[v]
        sw.q_av[i] += sw.q[i] * dt
    end
    for v in eachindex(sw.lake_index)
        i = sw.lake_index[v]

        q_in = get_inflow_waterbody(sw, links_at_node.src[i])
        update!(sw.lake, v, q_in + sw.inflow_wb[i], doy, dt)
        sw.q[i] = sw.lake.outflow[v]
        sw.q_av[i] += sw.q[i] * dt
    end
    if update_h
        @batch per = thread minbatch = 2000 for i in sw.active_n
            q_src = sum_at(sw.q, links_at_node.src[i])
            q_dst = sum_at(sw.q, links_at_node.dst[i])
            sw.volume[i] =
                sw.volume[i] + (q_src - q_dst + sw.inwater[i] - sw.abstraction[i]) * dt

            if sw.volume[i] < 0.0
                sw.error[i] = sw.error[i] + abs(sw.volume[i])
                sw.volume[i] = 0.0 # set volume to zero
            end
            sw.volume[i] = max(sw.volume[i] + sw.inflow[i] * dt, 0.0) # add external inflow

            if !isnothing(sw.floodplain)
                q_src = sum_at(sw.floodplain.q, links_at_node.src[i])
                q_dst = sum_at(sw.floodplain.q, links_at_node.dst[i])
                sw.floodplain.volume[i] = sw.floodplain.volume[i] + (q_src - q_dst) * dt
                # TODO check following approach:
                # if floodplain volume negative, extract from river volume first
                if sw.floodplain.volume[i] < 0.0
                    sw.floodplain.error[i] =
                        sw.floodplain.error[i] + abs(sw.floodplain.volume[i])
                    sw.floodplain.volume[i] = 0.0
                end
                volume_total = sw.volume[i] + sw.floodplain.volume[i]
                if volume_total > sw.bankfull_volume[i]
                    flood_volume = volume_total - sw.bankfull_volume[i]
                    h = flood_depth(sw.floodplain.profile, flood_volume, sw.dl[i], i)
                    sw.h[i] = sw.bankfull_depth[i] + h
                    sw.volume[i] = sw.h[i] * sw.width[i] * sw.dl[i]
                    sw.floodplain.volume[i] = max(volume_total - sw.volume[i], 0.0)
                    sw.floodplain.h[i] = sw.floodplain.volume[i] > 0.0 ? h : 0.0
                else
                    sw.h[i] = volume_total / (sw.dl[i] * sw.width[i])
                    sw.volume[i] = volume_total
                    sw.floodplain.h[i] = 0.0
                    sw.floodplain.volume[i] = 0.0
                end
                sw.floodplain.h_av[i] += sw.floodplain.h[i] * dt
            else
                sw.h[i] = sw.volume[i] / (sw.dl[i] * sw.width[i])
            end
            sw.h_av[i] += sw.h[i] * dt
        end
    end
    return nothing
end

function update!(sw::ShallowWaterRiver{T}, network, doy; update_h = true) where {T}
    if !isnothing(sw.reservoir)
        sw.reservoir.inflow .= 0.0
        sw.reservoir.totaloutflow .= 0.0
        sw.reservoir.actevap .= 0.0
    end
    if !isnothing(sw.lake)
        sw.lake.inflow .= 0.0
        sw.lake.totaloutflow .= 0.0
        sw.lake.actevap .= 0.0
    end
    if !isnothing(sw.floodplain)
        sw.floodplain.q_av .= 0.0
        sw.floodplain.h_av .= 0.0
    end
    sw.q_av .= 0.0
    sw.h_av .= 0.0

    t = T(0.0)
    while t < sw.dt
        dt = stable_timestep(sw)
        if t + dt > sw.dt
            dt = sw.dt - t
        end
        shallowwater_river_update!(sw, network, dt, doy, update_h)
        t = t + dt
    end
    sw.q_av ./= sw.dt
    sw.h_av ./= sw.dt

    if !isnothing(sw.floodplain)
        sw.floodplain.q_av ./= sw.dt
        sw.floodplain.h_av ./= sw.dt
        sw.q_channel_av .= sw.q_av
        sw.q_av .= sw.q_channel_av .+ sw.floodplain.q_av
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

@get_units @grid_loc @with_kw struct ShallowWaterLand{T}
    n::Int                                              # number of cells [-]
    xl::Vector{T} | "m"                                 # cell length x direction [m]
    yl::Vector{T} | "m"                                 # cell length y direction [m]
    xwidth::Vector{T} | "m" | "edge"                    # effective flow width x direction (floodplain) [m]
    ywidth::Vector{T} | "m" | "edge"                    # effective flow width y direction (floodplain) [m]
    g::T                                                # acceleration due to gravity [m s⁻²]
    theta::T                                            # weighting factor (de Almeida et al., 2012) [-]
    alpha::T                                            # stability coefficient (de Almeida et al., 2012) [-]
    h_thresh::T                                         # depth threshold for calculating flow [m]
    dt::T                                               # model time step [s]
    qy0::Vector{T} | "m3 s-1" | "edge"                  # flow in y direction at previous time step
    qx0::Vector{T} | "m3 s-1" | "edge"                  # flow in x direction at previous time step
    qx::Vector{T} | "m3 s-1" | "edge"                   # flow in x direction
    qy::Vector{T} | "m3 s-1" | "edge"                   # flow in y direction
    zx_max::Vector{T} | "m" | "edge"                    # maximum cell elevation (x direction)
    zy_max::Vector{T} | "m" | "edge"                    # maximum cell elevation (y direction)
    mannings_n_sq::Vector{T} | "(s m-1/3)2" | "edge"    # Manning's roughness squared
    volume::Vector{T} | "m3"                            # total volume of cell (including river volume for river cells)
    error::Vector{T} | "m3"                             # error volume
    runoff::Vector{T} | "m3 s-1"                        # runoff from hydrological model
    inflow_wb::Vector{T} | "m3 s-1"                     # inflow to water body from hydrological model
    h::Vector{T} | "m"                                  # water depth of cell (for river cells the reference is the river bed elevation `zb`)
    z::Vector{T} | "m"                                  # elevation of cell
    froude_limit::Bool                                  # if true a check is performed if froude number > 1.0 (algorithm is modified) [-]
    rivercells::Vector{Bool} | "-"                      # river cells
    h_av::Vector{T} | "m"                               # average water depth (for river cells the reference is the river bed elevation `zb`)
end

function initialize_shallowwater_land(
    nc,
    config,
    inds;
    modelsize_2d,
    indices_reverse, # maps from the 2D external domain to the 1D internal domain (Int for linear indexing).
    xlength,
    ylength,
    riverwidth,
    graph_riv,
    ldd_riv,
    inds_riv,
    river,
    waterbody,
    dt,
)
    froude_limit = get(config.model, "froude_limit", true)::Bool # limit flow to subcritical according to Froude number
    alpha = get(config.model, "inertial_flow_alpha", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    theta = get(config.model, "inertial_flow_theta", 0.8)::Float64 # weighting factor
    h_thresh = get(config.model, "h_thresh", 1.0e-03)::Float64 # depth threshold for flow at link

    @info "Local inertial approach is used for overlandflow." alpha theta h_thresh froude_limit

    n_land =
        ncread(nc, config, "lateral.land.n"; sel = inds, defaults = 0.072, type = Float)
    elevation_2d = ncread(
        nc,
        config,
        "lateral.land.elevation";
        optional = false,
        type = Float,
        fill = 0,
    )
    elevation = elevation_2d[inds]
    n = length(inds)

    # initialize links between cells in x and y direction.
    indices = Indices(; xu = zeros(n), xd = zeros(n), yu = zeros(n), yd = zeros(n))

    # links without neigbors are handled by an extra index (at n + 1, with n links), which
    # is set to a value of 0.0 m³ s⁻¹ for qx and qy fields at initialization.
    # links are defined as follows for the x and y direction, respectively:
    # node i => node xu (node i + CartesianIndex(1, 0))
    # node i => node yu (node i + CartesianIndex(0, 1))
    # where i is the index of inds
    nrow, ncol = modelsize_2d
    for (v, i) in enumerate(inds)
        for (m, neighbor) in enumerate(neighbors)
            j = i + neighbor
            dir = dirs[m]
            if (1 <= j[1] <= nrow) && (1 <= j[2] <= ncol) && (indices_reverse[j] != 0)
                getfield(indices, dir)[v] = indices_reverse[j]
            else
                getfield(indices, dir)[v] = n + 1
            end
        end
    end

    # determine z at links in x and y direction
    zx_max = fill(Float(0), n)
    zy_max = fill(Float(0), n)
    for i in 1:n
        xu = indices.xu[i]
        if xu <= n
            zx_max[i] = max(elevation[i], elevation[xu])
        end
        yu = indices.yu[i]
        if yu <= n
            zy_max[i] = max(elevation[i], elevation[yu])
        end
    end

    # set the effective flow width for river cells in the x and y direction at cell edges.
    # for waterbody cells (reservoir or lake), h is set to zero (fixed) and not updated, and
    # overland flow from a downstream cell is not possible (effective flowwidth is zero).
    we_x = copy(xlength)
    we_y = copy(ylength)
    set_effective_flowwidth!(
        we_x,
        we_y,
        indices,
        graph_riv,
        riverwidth,
        ldd_riv,
        waterbody,
        indices_reverse[inds_riv],
    )

    sw_land = ShallowWaterLand{Float}(;
        n = n,
        xl = xlength,
        yl = ylength,
        xwidth = we_x,
        ywidth = we_y,
        g = 9.80665,
        theta = theta,
        alpha = alpha,
        h_thresh = h_thresh,
        dt = tosecond(dt),
        qx0 = zeros(n + 1),
        qy0 = zeros(n + 1),
        qx = zeros(n + 1),
        qy = zeros(n + 1),
        zx_max = zx_max,
        zy_max = zy_max,
        mannings_n_sq = n_land .* n_land,
        volume = zeros(n),
        error = zeros(n),
        runoff = zeros(n),
        inflow_wb = zeros(n),
        h = zeros(n),
        h_av = zeros(n),
        z = elevation,
        froude_limit = froude_limit,
        rivercells = river,
    )

    return sw_land, indices
end

"""
    stable_timestep(sw::ShallowWaterRiver)
    stable_timestep(sw::ShallowWaterLand)

Compute a stable timestep size for the local inertial approach, based on Bates et al. (2010).

dt = alpha * (Δx / sqrt(g max(h))
"""
function stable_timestep(sw::ShallowWaterRiver{T})::T where {T}
    dt_min = T(Inf)
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(sw.n)
        @fastmath @inbounds dt = sw.alpha * sw.dl[i] / sqrt(sw.g * sw.h[i])
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? T(10.0) : dt_min
    return dt_min
end

function stable_timestep(sw::ShallowWaterLand{T})::T where {T}
    dt_min = T(Inf)
    @batch per = thread reduction = ((min, dt_min),) for i in 1:(sw.n)
        @fastmath @inbounds dt = if sw.rivercells[i] == 0
            sw.alpha * min(sw.xl[i], sw.yl[i]) / sqrt(sw.g * sw.h[i])
        else
            T(Inf)
        end
        dt_min = min(dt, dt_min)
    end
    dt_min = isinf(dt_min) ? T(10.0) : dt_min
    return dt_min
end

function update!(
    sw::ShallowWaterLand{T},
    swr::ShallowWaterRiver{T},
    network,
    doy;
    update_h = false,
) where {T}
    (; nodes_at_link, links_at_node) = network.river

    if !isnothing(swr.reservoir)
        swr.reservoir.inflow .= 0.0
        swr.reservoir.totaloutflow .= 0.0
        swr.reservoir.actevap .= 0.0
    end
    if !isnothing(swr.lake)
        swr.lake.inflow .= 0.0
        swr.lake.totaloutflow .= 0.0
        swr.lake.actevap .= 0.0
    end
    swr.q_av .= 0.0
    swr.h_av .= 0.0
    sw.h_av .= 0.0

    t = T(0.0)
    while t < swr.dt
        dt_river = stable_timestep(swr)
        dt_land = stable_timestep(sw)
        dt = min(dt_river, dt_land)
        if t + dt > swr.dt
            dt = swr.dt - t
        end
        shallowwater_river_update!(swr, network.river, dt, doy, update_h)
        shallowwater_update!(sw, swr, network, dt)
        t = t + dt
    end
    swr.q_av ./= swr.dt
    swr.h_av ./= swr.dt
    sw.h_av ./= sw.dt

    return nothing
end

function shallowwater_update!(
    sw::ShallowWaterLand{T},
    swr::ShallowWaterRiver{T},
    network,
    dt,
) where {T}
    indices = network.land.staggered_indices
    inds_riv = network.land.index_river

    (; links_at_node) = network.river

    sw.qx0 .= sw.qx
    sw.qy0 .= sw.qy

    # update qx
    @batch per = thread minbatch = 6000 for i in 1:(sw.n)
        yu = indices.yu[i]
        yd = indices.yd[i]
        xu = indices.xu[i]
        xd = indices.xd[i]

        # the effective flow width is zero when the river width exceeds the cell width (dy
        # for flow in x dir) and floodplain flow is not calculated.
        if xu <= sw.n && sw.ywidth[i] != T(0.0)
            zs_x = sw.z[i] + sw.h[i]
            zs_xu = sw.z[xu] + sw.h[xu]
            zs_max = max(zs_x, zs_xu)
            hf = (zs_max - sw.zx_max[i])

            if hf > sw.h_thresh
                length = T(0.5) * (sw.xl[i] + sw.xl[xu]) # can be precalculated
                sw.qx[i] = local_inertial_flow(
                    sw.theta,
                    sw.qx0[i],
                    sw.qx0[xd],
                    sw.qx0[xu],
                    zs_x,
                    zs_xu,
                    hf,
                    sw.ywidth[i],
                    length,
                    sw.mannings_n_sq[i],
                    sw.g,
                    sw.froude_limit,
                    dt,
                )
                # limit qx in case water is not available
                if sw.h[i] <= T(0.0)
                    sw.qx[i] = min(sw.qx[i], T(0.0))
                end
                if sw.h[xu] <= T(0.0)
                    sw.qx[i] = max(sw.qx[i], T(0.0))
                end
            else
                sw.qx[i] = T(0.0)
            end
        end

        # update qy

        # the effective flow width is zero when the river width exceeds the cell width (dx
        # for flow in y dir) and floodplain flow is not calculated.
        if yu <= sw.n && sw.xwidth[i] != T(0.0)
            zs_y = sw.z[i] + sw.h[i]
            zs_yu = sw.z[yu] + sw.h[yu]
            zs_max = max(zs_y, zs_yu)
            hf = (zs_max - sw.zy_max[i])

            if hf > sw.h_thresh
                length = T(0.5) * (sw.yl[i] + sw.yl[yu]) # can be precalculated
                sw.qy[i] = local_inertial_flow(
                    sw.theta,
                    sw.qy0[i],
                    sw.qy0[yd],
                    sw.qy0[yu],
                    zs_y,
                    zs_yu,
                    hf,
                    sw.xwidth[i],
                    length,
                    sw.mannings_n_sq[i],
                    sw.g,
                    sw.froude_limit,
                    dt,
                )
                # limit qy in case water is not available
                if sw.h[i] <= T(0.0)
                    sw.qy[i] = min(sw.qy[i], T(0.0))
                end
                if sw.h[yu] <= T(0.0)
                    sw.qy[i] = max(sw.qy[i], T(0.0))
                end
            else
                sw.qy[i] = T(0.0)
            end
        end
    end

    # change in volume and water levels based on horizontal fluxes for river and land cells
    @batch per = thread minbatch = 6000 for i in 1:(sw.n)
        yd = indices.yd[i]
        xd = indices.xd[i]

        if sw.rivercells[i]
            if swr.waterbody[inds_riv[i]]
                # for reservoir or lake set inflow from land part, these are boundary points
                # and update of volume and h is not required
                swr.inflow_wb[inds_riv[i]] =
                    sw.inflow_wb[i] +
                    sw.runoff[i] +
                    (sw.qx[xd] - sw.qx[i] + sw.qy[yd] - sw.qy[i])
            else
                sw.volume[i] +=
                    (
                        sum_at(swr.q, links_at_node.src[inds_riv[i]]) -
                        sum_at(swr.q, links_at_node.dst[inds_riv[i]]) + sw.qx[xd] -
                        sw.qx[i] + sw.qy[yd] - sw.qy[i] +
                        swr.inflow[inds_riv[i]] +
                        sw.runoff[i] - swr.abstraction[inds_riv[i]]
                    ) * dt
                if sw.volume[i] < T(0.0)
                    sw.error[i] = sw.error[i] + abs(sw.volume[i])
                    sw.volume[i] = T(0.0) # set volume to zero
                end
                if sw.volume[i] >= swr.bankfull_volume[inds_riv[i]]
                    swr.h[inds_riv[i]] =
                        swr.bankfull_depth[inds_riv[i]] +
                        (sw.volume[i] - swr.bankfull_volume[inds_riv[i]]) /
                        (sw.xl[i] * sw.yl[i])
                    sw.h[i] = swr.h[inds_riv[i]] - swr.bankfull_depth[inds_riv[i]]
                    swr.volume[inds_riv[i]] =
                        swr.h[inds_riv[i]] * swr.dl[inds_riv[i]] * swr.width[inds_riv[i]]
                else
                    swr.h[inds_riv[i]] =
                        sw.volume[i] / (swr.dl[inds_riv[i]] * swr.width[inds_riv[i]])
                    sw.h[i] = T(0.0)
                    swr.volume[inds_riv[i]] = sw.volume[i]
                end
                swr.h_av[inds_riv[i]] += swr.h[inds_riv[i]] * dt
            end
        else
            sw.volume[i] +=
                (sw.qx[xd] - sw.qx[i] + sw.qy[yd] - sw.qy[i] + sw.runoff[i]) * dt
            if sw.volume[i] < T(0.0)
                sw.error[i] = sw.error[i] + abs(sw.volume[i])
                sw.volume[i] = T(0.0) # set volume to zero
            end
            sw.h[i] = sw.volume[i] / (sw.xl[i] * sw.yl[i])
        end
        sw.h_av[i] += sw.h[i] * dt
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

@get_units @grid_loc @with_kw struct FloodPlain{T, P}
    profile::P                                         # floodplain profile
    mannings_n::Vector{T} | "s m-1/3"                   # manning's roughness
    mannings_n_sq::Vector{T} | "(s m-1/3)2" | "edge"    # manning's roughness squared
    volume::Vector{T} | "m3"                            # volume
    h::Vector{T} | "m"                                  # water depth
    h_av::Vector{T} | "m"                               # average water depth
    error::Vector{T} | "m3"                             # error volume
    a::Vector{T} | "m2" | "edge"                        # flow area
    r::Vector{T} | "m" | "edge"                         # hydraulic radius
    hf::Vector{T} | "m" | "edge"                        # water depth at edge/link
    zb_max::Vector{T} | "m" | "edge"                    # maximum bankfull elevation (edge/link)
    q0::Vector{T} | "m3 s-1" | "edge"                   # discharge at previous time step
    q::Vector{T} | "m3 s-1" | "edge"                    # discharge
    q_av::Vector{T} | "m" | "edge"                      # average river discharge
    hf_index::Vector{Int} | "-" | "edge"                # index with `hf` above depth threshold
end

"Determine the initial floodplain volume"
function initialize_volume!(river, nriv::Int)
    for i in 1:nriv
        i1, i2 =
            interpolation_indices(river.floodplain.h[i], river.floodplain.profile.depth)
        a = flow_area(
            river.floodplain.profile.width[i2, i],
            river.floodplain.profile.a[i1, i],
            river.floodplain.profile.depth[i1],
            river.floodplain.h[i],
        )
        a = max(a - (river.width[i] * river.floodplain.h[i]), 0.0)
        river.floodplain.volume[i] = river.dl[i] * a
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
function flood_depth(profile::FloodPlainProfile{T}, flood_volume, dl, i::Int)::T where {T}
    i1, i2 = interpolation_indices(flood_volume, @view profile.volume[:, i])
    ΔA = (flood_volume - profile.volume[i1, i]) / dl
    dh = ΔA / profile.width[i2, i]
    flood_depth = profile.depth[i1] + dh
    return flood_depth
end

"Initialize floodplain geometry and `FloodPlain` parameters"
function initialize_floodplain_1d(
    nc,
    config,
    inds,
    riverwidth,
    riverlength,
    zb,
    index_pit,
    n_edges,
    nodes_at_link,
)
    n_floodplain = ncread(
        nc,
        config,
        "lateral.river.floodplain.n";
        sel = inds,
        defaults = 0.072,
        type = Float,
    )
    volume = ncread(
        nc,
        config,
        "lateral.river.floodplain.volume";
        sel = inds,
        type = Float,
        dimname = :flood_depth,
    )
    n = length(inds)

    # for convenience (interpolation) flood depth 0.0 m is added, with associated area (a),
    # volume, width (river width) and wetted perimeter (p).
    volume = vcat(fill(Float(0), n)', volume)
    start_volume = volume
    flood_depths = Float.(nc["flood_depth"][:])
    pushfirst!(flood_depths, 0.0)
    n_depths = length(flood_depths)

    p = zeros(Float, n_depths, n)
    a = zeros(Float, n_depths, n)
    segment_volume = zeros(Float, n_depths, n)
    width = zeros(Float, n_depths, n)
    width[1, :] = riverwidth[1:n]

    # determine flow area (a), width and wetted perimeter (p)FloodPlain
    h = diff(flood_depths)
    incorrect_vol = 0
    riv_cells = 0
    error_vol = 0
    for i in 1:n
        riv_cell = 0
        diff_volume = diff(volume[:, i])

        for j in 1:(n_depths - 1)
            # assume rectangular shape of flood depth segment
            width[j + 1, i] = diff_volume[j] / (h[j] * riverlength[i])
            # check provided flood volume (floodplain width should be constant or increasing
            # as a function of flood depth)
            if width[j + 1, i] < width[j, i]
                # raise warning only if difference is larger than rounding error of 0.01 m³
                if ((width[j, i] - width[j + 1, i]) * h[j] * riverlength[i]) > 0.01
                    incorrect_vol += 1
                    riv_cell = 1
                    error_vol =
                        error_vol +
                        ((width[j, i] - width[j + 1, i]) * h[j] * riverlength[i])
                end
                width[j + 1, i] = width[j, i]
            end
            a[j + 1, i] = width[j + 1, i] * h[j]
            p[j + 1, i] = (width[j + 1, i] - width[j, i]) + 2.0 * h[j]
            segment_volume[j + 1, i] = a[j + 1, i] * riverlength[i]
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
    profile = FloodPlainProfile{Float, n_depths}(;
        volume = volume,
        width = width,
        depth = flood_depths,
        a = a,
        p = p,
    )

    # manning roughness at edges
    append!(n_floodplain, n_floodplain[index_pit]) # copy to ghost nodes
    mannings_n_sq = fill(Float(0), n_edges)
    zb_max = fill(Float(0), n_edges)
    for i in 1:n_edges
        src_node = nodes_at_link.src[i]
        dst_node = nodes_at_link.dst[i]
        mannings_n =
            (
                n_floodplain[dst_node] * riverlength[dst_node] +
                n_floodplain[src_node] * riverlength[src_node]
            ) / (riverlength[dst_node] + riverlength[src_node])
        mannings_n_sq[i] = mannings_n * mannings_n
        zb_max[i] = max(zb[src_node], zb[dst_node])
    end

    floodplain = FloodPlain(;
        profile = profile,
        mannings_n = n_floodplain,
        mannings_n_sq = mannings_n_sq,
        volume = zeros(n),
        error = zeros(n),
        h = zeros(n + length(index_pit)),
        h_av = zeros(n),
        a = zeros(n_edges),
        r = zeros(n_edges),
        hf = zeros(n_edges),
        zb_max = zb_max,
        q = zeros(n_edges),
        q_av = zeros(n_edges),
        q0 = zeros(n_edges),
        hf_index = zeros(Int, n_edges),
    )
    return floodplain
end

"""
    set_river_inwater!(model::Model, ssf_toriver)

Set `inwater` of the lateral river component for a `Model`. `ssf_toriver` is the subsurface
flow to the river.
"""
function set_river_inwater!(model::Model, ssf_toriver)
    (; lateral, vertical, network, config) = model
    (; net_runoff_river) = vertical.runoff.variables
    (; inwater) = lateral.river.boundary_conditions
    inds = network.index_river
    do_water_demand = haskey(config.model, "water_demand")
    if do_water_demand
        @. inwater = (
            ssf_toriver[inds] +
            lateral.land.to_river[inds] +
            # net_runoff_river
            (net_runoff_river[inds] * network.land.area[inds] * 0.001) / vertical.dt +
            (
                lateral.river.allocation.variables.nonirri_returnflow *
                0.001 *
                network.river.area
            ) / vertical.dt
        )
    else
        @. inwater = (
            ssf_toriver[inds] +
            lateral.land.variables.to_river[inds] +
            # net_runoff_river
            (net_runoff_river[inds] * network.land.area[inds] * 0.001) / vertical.dt
        )
    end
    return nothing
end

"""
    set_land_inwater!(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmGwfModel}

Set `inwater` of the lateral land component for the `SbmGwfModel` type.
"""
function set_land_inwater!(
    model::Model{N, L, V, R, W, T},
) where {N, L, V, R, W, T <: SbmGwfModel}
    (; lateral, vertical, network, config) = model
    (; net_runoff) = vertical.soil.variables
    (; inwater) = lateral.land.boundary_conditions
    do_drains = get(config.model, "drains", false)::Bool
    drainflux = zeros(length(net_runoff))
    do_water_demand = haskey(config.model, "water_demand")
    if do_drains
        drainflux[lateral.subsurface.drain.index] =
            -lateral.subsurface.drain.flux ./ tosecond(basetimestep)
    end
    if do_water_demand
        @. inwater =
            (net_runoff + vertical.allocation.variables.nonirri_returnflow) *
            network.land.area *
            0.001 / vertical.dt + drainflux
    else
        @. inwater = (net_runoff * network.land.area * 0.001) / vertical.dt + drainflux
    end
    return nothing
end

"""
    set_land_inwater!(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}

Set `inwater` of the lateral land component for the `SbmModel` type.
"""
function set_land_inwater!(
    model::Model{N, L, V, R, W, T},
) where {N, L, V, R, W, T <: SbmModel}
    (; lateral, vertical, network, config) = model
    (; net_runoff) = vertical.soil.variables
    (; inwater) = lateral.land.boundary_conditions
    do_water_demand = haskey(config.model, "water_demand")
    if do_water_demand
        @. inwater =
            (net_runoff + vertical.allocation.variables.nonirri_returnflow) *
            network.land.area *
            0.001 / vertical.dt
    else
        @. inwater = (net_runoff * network.land.area * 0.001) / vertical.dt
    end
    return nothing
end

# Computation of inflow from the lateral components `land` and `subsurface` to water bodies
# depends on the routing scheme (see different `get_inflow_waterbody` below). For the river
# kinematic wave, the variables `to_river` can be excluded, because this part is added to
# the river kinematic wave (kinematic wave is also solved for the water body cell). For
# local inertial river routing, `to_river` is included, because for the local inertial
# solution the water body cells are excluded (boundary condition). For `GroundwaterFlow`
# (Darcian flow in 4 directions), the lateral subsurface flow is excluded (for now) and
# inflow consists of overland flow.
"""
    set_inflow_waterbody!(
        model::Model{N,L,V,R,W,T},
    ) where {N,L<:NamedTuple{<:Any,<:Tuple{Any,SurfaceFlowLand,SurfaceFlowRiver}},V,R,W,T}

Set inflow from the subsurface and land components to a water body (reservoir or lake)
`inflow_wb` from a model type that contains the lateral components `SurfaceFlowLand` and
`SurfaceFlowRiver`.
"""
function set_inflow_waterbody!(
    model::Model{N, L, V, R, W, T},
) where {
    N,
    L <: NamedTuple{<:Any, <:Tuple{Any, SurfaceFlowLand, SurfaceFlowRiver}},
    V,
    R,
    W,
    T,
}
    (; lateral, network) = model
    (; subsurface, land, river) = lateral
    (; reservoir, lake, inflow_wb) = river.boundary_conditions
    inds = network.index_river

    if !isnothing(reservoir) || !isnothing(lake)
        if typeof(subsurface) <: LateralSSF || typeof(subsurface) <: GroundwaterExchange
            @. inflow_wb =
                subsurface.ssf[inds] / tosecond(basetimestep) + land.variables.q_av[inds]
        elseif typof(subsurface.flow) <: GroundwaterFlow || isnothing(subsurface)
            inflow_wb .= land.q_av[inds]
        end
    end
    return nothing
end

"""
    set_inflow_waterbody!(
        model::Model{N,L,V,R,W,T},
    ) where {N,L<:NamedTuple{<:Any,<:Tuple{Any,SurfaceFlowLand,ShallowWaterRiver}},V,R,W,T}

Set inflow from the subsurface and land components to a water body (reservoir or lake)
`inflow_wb` from a model type that contains the lateral components `SurfaceFlowLand` and
`ShallowWaterRiver`.
"""
function set_inflow_waterbody!(
    model::Model{N, L, V, R, W, T},
) where {
    N,
    L <: NamedTuple{<:Any, <:Tuple{Any, SurfaceFlowLand, ShallowWaterRiver}},
    V,
    R,
    W,
    T,
}
    (; lateral, network) = model
    (; subsurface, land, river) = lateral
    inds = network.index_river

    if !isnothing(lateral.river.reservoir) || !isnothing(lateral.river.lake)
        if typeof(subsurface) <: LateralSSF || typeof(subsurface) <: GroundwaterExchange
            @. river.inflow_wb =
                (subsurface.ssf[inds] + subsurface.to_river[inds]) /
                tosecond(basetimestep) +
                land.q_av[inds] +
                land.to_river[inds]
        elseif typeof(subsurface.flow) <: GroundwaterFlow || isnothing(subsurface)
            @. river.inflow_wb = lateral.land.q_av[inds] + lateral.land.to_river[inds]
        end
    end
    return nothing
end

"""
    set_inflow_waterbody!(
        model::Model{N,L,V,R,W,T},
    ) where {N,L<:NamedTuple{<:Any,<:Tuple{Any,ShallowWaterLand,ShallowWaterRiver}},V,R,W,T}

Set inflow from the subsurface and land components to a water body (reservoir or lake)
`inflow_wb` from a model type that contains the lateral components `ShallowWaterLand` and
`ShallowWaterRiver`.
"""
function set_inflow_waterbody!(
    model::Model{N, L, V, R, W, T},
) where {
    N,
    L <: NamedTuple{<:Any, <:Tuple{Any, ShallowWaterLand, ShallowWaterRiver}},
    V,
    R,
    W,
    T,
}
    (; lateral, network) = model
    (; subsurface, land, river) = lateral
    inds = network.index_river

    if !isnothing(lateral.river.reservoir) || !isnothing(lateral.river.lake)
        if typeof(subsurface) <: LateralSSF || typeof(subsurface) <: GroundwaterExchange
            @. land.inflow_wb[inds] =
                (subsurface.ssf[inds] + subsurface.to_river[inds]) / tosecond(basetimestep)
        end
    end
    return nothing
end

"""
    surface_routing!(model; ssf_toriver = 0.0)

Run surface routing (land and river). Kinematic wave for overland flow and kinematic wave or
local inertial model for river flow.
"""
function surface_routing!(model; ssf_toriver = 0.0)
    (; lateral, network, clock) = model

    # run kinematic wave for overland flow
    set_land_inwater!(model)
    update!(lateral.land, network.land, network.frac_toriver)

    # run river flow
    set_river_inwater!(model, ssf_toriver)
    set_inflow_waterbody!(model)
    update!(lateral.river, network.river, julian_day(clock.time - clock.dt))
    return nothing
end

"""
    surface_routing(
        model::Model{N,L,V,R,W,T};
        ssf_toriver = 0.0,
    ) where {N,L<:NamedTuple{<:Any,<:Tuple{Any,ShallowWaterLand,ShallowWaterRiver}},V,R,W,T}

Run surface routing (land and river) for a model type that contains the lateral components
`ShallowWaterLand` and `ShallowWaterRiver`.
"""
function surface_routing!(
    model::Model{N, L, V, R, W, T};
    ssf_toriver = 0.0,
) where {
    N,
    L <: NamedTuple{<:Any, <:Tuple{Any, ShallowWaterLand, ShallowWaterRiver}},
    V,
    R,
    W,
    T,
}
    (; lateral, vertical, network, clock) = model
    (; net_runoff) = vertical.soil.variables
    (; net_runoff_river) = vertical.runoff.variables

    @. lateral.land.runoff = (
        (net_runoff / 1000.0) * (network.land.area) / vertical.dt +
        ssf_toriver +
        # net_runoff_river
        ((net_runoff_river * network.land.area * 0.001) / vertical.dt)
    )
    set_inflow_waterbody!(model)
    update!(lateral.land, lateral.river, network, julian_day(clock.time - clock.dt))
    return nothing
end
