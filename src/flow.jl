
abstract type SurfaceFlow end

@get_units @exchange @grid_type @grid_location @with_kw struct SurfaceFlowRiver{T,R,L} <: SurfaceFlow
    β::T | "-" | _ | "scalar"                    # constant in Manning's equation
    sl::Vector{T} | "m m-1"                      # Slope [m m⁻¹]
    n::Vector{T} | "s m-1/3"                     # Manning's roughness [s m⁻⅓]
    dl::Vector{T} | "m"                          # Drain length [m]
    q::Vector{T} | "m3 s-1"                      # Discharge [m³ s⁻¹]
    qin::Vector{T} | "m3 s-1"                    # Inflow from upstream cells [m³ s⁻¹]
    q_av::Vector{T} | "m3 s-1"                   # Average discharge [m³ s⁻¹]
    qlat::Vector{T} | "m2 s-1"                   # Lateral inflow per unit length [m² s⁻¹]
    inwater::Vector{T} | "m3 s-1"                # Lateral inflow [m³ s⁻¹]
    inflow::Vector{T} | "m3 s-1"                 # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    volume::Vector{T} | "m3"                     # Kinematic wave volume [m³] (based on water level h)
    h::Vector{T} | "m"                           # Water level [m]
    h_av::Vector{T} | "m"                        # Average water level [m]
    bankfull_depth::Vector{T} | "m"              # Bankfull water level [m]
    Δt::T | "s" | 0 | "none" | "none"            # Model time step [s]
    its::Int | "-" | 0 | "none" | "none"         # Number of fixed iterations
    width::Vector{T} | "m"                       # Flow width [m]
    alpha_pow::T | "-" | _ | "scalar"            # Used in the power part of α
    alpha_term::Vector{T} | "-"                  # Term used in computation of α
    α::Vector{T} | "s3/5 m1/5"                   # Constant in momentum equation A = αQᵝ, based on Manning's equation
    cel::Vector{T} | "m s-1"                     # Celerity of the kinematic wave
    reservoir_index::Vector{Int} | "-"           # map cell to 0 (no reservoir) or i (pick reservoir i in reservoir field)
    lake_index::Vector{Int} | "-"                # map cell to 0 (no lake) or i (pick lake i in lake field)
    reservoir::R | "-" | 0                       # Reservoir model struct of arrays
    lake::L | "-" | 0                            # Lake model struct of arrays
    kinwave_it::Bool | "-" | 0 | "none" | "none" # Boolean for iterations kinematic wave

    # TODO unclear why this causes a MethodError
    # function SurfaceFlow{T,R,L}(args...) where {T,R,L}
    #     equal_size_vectors(args)
    #     return new(args...)
    # end
end

@get_units @exchange @grid_type @grid_location @with_kw struct SurfaceFlowLand{T} <: SurfaceFlow
    β::T | "-"  | _ | "scalar"                      # constant in Manning's equation
    sl::Vector{T} | "m m-1"                         # Slope [m m⁻¹]
    n::Vector{T} | "s m-1/3"                        # Manning's roughness [s m⁻⅓]
    dl::Vector{T} | "m"                             # Drain length [m]
    q::Vector{T} | "m3 s-1"                         # Discharge [m³ s⁻¹]
    qin::Vector{T} | "m3 s-1"                       # Inflow from upstream cells [m³ s⁻¹]
    q_av::Vector{T} | "m3 s-1"                      # Average discharge [m³ s⁻¹]
    qlat::Vector{T} | "m2 s-1"                      # Lateral inflow per unit length [m² s⁻¹]
    inwater::Vector{T} | "m3 s-1"                   # Lateral inflow [m³ s⁻¹]
    volume::Vector{T} | "m3"                        # Kinematic wave volume [m³] (based on water level h)
    h::Vector{T} | "m"                              # Water level [m]
    h_av::Vector{T} | "m"                           # Average water level [m]
    Δt::T | "s" | 0 | "none" | "none"               # Model time step [s]
    its::Int | "-" | 0 | "none" | "none"            # Number of fixed iterations
    width::Vector{T} | "m"                          # Flow width [m]
    alpha_pow::T | "-" | _ | "scalar"               # Used in the power part of α
    alpha_term::Vector{T} | "-"                     # Term used in computation of α
    α::Vector{T} | "s3/5 m1/5"                      # Constant in momentum equation A = αQᵝ, based on Manning's equation
    cel::Vector{T} | "m s-1"                        # Celerity of the kinematic wave
    to_river::Vector{T} | "m3 s-1"                  # Part of overland flow [m³ s⁻¹] that flows to the river
    kinwave_it::Bool | "-" | 0 |  "none" | "none"   # Boolean for iterations kinematic wave
end

function initialize_surfaceflow_land(nc, config, inds; sl, dl, width, iterate, tstep, Δt)
    @info "Kinematic wave approach is used for overland flow." iterate
    if tstep > 0
        @info "Using a fixed sub-timestep (seconds) $tstep for kinematic wave overland flow."
    end

    n_land =
        ncread(nc, config, "lateral.land.n"; sel = inds, defaults = 0.072, type = Float)
    n = length(inds)

    sf_land = SurfaceFlowLand(
        β = Float(0.6),
        sl = sl,
        n = n_land,
        dl = dl,
        q = zeros(Float, n),
        qin = zeros(Float, n),
        q_av = zeros(Float, n),
        qlat = zeros(Float, n),
        inwater = zeros(Float, n),
        volume = zeros(Float, n),
        h = zeros(Float, n),
        h_av = zeros(Float, n),
        Δt = Float(tosecond(Δt)),
        its = tstep > 0 ? Int(cld(tosecond(Δt), tstep)) : tstep,
        width = width,
        alpha_pow = Float((2.0 / 3.0) * 0.6),
        alpha_term = fill(mv, n),
        α = fill(mv, n),
        cel = zeros(Float, n),
        to_river = zeros(Float, n),
        kinwave_it = iterate,
    )

    return sf_land
end

function initialize_surfaceflow_river(
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
    Δt,
)
    @info "Kinematic wave approach is used for river flow." iterate
    if tstep > 0
        @info "Using a fixed sub-timestep (seconds) $tstep for kinematic wave river flow."
    end

    n_river =
        ncread(nc, config, "lateral.river.n"; sel = inds, defaults = 0.036, type = Float)
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
    sl = ncread(
        nc,
        config,
        "lateral.river.slope";
        optional = false,
        sel = inds,
        type = Float,
    )
    clamp!(sl, 0.00001, Inf)

    n = length(inds)

    sf_river = SurfaceFlowRiver(
        β = Float(0.6),
        sl = sl,
        n = n_river,
        dl = dl,
        q = zeros(Float, n),
        qin = zeros(Float, n),
        q_av = zeros(Float, n),
        qlat = zeros(Float, n),
        inwater = zeros(Float, n),
        inflow = zeros(Float, n),
        volume = zeros(Float, n),
        h = zeros(Float, n),
        h_av = zeros(Float, n),
        bankfull_depth = bankfull_depth,
        Δt = Float(tosecond(Δt)),
        its = tstep > 0 ? Int(cld(tosecond(Δt), tstep)) : tstep,
        width = width,
        alpha_pow = Float((2.0 / 3.0) * 0.6),
        alpha_term = fill(mv, n),
        α = fill(mv, n),
        cel = zeros(Float, n),
        reservoir_index = reservoir_index,
        lake_index = lake_index,
        reservoir = reservoir,
        lake = lake,
        kinwave_it = iterate,
    )

    return sf_river
end

statevars(::SurfaceFlowRiver) = (:q, :h, :h_av)
statevars(::SurfaceFlowLand) = (:q, :h, :h_av)

function update(sf::SurfaceFlowLand, network, frac_toriver)
    @unpack graph, subdomain_order, topo_subdomain, indices_subdomain, upstream_nodes =
        network

    ns = length(subdomain_order)

    @. sf.alpha_term = pow(sf.n / sqrt(sf.sl), sf.β)
    # use fixed alpha value based flow width
    @. sf.α = sf.alpha_term * pow(sf.width, sf.alpha_pow)
    @. sf.qlat = sf.inwater / sf.dl

    sf.q_av .= 0.0
    sf.h_av .= 0.0
    sf.to_river .= 0.0

    Δt, its = stable_timestep(sf)
    for _ = 1:its
        sf.qin .= 0.0
        for k = 1:ns
            threaded_foreach(eachindex(subdomain_order[k]), basesize = 1) do i
                m = subdomain_order[k][i]
                for (n, v) in zip(indices_subdomain[m], topo_subdomain[m])
                    # for a river cell without a reservoir or lake part of the upstream
                    # surface flow goes to the river (frac_toriver) and part goes to the
                    # surface flow reservoir (1.0 - frac_toriver), upstream nodes with a
                    # reservoir or lake are excluded
                    sf.to_river[v] += sum_at(
                        i -> sf.q[i] * frac_toriver[i],
                        upstream_nodes[n],
                        eltype(sf.to_river),
                    )
                    if sf.width[v] > 0.0
                        sf.qin[v] = sum_at(
                            i -> sf.q[i] * (1.0 - frac_toriver[i]),
                            upstream_nodes[n],
                            eltype(sf.q),
                        )
                    end

                    sf.q[v] = kinematic_wave(
                        sf.qin[v],
                        sf.q[v],
                        sf.qlat[v],
                        sf.α[v],
                        sf.β,
                        Δt,
                        sf.dl[v],
                    )

                    # update h, only if surface width > 0.0
                    if sf.width[v] > 0.0
                        crossarea = sf.α[v] * pow(sf.q[v], sf.β)
                        sf.h[v] = crossarea / sf.width[v]
                    end
                    sf.q_av[v] += sf.q[v]
                    sf.h_av[v] += sf.h[v]
                end
            end
        end
    end
    sf.q_av ./= its
    sf.h_av ./= its
    sf.to_river ./= its
    sf.volume .= sf.dl .* sf.width .* sf.h
end

function update(sf::SurfaceFlowRiver, network, inflow_wb, doy)

    @unpack graph, subdomain_order, topo_subdomain, indices_subdomain, upstream_nodes =
        network

    ns = length(subdomain_order)

    @. sf.alpha_term = pow(sf.n / sqrt(sf.sl), sf.β)
    # use fixed alpha value based on 0.5 * bankfull_depth
    @. sf.α = sf.alpha_term * pow(sf.width + sf.bankfull_depth, sf.alpha_pow)

    @. sf.qlat = sf.inwater / sf.dl

    sf.q_av .= 0.0
    sf.h_av .= 0.0
    # because of possible iterations set reservoir and lake inflow and total outflow at
    # start to zero, the total sum of inflow and outflow at each sub time step is calculated
    if !isnothing(sf.reservoir)
        sf.reservoir.inflow .= 0.0
        sf.reservoir.totaloutflow .= 0.0
    end
    if !isnothing(sf.lake)
        sf.lake.inflow .= 0.0
        sf.lake.totaloutflow .= 0.0
    end

    Δt, its = stable_timestep(sf)
    for _ = 1:its
        sf.qin .= 0.0
        for k = 1:ns
            threaded_foreach(eachindex(subdomain_order[k]), basesize = 1) do i
                m = subdomain_order[k][i]
                for (n, v) in zip(indices_subdomain[m], topo_subdomain[m])
                    # sf.qin by outflow from upstream reservoir or lake location is added
                    sf.qin[v] += sum_at(sf.q, upstream_nodes[n])

                    # Inflow supply/abstraction is added to qlat (divide by flow length)
                    # If inflow < 0, abstraction is limited
                    if sf.inflow[v] < 0.0
                        max_abstract = min(sf.qin[v] + sf.qlat[v] * sf.dl[v], -sf.inflow[v])
                        inflow = -max_abstract / sf.dl[v]
                    else
                        inflow = sf.inflow[v] / sf.dl[v]
                    end

                    sf.q[v] = kinematic_wave(
                        sf.qin[v],
                        sf.q[v],
                        sf.qlat[v] + inflow,
                        sf.α[v],
                        sf.β,
                        Δt,
                        sf.dl[v],
                    )

                    if !isnothing(sf.reservoir) && sf.reservoir_index[v] != 0
                        # run reservoir model and copy reservoir outflow to inflow (qin) of
                        # downstream river cell
                        i = sf.reservoir_index[v]
                        update(sf.reservoir, i, sf.q[v] + inflow_wb[v], Δt)

                        downstream_nodes = outneighbors(graph, v)
                        n_downstream = length(downstream_nodes)
                        if n_downstream == 1
                            j = only(downstream_nodes)
                            sf.qin[j] = sf.reservoir.outflow[i]
                        elseif n_downstream == 0
                            error(
                                """A reservoir without a downstream river node is not supported.
                                Add a downstream river node or move the reservoir to an upstream node (model schematization).
                                """,
                            )
                        else
                            error("bifurcations not supported")
                        end

                    elseif !isnothing(sf.lake) && sf.lake_index[v] != 0
                        # run lake model and copy lake outflow to inflow (qin) of downstream river
                        # cell
                        i = sf.lake_index[v]
                        update(sf.lake, i, sf.q[v] + inflow_wb[v], doy, Δt)

                        downstream_nodes = outneighbors(graph, v)
                        n_downstream = length(downstream_nodes)
                        if n_downstream == 1
                            j = only(downstream_nodes)
                            sf.qin[j] = sf.lake.outflow[i]
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
                    crossarea = sf.α[v] * pow(sf.q[v], sf.β)
                    sf.h[v] = crossarea / sf.width[v]

                    sf.q_av[v] += sf.q[v]
                    sf.h_av[v] += sf.h[v]
                end
            end
        end
    end
    sf.q_av ./= its
    sf.h_av ./= its
    sf.volume .= sf.dl .* sf.width .* sf.h
end

function stable_timestep(sf::S) where {S<:SurfaceFlow}
    n = length(sf.q)
    # two options for iteration, fixed or based on courant number.
    if sf.kinwave_it
        if sf.its > 0
            its = sf.its
        else
            # calculate celerity
            courant = zeros(n)
            for v = 1:n
                if sf.q[v] > 0.0
                    sf.cel[v] = 1.0 / (sf.α[v] * sf.β * pow(sf.q[v], (sf.β - 1.0)))
                    courant[v] = (sf.Δt / sf.dl[v]) * sf.cel[v]
                end
            end
            filter!(x -> x ≠ 0.0, courant)
            its = isempty(courant) ? 1 : ceil(Int, (1.25 * quantile!(courant, 0.95)))
        end
    else
        its = 1
    end

    # sub time step
    Δt = sf.Δt / its
    return Δt, its
end

@get_units @exchange @grid_type  @grid_location @with_kw struct LateralSSF{T}
    kh₀::Vector{T} | "m d-1"               # Horizontal hydraulic conductivity at soil surface [m d⁻¹]
    f::Vector{T} | "m-1"                   # A scaling parameter [m⁻¹] (controls exponential decline of kh₀)
    soilthickness::Vector{T} | "m"         # Soil thickness [m]
    θₛ::Vector{T} | "-"                     # Saturated water content (porosity) [-]
    θᵣ::Vector{T} | "-"                    # Residual water content [-]
    Δt::T | "d" | 0 | "none" | "none"      # model time step [d]
    βₗ::Vector{T} | "m m-1"                 # Slope [m m⁻¹]
    dl::Vector{T} | "m"                    # Drain length [m]
    dw::Vector{T} | "m"                    # Flow width [m]
    zi::Vector{T} | "m"                    # Pseudo-water table depth [m] (top of the saturated zone)
    exfiltwater::Vector{T} | "m Δt-1"      # Exfiltration [m Δt⁻¹] (groundwater above surface level, saturated excess conditions)
    recharge::Vector{T} | "m2 Δt-1"        # Net recharge to saturated store [m² Δt⁻¹]
    ssf::Vector{T} | "m3 d-1"              # Subsurface flow [m³ d⁻¹]
    ssfin::Vector{T} | "m3 d-1"            # Inflow from upstream cells [m³ d⁻¹]
    ssfmax::Vector{T} | "m2 d-1"           # Maximum subsurface flow [m² d⁻¹]
    to_river::Vector{T} | "m3 d-1"         # Part of subsurface flow [m³ d⁻¹] that flows to the river

    function LateralSSF{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::LateralSSF) = (:ssf,)

function update(ssf::LateralSSF, network, frac_toriver)
    @unpack subdomain_order, topo_subdomain, indices_subdomain, upstream_nodes = network

    ns = length(subdomain_order)
    for k = 1:ns
        threaded_foreach(eachindex(subdomain_order[k]), basesize = 1) do i
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
                    ssf.kh₀[v],
                    ssf.βₗ[v],
                    ssf.θₛ[v] - ssf.θᵣ[v],
                    ssf.f[v],
                    ssf.soilthickness[v],
                    ssf.Δt,
                    ssf.dl[v],
                    ssf.dw[v],
                    ssf.ssfmax[v],
                )
            end
        end
    end
end

@get_units @exchange @grid_type @grid_location @with_kw struct GroundwaterExchange{T}
    Δt::T | "d" | 0 | "none" | "none"   # model time step [d]
    exfiltwater::Vector{T} | "m Δt-1"   # Exfiltration [m Δt⁻¹]  (groundwater above surface level, saturated excess conditions)
    zi::Vector{T} | "m"                 # Pseudo-water table depth [m] (top of the saturated zone)
    to_river::Vector{T} | "m3 d-1"      # Part of subsurface flow [m³ d⁻¹] that flows to the river
    ssf::Vector{T} | "m3 d-1"           # Subsurface flow [m³ d⁻¹]
end

@get_units @exchange @grid_type @grid_location @with_kw struct ShallowWaterRiver{T,R,L,F}
    n::Int | "-" | 0 | "none" | "none"                      # number of cells
    ne::Int | "-" | 0 | "none" | "none"                     # number of edges/links
    active_n::Vector{Int} | "-"                             # active nodes
    active_e::Vector{Int} | "-" | _ | "edge"                # active edges/links
    g::T | "m s-2" | _ | "scalar"                           # acceleration due to gravity
    α::T | "-" | _ | "scalar"                               # stability coefficient (Bates et al., 2010)
    h_thresh::T | "m" | _ | "scalar"                        # depth threshold for calculating flow
    Δt::T | "s" | 0 | "none" | "none"                       # model time step [s]
    q::Vector{T} | "m3 s-1" | _ | "edge"                    # river discharge (subgrid channel)
    q0::Vector{T} | "m3 s-1"| _ | "edge"                    # river discharge (subgrid channel) at previous time step
    q_av::Vector{T} | "m3 s-1" | _ | "edge"                 # average river channel (+ floodplain) discharge [m³ s⁻¹]
    q_channel_av::Vector{T} | "m3 s-1"                      # average river channel discharge [m³ s⁻¹]
    zb_max::Vector{T} | "m"                                 # maximum channel bed elevation
    mannings_n_sq::Vector{T} | "(s m-1/3)2" | _ | "edge"    # Manning's roughness squared at edge/link
    mannings_n::Vector{T} | "s m-1/3"                       # Manning's roughness at node
    h::Vector{T} | "m"                                      # water depth
    η_max::Vector{T} | "m" | _ | "edge"                     # maximum water elevation at edge
    η_src::Vector{T} | "m"                                  # water elevation of source node of edge
    η_dst::Vector{T} | "m"                                  # water elevation of downstream node of edge
    hf::Vector{T} | "m" | _ | "edge"                        # water depth at edge/link
    h_av::Vector{T} | "m"                                   # average water depth
    dl::Vector{T} | "m"                                     # river length
    dl_at_link::Vector{T} | "m" | _ | "edge"                # river length at edge/link
    width::Vector{T} | "m"                      	        # river width
    width_at_link::Vector{T} | "m" | _ | "edge"             # river width at edge/link
    a::Vector{T} | "m2" | _ | "edge"                        # flow area at edge/link
    r::Vector{T} | "m" | _ | "edge"                         # wetted perimeter at edge/link
    volume::Vector{T} | "m3"                                # river volume
    error::Vector{T} | "m3"                                 # error volume
    inwater::Vector{T} | "m3 s-1"                           # lateral inflow [m³ s⁻¹]
    inflow::Vector{T} | "m3 s-1"                            # external inflow (abstraction/supply/demand) [m³ s⁻¹]
    inflow_wb::Vector{T} | "m3 s-1"                         # inflow waterbody (lake or reservoir model) from land part [m³ s⁻¹]
    bankfull_volume::Vector{T} | "m3"                       # bankfull volume
    bankfull_depth::Vector{T} | "m"                         # bankfull depth
    zb::Vector{T} | "m"                                     # river bed elevation
    froude_limit::Bool | "-" | 0 | "none" | "none"          # if true a check is performed if froude number > 1.0 (algorithm is modified)
    reservoir_index::Vector{Int} | "-"                      # river cell index with a reservoir (each index of reservoir_index maps to reservoir i in reservoir field)
    lake_index::Vector{Int} | "-"                           # river cell index with a lake (each index of lake_index maps to lake i in lake field)
    waterbody::Vector{Bool} | "-"                           # water body cells (reservoir or lake)
    reservoir::R | "-" | 0                                  # Reservoir model struct of arrays
    lake::L | "-" | 0                                       # Lake model struct of arrays
    floodplain::F | "-" | 0                                 # Floodplain (1D) schematization
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
    Δt,
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
    riverlength_bc = get(config.model, "riverlength_bc", 1.0e04)::Float64 # river length at boundary point (ghost point)
    alpha = get(config.model, "inertial_flow_alpha", 0.7)::Float64 # stability coefficient for model time step (0.2-0.7)
    h_thresh = get(config.model, "h_thresh", 1.0e-03)::Float64 # depth threshold for flow at link
    froude_limit = get(config.model, "froude_limit", true)::Bool # limit flow to subcritical according to Froude number
    floodplain_1d = floodplain

    @info "Local inertial approach is used for river flow." alpha h_thresh froude_limit riverlength_bc floodplain_1d

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
    # set ghost points for boundary condition (downstream river outlet): river width, bed
    # elevation, manning n is copied from the upstream cell, river depth h is set at 0.0
    # (fixed). river length at boundary point is by default 1.0e4 m (riverlength_bc).
    index_pit = findall(x -> x == 5, ldd)
    npits = length(index_pit)
    add_vertex_edge_graph!(graph, index_pit)
    append!(dl, fill(riverlength_bc, npits))
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
    for i = 1:_ne
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

    # set depth h to zero (including reservoir and lake locations)
    h = fill(0.0, n + length(index_pit))
    q_av = zeros(_ne)

    waterbody = !=(0).(reservoir_index .+ lake_index)
    active_index = findall(x -> x == 0, waterbody)

    sw_river = ShallowWaterRiver(
        n = n,
        ne = _ne,
        active_n = active_index,
        active_e = active_index,
        g = 9.80665,
        α = alpha,
        h_thresh = h_thresh,
        Δt = tosecond(Δt),
        q = zeros(_ne),
        q0 = zeros(_ne),
        q_av = q_av,
        q_channel_av = isnothing(floodplain) ? q_av : zeros(_ne),
        zb_max = zb_max,
        mannings_n_sq = mannings_n_sq,
        mannings_n = n_river,
        h = h,
        η_max = zeros(_ne),
        η_src = zeros(_ne),
        η_dst = zeros(_ne),
        hf = zeros(_ne),
        h_av = zeros(n),
        width = width,
        width_at_link = width_at_link,
        a = zeros(_ne),
        r = zeros(_ne),
        volume = fill(0.0, n),
        error = zeros(n),
        inflow = zeros(n),
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
    )
    return sw_river, nodes_at_link
end

function shallowwater_river_update(
    sw::ShallowWaterRiver,
    network,
    Δt,
    inflow_wb,
    doy,
    update_h,
)

    @unpack nodes_at_link, links_at_node = network

    sw.q0 .= sw.q
    if !isnothing(sw.floodplain)
        sw.floodplain.q0 .= sw.floodplain.q
    end
    # For reservoir and lake locations the local inertial solution is replaced by the
    # reservoir or lake model. These locations are handled as boundary conditions in the
    # local inertial model (fixed h).
    for v in eachindex(sw.reservoir_index)
        i = sw.reservoir_index[v]
        update(
            sw.reservoir,
            v,
            sum_at(sw.q0, links_at_node.src[i]) + inflow_wb[i] + sw.inflow_wb[i],
            Δt,
        )
        sw.q[i] = sw.reservoir.outflow[v]
        sw.q_av[i] += sw.q[i] * Δt
    end
    for v in eachindex(sw.lake_index)
        i = sw.lake_index[v]
        update(
            sw.lake,
            v,
            sum_at(sw.q0, links_at_node.src[i]) + inflow_wb[i] + sw.inflow_wb[i],
            doy,
            Δt,
        )
        sw.q[i] = sw.lake.outflow[v]
        sw.q_av[i] += sw.q[i] * Δt
    end
    @tturbo for j in eachindex(sw.active_e)
        i = sw.active_e[j]
        i_src = nodes_at_link.src[i]
        i_dst = nodes_at_link.dst[i]
        sw.η_src[i] = sw.zb[i_src] + sw.h[i_src]
        sw.η_dst[i] = sw.zb[i_dst] + sw.h[i_dst]

        sw.η_max[i] = max(sw.η_src[i], sw.η_dst[i])
        sw.hf[i] = (sw.η_max[i] - sw.zb_max[i])

        sw.a[i] = sw.width_at_link[i] * sw.hf[i] # flow area (rectangular channel)
        sw.r[i] = sw.a[i] / (sw.width_at_link[i] + 2.0 * sw.hf[i]) # hydraulic radius (rectangular channel)

        sw.q[i] = IfElse.ifelse(
            sw.hf[i] > sw.h_thresh,
            local_inertial_flow(
                sw.q0[i],
                sw.η_src[i],
                sw.η_dst[i],
                sw.hf[i],
                sw.a[i],
                sw.r[i],
                sw.dl_at_link[i],
                sw.mannings_n_sq[i],
                sw.g,
                sw.froude_limit,
                Δt,
            ),
            0.0,
        )

        # limit q in case water is not available
        sw.q[i] = IfElse.ifelse(sw.h[i_src] <= 0.0, min(sw.q[i], 0.0), sw.q[i])
        sw.q[i] = IfElse.ifelse(sw.h[i_dst] <= 0.0, max(sw.q[i], 0.0), sw.q[i])

        sw.q_av[i] += sw.q[i] * Δt
    end
    if !isnothing(sw.floodplain)
        @tturbo @. sw.floodplain.hf = max(sw.η_max - sw.floodplain.zb_max, 0.0)

        n = 0
        @inbounds for i in sw.active_e
            @inbounds if sw.floodplain.hf[i] > sw.h_thresh
                n += 1
                sw.floodplain.hf_index[n] = i
            else
                sw.floodplain.q[i] = 0.0
            end
        end

        @tturbo for j = 1:n
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
                    sw.η_src[i],
                    sw.η_dst[i],
                    sw.floodplain.hf[i],
                    sw.floodplain.a[i],
                    sw.floodplain.r[i],
                    sw.dl_at_link[i],
                    sw.floodplain.mannings_n_sq[i],
                    sw.g,
                    sw.froude_limit,
                    Δt,
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
            sw.floodplain.q_av[i] += sw.floodplain.q[i] * Δt
        end
    end
    if update_h
        @batch per = thread minbatch = 2000 for i in sw.active_n

            q_src = sum_at(sw.q, links_at_node.src[i])
            q_dst = sum_at(sw.q, links_at_node.dst[i])
            sw.volume[i] = sw.volume[i] + (q_src - q_dst + sw.inwater[i]) * Δt

            if sw.volume[i] < 0.0
                sw.error[i] = sw.error[i] + abs(sw.volume[i])
                sw.volume[i] = 0.0 # set volume to zero
            end
            sw.volume[i] = max(sw.volume[i] + sw.inflow[i] * Δt, 0.0) # add external inflow

            if !isnothing(sw.floodplain)
                q_src = sum_at(sw.floodplain.q, links_at_node.src[i])
                q_dst = sum_at(sw.floodplain.q, links_at_node.dst[i])
                sw.floodplain.volume[i] = sw.floodplain.volume[i] + (q_src - q_dst) * Δt
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
                sw.floodplain.h_av[i] += sw.floodplain.h[i] * Δt
            else
                sw.h[i] = sw.volume[i] / (sw.dl[i] * sw.width[i])
            end
            sw.h_av[i] += sw.h[i] * Δt
        end
    end
end

function update(
    sw::ShallowWaterRiver{T},
    network,
    inflow_wb,
    doy;
    update_h = true,
) where {T}
    @unpack nodes_at_link, links_at_node = network

    if !isnothing(sw.reservoir)
        sw.reservoir.inflow .= 0.0
        sw.reservoir.totaloutflow .= 0.0
    end
    if !isnothing(sw.lake)
        sw.lake.inflow .= 0.0
        sw.lake.totaloutflow .= 0.0
    end
    if !isnothing(sw.floodplain)
        sw.floodplain.q_av .= 0.0
        sw.floodplain.h_av .= 0.0
    end
    sw.q_av .= 0.0
    sw.h_av .= 0.0

    t = T(0.0)
    while t < sw.Δt
        Δt = stable_timestep(sw)
        if t + Δt > sw.Δt
            Δt = sw.Δt - t
        end
        shallowwater_river_update(sw, network, Δt, inflow_wb, doy, update_h)
        t = t + Δt
    end
    sw.q_av ./= sw.Δt
    sw.h_av ./= sw.Δt

    if !isnothing(sw.floodplain)
        sw.floodplain.q_av ./= sw.Δt
        sw.floodplain.h_av ./= sw.Δt
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

@get_units @exchange @grid_type @grid_location @with_kw struct ShallowWaterLand{T}
    n::Int | "-" | 0 | "none" | "none"                      # number of cells
    xl::Vector{T} | "m"                                     # cell length x direction
    yl::Vector{T} | "m"                                     # cell length y direction
    xwidth::Vector{T} | "m" | _ | "edge"                    # effective flow width x direction (floodplain)
    ywidth::Vector{T} | "m" | _ | "edge"                    # effective flow width y direction (floodplain)
    g::T | "m2 s-1" | _ | "scalar"                          # acceleration due to gravity
    θ::T | "-" | _ | "scalar"                               # weighting factor (de Almeida et al., 2012)
    α::T | "-" | _ | "scalar"                               # stability coefficient (de Almeida et al., 2012)
    h_thresh::T | "m" | _ | "scalar"                        # depth threshold for calculating flow
    Δt::T | "s" | _ | "scalar"                              # model time step [s]
    qy0::Vector{T} | "m3 s-1" | _ | "edge"                  # flow in y direction at previous time step
    qx0::Vector{T} | "m3 s-1" | _ | "edge"                  # flow in x direction at previous time step
    qx::Vector{T} | "m3 s-1" | _ | "edge"                   # flow in x direction
    qy::Vector{T} | "m3 s-1" | _ | "edge"                   # flow in y direction
    zx_max::Vector{T} | "m" | _ | "edge"                    # maximum cell elevation (x direction)
    zy_max::Vector{T} | "m" | _ | "edge"                    # maximum cell elevation (y direction)
    mannings_n_sq::Vector{T} | "(s m-1/3)2" | _ | "edge"    # Manning's roughness squared
    volume::Vector{T} | "m3"                                # total volume of cell (including river volume for river cells)
    error::Vector{T} | "m3"                                 # error volume
    runoff::Vector{T} | "m3 s-1"                            # runoff from hydrological model
    h::Vector{T} | "m"                                      # water depth of cell (for river cells the reference is the river bed elevation `zb`)
    z::Vector{T} | "m"                                      # elevation of cell
    froude_limit::Bool | "-" | 0  |"none" | "none"          # if true a check is performed if froude number > 1.0 (algorithm is modified)
    rivercells::Vector{Bool} | "-"                          # river cells
    h_av::Vector{T} | "m"                                   # average water depth (for river cells the reference is the river bed elevation `zb`)
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
    Δt,
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
    indices = Indices(xu = zeros(n), xd = zeros(n), yu = zeros(n), yd = zeros(n))

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
    for i = 1:n
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

    sw_land = ShallowWaterLand{Float}(
        n = n,
        xl = xlength,
        yl = ylength,
        xwidth = we_x,
        ywidth = we_y,
        g = 9.80665,
        θ = theta,
        α = alpha,
        h_thresh = h_thresh,
        Δt = tosecond(Δt),
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

Δt = α * (Δx / sqrt(g max(h))
"""
function stable_timestep(sw::ShallowWaterRiver{T})::T where {T}
    Δtₘᵢₙ = T(Inf)
    @tturbo for i = 1:sw.n
        Δt = sw.α * sw.dl[i] / sqrt(sw.g * sw.h[i])
        Δtₘᵢₙ = Δt < Δtₘᵢₙ ? Δt : Δtₘᵢₙ
    end
    Δtₘᵢₙ = isinf(Δtₘᵢₙ) ? T(10.0) : Δtₘᵢₙ
    return Δtₘᵢₙ
end

function stable_timestep(sw::ShallowWaterLand{T})::T where {T}
    Δtₘᵢₙ = T(Inf)
    @tturbo for i = 1:sw.n
        Δt = IfElse.ifelse(
            sw.rivercells[i] == 0,
            sw.α * min(sw.xl[i], sw.yl[i]) / sqrt(sw.g * sw.h[i]),
            T(Inf),
        )
        Δtₘᵢₙ = Δt < Δtₘᵢₙ ? Δt : Δtₘᵢₙ
    end
    Δtₘᵢₙ = isinf(Δtₘᵢₙ) ? T(10.0) : Δtₘᵢₙ
    return Δtₘᵢₙ
end

function update(
    sw::ShallowWaterLand{T},
    swr::ShallowWaterRiver{T},
    network,
    inflow_wb,
    doy;
    update_h = false,
) where {T}

    @unpack nodes_at_link, links_at_node = network.river

    if !isnothing(swr.reservoir)
        swr.reservoir.inflow .= 0.0
        swr.reservoir.totaloutflow .= 0.0
    end
    if !isnothing(swr.lake)
        swr.lake.inflow .= 0.0
        swr.lake.totaloutflow .= 0.0
    end
    swr.q_av .= 0.0
    swr.h_av .= 0.0
    sw.h_av .= 0.0

    t = T(0.0)
    while t < swr.Δt
        Δt_river = stable_timestep(swr)
        Δt_land = stable_timestep(sw)
        Δt = min(Δt_river, Δt_land)
        if t + Δt > swr.Δt
            Δt = swr.Δt - t
        end
        shallowwater_river_update(swr, network.river, Δt, inflow_wb, doy, update_h)
        update(sw, swr, network, Δt)
        t = t + Δt
    end
    swr.q_av ./= swr.Δt
    swr.h_av ./= swr.Δt
    sw.h_av ./= sw.Δt
end

function update(sw::ShallowWaterLand{T}, swr::ShallowWaterRiver{T}, network, Δt) where {T}

    indices = network.land.staggered_indices
    inds_riv = network.land.index_river

    @unpack nodes_at_link, links_at_node = network.river

    sw.qx0 .= sw.qx
    sw.qy0 .= sw.qy

    # update qx
    @batch per = thread minbatch = 6000 for i = 1:sw.n
        yu = indices.yu[i]
        yd = indices.yd[i]
        xu = indices.xu[i]
        xd = indices.xd[i]

        # the effective flow width is zero when the river width exceeds the cell width (dy
        # for flow in x dir) and floodplain flow is not calculated.
        if xu <= sw.n && sw.ywidth[i] != T(0.0)

            η_x = sw.z[i] + sw.h[i]
            η_xu = sw.z[xu] + sw.h[xu]
            η_max = max(η_x, η_xu)
            hf = (η_max - sw.zx_max[i])

            if hf > sw.h_thresh
                length = T(0.5) * (sw.xl[i] + sw.xl[xu]) # can be precalculated
                sw.qx[i] = local_inertial_flow(
                    sw.θ,
                    sw.qx0[i],
                    sw.qx0[xd],
                    sw.qx0[xu],
                    η_x,
                    η_xu,
                    hf,
                    sw.ywidth[i],
                    length,
                    sw.mannings_n_sq[i],
                    sw.g,
                    sw.froude_limit,
                    Δt,
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

            η_y = sw.z[i] + sw.h[i]
            η_yu = sw.z[yu] + sw.h[yu]
            η_max = max(η_y, η_yu)
            hf = (η_max - sw.zy_max[i])

            if hf > sw.h_thresh
                length = T(0.5) * (sw.yl[i] + sw.yl[yu]) # can be precalculated
                sw.qy[i] = local_inertial_flow(
                    sw.θ,
                    sw.qy0[i],
                    sw.qy0[yd],
                    sw.qy0[yu],
                    η_y,
                    η_yu,
                    hf,
                    sw.xwidth[i],
                    length,
                    sw.mannings_n_sq[i],
                    sw.g,
                    sw.froude_limit,
                    Δt,
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
    @batch per = thread minbatch = 6000 for i = 1:sw.n
        yd = indices.yd[i]
        xd = indices.xd[i]

        if sw.rivercells[i]
            if swr.waterbody[inds_riv[i]]
                # for reservoir or lake set inflow from land part, these are boundary points
                # and update of volume and h is not required
                swr.inflow_wb[inds_riv[i]] =
                    sw.runoff[i] + (sw.qx[xd] - sw.qx[i] + sw.qy[yd] - sw.qy[i])
            else
                sw.volume[i] +=
                    (
                        sum_at(swr.q, links_at_node.src[inds_riv[i]]) -
                        sum_at(swr.q, links_at_node.dst[inds_riv[i]]) + sw.qx[xd] -
                        sw.qx[i] + sw.qy[yd] - sw.qy[i] +
                        swr.inflow[inds_riv[i]] +
                        sw.runoff[i]
                    ) * Δt
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
                swr.h_av[inds_riv[i]] += swr.h[inds_riv[i]] * Δt
            end
        else
            sw.volume[i] +=
                (sw.qx[xd] - sw.qx[i] + sw.qy[yd] - sw.qy[i] + sw.runoff[i]) * Δt
            if sw.volume[i] < T(0.0)
                sw.error[i] = sw.error[i] + abs(sw.volume[i])
                sw.volume[i] = T(0.0) # set volume to zero
            end
            sw.h[i] = sw.volume[i] / (sw.xl[i] * sw.yl[i])
        end
        sw.h_av[i] += sw.h[i] * Δt
    end
end

""" 
    FloodPlainProfile

Floodplain `volume` is a function of `depth` (flood depth intervals). Based on the
cumulative floodplain `volume` a floodplain profile as a function of `flood_depth` is
derived with floodplain area `a` (cumulative) and wetted perimeter radius `p` (cumulative).
"""
@get_units @exchange @grid_type @grid_location @with_kw struct FloodPlainProfile{T,N}
    depth::Vector{T} | "m" | 0                # Flood depth
    volume::Array{T,2} | "m3"                 # Flood volume (cumulative)
    width::Array{T,2} | "m"                   # Flood width
    a::Array{T,2} | "m2"                      # Flow area (cumulative)
    p::Array{T,2} | "m"                       # Wetted perimeter (cumulative)
end

@get_units @exchange @grid_type @grid_location @with_kw struct FloodPlain{T,P}
    profile::P | "-" | 0                                    # floodplain profile
    mannings_n::Vector{T} | "s m-1/3"                       # manning's roughness
    mannings_n_sq::Vector{T} | "(s m-1/3)2" | _ | "edge"    # manning's roughness squared
    volume::Vector{T} | "m3"                                # volume
    h::Vector{T} | "m"                                      # water depth
    h_av::Vector{T} | "m"                                   # average water depth
    error::Vector{T} | "m3"                                 # error volume
    a::Vector{T} | "m2" | _ | "edge"                        # flow area
    r::Vector{T} | "m" | _ | "edge"                         # hydraulic radius
    hf::Vector{T} | "m" | _ | "edge"                        # water depth at edge/link
    zb_max::Vector{T} | "m" | _ | "edge"                    # maximum bankfull elevation (edge/link)
    q0::Vector{T} | "m3 s-1" | _ | "edge"                   # discharge at previous time step
    q::Vector{T} | "m3 s-1" | _ | "edge"                    # discharge
    q_av::Vector{T} | "m" | _ | "edge"                      # average river discharge
    hf_index::Vector{Int} | "-" | _ | "edge"                # index with `hf` above depth threshold
end

"Determine the initial floodplain volume"
function initialize_volume!(river, nriv::Int)
    for i = 1:nriv
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
    Δh = h - depth  # depth at i1
    area = area + (width * Δh) # area at i1, width at i2
    return area
end

"""
    function wetted_perimeter(p, depth, h)

Compute floodplain wetted perimeter based on flow depth `h` and floodplain `depth` and
wetted perimeter `p` of a floodplain profile.
"""
function wetted_perimeter(p, depth, h)
    Δh = h - depth # depth at i1
    p = p + (2.0 * Δh) # p at i1
    return p
end

"Compute flood depth by interpolating flood volume `flood_volume` using flood depth intervals."
function flood_depth(profile::FloodPlainProfile{T}, flood_volume, dl, i::Int)::T where {T}
    i1, i2 = interpolation_indices(flood_volume, @view profile.volume[:, i])
    ΔA = (flood_volume - profile.volume[i1, i]) / dl
    Δh = ΔA / profile.width[i2, i]
    flood_depth = profile.depth[i1] + Δh
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
    for i = 1:n
        riv_cell = 0
        diff_volume = diff(volume[:, i])

        for j = 1:(n_depths-1)
            # assume rectangular shape of flood depth segment
            width[j+1, i] = diff_volume[j] / (h[j] * riverlength[i])
            # check provided flood volume (floodplain width should be constant or increasing
            # as a function of flood depth)
            if width[j+1, i] < width[j, i]
                # raise warning only if difference is larger than rounding error of 0.01 m³
                if ((width[j, i] - width[j+1, i]) * h[j] * riverlength[i]) > 0.01
                    incorrect_vol += 1
                    riv_cell = 1
                    error_vol =
                        error_vol + ((width[j, i] - width[j+1, i]) * h[j] * riverlength[i])
                end
                width[j+1, i] = width[j, i]
            end
            a[j+1, i] = width[j+1, i] * h[j]
            p[j+1, i] = (width[j+1, i] - width[j, i]) + 2.0 * h[j]
            segment_volume[j+1, i] = a[j+1, i] * riverlength[i]
            if j == 1
                # for interpolation wetted perimeter at flood depth 0.0 is required
                p[j, i] = p[j+1, i] - 2.0 * h[j]
            end
        end

        p[2:end, i] = cumsum(p[2:end, i])
        a[:, i] = cumsum(a[:, i])
        volume[:, i] = cumsum(segment_volume[:, i])

        riv_cells += riv_cell
    end

    if incorrect_vol > 0
        perc_riv_cells = round(100.0 * (riv_cells / n), digits = 2)
        perc_error_vol = round(100.0 * (error_vol / sum(start_volume[end, :])), digits = 2)
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
    profile = FloodPlainProfile{Float,n_depths}(
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
    for i = 1:n_edges
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

    floodplain = FloodPlain(
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
    set_river_inwater(model::Model{N,L,V,R,W,T}, ssf_toriver) where {N,L,V<:SBM,R,W,T}

Set `inwater` of the river component for a `Model` with vertical `SBM` concept.
`ssf_toriver` is the subsurface flow to the river.
"""
function set_river_inwater(model::Model{N,L,V,R,W,T}, ssf_toriver) where {N,L,V<:SBM,R,W,T}
    @unpack lateral, vertical, network = model
    inds = network.index_river

    @. lateral.river.inwater = (
        ssf_toriver[inds] +
        lateral.land.to_river[inds] +
        # net_runoff_river
        (
            (
                vertical.net_runoff_river[inds] *
                network.land.xl[inds] *
                network.land.yl[inds] *
                0.001
            ) / vertical.Δt
        )
    )
end

"""
    set_river_inwater(model, ssf_toriver)

Set `inwater` of the river component (based on overland flow).
"""
function set_river_inwater(model, ssf_toriver)
    @unpack lateral, network = model
    inds = network.index_river
    lateral.river.inwater .= lateral.land.to_river[inds]
end

"""
    set_land_inwater(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmGwfModel}

Set `inwater` of the land component for the `SbmGwgModel` type.
"""
function set_land_inwater(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmGwfModel}
    @unpack lateral, vertical, network, config = model

    do_drains = get(config.model, "drains", false)::Bool
    drainflux = zeros(vertical.n)
    if do_drains
        drainflux[lateral.subsurface.drain.index] =
            -lateral.subsurface.drain.flux ./ tosecond(basetimestep)
    end

    lateral.land.inwater .=
        (vertical.runoff .* network.land.xl .* network.land.yl .* 0.001) ./
        lateral.land.Δt .+ drainflux
end

"""
    set_land_inwater(model)

Set `inwater` of the land component, based on `runoff` of the `vertical` concept.
"""
function set_land_inwater(model)
    @unpack lateral, vertical, network = model
    lateral.land.inwater .=
        (vertical.runoff .* network.land.xl .* network.land.yl .* 0.001) ./ lateral.land.Δt
end

"""
    get_inflow_waterbody(model)

Get inflow to a water body (reservoir or lake) `inflow_wb` based on overland flow.
"""
function get_inflow_waterbody(model)
    @unpack lateral, network = model
    inds = network.index_river
    if !isnothing(lateral.river.reservoir) || !isnothing(lateral.river.lake)
        inflow_wb = lateral.land.q_av[inds]
    else
        inflow_wb = nothing
    end
    return inflow_wb
end

"""
    get_inflow_waterbody(
        model::Model{N,L,V,R,W,T},
    ) where {N,L<:NamedTuple{<:Any,<:Tuple{LateralSSF,SurfaceFlow,Any}},V,R,W,T}

Get inflow to a water body (reservoir or lake) `inflow_wb` from a model type that contains
the lateral components `LateralSSF` and `SurfaceFlow`.
"""
function get_inflow_waterbody(
    model::Model{N,L,V,R,W,T},
) where {N,L<:NamedTuple{<:Any,<:Tuple{LateralSSF,SurfaceFlow,Any}},V,R,W,T}
    @unpack lateral, network = model

    inds = network.index_river
    if !isnothing(lateral.river.reservoir) || !isnothing(lateral.river.lake)
        inflow_wb =
            lateral.subsurface.ssf[inds] ./ tosecond(basetimestep) .+
            lateral.land.q_av[inds]
    else
        inflow_wb = nothing
    end
    return inflow_wb
end

"""
    get_inflow_waterbody(
        model::Model{N,L,V,R,W,T},
    ) where {N,L<:NamedTuple{<:Any,<:Tuple{LateralSSF,ShallowWaterLand,Any}},V,R,W,T}

Get inflow to a water body (reservoir or lake) `inflow_wb` from a model type that contains
the lateral components `LateralSSF` and `ShallowWaterLand`.
"""
function get_inflow_waterbody(
    model::Model{N,L,V,R,W,T},
) where {N,L<:NamedTuple{<:Any,<:Tuple{LateralSSF,ShallowWaterLand,Any}},V,R,W,T}
    @unpack lateral, network = model

    inds = network.index_river
    if !isnothing(lateral.river.reservoir) || !isnothing(lateral.river.lake)
        inflow_wb = lateral.subsurface.ssf[inds] ./ tosecond(basetimestep)
    else
        inflow_wb = nothing
    end
    return inflow_wb
end

"""
    surface_routing(model; ssf_toriver = 0.0)

Run surface routing (land and river). Kinematic wave for overland flow and kinematic wave or
local inertial model for river flow.
"""
function surface_routing(model; ssf_toriver = 0.0)
    @unpack lateral, network, clock = model

    # run kinematic wave for overland flow
    set_land_inwater(model)
    update(lateral.land, network.land, network.frac_toriver)

    # run river flow
    set_river_inwater(model, ssf_toriver)
    update(
        lateral.river,
        network.river,
        get_inflow_waterbody(model),
        julian_day(clock.time - clock.Δt),
    )
end

"""
    surface_routing(
        model::Model{N,L,V,R,W,T};
        ssf_toriver = 0.0,
    ) where {N,L<:NamedTuple{<:Any,<:Tuple{Any,ShallowWaterLand,ShallowWaterRiver}},V,R,W,T}

Run surface routing (land and river) for a model type that contains the lateral components
`ShallowWaterLand` and `ShallowWaterRiver`.
"""
function surface_routing(
    model::Model{N,L,V,R,W,T};
    ssf_toriver = 0.0,
) where {N,L<:NamedTuple{<:Any,<:Tuple{Any,ShallowWaterLand,ShallowWaterRiver}},V,R,W,T}

    @unpack lateral, vertical, network, clock = model

    @. lateral.land.runoff = (
        (vertical.runoff / 1000.0) * (network.land.xl * network.land.yl) / vertical.Δt +
        ssf_toriver +
        # net_runoff_river
        (
            (vertical.net_runoff_river * network.land.xl * network.land.yl * 0.001) /
            vertical.Δt
        )
    )

    update(
        lateral.land,
        lateral.river,
        network,
        get_inflow_waterbody(model),
        julian_day(clock.time - clock.Δt),
    )
end
