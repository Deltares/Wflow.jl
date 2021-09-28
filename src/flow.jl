@get_units @with_kw struct SurfaceFlow{T,R,L}
    β::T | "-"                              # constant in Manning's equation
    sl::Vector{T} | "m m-1"                 # Slope [m m⁻¹]
    n::Vector{T} | "s m-1/3"                # Manning's roughness [s m⁻⅓]
    dl::Vector{T} | "m"                     # Drain length [m]
    q::Vector{T} | "m3 s-1"                 # Discharge [m³ s⁻¹]
    qin::Vector{T} | "m3 s-1"               # Inflow from upstream cells [m³ s⁻¹]
    q_av::Vector{T} | "m3 s-1"              # Average discharge [m³ s⁻¹]
    qlat::Vector{T} | "m2 s-1"              # Lateral inflow per unit length [m² s⁻¹]
    inwater::Vector{T} | "m3 s-1"           # Lateral inflow [m³ s⁻¹]
    inflow::Vector{T} | "m3 s-1"            # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    volume::Vector{T} | "m3"                # Kinematic wave volume [m³] (based on water level h)
    h::Vector{T} | "m"                      # Water level [m]
    h_av::Vector{T} | "m"                   # Average water level [m]
    h_bankfull::Vector{T} | "m"             # Bankfull water level [m]
    Δt::T | "s"                             # Model time step [s]
    its::Int | "-"                          # Number of fixed iterations
    width::Vector{T} | "m"                  # Flow width [m]
    alpha_pow::T | "-"                      # Used in the power part of α
    alpha_term::Vector{T} | "-"             # Term used in computation of α
    α::Vector{T} | "s3/5 m1/5"              # Constant in momentum equation A = αQᵝ, based on Manning's equation
    cel::Vector{T} | "m s-1"                # Celerity of the kinematic wave
    to_river::Vector{T} | "m3 s-1"          # Part of overland flow [m³ s⁻¹] that flows to the river
    rivercells::Vector{Bool} | "-"          # Location of river cells (0 or 1)
    wb_pit::Vector{Bool} | "-"              # Boolean location (0 or 1) of a waterbody (wb, reservoir or lake).
    reservoir_index::Vector{Int} | "-"      # map cell to 0 (no reservoir) or i (pick reservoir i in reservoir field)
    lake_index::Vector{Int} | "-"           # map cell to 0 (no lake) or i (pick lake i in lake field)
    reservoir::R                            # Reservoir model struct of arrays
    lake::L                                 # Lake model struct of arrays
    kinwave_it::Bool                        # Boolean for iterations kinematic wave

    # TODO unclear why this causes a MethodError
    # function SurfaceFlow{T,R,L}(args...) where {T,R,L}
    #     equal_size_vectors(args)
    #     return new(args...)
    # end
end

statevars(::SurfaceFlow) = (:q, :h, :h_av)

function update(
    sf::SurfaceFlow,
    network;
    frac_toriver = nothing,
    inflow_wb = nothing,
    doy = 0,
)
    @unpack graph,
    order,
    subdomain_order,
    topo_subdomain,
    indices_subdomain,
    upstream_nodes = network

    n = length(order)
    ns = length(subdomain_order)

    @. sf.alpha_term = pow(sf.n / sqrt(sf.sl), sf.β)
    # use fixed alpha value based on 0.5 * h_bankfull
    @. sf.α = sf.alpha_term * pow(sf.width + sf.h_bankfull, sf.alpha_pow)

    # two options for iteration, fixed or based on courant number.
    if sf.kinwave_it
        if sf.its > 0
            its = sf.its
        else
            # calculate celerity
            courant = zeros(n)
            for k = 1:ns
                for m in subdomain_order[k]
                    for v in topo_subdomain[m]
                        if sf.q[v] > 0.0
                            sf.cel[v] = 1.0 / (sf.α[v] * sf.β * pow(sf.q[v], (sf.β - 1.0)))
                            courant[v] = (sf.Δt / sf.dl[v]) * sf.cel[v]
                        end
                    end
                end
            end
            filter!(x -> x ≠ 0.0, courant)
            its = isempty(courant) ? 1 : ceil(Int, (1.25 * quantile!(courant, 0.95)))
        end
    else
        its = 1
    end

    # sub time step
    adt = sf.Δt / its

    sf.q_av .= 0.0
    sf.h_av .= 0.0
    sf.to_river .= 0.0
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

    for _ = 1:its
        sf.qin .= 0.0
        for k = 1:ns
            @threads for m in subdomain_order[k]
                for (n, v) in zip(indices_subdomain[m], topo_subdomain[m])

                    # for overland flow frac_toriver needs to be defined
                    if frac_toriver !== nothing # run kinematic wave for land domain
                        # for a river cell without a reservoir or lake (wb_pit is false) part of the
                        # upstream surface flow goes to the river (frac_toriver) and part goes to
                        # the surface flow reservoir (1.0 - frac_toriver), upstream nodes with a
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
                    else # run kinematic wave for river domain (including reservoirs and lakes)
                        # sf.qin by outflow from upstream reservoir or lake location is added
                        sf.qin[v] += sum_at(sf.q, upstream_nodes[n])
                    end

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
                        adt,
                        sf.dl[v],
                    )

                    if !isnothing(sf.reservoir) && sf.reservoir_index[v] != 0
                        # run reservoir model and copy reservoir outflow to inflow (qin) of
                        # downstream river cell
                        i = sf.reservoir_index[v]
                        update(sf.reservoir, i, sf.q[v] + inflow_wb[v], adt)

                        downstream_nodes = outneighbors(graph, v)
                        n_downstream = length(downstream_nodes)
                        if n_downstream == 1
                            j = only(downstream_nodes)
                            sf.qin[j] = sf.reservoir.outflow[i]
                        elseif n_downstream == 0
                            error(
                                """A reservoir without a downstream river node is not supported. 
                                Add a downstream river node or move the reservoir to an upstream node (model schematization).
                                """)
                        else
                            error("bifurcations not supported")
                        end

                    elseif !isnothing(sf.lake) && sf.lake_index[v] != 0
                        # run lake model and copy lake outflow to inflow (qin) of downstream river
                        # cell
                        i = sf.lake_index[v]
                        update(sf.lake, i, sf.q[v] + inflow_wb[v], doy, adt)

                        downstream_nodes = outneighbors(graph, v)
                        n_downstream = length(downstream_nodes)
                        if n_downstream == 1
                            j = only(downstream_nodes)
                            sf.qin[j] = sf.lake.outflow[i]
                        elseif n_downstream == 0
                            error(
                                """A lake without a downstream river node is not supported. 
                                Add a downstream river node or move the lake to an upstream node (model schematization).
                                """)
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
    sf.to_river ./= its
    sf.volume .= sf.dl .* sf.width .* sf.h
end

@get_units @with_kw struct LateralSSF{T}
    kh₀::Vector{T} | "m Δt-1"              # Horizontal hydraulic conductivity at soil surface [m Δt⁻¹]
    f::Vector{T} | "m-1"                   # A scaling parameter [m⁻¹] (controls exponential decline of kh₀)
    soilthickness::Vector{T} | "m"         # Soil thickness [m]
    θₛ::Vector{T} | "-"                     # Saturated water content (porosity) [-]
    θᵣ::Vector{T} | "-"                    # Residual water content [-]
    t::T | "Δt s"                          # time step [Δt s]
    Δt::T | "s"                            # model time step [s]
    βₗ::Vector{T} | "m m-1"                 # Slope [m m⁻¹]
    dl::Vector{T} | "m"                    # Drain length [m]
    dw::Vector{T} | "m"                    # Flow width [m]
    zi::Vector{T} | "m"                    # Pseudo-water table depth [m] (top of the saturated zone)
    exfiltwater::Vector{T} | "m Δt-1"      # Exfiltration [m Δt⁻¹] (groundwater above surface level, saturated excess conditions)
    recharge::Vector{T} | "m Δt-1"         # Net recharge to saturated store [m Δt⁻¹]
    ssf::Vector{T} | "m3 Δt-1"             # Subsurface flow [m³ Δt⁻¹]
    ssfin::Vector{T} | "m3 Δt-1"           # Inflow from upstream cells [m³ Δt⁻¹]
    ssfmax::Vector{T} | "m2 Δt-1"          # Maximum subsurface flow [m² Δt⁻¹]
    to_river::Vector{T} | "m3 Δt-1"        # Part of subsurface flow [m³ Δt⁻¹] that flows to the river
    wb_pit::Vector{Bool} | "-"             # Boolean location (0 or 1) of a waterbody (wb, reservoir or lake).

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
        @threads for m in subdomain_order[k]
            for (n, v) in zip(indices_subdomain[m], topo_subdomain[m])
                # for a river cell without a reservoir or lake (wb_pit is false) part of the
                # upstream subsurface flow goes to the river (frac_toriver) and part goes to the
                # subsurface flow reservoir (1.0 - frac_toriver) upstream nodes with a reservoir or
                # lake are excluded
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
                    ssf.t,
                    ssf.dl[v],
                    ssf.dw[v],
                    ssf.ssfmax[v],
                )
            end
        end
    end
end

@get_units @with_kw struct GroundwaterExchange{T}
    Δt::T | "s"                         # model time step [s]
    exfiltwater::Vector{T} | "m Δt-1"   # Exfiltration [m Δt⁻¹]  (groundwater above surface level, saturated excess conditions)
    zi::Vector{T} | "m"                 # Pseudo-water table depth [m] (top of the saturated zone)
    to_river::Vector{T} | "m3 Δt-1"     # Part of subsurface flow [m³ Δt⁻¹] that flows to the river
    ssf::Vector{T} | "m3 Δt-1"          # Subsurface flow [m³ Δt⁻¹]
end

@get_units @with_kw struct ShallowWaterRiver{T,R,L}
    n::Int | "-"                        # number of cells
    ne::Int | "-"                       # number of edges/links
    g::T | "m2 s-1"                     # acceleration due to gravity
    α::T | "-"                          # stability coefficient (de Almeida et al., 2012.)
    h_thresh::T | "m"                   # stability coefficient (de Almeida et al., 2012.)
    Δt::T | "s"                         # model time step [s]
    q::Vector{T} | "m3 s-1"             # river discharge (subgrid channel)
    q_av::Vector{T} | "m3 s-1"          # average discharge [m³ s⁻¹]
    zmax::Vector{T} | "m"               # maximum channel bed elevation
    mannings_n::Vector{T} | "s m-1/3"   # Manning's roughness at edge/link
    h::Vector{T} | "m"                  # water depth
    η_max::Vector{T} | "m"              # maximum water elevation
    hf::Vector{T} | "m"                 # water depth at edge/link
    h_av::Vector{T} | "m"               # average water depth
    length::Vector{T} | "m"             # river length
    length_at_link::Vector{T} | "m"     # river length at edge/link
    width::Vector{T} | "m"              # river width
    width_at_link::Vector{T} | "m"      # river width at edge/link
    a::Vector{T} | "m2"                 # flow area at edge/link
    r::Vector{T} | "m"                  # wetted perimeter at edge/link
    volume::Vector{T} | "m3"            # river volume
    error::Vector{T} | "m3"             # error volume
    inwater::Vector{T} | "m3 s-1"       # Lateral inflow [m³ s⁻¹]
    bankvolume::Vector{T} | "m3"        # bank volume
    bankheight::Vector{T} | "m"         # bank height
    slp::Vector{T} | "-"                # river slope between river cells
    z::Vector{T} | "m"                  # river bed elevation
    froude::Bool | "-"                  # if true a check is performed if froude number > 1.0 (algorithm is modified)
    reservoir_index::Vector{Int} | "-"  # map cell to 0 (no reservoir) or i (pick reservoir i in reservoir field)
    lake_index::Vector{Int} | "-"       # map cell to 0 (no lake) or i (pick lake i in lake field)
    reservoir::R                        # Reservoir model struct of arrays
    lake::L                             # Lake model struct of arrays
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

    @threads for i = 1:sw.ne
        ηsrc = sw.z[nodes_at_link.src[i]] + sw.h[nodes_at_link.src[i]]
        ηdst = sw.z[nodes_at_link.dst[i]] + sw.h[nodes_at_link.dst[i]]

        sw.η_max[i] = max(ηsrc, ηdst)
        sw.hf[i] = (sw.η_max[i] - sw.zmax[i])

        if sw.hf[i] > sw.h_thresh
            sw.a[i] = sw.width_at_link[i] * sw.hf[i] # cross area (rectangular channel)
            sw.r[i] = sw.a[i] / (sw.width_at_link[i] + 2.0 * sw.hf[i]) # wetted perimeter (rectangular channel)
            sw.q[i] = local_inertial_riverflow(
                sw.q[i],
                ηsrc,
                ηdst,
                sw.hf[i],
                sw.a[i],
                sw.r[i],
                sw.length_at_link[i],
                sw.mannings_n[i],
                sw.g,
                sw.froude,
                Δt,
            )
        else
            sw.q[i] = 0.0
        end

        if !isnothing(sw.reservoir) && sw.reservoir_index[i] != 0
            sw.q[i] = max(sw.q[i], 0.0)
            v = sw.reservoir_index[i]
            update(sw.reservoir, v, sw.q[i] + inflow_wb[i], Δt)
            sw.inwater[nodes_at_link.dst[i]] = sw.reservoir.outflow[v]
        elseif !isnothing(sw.lake) && sw.lake_index[i] != 0
            sw.q[i] = max(sw.q[i], 0.0)
            v = sw.lake_index[i]
            update(sw.lake, v, sw.q[i] + inflow_wb[i], doy, Δt)
            sw.inwater[nodes_at_link.dst[i]] = sw.lake.outflow[v]
        end
        sw.q_av[i] += sw.q[i] * Δt
    end
    if update_h
        @threads for i = 1:sw.n
            sw.volume[i] =
                sw.volume[i] +
                (
                    sum_at(sw.q, links_at_node.src[i]) -
                    sum_at(sw.q, links_at_node.dst[i]) + sw.inwater[i]
                ) * Δt
            if sw.volume[i] < 0.0
                sw.error[i] = sw.error[i] + abs(sw.volume[i])
            end
            sw.volume[i] = max(sw.volume[i], 0.0) # set volume to zero if negative
            sw.h[i] = sw.volume[i] / (sw.length[i] * sw.width[i])
            sw.h_av[i] += sw.h[i] * Δt
        end
    end
end

function update(
    sw::ShallowWaterRiver,
    network;
    inflow_wb = nothing,
    doy = 0,
    update_h = true,
)
    @unpack nodes_at_link, links_at_node = network

    if !isnothing(sw.reservoir)
        sw.reservoir.inflow .= 0.0
        sw.reservoir.totaloutflow .= 0.0
    end
    if !isnothing(sw.lake)
        sw.lake.inflow .= 0.0
        sw.lake.totaloutflow .= 0.0
    end
    sw.q_av .= 0.0
    sw.h_av .= 0.0

    t = 0.0
    while t < sw.Δt
        Δt = stable_timestep(sw)
        shallowwater_river_update(sw, network, Δt, inflow_wb, doy, update_h)
        if t + Δt > sw.Δt
            Δt = sw.Δt - t
        end
        t = t + Δt
    end
    sw.q_av ./= sw.Δt
    sw.h_av ./= sw.Δt
end

function stable_timestep(sw::ShallowWaterRiver)
    Δtₘᵢₙ = Inf
    for i = 1:sw.n
        Δt = sw.α * sw.length[i] / sqrt.(sw.g .* sw.h[i])
        Δtₘᵢₙ = Δt < Δtₘᵢₙ ? Δt : Δtₘᵢₙ
    end
    Δtₘᵢₙ = isinf(Δtₘᵢₙ) ? 10.0 : Δtₘᵢₙ
    return Δtₘᵢₙ
end
