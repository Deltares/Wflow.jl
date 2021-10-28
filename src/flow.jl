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

function initialize_surfaceflow_land(
    nc,
    config,
    inds;
    sl,
    dl,
    width,
    wb_pit,
    iterate,
    tstep,
    Δt,
)

    n_land = ncread(
        nc,
        param(config, "input.lateral.land.n", nothing);
        sel = inds,
        defaults = 0.072,
        type = Float,
    )
    n = length(inds)

    sf_land = SurfaceFlow(
        β = Float(0.6),
        sl = sl,
        n = n_land,
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
        h_bankfull = zeros(Float, n),
        Δt = Float(tosecond(Δt)),
        its = tstep > 0 ? Int(cld(tosecond(Δt), tstep)) : tstep,
        width = width,
        wb_pit = wb_pit,
        alpha_pow = Float((2.0 / 3.0) * 0.6),
        alpha_term = fill(mv, n),
        α = fill(mv, n),
        cel = zeros(Float, n),
        to_river = zeros(Float, n),
        rivercells = fill(false, n),
        reservoir_index = fill(0, n),
        lake_index = fill(0, n),
        reservoir = nothing,
        lake = nothing,
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
    wb_pit,
    reservoir_index,
    reservoir,
    lake_index,
    lake,
    river,
    iterate,
    tstep,
    Δt,
)

    n_river = ncread(
        nc,
        param(config, "input.lateral.river.n", nothing);
        sel = inds,
        defaults = 0.036,
        type = Float,
    )
    h_bankfull = ncread(
        nc,
        param(config, "input.lateral.river.h_bankfull", nothing);
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    sl = ncread(nc, param(config, "input.lateral.river.slope"); sel = inds, type = Float)
    clamp!(sl, 0.00001, Inf)

    n = length(inds)

    sf_river = SurfaceFlow(
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
        h_bankfull = h_bankfull,
        Δt = Float(tosecond(Δt)),
        its = tstep > 0 ? ceil(Int(tosecond(Δt) / tstep)) : tstep,
        width = width,
        wb_pit = wb_pit,
        alpha_pow = Float((2.0 / 3.0) * 0.6),
        alpha_term = fill(mv, n),
        α = fill(mv, n),
        cel = zeros(Float, n),
        to_river = zeros(Float, n),
        reservoir_index = reservoir_index,
        lake_index = lake_index,
        reservoir = reservoir,
        lake = lake,
        rivercells = river,
        kinwave_it = iterate,
    )

    return sf_river
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

    @. sf.qlat .= sf.inwater ./ sf.dl

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
                                """,
                            )
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
    g::T | "m s-2"                      # acceleration due to gravity
    α::T | "-"                          # stability coefficient (Bates et al., 2010)
    h_thresh::T | "m"                   # depth threshold for calculating flow
    Δt::T | "s"                         # model time step [s]
    q::Vector{T} | "m3 s-1"             # river discharge (subgrid channel)
    q_av::Vector{T} | "m3 s-1"          # average river discharge [m³ s⁻¹]
    zmax::Vector{T} | "m"               # maximum channel bed elevation
    mannings_n::Vector{T} | "s m-1/3"   # Manning's roughness at edge/link
    h::Vector{T} | "m"                  # water depth
    η_max::Vector{T} | "m"              # maximum water elevation
    hf::Vector{T} | "m"                 # water depth at edge/link
    h_av::Vector{T} | "m"               # average water depth
    dl::Vector{T} | "m"                 # river length
    dl_at_link::Vector{T} | "m"         # river length at edge/link
    width::Vector{T} | "m"              # river width
    width_at_link::Vector{T} | "m"      # river width at edge/link
    a::Vector{T} | "m2"                 # flow area at edge/link
    r::Vector{T} | "m"                  # wetted perimeter at edge/link
    volume::Vector{T} | "m3"            # river volume
    error::Vector{T} | "m3"             # error volume
    inwater::Vector{T} | "m3 s-1"       # lateral inflow [m³ s⁻¹]
    inwater0::Vector{T} | "m3 s-1"      # lateral inflow at previous time step [m³ s⁻¹]
    inflow::Vector{T} | "m3 s-1"        # External inflow (abstraction/supply/demand) [m³ s⁻¹]
    bankvolume::Vector{T} | "m3"        # bank volume
    bankheight::Vector{T} | "m"         # bank height
    z::Vector{T} | "m"                  # river bed elevation
    froude_limit::Bool | "-"            # if true a check is performed if froude number > 1.0 (algorithm is modified)
    reservoir_index::Vector{Int} | "-"  # map cell to 0 (no reservoir) or i (pick reservoir i in reservoir field)
    lake_index::Vector{Int} | "-"       # map cell to 0 (no lake) or i (pick lake i in lake field)
    reservoir::R                        # Reservoir model struct of arrays
    lake::L                             # Lake model struct of arrays
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
    riverlength_bc = get(config.model, "riverlength_bc", 1.0e5) # river length at boundary point (ghost point)
    alpha = get(config.model, "inertial_flow_alpha", 0.7) # stability coefficient for model time step (0.2-0.7)
    h_thresh = get(config.model, "h_thresh", 1.0e-03) # depth threshold for flow at link
    froude_limit = get(config.model, "froude_limit", true) # limit flow to subcritical according to Froude number

    river_elevation_2d =
        ncread(nc, param(config, "input.lateral.river.elevation"); type = Float, fill = 0)
    river_elevation = river_elevation_2d[inds]
    n_river = ncread(
        nc,
        param(config, "input.lateral.river.n", nothing);
        sel = inds,
        defaults = 0.036,
        type = Float,
    )

    n = length(inds)

    # set ghost points for boundary condition (downstream river outlet): river width and
    # manning n is copied from the upstream cell, river elevation and h are set at 0.0 (sea
    # level). river length at boundary point is by default 1.0e5 m (riverlength_bc).
    index_pit = findall(x -> x == 5, ldd)
    npits = length(index_pit)
    add_vertex_edge_graph!(graph, index_pit)
    append!(river_elevation, fill(Float(0), npits))
    append!(dl, fill(riverlength_bc, npits))
    append!(width, width[index_pit])
    append!(n_river, n_river[index_pit])

    # for each link the src and dst node is required
    nodes_at_link = adjacent_nodes_at_link(graph)
    _ne = ne(graph)

    # determine z, width, length and manning's n at links
    zmax = fill(Float(0), _ne)
    width_at_link = fill(Float(0), _ne)
    length_at_link = fill(Float(0), _ne)
    mannings_n = fill(Float(0), _ne)
    for i = 1:_ne
        zmax[i] = max(
            river_elevation[nodes_at_link.src[i]],
            river_elevation[nodes_at_link.dst[i]],
        )
        width_at_link[i] = min(width[nodes_at_link.dst[i]], width[nodes_at_link.src[i]])
        length_at_link[i] = 0.5 * (dl[nodes_at_link.dst[i]] + dl[nodes_at_link.src[i]])
        mannings_n[i] =
            (
                n_river[nodes_at_link.dst[i]] * dl[nodes_at_link.dst[i]] +
                n_river[nodes_at_link.src[i]] * dl[nodes_at_link.src[i]]
            ) / (dl[nodes_at_link.dst[i]] + dl[nodes_at_link.src[i]])
    end

    sw_river = ShallowWaterRiver(
        n = n,
        ne = _ne,
        g = 9.80665,
        α = alpha,
        h_thresh = h_thresh,
        Δt = tosecond(Δt),
        q = zeros(_ne),
        q_av = zeros(_ne),
        zmax = zmax,
        mannings_n = mannings_n,
        h = fill(0.0, n + length(index_pit)),
        η_max = zeros(_ne),
        hf = zeros(_ne),
        h_av = zeros(n),
        width = width,
        width_at_link = width_at_link,
        a = zeros(_ne),
        r = zeros(_ne),
        volume = fill(0.0, n),
        error = zeros(Float, n),
        inflow = zeros(n),
        inwater = zeros(n),
        inwater0 = fill(mv, n),
        dl = dl,
        dl_at_link = length_at_link,
        bankvolume = fill(mv, n),
        bankheight = fill(mv, n),
        z = river_elevation,
        froude_limit = froude_limit,
        reservoir_index = reservoir_index,
        lake_index = lake_index,
        reservoir = reservoir,
        lake = lake,
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
                sw.dl_at_link[i],
                sw.mannings_n[i],
                sw.g,
                sw.froude_limit,
                Δt,
            )
        else
            sw.q[i] = 0.0
        end

        if !isnothing(sw.reservoir) && sw.reservoir_index[i] != 0
            sw.q[i] = max(sw.q[i], 0.0)
            v = sw.reservoir_index[i]
            update(sw.reservoir, v, sw.q[i] + inflow_wb[i], Δt)
            # add lake outflow to inwater of destination node
            sw.inwater[nodes_at_link.dst[i]] =
                sw.inwater0[nodes_at_link.dst[i]] + sw.reservoir.outflow[v]
        elseif !isnothing(sw.lake) && sw.lake_index[i] != 0
            sw.q[i] = max(sw.q[i], 0.0)
            v = sw.lake_index[i]
            update(sw.lake, v, sw.q[i] + inflow_wb[i], doy, Δt)
            # add reservoir outflow to inwater of destination node
            sw.inwater[nodes_at_link.dst[i]] =
                sw.inwater0[nodes_at_link.dst[i]] + sw.lake.outflow[v]
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
                sw.volume[i] = 0.0 # set volume to zero
            end
            sw.volume[i] = max(sw.volume[i] + sw.inflow[i] * Δt, 0.0) # add external inflow
            sw.h[i] = sw.volume[i] / (sw.dl[i] * sw.width[i])
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
    if !isnothing(sw.reservoir) || !isnothing(sw.lake)
        sw.inwater0 .= sw.inwater
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


"""
    stable_timestep(sw::ShallowWaterRiver)

Compute a stable timestep size for the local inertial approach for river flow, based on
Bates et al. (2010).

Δt = α * (Δx / sqrt(g max(h))
"""

function stable_timestep(sw::ShallowWaterRiver)
    Δtₘᵢₙ = Inf
    for i = 1:sw.n
        Δt = sw.α * sw.dl[i] / sqrt.(sw.g .* sw.h[i])
        Δtₘᵢₙ = Δt < Δtₘᵢₙ ? Δt : Δtₘᵢₙ
    end
    Δtₘᵢₙ = isinf(Δtₘᵢₙ) ? 10.0 : Δtₘᵢₙ
    return Δtₘᵢₙ
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
        ssf_toriver .+ lateral.land.to_river[inds] +
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
    get_inflow_waterbody(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}

Get inflow to a water body (reservoir or lake) `inflow_wb` based on overland and lateral
subsurface flow for the `SbmModel` model type.
"""
function get_inflow_waterbody(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:SbmModel}
    @unpack lateral, network = model

    inds = network.index_river
    if !isnothing(lateral.river.reservoir) || !isnothing(lateral.river.lake)
        inflow_wb =
            lateral.subsurface.ssf[inds] ./ lateral.subsurface.Δt .+ lateral.land.q_av[inds]
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
    update(lateral.land, network.land, frac_toriver = network.frac_toriver)

    # run river flow
    set_river_inwater(model, ssf_toriver)
    update(
        lateral.river,
        network.river,
        inflow_wb = get_inflow_waterbody(model),
        doy = dayofyear(clock.time),
    )
end
