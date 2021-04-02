@get_units @with_kw struct SurfaceFlow{T,R,L}
    β::T | "-"                              # constant in Manning's equation
    sl::Vector{T} | "m m-1"                 # Slope [m m⁻¹]
    n::Vector{T} | "s m-1/3"                # Manning's roughness [s m⁻⅓]
    dl::Vector{T} | "m"                     # Drain length [m]
    q::Vector{T} | "m3 s-1"                 # Discharge [m³ s⁻¹]
    qin::Vector{T} | "m3 s-1"               # Inflow from upstream cells [m³ s⁻¹]
    q_av::Vector{T} | "m3 s-1"              # Average discharge [m³ s⁻¹]
    qlat::Vector{T} | "m3 s-1"              # Lateral discharge [m³ s⁻¹]
    h::Vector{T} | "m"                      # Water level [m]
    h_av::Vector{T} | "m"                   # Average water level [m]
    Δt::T | "s"                             # Model time step [s]
    its::Int | "-"                          # Number of fixed iterations
    width::Vector{T} | "m"                  # Flow width [m]
    alpha_pow::T | "-"                      # Used in the power part of α
    α::Vector{T} | "s3/5 m1/5"              # Constant in momentum equation A = αQᵝ, based on Manning's equation
    eps::T | "s3/5 m1/5"                    # Maximum allowed change in α, if exceeded cross sectional area and h is recalculated
    cel::Vector{T} | "m s-1"                # Celerity of the kinematic wave
    to_river::Vector{T} | "m3 s-1"          # Part of overland flow [m³ s⁻¹] that flows to the river
    rivercells::Vector{Bool} | "-"          # Location of river cells (0 or 1)
    wb_pit::Vector{Bool} | "-"              # Boolean location (0 or 1) of a waterbody (wb, reservoir or lake).
    reservoir_index::Vector{Int} | "-"      # map cell to 0 (no reservoir) or i (pick reservoir i in reservoir field)
    lake_index::Vector{Int} | "-"           # map cell to 0 (no lake) or i (pick lake i in lake field)
    reservoir::R                            # Reservoir model struct of arrays
    lake::L                                 # Lake model struct of arrays

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
    do_iter = false,
    doy = 0,
)
    @unpack graph, order, upstream_nodes = network


    n = length(order)
    # two options for iteration, fixed or based on courant number.
    if do_iter
        if sf.its > 0
            its = sf.its
        else
            # calculate celerity
            courant = zeros(n)
            for v in order
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
    adt = sf.Δt / its

    alpha_term = pow.(sf.n ./ sqrt.(sf.sl),sf.β)
    sf.α .= alpha_term .* pow.(sf.width .+ 2.0 .* sf.h, sf.alpha_pow)

    q_sum = zeros(n)
    h_sum = zeros(n)
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
        for (n, v) in enumerate(order)
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


                if !isnothing(sf.reservoir) && sf.reservoir_index[v] != 0
                    # run reservoir model and copy reservoir outflow to inflow (qin) of
                    # downstream river cell
                    i = sf.reservoir_index[v]
                    update(sf.reservoir, i, sf.q[v] + inflow_wb[v], adt)
                    j = only(outneighbors(graph, v))
                    sf.qin[j] = sf.reservoir.outflow[i]

                elseif !isnothing(sf.lake) && sf.lake_index[v] != 0
                    # run lake model and copy lake outflow to inflow (qin) of downstream river
                    # cell
                    i = sf.lake_index[v]
                    update(sf.lake, i, sf.q[v] + inflow_wb[v], doy, adt)
                    j = only(outneighbors(graph, v))
                    sf.qin[j] = sf.lake.outflow[i]
                end
            end

            sf.q[v] =
                kinematic_wave(sf.qin[v], sf.q[v], sf.qlat[v], sf.α[v], sf.β, adt, sf.dl[v])

            # update alpha
            crossarea = sf.α[v] * pow(sf.q[v], sf.β)
            sf.h[v] = crossarea / sf.width[v]
            wetper = sf.width[v] + (2.0 * sf.h[v]) # wetted perimeter
            α = sf.α[v]
            sf.α[v] = alpha_term[v] * pow(wetper, sf.alpha_pow)

            if abs(α - sf.α[v]) > sf.eps
                crossarea = sf.α[v] * pow(sf.q[v], sf.β)
                sf.h[v] = crossarea / sf.width[v]
            end

            q_sum[v] += sf.q[v]
            h_sum[v] += sf.h[v]
        end

    end
    sf.q_av .= q_sum ./ its
    sf.h_av .= h_sum ./ its
    sf.to_river .= sf.to_river ./ its
end

@get_units @with_kw struct LateralSSF{T}
    kh₀::Vector{T} | "mm Δt-1"              # Horizontal hydraulic conductivity at soil surface [mm Δt⁻¹]
    f::Vector{T} | "mm-1"                   # A scaling parameter [mm⁻¹] (controls exponential decline of kh₀)
    soilthickness::Vector{T} | "mm"         # Soil thickness [mm]
    θₛ::Vector{T} | "mm mm-1"               # Saturated water content (porosity) [-]
    θᵣ::Vector{T} | "mm mm-1"               # Residual water content [-]
    t::T | "Δt s"                           # time step [Δt s]
    Δt::T | "s"                             # model time step [s]
    βₗ::Vector{T} | "m m-1"                  # Slope [m m⁻¹]
    dl::Vector{T}                           # Drain length [mm]
    dw::Vector{T}                           # Flow width [mm]
    zi::Vector{T}                           # Pseudo-water table depth [mm] (top of the saturated zone)
    exfiltwater::Vector{T}                  # Exfiltration [mm]  (groundwater above surface level, saturated excess conditions)
    recharge::Vector{T}                     # Net recharge to saturated store [mm]
    ssf::Vector{T} | "mm3 Δt-1"             # Subsurface flow [mm³ Δt⁻¹]
    ssfin::Vector{T} | "mm3 Δt-1"           # Inflow from upstream cells [mm³ Δt⁻¹]
    ssfmax::Vector{T} | "mm2 Δt-1"          # Maximum subsurface flow [mm² Δt⁻¹]
    to_river::Vector{T} | "mm3 Δt-1"        # Part of subsurface flow [mm³ Δt⁻¹] that flows to the river
    wb_pit::Vector{Bool} | "-"              # Boolean location (0 or 1) of a waterbody (wb, reservoir or lake).

    function LateralSSF{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::LateralSSF) = (:ssf,)

function update(ssf::LateralSSF, network, frac_toriver)
    @unpack order, upstream_nodes = network
    for (n, v) in enumerate(order)

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

@get_units @with_kw struct GroundwaterExchange{T}
    Δt::T | "s"                         # model time step [s]
    exfiltwater::Vector{T}              # Exfiltration [mm]  (groundwater above surface level, saturated excess conditions)
    zi::Vector{T}                       # Pseudo-water table depth [mm] (top of the saturated zone)
    to_river::Vector{T} | "m3 s-1"      # Part of subsurface flow [mm³ Δt⁻¹] that flows to the river
    ssf::Vector{T} | "mm3 Δt-1"         # Subsurface flow [mm³ Δt⁻¹]
end