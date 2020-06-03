Base.@kwdef struct SurfaceFlow{T}
    β::T = 0.6                              # constant in Manning's equation
    sl::Vector{T}                           # Slope [m m⁻¹]
    n::Vector{T}                            # Manning's roughness [sl m⁻⅓]
    dcl::Vector{T}                          # Drain length [m]
    q::Vector{T} = fill(0.0, length(sl))    # Discharge [m³ s⁻¹]
    q_av::Vector{T} = fill(0.0, length(sl)) # Average discharge [m³ s⁻¹]
    qlat::Vector{T} = fill(0.0, length(sl)) # Lateral discharge [m³ s⁻¹]
    h::Vector{T} = fill(0.0, length(sl))    # Water level [m]
    h_av::Vector{T} = fill(0.0, length(sl)) # Average water level [m]
    Δt::T                                   # Model time step [s]
    riverwidth::Vector{T}                   # River width [m]
    alpha_term::Vector{T} = pow.(n ./ sqrt.(sl), β)
    alpha_pow::T = (2.0 / 3.0) * β
    α::Vector{T} = alpha_term .* pow.(riverwidth .+ 2.0 .* h, alpha_pow) # constant in Manning's equation
    eps::T = 1e-03
    cel::Vector{T} = fill(0.0, length(sl))
end


"""
    statenames(::Type{SurfaceFlow})

Returns Array{Symbol,1} for extracting model state fields.
"""
function statenames(::Type{SurfaceFlow})

    states = [:q, :h]
    # TODO: (warm) states read from netcdf file or cold state (reinit=1, setting in ini file)

end

function update(
    sf::SurfaceFlow,
    dag,
    toposort,
    n;
    do_iter = false,
    do_tstep = false,
    ts = 1,
)

    if do_iter
        if do_tstep
            ts = ceil(Int(sf.Δt / tstep))
        else
            sf.cel[sf.q>0.0] =
                1.0 ./
                (sf.α[sf.q>0.0] .* sf.β * pow.(sf.q[sf.q>0.0], (sf.β - 1.0)))
            courant = (sf.Δt ./ sf.dcl) .* sf.cel
            ts = ceil(Int, (1.25 * quantile!(courant, 0.95)))
        end
    end

    adt = sf.Δt / ts

    q_sum = zeros(n)
    h_sum = zeros(n)
    for _ = 1:ts
        for v in toposort
            upstream_nodes = inneighbors(dag, v)
            qin = isempty(upstream_nodes) ? 0.0 :
                sum(sf.q[i] for i in upstream_nodes)
            sf.q[v] = kinematic_wave(
                qin,
                sf.q[v],
                sf.qlat[v],
                sf.α[v],
                sf.β,
                adt,
                sf.dcl[v],
            )

            # update alpha
            crossarea = sf.α[v] * pow(sf.q[v], sf.β)
            sf.h[v] = crossarea / sf.riverwidth[v]
            wetper = sf.riverwidth[v] + (2.0 * sf.h[v]) # wetted perimeter
            α = sf.α[v]
            sf.α[v] = sf.alpha_term[v] * pow(wetper, sf.alpha_pow)

            if abs(α - sf.α[v]) > sf.eps
                crossarea = sf.α[v] * pow(sf.q[v], sf.β)
                sf.h[v] = crossarea / sf.riverwidth[v]
            end

            q_sum[v] += sf.h[v]
            h_sum[v] += sf.q[v]

        end

        sf.q_av = q_sum ./ ts
        sf.h_av = h_sum ./ ts

    end

end
