Base.@kwdef struct SurfaceFlow{T}
    β::T = 0.6                              # constant in Manning's equation
    sl::Vector{T}                           # Slope [m m⁻¹]
    n::Vector{T}                            # Manning's roughness [sl m⁻⅓]
    dl::Vector{T}                           # Drain length [m]
    q::Vector{T} = fill(0.0, length(sl))    # Discharge [m³ s⁻¹]
    q_av::Vector{T} = fill(0.0, length(sl)) # Average discharge [m³ s⁻¹]
    qlat::Vector{T} = fill(0.0, length(sl)) # Lateral discharge [m³ s⁻¹]
    h::Vector{T} = fill(0.0, length(sl))    # Water level [m]
    h_av::Vector{T} = fill(0.0, length(sl)) # Average water level [m]
    Δt::T                                   # Model time step [s]
    width::Vector{T}                        # Flow width [m]
    alpha_term::Vector{T} = pow.(n ./ sqrt.(sl), β)
    alpha_pow::T = (2.0 / 3.0) * β
    α::Vector{T} = alpha_term .* pow.(width .+ 2.0 .* h, alpha_pow) # constant in Manning's equation
    eps::T = 1e-03
    cel::Vector{T} = fill(0.0, length(sl))
    to_river::Vector{T} = fill(0.0, length(sl))
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
    frac_toriver = nothing,
    river = nothing,
    do_iter = false,
    do_tstep = false,
    tstep = 0.0,
)

    if do_iter
        if do_tstep
            ts = ceil(Int(sf.Δt / tstep))
        else
            for v in toposort
                if sf.q[v] > 0.0
                    sf.cel[v] = 1.0 / (sf.α[v] * sf.β * pow(sf.q[v], (sf.β - 1.0)))
                else
                    sf.cel[v] = 0.0
                end
            end
            courant = (sf.Δt ./ sf.dl) .* sf.cel
            ts = max(ceil(Int, (1.25 * quantile!(courant, 0.95))), 1)
        end
    end

    adt = sf.Δt / ts

    q_sum = zeros(n)
    h_sum = zeros(n)
    for _ = 1:ts
        for v in toposort
            upstream_nodes = inneighbors(dag, v)
            if (frac_toriver != nothing) & (river != nothing)
                if Bool(river[v])
                    qin = isempty(upstream_nodes) ? 0.0 :
                        sum(sf.q[i] * (1.0 - frac_toriver[i]) for i in upstream_nodes)
                    sf.to_river[v] = isempty(upstream_nodes) ? 0.0 :
                        sum(sf.q[i] * frac_toriver[i] for i in upstream_nodes)
                else
                    qin =
                        isempty(upstream_nodes) ? 0.0 : sum(sf.q[i] for i in upstream_nodes)
                end
            else
                qin = isempty(upstream_nodes) ? 0.0 : sum(sf.q[i] for i in upstream_nodes)
            end
            sf.q[v] = kinematic_wave(qin, sf.q[v], sf.qlat[v], sf.α[v], sf.β, adt, sf.dl[v])

            # update alpha
            crossarea = sf.α[v] * pow(sf.q[v], sf.β)
            sf.h[v] = crossarea / sf.width[v]
            wetper = sf.width[v] + (2.0 * sf.h[v]) # wetted perimeter
            α = sf.α[v]
            sf.α[v] = sf.alpha_term[v] * pow(wetper, sf.alpha_pow)

            if abs(α - sf.α[v]) > sf.eps
                crossarea = sf.α[v] * pow(sf.q[v], sf.β)
                sf.h[v] = crossarea / sf.width[v]
            end

            q_sum[v] += sf.h[v]
            h_sum[v] += sf.q[v]

        end

        sf.q_av[:] = q_sum ./ ts
        sf.h_av[:] = h_sum ./ ts

    end

end
