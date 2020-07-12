Base.@kwdef struct SurfaceFlow{T,R,L}
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
    alpha_term::Vector{T} = pow.(n ./ sqrt.(sl), β)  # Constant part of α
    alpha_pow::T = (2.0 / 3.0) * β          # Used in the power part of α
    α::Vector{T} = alpha_term .* pow.(width .+ 2.0 .* h, alpha_pow) # Constant in momentum equation A = αQᵝ, based on Manning's equation
    eps::T = 1e-03                          # Maximum allowed change in α, if exceeded cross sectional area and h is recalculated
    cel::Vector{T} = fill(0.0, length(sl))  # Celerity of the kinematic wave
    to_river::Vector{T} = fill(0.0, length(sl)) # Part of overland flow [m³ s⁻¹] that flows to the river
    rivercells::Vector{Bool} = fill(false, length(sl)) # Location of river cells (0 or 1)
    wb_pit::Vector{Bool} = fill(false, length(sl)) # Boolean location (0 or 1) of a waterbody (wb, reservoir or lake).
    reservoir::R = Fill(missing, length(sl)) # Reservoir model, only located in river cells
    lake::L = Fill(missing, length(sl))                 # Lake model, only located in river cells
end

statenames(::SurfaceFlow) = (:q, :h)

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
    doy = 0,
)
    # two options for iteration, a fixed sub time step or based on courant number.
    if do_iter
        if do_tstep
            # better to set this during initialization
            ts = ceil(Int(sf.Δt / tstep))
        else
            # calculate celerity
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

    # sub time step
    adt = sf.Δt / ts
    p = 3.0 / ts # dummy precipitation
    pet = 4.0 / ts # dummy evaporation

    q_sum = zeros(n)
    h_sum = zeros(n)
    for _ = 1:ts
        for v in toposort
            upstream_nodes = inneighbors(dag, v)
            # for overland flow frac_toriver and river cells need to be defined
            if (frac_toriver != nothing) && (river != nothing)
                # for a river cell without a reservoir or lake (wb_pit = 0) part of the upstream surface flow
                # goes to the river (frac_toriver) and part goes to the surface flow reservoir (1.0 - frac_toriver)
                # upstream nodes with a reservoir or lake are excluded
                if river[v] && !sf.wb_pit[v]
                    qin = isempty(upstream_nodes) ? 0.0 :
                        sum_empty([
                        sf.q[i] * (1.0 - frac_toriver[i])
                        for i in upstream_nodes if sf.wb_pit[i] == 0
                    ])
                    sf.to_river[v] = isempty(upstream_nodes) ? 0.0 :
                        sum_empty([
                        sf.q[i] * frac_toriver[i]
                        for i in upstream_nodes if sf.wb_pit[i] == 0
                    ])
                    # for a river cell with a reservoir or lake (wb_pit = 1) all upstream surface flow goes
                    # to the river.
                elseif river[v] && sf.wb_pit[v]
                    sf.to_river[v] = sum_at(sf.q, upstream_nodes)
                    qin = 0.0
                else
                    qin = sum_at(sf.q, upstream_nodes)
                end
                # for all the other cells all upstream surface flow goes to the surface flow reservoir.
            else
                qin = sum_at(sf.q, upstream_nodes)
            end
            # run reservoir model and copy reservoir outflow to river cell
            # dummy values now for reservoir precipitation and evaporation (3.0 and 4.0)
            if !ismissing(sf.reservoir[v])
                sf.reservoir[v] = update(sf.reservoir[v], qin, p, pet, adt)
                sf.q[v] = sf.reservoir[v].outflow
                # run lake model and copy lake outflow to river cell
                # dummy values now for lake precipitation and evaporation (3.0 and 4.0)
            elseif !ismissing(sf.lake[v])
                if sf.lake[v].lowerlake_ind != 0
                    lower_lake = sf.lake[sf.lake[v].lowerlake_ind]
                    sf.lake[v] = update(sf.lake[v], qin, p, pet, doy, lower_lake)
                    sf.q[v] = sf.lake[v].outflow
                else
                    sf.lake[v] = update(sf.lake[v], qin, p, pet, doy)
                    sf.q[v] = sf.lake[v].outflow
                end

            else
                sf.q[v] =
                    kinematic_wave(qin, sf.q[v], sf.qlat[v], sf.α[v], sf.β, adt, sf.dl[v])
            end

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

        sf.q_av .= q_sum ./ ts
        sf.h_av .= h_sum ./ ts

    end

end
