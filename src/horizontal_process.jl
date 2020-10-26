"Convert a gridded drainage direction to a directed graph"
function flowgraph(ldd::AbstractVector, inds::AbstractVector, pcrdir::AbstractVector)
    # prepare a directed graph to be filled
    n = length(inds)
    graph = DiGraph(n)

    # loop over ldd, adding the edge to the downstream node
    for (from_node, from_index) in enumerate(inds)
        ldd_val = ldd[from_node]
        # skip pits to prevent cycles
        ldd_val == 5 && continue
        to_index = from_index + pcrdir[ldd_val]
        # find the node id of the downstream cell
        to_node = searchsortedfirst(inds, to_index)
        add_edge!(graph, from_node, to_node)
    end
    @assert is_directed(graph)
    @assert !is_cyclic(graph)
    return graph
end

"Kinematic wave flow rate for a single cell and timestep"
function kinematic_wave(Qin, Qold, q, α, β, Δt, Δx)
    ϵ = 1.0e-12
    max_iters = 3000

    if Qin + Qold + q ≈ 0.0
        return 0.0
    else
        # common terms
        ab_pQ = α * β * pow(((Qold + Qin) / 2.0), (β - 1.0))
        Δtx = Δt / Δx
        C = Δtx * Qin + α * pow(Qold, β) + Δt * q

        Qkx = (Δtx * Qin + Qold * ab_pQ + Δt * q) / (Δtx + ab_pQ)
        if isnan(Qkx)
            Qkx = 0.0
        end
        Qkx = max(Qkx, 1.0e-30)
        count = 1

        while true
            fQkx = Δtx * Qkx + α * pow(Qkx, β) - C
            dfQkx = Δtx + α * β * pow(Qkx, (β - 1.0))
            Qkx = Qkx - fQkx / dfQkx
            Qkx = max(Qkx, 1.0e-30)
            if (abs(fQkx) <= ϵ) || (count >= max_iters)
                break
            end
            count += 1
        end

        return Qkx
    end
end

"Kinematic wave flow rate over the whole network for a single timestep"
function kin_wave!(Q, graph, toposort, Qold, q, α, β, DCL, Δt)
    for v in toposort
        upstream_nodes = inneighbors(graph, v)
        Qin = sum_at(Q, upstream_nodes)
        Q[v] = kinematic_wave(Qin, Qold[v], q, α[v], β, Δt, DCL[v])
    end
    return Q
end

"Kinematic wave for lateral subsurface flow for a single cell and timestep"
function kinematic_wave_ssf(ssfin, ssfₜ₋₁, ziₜ₋₁, r, kh₀, β, θₑ, f, d, Δt, Δx, dw, ssfmax)

    ϵ = 1.0e-3
    max_iters = 3000

    if ssfin + ssfₜ₋₁ ≈ 0.0 && r <= 0.0
        return 0.0, d, 0.0
    else
        # initial estimate
        ssf = (ssfₜ₋₁ + ssfin) / 2.0
        count = 1

        # Estimate zi on the basis of the relation between subsurfacel flow and zi
        zi = log((f * ssf) / (dw * kh₀ * β) + exp(-f * d)) / -f
        # Reciprocal of derivative delta Q/ delta z_i, constrained w.r.t. neff on the basis of the continuity equation)
        Cn = (kh₀ * exp(-f * zi) * β) / θₑ
        # Term of the continuity equation for Newton-Raphson iteration for iteration 1
        # because celerity Cn is depending on zi, the increase or decrease of zi is moved to the recharge term of the continuity equation
        # then (1./Cn)*ssfₜ₋₁ can be replaced with (1./Cn)*ssf, and thus celerity and lateral flow rate ssf are then in line
        c = (Δt / Δx) * ssfin + (1.0 / Cn) * ssf + Δt * (r - (ziₜ₋₁ - zi) * θₑ * dw)

        # Continuity equation of which solution should be zero
        fQ = (Δt / Δx) * ssf + (1.0 / Cn) * ssf - c
        # Derivative of the continuity equation w.r.t. Q_out for iteration 1
        dfQ = (Δt / Δx) + 1.0 / Cn
        # Update lateral outflow estimate ssf (Q_out) for iteration 1
        ssf = ssf - (fQ / dfQ)
        if isnan(ssf)
            ssf = 0.0
        end
        ssf = max(ssf, 1.0e-30)

        # Start while loop of Newton-Raphson iteration m until continuity equation approaches zero
        while true
            # Estimate zi on the basis of the relation between lateral flow rate and groundwater level
            zi = log((f * ssf) / (dw * kh₀ * β) + exp(-f * d)) / -f
            # Reciprocal of derivative delta Q/ delta z_i, constrained w.r.t. neff on the basis of the continuity equation
            Cn = (kh₀ * exp(-f * zi) * β) / θₑ

            # Term of the continuity equation for given Newton-Raphson iteration m
            # because celerity Cn is depending on zi, the increase or decrease of zi is moved to the recharge term of the continuity equation
            # then (1./Cn)*ssfₜ₋₁ can be replaced with (1./Cn)*ssf, and thus celerity and lateral flow rate ssf are then in line
            c = (Δt / Δx) * ssfin + (1.0 / Cn) * ssf + Δt * (r - (ziₜ₋₁ - zi) * θₑ * dw)

            # Continuity equation of which solution should be zero
            fQ = (Δt / Δx) * ssf + (1.0 / Cn) * ssf - c
            # Derivative of the continuity equation w.r.t. Q_out for iteration m+1
            dfQ = (Δt / Δx) + 1.0 / Cn
            # Update lateral outflow estimate ssf (Q_out) for iteration m+1
            ssf = ssf - (fQ / dfQ)
            if isnan(ssf)
                ssf = 0.0
            end
            ssf = max(ssf, 1.0e-30)
            if (abs(fQ) <= ϵ) || (count >= max_iters)
                break
            end
            count += 1
        end

        # Constrain the lateral flow rate ssf
        ssf = min(ssf, (ssfmax * dw))
        # On the basis of the lateral flow rate, estimate the amount of groundwater level above surface (saturation excess conditions), then rest = negative
        rest = ziₜ₋₁ - (ssfin + r * Δx - ssf) / (dw * Δx) / θₑ
        # In case the groundwater level lies above surface (saturation excess conditions, rest = negative), calculate the exfiltration rate and set groundwater back to zero.
        exfilt = min(rest, 0.0) * -θₑ
        zi = clamp(zi, 0.0, d)

        return ssf, zi, exfilt

    end
end

"Transport of material downstream with a limited transport capacity over a directed graph"
function accucapacityflux(network, material, capacity)
    @unpack graph, order = network
    for v in order
        upstream_nodes = inneighbors(graph, v)
        if !isempty(upstream_nodes)
            flux = sum(min(material[i], capacity[i]) for i in upstream_nodes)
            material[v] += flux
        end
    end
    return material
end

function accucapacitystate(network, material, capacity)
    @unpack graph, order = network
    for v in order
        upstream_nodes = inneighbors(graph, v)
        if !isempty(upstream_nodes)
            state = sum(max(material[i] - capacity[i], 0.0) for i in upstream_nodes)
            material[v] += state
        end
    end
    return material
end