"Convert a gridded drainage direction to a directed graph"
function flowgraph(ldd::AbstractVector, inds::AbstractVector, pcr_dir::AbstractVector)
    # prepare a directed graph to be filled
    n = length(inds)
    graph = DiGraph(n)

    # loop over ldd, adding the edge to the downstream node
    for (from_node, from_index) in enumerate(inds)
        ldd_val = ldd[from_node]
        # skip pits to prevent cycles
        ldd_val == 5 && continue
        to_index = from_index + pcr_dir[ldd_val]
        # find the node id of the downstream cell
        to_node = searchsortedfirst(inds, to_index)
        add_edge!(graph, from_node, to_node)
    end
    if is_cyclic(graph)
        error("""One or more cycles detected in flow graph.
            The provided local drainage direction map may be unsound.
            Verify that each active flow cell flows towards a pit.
            """)
    end
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

    ϵ = 1.0e-12
    max_iters = 3000

    if ssfin + ssfₜ₋₁ ≈ 0.0 && r <= 0.0
        return 0.0, d, 0.0
    else
        # initial estimate
        ssf = (ssfₜ₋₁ + ssfin) / 2.0
        count = 1

        # celerity (Cn) based on average hydraulic conductivity (ksat)
        ksat = (kh₀ / f * (exp(-f * ziₜ₋₁) - exp(-f * d))) / (d - ziₜ₋₁)
        Cn = (β * ksat) / θₑ

        # constant term of the continuity equation for Newton-Raphson
        c = (Δt / Δx) * ssfin + (1.0 / Cn) * ssfₜ₋₁ + r

        # continuity equation of which solution should be zero
        fQ = (Δt / Δx) * ssf + (1.0 / Cn) * ssf - c
        # derivative of the continuity equation
        dfQ = (Δt / Δx) + 1.0 / Cn
        # update lateral outflow estimate ssf
        ssf = ssf - (fQ / dfQ)
        if isnan(ssf)
            ssf = 0.0
        end
        ssf = max(ssf, 1.0e-30)

        # start Newton-Raphson iterations
        while true
            fQ = (Δt / Δx) * ssf + (1.0 / Cn) * ssf - c
            dfQ = (Δt / Δx) + 1.0 / Cn
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

        # constrain the lateral flow rate ssf
        ssf = min(ssf, (ssfmax * dw))
        # estimate zi (water balance)
        zi = ziₜ₋₁ - (ssfin * Δt + r * Δx - ssf * Δt) / (dw * Δx) / θₑ
        # if zi is above surface (saturation excess conditions, zi = negative), compute the
        # exfiltration rate and set zi to zero.
        exfilt = min(zi, 0.0) * -θₑ
        zi = clamp(zi, 0.0, d)

        return ssf, zi, exfilt

    end
end

"""
    accucapacitystate!(material, network, capacity) -> material

Transport of material downstream with a limited transport capacity over a directed graph.
Mutates the material input. The network is expected to hold a graph and order field, where
the graph implements the Graphs interface, and the order is a valid topological ordering
such as that returned by `Graphs.topological_sort_by_dfs`.

Returns the material state after transport.
"""
function accucapacitystate!(material, network, capacity)
    @unpack graph, order = network
    for v in order
        downstream_nodes = outneighbors(graph, v)
        n = length(downstream_nodes)
        flux_val = min(material[v], capacity[v])
        material[v] -= flux_val
        if n == 0
            # pit: material is transported out of the map if a capacity is set,
            # cannot add the material anywhere
        elseif n == 1
            material[only(downstream_nodes)] += flux_val
        else
            error("bifurcations not supported")
        end
    end
    return material
end

"""
    accucapacitystate!(material, network, capacity) -> material

Non mutating version of `accucapacitystate!`.
"""
function accucapacitystate(material, network, capacity)
    accucapacitystate!(copy(material), network, capacity)
end

"""
    accucapacityflux!(flux, material, network, capacity) -> flux, material

Transport of material downstream with a limited transport capacity over a directed graph.
Updates the material input, and overwrites the flux input, not using existing values. The
network is expected to hold a graph and order field, where the graph implements the Graphs
interface, and the order is a valid topological ordering such as that returned by
`Graphs.topological_sort_by_dfs`.

Returns the flux (material leaving each cell), and material (left after transport).
"""
function accucapacityflux!(flux, material, network, capacity)
    @unpack graph, order = network
    for v in order
        downstream_nodes = outneighbors(graph, v)
        n = length(downstream_nodes)
        flux_val = min(material[v], capacity[v])
        material[v] -= flux_val
        flux[v] = flux_val
        if n == 0
            # pit: material is transported out of the map if a capacity is set,
            # cannot add the material anywhere
        elseif n == 1
            material[only(downstream_nodes)] += flux_val
        else
            error("bifurcations not supported")
        end
    end
    return flux, material
end

"""
    accucapacityflux!(material, network, capacity) -> flux, material

Non mutating version of `accucapacityflux!`.
"""
function accucapacityflux(material, network, capacity)
    accucapacityflux!(zero(material), copy(material), network, capacity)
end

"""
    lateral_snow_transport!(snow, snowwater, slope, network)

Lateral snow transport. Transports snow downhill. Mutates `snow` and `snowwater`.
"""
function lateral_snow_transport!(snow, snowwater, slope, network)
    snowflux_frac = min.(0.5, slope ./ 5.67) .* min.(1.0, snow ./ 10000.0)
    maxflux = snowflux_frac .* snow
    snow = accucapacitystate!(snow, network, maxflux)
    snowwater = accucapacitystate!(snowwater, network, snowwater .* snowflux_frac)
    return snow, snowwater
end

"""
    local_inertial_flow(q0, η0, η1, hf, A, R, length, mannings_n, g, froude_limit, Δt)

Local inertial approach for flow through area `A`. Returns the flow `q` between two adjacent
river cells (nodes) for a single timestep.
"""
function local_inertial_flow(
    q0,
    η0,
    η1,
    hf,
    A,
    R,
    length,
    mannings_n_sq,
    g,
    froude_limit,
    Δt,
)

    slope = (η1 - η0) / length
    pow_R = cbrt(R * R * R * R)
    unit = one(hf)
    q = (
        (q0 - g * A * Δt * slope) / (unit + g * Δt * mannings_n_sq * abs(q0) / (pow_R * A))
    )

    # if froude number > 1.0, limit flow
    fr = ((q / A) / sqrt(g * hf)) * froude_limit
    q = IfElse.ifelse((abs(fr) > 1.0) * (q > 0.0), sqrt(g * hf) * A, q)
    q = IfElse.ifelse((abs(fr) > 1.0) * (q < 0.0), -sqrt(g * hf) * A, q)

    return q
end

"""
    local_inertial_flow(θ, q0, qd, qu, η0, η1, hf, width, length, mannings_n, g, froude_limit, Δt)

Local inertial approach for flow through a rectangular area. Returns the flow `q` between
two adjacent cells (nodes) for a single timestep. Algorithm is based on de Almeida et al.
(2012).
"""
function local_inertial_flow(
    θ,
    q0,
    qd,
    qu,
    η0,
    η1,
    hf,
    width,
    length,
    mannings_n_sq,
    g,
    froude_limit,
    Δt,
)

    slope = (η1 - η0) / length
    unit = one(θ)
    half = oftype(θ, 0.5)
    pow_hf = cbrt(hf * hf * hf * hf * hf * hf * hf)

    q = (
        ((θ * q0 + half * (unit - θ) * (qu + qd)) - g * hf * width * Δt * slope) /
        (unit + g * Δt * mannings_n_sq * abs(q0) / (pow_hf * width))
    )
    # if froude number > 1.0, limit flow
    if froude_limit
        fr = (q / width / hf) / sqrt(g * hf)
        if abs(fr) > 1.0 && q > 0.0
            q = hf * sqrt(g * hf) * width
        elseif abs(fr) > 1.0 && q < 0.0
            q = -hf * sqrt(g * hf) * width
        end
    end

    return q
end
