"Convert a gridded drainage direction to a directed graph"
function flowgraph(ldd::AbstractVector, indices::AbstractVector, PCR_DIR::AbstractVector)
    # prepare a directed graph to be filled
    n = length(indices)
    graph = DiGraph{Int}(n)

    # loop over ldd, adding the edge to the downstream node
    for (from_node, from_index) in enumerate(indices)
        ldd_val = ldd[from_node]
        # skip pits to prevent cycles
        ldd_val == 5 && continue
        to_index = from_index + PCR_DIR[ldd_val]
        # find the node id of the downstream cell
        to_node = searchsortedfirst(indices, to_index)
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
function kinematic_wave(Qin, Qold, q, alpha, beta, dt, dx)
    epsilon = 1.0e-12
    max_iters = 3000

    if Qin + Qold + q ≈ 0.0
        return 0.0
    else
        # common terms
        ab_pQ = alpha * beta * pow(((Qold + Qin) / 2.0), (beta - 1.0))
        dtx = dt / dx
        C = dtx * Qin + alpha * pow(Qold, beta) + dt * q

        Qkx = (dtx * Qin + Qold * ab_pQ + dt * q) / (dtx + ab_pQ)
        if isnan(Qkx)
            Qkx = 0.0
        end
        Qkx = max(Qkx, 1.0e-30)
        count = 1

        while true
            fQkx = dtx * Qkx + alpha * pow(Qkx, beta) - C
            dfQkx = dtx + alpha * beta * pow(Qkx, (beta - 1.0))
            Qkx = Qkx - fQkx / dfQkx
            Qkx = max(Qkx, 1.0e-30)
            if (abs(fQkx) <= epsilon) || (count >= max_iters)
                break
            end
            count += 1
        end

        return Qkx
    end
end

"Kinematic wave flow rate over the whole network for a single timestep"
function kin_wave!(Q, graph, toposort, Qold, q, alpha, beta, DCL, dt)
    for v in toposort
        upstream_nodes = inneighbors(graph, v)
        Qin = sum_at(Q, upstream_nodes)
        Q[v] = kinematic_wave(Qin, Qold[v], q, alpha[v], beta, dt, DCL[v])
    end
    return Q
end

"Returns water table depth `zi` based on lateral subsurface flow `ssf` and hydraulic conductivity profile `KhExponential`"
function ssf_water_table_depth(ssf, slope, d, dw, kh_profile::KhExponential, i)
    (; f, kh_0) = kh_profile
    zi = log((f[i] * ssf) / (dw * kh_0[i] * slope) + exp(-f[i] * d)) / -f[i]
    return zi
end

"Returns water table depth `zi` based on lateral subsurface flow `ssf` and hydraulic conductivity profile `KhExponentialConstant`"
function ssf_water_table_depth(ssf, slope, d, dw, kh_profile::KhExponentialConstant, i)
    (; z_exp) = kh_profile
    (; kh_0, f) = kh_profile.exponential
    ssf_constant = kh_0[i] * slope * exp(-f[i] * z_exp[i]) * (d - z_exp[i]) * dw
    if ssf > ssf_constant
        zi =
            log(
                (f[i] * (ssf - ssf_constant)) / (dw * kh_0[i] * slope) +
                exp(-f[i] * z_exp[i]),
            ) / -f[i]
    else
        zi = d - ssf / (dw * kh_0[i] * slope * exp(-f[i] * z_exp[i]))
    end
    return zi
end

"Returns kinematic wave celecity `Cn` of lateral subsurface flow based on hydraulic conductivity profile `KhExponential`"
function ssf_celerity(zi, slope, theta_e, kh_profile::KhExponential, i)
    (; kh_0, f) = kh_profile
    Cn = (kh_0[i] * exp(-f[i] * zi) * slope) / theta_e
    return Cn
end

"Returns kinematic wave celecity `Cn` of lateral subsurface flow based on hydraulic conductivity profile `KhExponentialConstant`"
function ssf_celerity(zi, slope, theta_e, kh_profile::KhExponentialConstant, i)
    (; z_exp) = kh_profile
    (; kh_0, f) = kh_profile.exponential
    Cn_const = (kh_0[i] * exp(-f[i] * z_exp[i]) * slope) / theta_e
    return if zi < z_exp[i]
        (kh_0[i] * exp(-f[i] * zi) * slope) / theta_e + Cn_const
    else
        Cn_const
    end
end

"""
    kinematic_wave_ssf(ssfin, ssf_prev, zi_prev, r, kh_0, slope, theta_e, f, d, dt, dx, dw, ssfmax, z_exp, ksat_profile)

Kinematic wave for lateral subsurface flow for a single cell and timestep. An exponential
decline of hydraulic conductivity at the soil surface `kh_0`, controllled by parameter `f`,
is assumed. The hydraulic conductivity profile `ksat_profile` is either `exponential` or
`exponential_constant`, with `z_exp` the depth from the soil surface for which the
exponential decline of `kh_0` is valid.

Returns lateral subsurface flow `ssf`, water table depth `zi` and exfiltration rate
`exfilt`.
"""
function kinematic_wave_ssf(
    ssfin,
    ssf_prev,
    zi_prev,
    r,
    slope,
    theta_e,
    d,
    dt,
    dx,
    dw,
    ssfmax,
    kh_profile::Union{KhExponential, KhExponentialConstant},
    i,
)
    epsilon = 1.0e-12
    max_iters = 3000

    if ssfin + ssf_prev ≈ 0.0 && r <= 0.0
        return 0.0, d, 0.0
    else
        # initial estimate
        ssf = (ssf_prev + ssfin) / 2.0
        count = 1

        # Estimate zi on the basis of the relation between subsurface flow and zi
        zi = ssf_water_table_depth(ssf, slope, d, dw, kh_profile, i)
        # Reciprocal of derivative delta Q/ delta z_i, constrained w.r.t. neff on the basis of the continuity equation)
        Cn = ssf_celerity(zi, slope, theta_e, kh_profile, i)
        # Term of the continuity equation for Newton-Raphson iteration for iteration 1
        # because celerity Cn is depending on zi, the increase or decrease of zi is moved to the recharge term of the continuity equation
        # then (1./Cn)*ssf_prev can be replaced with (1./Cn)*ssf, and thus celerity and lateral flow rate ssf are then in line
        c = (dt / dx) * ssfin + (1.0 / Cn) * ssf + (r - (zi_prev - zi) * theta_e * dw)

        # Continuity equation of which solution should be zero
        fQ = (dt / dx) * ssf + (1.0 / Cn) * ssf - c
        # Derivative of the continuity equation w.r.t. Q_out for iteration 1
        dfQ = (dt / dx) + 1.0 / Cn
        # Update lateral outflow estimate ssf (Q_out) for iteration 1
        ssf = ssf - (fQ / dfQ)
        if isnan(ssf)
            ssf = 0.0
        end
        ssf = max(ssf, 1.0e-30)

        # Start while loop of Newton-Raphson iteration m until continuity equation approaches zero
        while true
            # Estimate zi on the basis of the relation between lateral flow rate and groundwater level
            zi = ssf_water_table_depth(ssf, slope, d, dw, kh_profile, i)
            # Reciprocal of derivative delta Q/ delta z_i, constrained w.r.t. neff on the basis of the continuity equation
            Cn = ssf_celerity(zi, slope, theta_e, kh_profile, i)
            # Term of the continuity equation for given Newton-Raphson iteration m
            # because celerity Cn is depending on zi, the increase or decrease of zi is moved to the recharge term of the continuity equation
            # then (1./Cn)*ssf_prev can be replaced with (1./Cn)*ssf, and thus celerity and lateral flow rate ssf are then in line
            c = (dt / dx) * ssfin + (1.0 / Cn) * ssf + (r - (zi_prev - zi) * theta_e * dw)

            # Continuity equation of which solution should be zero
            fQ = (dt / dx) * ssf + (1.0 / Cn) * ssf - c
            # Derivative of the continuity equation w.r.t. Q_out for iteration m+1
            dfQ = (dt / dx) + 1.0 / Cn
            # Update lateral outflow estimate ssf (Q_out) for iteration m+1
            ssf = ssf - (fQ / dfQ)
            if isnan(ssf)
                ssf = 0.0
            end
            ssf = max(ssf, 1.0e-30)
            if (abs(fQ) <= epsilon) || (count >= max_iters)
                break
            end
            count += 1
        end

        # Constrain the lateral flow rate ssf
        ssf = min(ssf, (ssfmax * dw))
        # On the basis of the lateral flow rate, estimate the amount of groundwater level above surface (saturation excess conditions), then rest = negative
        rest = zi_prev - (ssfin * dt + r * dx - ssf * dt) / (dw * dx) / theta_e
        # In case the groundwater level lies above surface (saturation excess conditions, rest = negative), calculate the exfiltration rate and set groundwater back to zero.
        exfilt = min(rest, 0.0) * -theta_e
        zi = clamp(zi, 0.0, d)

        return ssf, zi, exfilt
    end
end

"""
    kinematic_wave_ssf(ssfin, ssf_prev, zi_prev, r, kh, slope, theta_e, d, dt, dx, dw, ssfmax)

Kinematic wave for lateral subsurface flow for a single cell and timestep, based on
(average) hydraulic conductivity `kh`.

Returns lateral subsurface flow `ssf`, water table depth `zi` and exfiltration rate
`exfilt`.
"""
function kinematic_wave_ssf(
    ssfin,
    ssf_prev,
    zi_prev,
    r,
    slope,
    theta_e,
    d,
    dt,
    dx,
    dw,
    ssfmax,
    kh_profile::KhLayered,
    i,
)
    epsilon = 1.0e-12
    max_iters = 3000

    if ssfin + ssf_prev ≈ 0.0 && r <= 0.0
        return 0.0, d, 0.0
    else
        # initial estimate
        ssf = (ssf_prev + ssfin) / 2.0
        count = 1
        # celerity (Cn)
        Cn = (slope * kh_profile.kh[i]) / theta_e
        # constant term of the continuity equation for Newton-Raphson
        c = (dt / dx) * ssfin + (1.0 / Cn) * ssf_prev + r
        # continuity equation of which solution should be zero
        fQ = (dt / dx) * ssf + (1.0 / Cn) * ssf - c
        # Derivative of the continuity equation
        dfQ = (dt / dx) + 1.0 / Cn
        # Update lateral subsurface flow estimate ssf
        ssf = ssf - (fQ / dfQ)
        if isnan(ssf)
            ssf = 0.0
        end
        ssf = max(ssf, 1.0e-30)

        # Start while loop of Newton-Raphson iteration
        while true
            fQ = (dt / dx) * ssf + (1.0 / Cn) * ssf - c
            dfQ = (dt / dx) + 1.0 / Cn
            ssf = ssf - (fQ / dfQ)
            if isnan(ssf)
                ssf = 0.0
            end
            ssf = max(ssf, 1.0e-30)
            if (abs(fQ) <= epsilon) || (count >= max_iters)
                break
            end
            count += 1
        end

        # Constrain the lateral subsurface flow rate ssf
        ssf = min(ssf, (ssfmax * dw))
        # On the basis of the lateral flow rate, estimate the amount of groundwater level above surface (saturation excess conditions), then rest = negative
        zi = zi_prev - (ssfin * dt + r * dx - ssf * dt) / (dw * dx) / theta_e
        if zi > d
            ssf = max(ssf - (dw * dx) * theta_e * (zi - d), 1.0e-30)
        end
        # Exfiltration rate and set zi to zero.
        exfilt = min(zi, 0.0) * -theta_e
        zi = clamp(zi, 0.0, d)

        return ssf, zi, exfilt
    end
end

"""
    accucapacitystate!(material, network, capacity)

Transport of material downstream with a limited transport capacity over a directed graph.
Mutates the material input. The network is expected to hold a graph and order field, where
the graph implements the Graphs interface, and the order is a valid topological ordering
such as that returned by `Graphs.topological_sort_by_dfs`.
"""
function accucapacitystate!(material, network, capacity)
    (; graph, order) = network
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
    return nothing
end

"""
    accucapacitystate!(material, network, capacity) -> material

Non mutating version of `accucapacitystate!`.
"""
function accucapacitystate(material, network, capacity)
    material = copy(material)
    accucapacitystate!(material, network, capacity)
    return material
end

"""
    accucapacityflux!(flux, material, network, capacity)

Transport of material downstream with a limited transport capacity over a directed graph.
Updates the material input, and overwrites the flux input, not using existing values. The
network is expected to hold a graph and order field, where the graph implements the Graphs
interface, and the order is a valid topological ordering such as that returned by
`Graphs.topological_sort_by_dfs`.
"""
function accucapacityflux!(flux, material, network, capacity)
    (; graph, order) = network
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
    return nothing
end

"""
    accucapacityflux!(material, network, capacity) -> flux, material

Non mutating version of `accucapacityflux!`.
"""
function accucapacityflux(material, network, capacity)
    flux = zero(material)
    material = copy(material)
    accucapacityflux!(flux, material, network, capacity)
    return flux, material
end

"""
    lateral_snow_transport!(snow, snowwater, slope, network)

Lateral snow transport. Transports snow downhill. Mutates `snow` and `snowwater`.
"""
function lateral_snow_transport!(snow, snowwater, slope, network)
    snowflux_frac = min.(0.5, slope ./ 5.67) .* min.(1.0, snow ./ 10000.0)
    maxflux = snowflux_frac .* snow
    accucapacitystate!(snow, network, maxflux)
    accucapacitystate!(snowwater, network, snowwater .* snowflux_frac)
    return nothing
end

"""
    local_inertial_flow(q0, zs0, zs1, hf, A, R, length, mannings_n, g, froude_limit, dt)

Local inertial approach for flow through area `A`. Returns the flow `q` between two adjacent
river cells (nodes) for a single timestep.
"""
function local_inertial_flow(
    q0,
    zs0,
    zs1,
    hf,
    A,
    R,
    length,
    mannings_n_sq,
    g,
    froude_limit,
    dt,
)
    slope = (zs1 - zs0) / length
    pow_R = R^(4.0f0 / 3.0f0)
    unit = one(hf)
    q = (
        (q0 - g * A * dt * slope) / (unit + g * dt * mannings_n_sq * abs(q0) / (pow_R * A))
    )

    # if froude number > 1.0, limit flow
    fr = ((q / A) / sqrt(g * hf)) * froude_limit
    q = ifelse((abs(fr) > Float(1.0)) * (q > Float(0.0)), sqrt(g * hf) * A, q)
    q = ifelse((abs(fr) > Float(1.0)) * (q < Float(0.0)), -sqrt(g * hf) * A, q)

    return q
end

"""
    local_inertial_flow(theta, q0, qd, qu, zs0, zs1, hf, width, length, mannings_n, g, froude_limit, dt)

Local inertial approach for flow through a rectangular area. Returns the flow `q` between
two adjacent cells (nodes) for a single timestep. Algorithm is based on de Almeida et al.
(2012).
"""
function local_inertial_flow(
    theta,
    q0,
    qd,
    qu,
    zs0,
    zs1,
    hf,
    width,
    length,
    mannings_n_sq,
    g,
    froude_limit,
    dt,
)
    slope = (zs1 - zs0) / length
    unit = one(theta)
    half = oftype(theta, 0.5)
    pow_hf = cbrt(hf * hf * hf * hf * hf * hf * hf)

    q = (
        ((theta * q0 + half * (unit - theta) * (qu + qd)) - g * hf * width * dt * slope) / (unit + g * dt * mannings_n_sq * abs(q0) / (pow_hf * width))
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
