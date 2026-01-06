"Convert a gridded drainage direction to a directed graph"
function flowgraph(ldd::AbstractVector, indices::AbstractVector, PCR_DIR::AbstractVector)
    # prepare a directed graph to be filled
    n = length(indices)
    graph = DiGraph(n)

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

const KIN_WAVE_MIN_FLOW = 1e-30 # [m³ s⁻¹]

"Kinematic wave surface flow rate for a single cell and timestep"
function kinematic_wave(q_in, q_prev, q_lat, alpha, beta, dt, dx)
    if q_in + q_prev + q_lat ≈ 0.0
        return 0.0
    else
        dt_dx = dt / dx
        exponent = beta - 1.0
        # initial estimate using linear scheme
        alpha_beta = alpha * beta
        ab_pq = alpha_beta * pow(((q_prev + q_in) / 2.0), exponent)
        q = (dt_dx * q_in + q_prev * ab_pq + dt * q_lat) / (dt_dx + ab_pq)
        if isnan(q)
            q = 0.0
        end
        q = max(q, KIN_WAVE_MIN_FLOW)
        # newton-raphson
        max_iters = 3000
        epsilon = 1.0e-12
        count = 0
        constant_term = dt_dx * q_in + alpha * pow(q_prev, beta) + dt * q_lat
        while true
            f_q = dt_dx * q + alpha * pow(q, beta) - constant_term
            df_q = dt_dx + alpha * beta * pow(q, exponent)
            q -= (f_q / df_q)
            if isnan(q)
                q = 0.0
            end
            q = max(q, KIN_WAVE_MIN_FLOW)
            if (abs(f_q) <= epsilon) || (count >= max_iters)
                break
            end
            count += 1
        end
        return q
    end
end

"Kinematic wave surface flow rate over the whole network for a single timestep"
function kin_wave!(q, graph, toposort, q_prev, q_lat, alpha, beta, flow_length, dt)
    for v in toposort
        upstream_nodes = inneighbors(graph, v)
        q_in = sum_at(q, upstream_nodes)
        q[v] = kinematic_wave(q_in, q_prev[v], q_lat, alpha[v], beta, dt, flow_length[v])
    end
    return q
end

"Return kinematic wave `celerity` of lateral subsurface flow based on hydraulic conductivity profile `KhExponential`"
function ssf_celerity(zi, slope, theta_e, kh_profile::KhExponential, i)
    (; kh_0, f) = kh_profile
    celerity = (kh_0[i] * exp(-f[i] * zi) * slope) / theta_e
    return celerity
end

"Return kinematic wave `celerity` of lateral subsurface flow based on hydraulic conductivity profile `KhExponentialConstant`"
function ssf_celerity(zi, slope, theta_e, kh_profile::KhExponentialConstant, i)
    (; z_exp) = kh_profile
    (; kh_0, f) = kh_profile.exponential
    z = zi < z_exp[i] ? zi : z_exp[i]
    celerity = (kh_0[i] * exp(-f[i] * z) * slope) / theta_e
    return celerity
end

"""
Return kinematic wave subsurface flow `ssf` for a single cell and timestep using the Newton-
Raphson method.
"""
function kw_ssf_newton_raphson(ssf, constant_term, celerity, dt, dx)
    epsilon = 1.0e-12
    max_iters = 3000
    count = 0
    dt_dx = dt / dx
    celerity_inv = inv(celerity)
    while true
        f = dt_dx * ssf + celerity_inv * ssf - constant_term
        df = dt_dx + celerity_inv
        ssf -= (f / df)
        if isnan(ssf)
            ssf = 0.0
        end
        ssf = max(ssf, KIN_WAVE_MIN_FLOW)
        if (abs(f) <= epsilon) || (count >= max_iters)
            break
        end
        count += 1
    end
    return ssf
end

"""
    kinematic_wave_ssf(ssfin, ssf_prev, zi_prev, r, slope, theta_e, d, dt, dx, dw, ssfmax, kh_profile, i)

Kinematic wave for lateral subsurface flow for a single cell and timestep. The hydraulic
conductivity profile `kh_profile` is either `KhExponential` or `KhExponentialConstant`.

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
    if ssfin + ssf_prev ≈ 0.0 && r <= 0.0
        return 0.0, d, 0.0
    else
        # initial estimate
        ssf = (ssf_prev + ssfin) / 2.0
        # newton-raphson
        celerity = ssf_celerity(zi_prev, slope, theta_e, kh_profile, i)
        constant_term = (dt / dx) * ssfin + (1.0 / celerity) * ssf_prev + r * dt
        ssf = kw_ssf_newton_raphson(ssf, constant_term, celerity, dt, dx)

        # constrain maximum lateral subsurface flow rate ssf
        ssf = min(ssf, (ssfmax * dw))
        # estimate water table depth zi, exfiltration rate and constrain zi and
        # lower boundary ssf
        zi = zi_prev - (ssfin * dt + r * dt * dx - ssf * dt) / (dw * dx) / theta_e
        if zi > d
            ssf = max(ssf - (dw * dx) * theta_e * (zi - d), KIN_WAVE_MIN_FLOW)
        end
        exfilt = min(zi, 0.0) * -theta_e
        zi = clamp(zi, 0.0, d)

        # constrain water table depth change to 0.1 m per (sub) timestep based on first `zi`
        # computation
        max_delta_zi = 0.1
        its = Int(cld(abs(max(zi, 0.0) - zi_prev), max_delta_zi))
        if its > 1
            dt_s = dt / its
            ssf_sum = 0.0
            exfilt_sum = 0.0
            for _ in 1:its
                celerity = ssf_celerity(zi_prev, slope, theta_e, kh_profile, i)
                constant_term = (dt_s / dx) * ssfin + ssf_prev / celerity + r * dt_s
                ssf = kw_ssf_newton_raphson(ssf_prev, constant_term, celerity, dt_s, dx)
                # constrain maximum lateral subsurface flow rate ssf
                ssf = min(ssf, (ssfmax * dw))
                # estimate water table depth zi, exfiltration rate and constrain zi and
                # lower boundary ssf
                zi =
                    zi_prev -
                    (ssfin * dt_s + r * dt_s * dx - ssf * dt_s) / (dw * dx) / theta_e
                if zi > d
                    ssf = max(ssf - (dw * dx) * theta_e * (zi - d), KIN_WAVE_MIN_FLOW)
                end
                exfilt_sum += min(zi, 0.0) * -theta_e
                zi = clamp(zi, 0.0, d)
                ssf_sum += ssf
                ssf_prev = ssf
                zi_prev = zi
            end
            ssf = ssf_sum / its
            exfilt = exfilt_sum
        end

        return ssf, zi, exfilt
    end
end

"""
    kinematic_wave_ssf(ssfin, ssf_prev, zi_prev, r, slope, theta_e, d, dt, dx, dw, ssfmax, kh_profile, i)

Kinematic wave for lateral subsurface flow for a single cell and timestep with a `KhLayered`
conductivity profile, using (average) hydraulic conductivity `kh`.

Return lateral subsurface flow `ssf`, water table depth `zi` and exfiltration rate `exfilt`.
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
    if ssfin + ssf_prev ≈ 0.0 && r <= 0.0
        return 0.0, d, 0.0
    else
        # initial estimate
        ssf_ini = (ssf_prev + ssfin) / 2.0
        # newton-raphson
        celerity = (slope * kh_profile.kh[i]) / theta_e
        constant_term = (dt / dx) * ssfin + ssf_prev / celerity + r * dt
        ssf = kw_ssf_newton_raphson(ssf_ini, constant_term, celerity, dt, dx)
        # constrain maximum lateral subsurface flow rate ssf
        ssf = min(ssf, (ssfmax * dw))
        # estimate water table depth zi, exfiltration rate and constrain zi and lower
        # boundary ssf
        zi = zi_prev - (ssfin * dt + r * dt * dx - ssf * dt) / (dw * dx) / theta_e
        if zi > d
            ssf = max(ssf - (dw * dx) * theta_e * (zi - d), KIN_WAVE_MIN_FLOW)
        end
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
    accucapacityflux(material, network, capacity) -> flux

Non mutating version of `accucapacityflux!`.
"""
function accucapacityflux(material, network, capacity)
    flux = zero(material)
    accucapacityflux!(flux, material, network, capacity)
    return flux
end

"""
    accucapacityflux_state(material, network, capacity) -> flux, material

Non mutating version of combined `accucapacityflux!` and `accucapacitystate!`.
"""
function accucapacityflux_state(material, network, capacity)
    flux = zero(material)
    material = copy(material)
    accucapacityflux!(flux, material, network, capacity)
    return flux, material
end

function flux_in!(flux_in, flux, network)
    (; upstream_nodes, order) = network
    for (i, v) in enumerate(order)
        flux_in[v] = sum_at(flux, upstream_nodes[i])
    end
    return nothing
end

"""
    lateral_snow_transport!(snow, slope, network)

Lateral snow transport. Transports snow downhill. Mutates `snow_storage` and `snow_water` of
a `snow` model.
"""
function lateral_snow_transport!(snow::AbstractSnowModel, domain::DomainLand)
    (; snow_storage, snow_water, snow_in, snow_out) = snow.variables
    (; slope) = domain.parameters
    snowflux_frac = min.(0.5, slope ./ 5.67) .* min.(1.0, snow_storage ./ 10000.0)
    maxflux = snowflux_frac .* snow_storage
    snow_out .= accucapacityflux(snow_storage, domain.network, maxflux)
    snow_out .+= accucapacityflux(snow_water, domain.network, snow_water .* snowflux_frac)
    flux_in!(snow_in, snow_out, domain.network)
end

lateral_snow_transport!(snow::NoSnowModel, domain::DomainLand) = nothing

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
    froude_limit,
    dt,
)
    slope = (zs1 - zs0) / length
    pow_R = cbrt(R * R * R * R)
    unit = one(hf)
    q = (
        (q0 - GRAVITATIONAL_ACCELERATION * A * dt * slope) / (
            unit + GRAVITATIONAL_ACCELERATION * dt * mannings_n_sq * abs(q0) / (pow_R * A)
        )
    )

    # if froude number > 1.0, limit flow
    fr = ((q / A) / sqrt(GRAVITATIONAL_ACCELERATION * hf)) * froude_limit
    q = ifelse((abs(fr) > 1.0) * (q > 0.0), sqrt(GRAVITATIONAL_ACCELERATION * hf) * A, q)
    q = ifelse((abs(fr) > 1.0) * (q < 0.0), -sqrt(GRAVITATIONAL_ACCELERATION * hf) * A, q)

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
    froude_limit,
    dt,
)
    slope = (zs1 - zs0) / length
    unit = one(theta)
    half = oftype(theta, 0.5)
    pow_hf = cbrt(hf * hf * hf * hf * hf * hf * hf)

    q = (
        (
            (theta * q0 + half * (unit - theta) * (qu + qd)) -
            GRAVITATIONAL_ACCELERATION * hf * width * dt * slope
        ) / (
            unit +
            GRAVITATIONAL_ACCELERATION * dt * mannings_n_sq * abs(q0) / (pow_hf * width)
        )
    )
    # if froude number > 1.0, limit flow
    if froude_limit
        fr = (q / width / hf) / sqrt(GRAVITATIONAL_ACCELERATION * hf)
        if abs(fr) > 1.0 && q > 0.0
            q = hf * sqrt(GRAVITATIONAL_ACCELERATION * hf) * width
        elseif abs(fr) > 1.0 && q < 0.0
            q = -hf * sqrt(GRAVITATIONAL_ACCELERATION * hf) * width
        end
    end

    return q
end
