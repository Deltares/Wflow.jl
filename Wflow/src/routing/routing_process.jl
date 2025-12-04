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

"Kinematic wave surface flow rate for a single cell and timestep"
function kinematic_wave(q_in, q_prev, q_lat, alpha, beta, dt, dx)
    if q_in + q_prev + q_lat ≈ 0.0
        return 0.0
    else
        C = dt / dx
        exponent = beta - 1.0
        product = alpha * beta
        # initial estimate using linear scheme
        ab_pq = product * pow((q_prev + q_in) / 2, exponent)
        q = (C * q_in + q_prev * ab_pq + dt * q_lat) / (C + ab_pq)
        if isnan(q)
            q = 1e-30
        else
            q = max(q, 1e-30)
        end
        # newton-raphson
        max_iters = 3000
        epsilon = 1.0e-12
        count = 0
        constant_term = C * q_in + alpha * pow(q_prev, beta) + dt * q_lat
        while true
            f_q = C * q + alpha * pow(q, beta) - constant_term
            df_q = C + product * pow(q, exponent)
            q -= (f_q / df_q)
            if isnan(q)
                q = 1e-30
            else
                q = max(q, 1.0e-30)
            end
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
    C = dt / dx
    celerity_inv = inv(celerity)
    while true
        f = C * ssf + ssf / celerity - constant_term
        df = C + celerity_inv
        ssf -= (f / df)
        if isnan(ssf)
            ssf = 0.0
        end
        ssf = max(ssf, 1.0e-30)
        if (abs(f) <= epsilon) || (count >= max_iters)
            break
        end
        count += 1
    end
    return ssf
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
    if ssfin + ssf_prev ≈ 0.0 && r <= 0.0
        return 0.0, d, 0.0
    else
        # initial estimate
        ssf = (ssf_prev + ssfin) / 2.0
        # newton-raphson
        celerity = ssf_celerity(zi_prev, slope, theta_e, kh_profile, i)
        constant_term = (dt / dx) * ssfin + ssf_prev / celerity + r * dt
        ssf = kw_ssf_newton_raphson(ssf, constant_term, celerity, dt, dx)

        # constrain maximum lateral subsurface flow rate ssf
        ssf = min(ssf, (ssfmax * dw))
        # estimate water table depth zi, exfiltration rate and constrain zi and
        # lower boundary ssf
        zi = zi_prev - (ssfin * dt + r * dt * dx - ssf * dt) / (dw * dx) / theta_e
        if zi > d
            ssf = max(ssf - (dw * dx) * theta_e * (zi - d), 1.0e-30)
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
                ssf = min(ssf, ssfmax * dw)
                # estimate water table depth zi, exfiltration rate and constrain zi and
                # lower boundary ssf
                zi =
                    zi_prev -
                    (ssfin * dt_s + r * dt_s * dx - ssf * dt_s) / (dw * dx) / theta_e
                if zi > d
                    ssf = max(ssf - (dw * dx) * theta_e * (zi - d), 1.0e-30)
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
        ssf = min(ssf, ssfmax * dw)
        # estimate water table depth zi, exfiltration rate and constrain zi and lower
        # boundary ssf
        zi = zi_prev - (ssfin * dt + r * dt * dx - ssf * dt) / (dw * dx) / theta_e
        if zi > d
            ssf = max(ssf - (dw * dx) * theta_e * (zi - d), 1.0e-30)
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
function accucapacityflux!(flux, material, network, capacity, dt)
    (; graph, order) = network
    for v in order
        downstream_nodes = outneighbors(graph, v)
        n = length(downstream_nodes)
        flux_val = min(material[v] / dt, capacity[v])
        material[v] -= flux_val * dt
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
function accucapacityflux(material, network, capacity, dt)
    flux = zero(material)
    accucapacityflux!(flux, material, network, capacity, dt)
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

const tan80 = 5.67

"""
    lateral_snow_transport!(snow, domain, dt)

Lateral snow transport. Transports snow downhill. Mutates `snow_storage` and `snow_water` of
a `snow` model.
"""
function lateral_snow_transport!(snow::AbstractSnowModel, domain::DomainLand, dt::Number)
    (; snow_storage, snow_water, snow_in, snow_out) = snow.variables
    (; slope) = domain.parameters
    # [m]
    snow_storage_max = 10.0
    # [-] = min([-], [-]) * min([-], [m] / [m])
    snowflux_frac = @. min(0.5, slope / tan80) * min(1.0, snow_storage / snow_storage_max)
    # [m s⁻¹] = [-] * [m] / [s]
    maxflux = snowflux_frac .* snow_storage / dt
    # [m s⁻¹]
    snow_out .= accucapacityflux(snow_storage, domain.network, maxflux, dt)
    snow_out .+=
        accucapacityflux(snow_water, domain.network, snow_water .* snowflux_frac, dt)
    flux_in!(snow_in, snow_out, domain.network)
end

lateral_snow_transport!(snow::NoSnowModel, domain::DomainLand, dt::Number) = nothing

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
    # [-] = ([m] - [m]) / [m]
    slope = (zs1 - zs0) / length
    # [m^4/3]
    pow_R = cbrt(R^4)
    unit = one(hf)
    # [m³ s⁻¹] = ([m³ s⁻¹] - [m s⁻²] * [m²] * [s] * [-]) / ([-] + [m s⁻²] * [s] * [(s m-1/3)²] * [m³ s⁻¹] / ([m^4/3] * [m²]))
    q = (
        (q0 - g * A * dt * slope) / (unit + g * dt * mannings_n_sq * abs(q0) / (pow_R * A))
    )

    # if froude number > 1.0, limit flow
    # [-] = (([m³ s⁻¹] / [m²]) / ([m s⁻²] * [m])^1/2) * [-]
    fr = ((q / A) / sqrt(g * hf)) * froude_limit
    q = ifelse((abs(fr) > 1.0) * (q > 0.0), sqrt(g * hf) * A, q)
    q = ifelse((abs(fr) > 1.0) * (q < 0.0), -sqrt(g * hf) * A, q)

    return q
end

"""
    local_inertial_flow(theta, q0, qd, qu, zs0, zs1, hf, width, length, mannings_n, g, froude_limit, dt)

Local inertial approach for flow through a rectangular area. Returns the flow `q` between
two adjacent cells (nodes) for a single timestep. Algorithm is based on de Almeida et al.
(2012).
"""
function local_inertial_flow(
    theta::T,
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
) where {T}
    # [-] = ([m] - [m]) / [m]
    slope = (zs1 - zs0) / length

    unit = one(T)
    half = T(0.5)
    # [m^7/3]
    pow_hf = cbrt(hf^7)

    # [m³ s⁻¹] = (([-] * [m³ s⁻¹] + [-] * ([-] - [-]) * ([m³ s⁻¹] + [m³ s⁻¹])) - [m s⁻²] * [m] * [m] * [s] * [-]) / ([-] + [m s⁻²] * [s] * [(s m-1/3)²] * [m³ s⁻¹] / ([m^7/3] * [m]))
    q = (
        ((theta * q0 + half * (unit - theta) * (qu + qd)) - g * hf * width * dt * slope) / (unit + g * dt * mannings_n_sq * abs(q0) / (pow_hf * width))
    )
    # if froude number > 1.0, limit flow
    if froude_limit
        fr = (q / (width * hf)) / sqrt(g * hf)
        if abs(fr) > 1.0 && q > 0.0
            q = hf * sqrt(g * hf) * width
        elseif abs(fr) > 1.0 && q < 0.0
            q = -hf * sqrt(g * hf) * width
        end
    end

    return q
end
