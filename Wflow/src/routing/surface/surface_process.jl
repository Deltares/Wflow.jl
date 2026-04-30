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
    kinematic_wave(q_in, q_prev, q_lat, alpha, beta, dt, dx)

Solve the kinematic wave equation for discharge `q` [m³ s⁻¹] at a single cell and timestep.

Finds the root of the implicit finite-difference discretization:

    (dt/dx) * q + alpha * q^beta - C = 0

where C = (dt/dx) * q_in + alpha * q_prev^beta + dt * q_lat.

With beta = 3/5, the substitution u = q^(1/5) transforms this into a sparse quintic
polynomial in u:

    p(u) = a * u^5 + b * u^3 - C = 0,  with a = dt/dx, b = alpha

This avoids expensive `pow` calls in the Newton-Raphson loop. The polynomial is convex for
u > 0 (p''(u) = 20a*u³ + 6b*u > 0), so Newton-Raphson from an overestimate converges
monotonically without risk of divergence.
"""
function kinematic_wave(q_in, q_prev, q_lat, alpha, beta, dt, dx)
    if q_in + q_prev + q_lat ≈ 0.0
        return 0.0
    end

    a = dt / dx
    b = alpha
    C = a * q_in + b * pow(q_prev, beta) + dt * q_lat

    # Initial overestimate: u₀ = cbrt(C / b) satisfies p(u₀) > 0
    u = cbrt(C / b)

    # Newton-Raphson in u-space with equivalent convergence criterion |p(u)| <= epsilon
    epsilon = 1.0e-12
    max_iters = 3000
    for _ in 1:max_iters
        u2 = u * u
        u3 = u2 * u
        u4 = u2 * u2
        u5 = u4 * u
        p = a * u5 + b * u3 - C
        abs(p) <= epsilon && break
        dp = 5.0 * a * u4 + 3.0 * b * u2
        u -= p / dp
    end

    q = u * u * u * u * u  # u^5
    return max(q, KIN_WAVE_MIN_FLOW)
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
