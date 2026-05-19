"""
    lateral_snow_transport!(snow, slope, network)

Lateral snow transport. Transports snow downhill. Mutates `snow_storage` and `snow_water` of
a `snow` model.
"""
function lateral_snow_transport!(snow::AbstractSnowModel, domain::DomainLand)
    (; snow_storage, snow_water, snow_in, snow_out) = snow.variables
    (; slope) = domain.parameters
    snowflux_frac = min.(0.5, slope ./ 5.67) .* min.(ONE, snow_storage ./ 10000.0)
    maxflux = snowflux_frac .* snow_storage
    snow_out .= accucapacityflux(snow_storage, domain.network, maxflux)
    snow_out .+= accucapacityflux(snow_water, domain.network, snow_water .* snowflux_frac)
    flux_in!(snow_in, snow_out, domain.network)
end

lateral_snow_transport!(snow::NoSnowModel, domain::DomainLand) = nothing

"Kinematic wave surface flow rate for a single cell and timestep"
function kinematic_wave(q_in, q_prev, q_lat, alpha, dt, dx)
    if q_in + q_prev + q_lat ≈ ZERO
        return ZERO, ZERO
    else
        dt_dx = dt / dx
        # constant_term = (dt/dx)*q_in + alpha*q_prev^(3/5) + dt*q_lat
        # Use q_prev^(3/5) = (q_prev^(1/5))^3
        u_prev = q_prev >= ZERO ? Wflow.pow(q_prev, 1 // 5) : ZERO
        constant_term = dt_dx * q_in + alpha * u_prev * u_prev * u_prev + dt * q_lat

        # Initial estimate: use u_prev as starting point (cheap, no pow calls)
        # Since q changes smoothly between timesteps, u_prev is already close to the root.
        # For the case q_prev == 0, use C/alpha as a rough estimate for u^3, so u ≈ cbrt(C/alpha)
        u = if u_prev > ZERO
            u_prev
        else
            cbrt(constant_term / alpha)
        end

        # Halley's method on f(u) = dt_dx * u^5 + alpha * u^3 - C = 0
        # f'(u)  = 5*dt_dx * u^4 + 3*alpha * u^2
        # f''(u) = 20*dt_dx * u^3 + 6*alpha * u
        max_iters = 3000
        epsilon = to_precision(1e-12)

        const_1 = 5 * dt_dx
        const_2 = 3 * alpha
        const_3 = 20 * dt_dx
        const_4 = 6 * alpha

        for _ in 1:max_iters
            u2 = u * u
            u3 = u2 * u
            f_u = u3 * (dt_dx * u2 + alpha) - constant_term
            df_u = u2 * (const_1 * u2 + const_2)
            d2f_u = u * (const_3 * u2 + const_4)
            # Halley's correction: u -= f*f' / (f'^2 - f*f''/2)
            u -= (f_u * df_u) / (df_u * df_u - f_u * d2f_u / 2)
            if isnan(u) || u <= ZERO
                u = KIN_WAVE_MIN_FLOW_QROOT
            end
            if (abs(f_u) <= epsilon)
                break
            end
        end
        u = max(u, KIN_WAVE_MIN_FLOW_QROOT)
        u3 = u * u * u
        crossarea = alpha * u3
        q = u3 * u * u
        return q, crossarea
    end
end

"Kinematic wave surface flow rate over the whole network for a single timestep"
function kin_wave!(q, graph, toposort, q_prev, q_lat, alpha, flow_length, dt)
    for v in toposort
        upstream_nodes = inneighbors(graph, v)
        q_in = sum_at(q, upstream_nodes)
        q[v], _ = kinematic_wave(q_in, q_prev[v], q_lat, alpha[v], dt, flow_length[v])
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

    # if froude number > ONE, limit flow
    fr = ((q / A) / sqrt(GRAVITATIONAL_ACCELERATION * hf)) * froude_limit
    q = ifelse((abs(fr) > ONE) * (q > ZERO), sqrt(GRAVITATIONAL_ACCELERATION * hf) * A, q)
    q = ifelse((abs(fr) > ONE) * (q < ZERO), -sqrt(GRAVITATIONAL_ACCELERATION * hf) * A, q)

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
    # if froude number > ONE, limit flow
    if froude_limit
        fr = (q / width / hf) / sqrt(GRAVITATIONAL_ACCELERATION * hf)
        if abs(fr) > ONE && q > ZERO
            q = hf * sqrt(GRAVITATIONAL_ACCELERATION * hf) * width
        elseif abs(fr) > ONE && q < ZERO
            q = -hf * sqrt(GRAVITATIONAL_ACCELERATION * hf) * width
        end
    end

    return q
end
