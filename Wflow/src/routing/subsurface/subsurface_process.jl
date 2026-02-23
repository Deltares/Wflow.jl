
"""
Return kinematic wave `celerity` of lateral subsurface flow based on hydraulic conductivity
profile `KhExponential`
"""
function ssf_celerity(zi, slope, specific_yield, kh_profile::KhExponential, i)
    (; kh_0, f) = kh_profile
    celerity = (kh_0[i] * exp(-f[i] * zi) * slope) / specific_yield
    return celerity
end

"""
Return kinematic wave `celerity` of lateral subsurface flow based on hydraulic conductivity
profile `KhExponentialConstant`
"""
function ssf_celerity(zi, slope, specific_yield, kh_profile::KhExponentialConstant, i)
    (; z_exp) = kh_profile
    (; kh_0, f) = kh_profile.exponential
    z = zi < z_exp[i] ? zi : z_exp[i]
    celerity = (kh_0[i] * exp(-f[i] * z) * slope) / specific_yield
    return celerity
end

"""
Return kinematic wave `celerity` of lateral subsurface flow based on hydraulic conductivity
profile `KhLayered`
"""
function ssf_celerity(zi, slope, specific_yield, kh_profile::KhLayered, i)
    (; kh) = kh_profile
    celerity = (slope * kh[i]) / specific_yield
    return celerity
end

"""
Return kinematic wave subsurface flow `q` for a single cell and timestep using the Newton-
Raphson method.
"""
function kw_ssf_newton_raphson(q, constant_term, celerity, dt, dx)
    epsilon = 1.0e-12
    max_iters = 3000
    count = 0
    dt_dx = dt / dx
    celerity_inv = inv(celerity)
    while true
        f = dt_dx * q + celerity_inv * q - constant_term
        df = dt_dx + celerity_inv
        q -= (f / df)
        if isnan(q)
            q = 0.0
        end
        q = max(q, KIN_WAVE_MIN_FLOW)
        if (abs(f) <= epsilon) || (count >= max_iters)
            break
        end
        count += 1
    end
    return q
end

"""
    kinematic_wave_ssf(q_in, q_prev, zi_prev, q_net_bnds, slope, sy, d, dt, dx, dw, q_max, kh_profile, soil, i)

Kinematic wave for lateral subsurface flow for a single cell and timestep. The hydraulic
conductivity profile `kh_profile` is either `KhExponential` or `KhExponentialConstant`.

Returns lateral subsurface flow `q`, water table depth `zi`, exfiltration rate `exfilt`
and dynamic specific yield `sy_d`.
"""
function kinematic_wave_ssf(
    q_in,
    q_prev,
    zi_prev,
    q_net_bnds,
    slope,
    sy,
    d,
    dt,
    dx,
    dw,
    q_max,
    kh_profile::Union{KhExponential, KhExponentialConstant},
    soil::SbmSoilModel,
    i,
)
    if q_in + q_prev ≈ 0.0 && r <= 0.0
        return 0.0, d, 0.0
    else
        # initial estimate
        q = (q_prev + q_in) / 2.0
        # newton-raphson
        celerity = ssf_celerity(zi_prev, slope, sy, kh_profile, i)
        constant_term = (dt / dx) * q_in + (1.0 / celerity) * q_prev + q_net_bnds * (dt / dx)
        q = kw_ssf_newton_raphson(q, constant_term, celerity, dt, dx)

        # constrain maximum lateral subsurface flow rate q
        q = min(q, (q_max * dw))
        # estimate water table depth zi, exfiltration rate and constrain zi and
        # lower boundary q
        net_flux = (q_in * dt + q_net_bnds * dt - q * dt) / (dw * dx)
        dh, exfilt = water_table_change(soil, net_flux, sy, i)
        zi = zi_prev - dh
        if zi > d
            q_excess = (dw * dx) * sy * (zi - d) / dt
            q = max(q - q_excess, KIN_WAVE_MIN_FLOW)
        end
        zi = clamp(zi, 0.0, d)
        # constrain water table depth change to 0.1 m per (sub) timestep based on first `zi`
        # computation
        max_delta_zi = 0.1
        its = Int(cld(abs(max(zi, 0.0) - zi_prev), max_delta_zi))
        if its > 1
            dt_s = dt / its
            q_sum = 0.0
            exfilt_sum = 0.0
            net_flux_sum = 0.0
            zi_start = zi_prev
            for _ in 1:its
                celerity = ssf_celerity(zi_prev, slope, sy, kh_profile, i)
                constant_term = (dt_s / dx) * q_in + q_prev / celerity + q_net_bnds * (dt_s / dx)
                q = kw_ssf_newton_raphson(q_prev, constant_term, celerity, dt_s, dx)
                # constrain maximum lateral subsurface flow rate q
                q = min(q, (q_max * dw))
                # estimate water table depth zi, exfiltration rate and constrain zi and
                # lower boundary q
                net_flux = (q_in * dt_s + q_net_bnds * dt_s - q * dt_s) / (dw * dx)
                dh, exfilt = water_table_change(soil, net_flux, sy, i)
                zi = zi_prev - dh
                if zi > d
                    q_excess = (dw * dx) * sy * (zi - d) / dt_s
                    q = max(q - q_excess, KIN_WAVE_MIN_FLOW)
                end
                zi = clamp(zi, 0.0, d)
                # update unsaturated zone
                zi_prev_mm = zi_prev * 1000.0
                zi_mm = zi * 1000.0
                update_ustorelayerdepth!(soil, zi_prev_mm, zi_mm, i)
                exfilt_sum += exfilt
                net_flux_sum += net_flux
                q_sum += q
                q_prev = q
                zi_prev = zi
            end
            q = q_sum / its
            exfilt = exfilt_sum
            net_flux = net_flux_sum
            dh = zi_start - zi
        else
            zi_prev_mm = zi_prev * 1000.0
            zi_mm = zi * 1000.0
            update_ustorelayerdepth!(soil, zi_prev_mm, zi_mm, i)
        end
        return q, zi, exfilt, net_flux
    end
end

"""
    kinematic_wave_ssf(q_in, q_prev, zi_prev, q_net_bnds, slope, sy, d, dt, dx, dw, q_max, kh_profile, soil, i)

Kinematic wave for lateral subsurface flow for a single cell and timestep with a `KhLayered`
conductivity profile, using (average) hydraulic conductivity `kh`.

Return lateral subsurface flow `q`, water table depth `zi`, exfiltration rate `exfilt` and
dynamic specific yield `sy_d`.
"""
function kinematic_wave_ssf(
    q_in,
    q_prev,
    zi_prev,
    q_net_bnds,
    slope,
    sy,
    d,
    dt,
    dx,
    dw,
    q_max,
    kh_profile::KhLayered,
    soil::SbmSoilModel,
    i,
)
    if q_in + q_prev ≈ 0.0 && r <= 0.0
        return 0.0, d, 0.0
    else
        # initial estimate
        q_ini = (q_prev + q_in) / 2.0
        # newton-raphson
        celerity = ssf_celerity(zi_prev, slope, sy, kh_profile, i)
        constant_term = (dt / dx) * q_in + q_prev / celerity + q_net_bnds * (dt / dx)
        q = kw_ssf_newton_raphson(q_ini, constant_term, celerity, dt, dx)
        # constrain maximum lateral subsurface flow rate q
        q = min(q, (q_max * dw))

        # estimate water table depth zi, exfiltration rate and constrain zi and lower
        # boundary q
        net_flux = (q_in * dt + q_net_bnds * dt - q * dt) / (dw * dx)
        dh, exfilt = water_table_change(soil, net_flux, sy, i)
        zi = zi_prev - dh
        sy_d = dh > 0.0 ? (net_flux - exfilt) / dh : sy
        if zi > d
            q_excess = (dw * dx) * sy_d * (zi - d) / dt
            q = max(q - q_excess, KIN_WAVE_MIN_FLOW)
        end
        zi = clamp(zi, 0.0, d)

        # update unsaturated zone
        zi_prev_mm = zi_prev * 1000.0
        zi_mm = zi * 1000.0
        update_ustorelayerdepth!(soil, zi_prev_mm, zi_mm, i)

        return q, zi, exfilt, net_flux
    end
end
