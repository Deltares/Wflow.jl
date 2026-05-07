
"""
Return kinematic wave `celerity` of lateral subsurface flow based on hydraulic conductivity
profile `KhExponential`
"""
function ssf_celerity(zi, slope, specific_yield, kh_profile::KhExponential, i)
    (; kh_0, f) = kh_profile
    # [m s⁻¹] = ([m s⁻¹] * exp(- [m⁻¹] * [m]) * [-]) / [-]
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
    # [m s⁻¹] = ([m s⁻¹] * exp(- [m⁻¹] * [m]) * [-]) / [-]
    celerity = (kh_0[i] * exp(-f[i] * z) * slope) / specific_yield
    return celerity
end

"""
Return kinematic wave `celerity` of lateral subsurface flow based on hydraulic conductivity
profile `KhLayered`
"""
function ssf_celerity(zi, slope, specific_yield, kh_profile::KhLayered, i)
    (; kh) = kh_profile
    # [m s⁻¹] = ([-] * [m s⁻¹]) / [-]
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
    # [s m⁻¹] = [s] / [m]
    dt_dx = dt / dx
    # [s m⁻¹] = inv([m s⁻¹])
    celerity_inv = inv(celerity)
    # [s m⁻¹] = [s m⁻¹] + [s m⁻¹]
    df = dt_dx + celerity_inv
    while true
        # [m²] = [m³ s⁻¹] * ([s m⁻¹] + [s m⁻¹]) - [m²]
        f = dt_dx * q + celerity_inv * q - constant_term
        # [m³ s⁻¹] -= [m²] / [s m⁻¹]
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

Returns lateral subsurface flow `q`, water table depth `zi`, exfiltration rate `exfilt` and
net flux `net_flux`.
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
    if q_in + q_prev ≈ 0.0
        return 0.0, d, 0.0, 0.0
    else
        # initial estimate
        # [m³ s⁻¹] = ([m³ s⁻¹] + [m³ s⁻¹]) / [-]
        q = (q_prev + q_in) / 2.0
        # newton-raphson
        # [m s⁻¹]
        celerity = ssf_celerity(zi_prev, slope, sy, kh_profile, i)
        # [m²] = ([s] / [m]) * ([m³ s⁻¹] * [m³ s⁻¹]) + [m³ s⁻¹] / [m s⁻¹]
        constant_term = (dt / dx) * (q_in + q_net_bnds) + q_prev / celerity
        # [m³ s⁻¹]
        q = kw_ssf_newton_raphson(q, constant_term, celerity, dt, dx)

        # constrain maximum lateral subsurface flow rate q
        # [m³ s⁻¹] = min([m³ s⁻¹], ([m² s⁻¹] * [m]))
        q = min(q, (q_max * dw))
        # estimate water table depth zi, exfiltration rate and constrain zi and
        # lower boundary q
        # [m s⁻¹] = ([m³ s⁻¹] + [m³ s⁻¹] - [m³ s⁻¹]) / ([m] * [m])
        net_flux = (q_in + q_net_bnds - q) / (dw * dx)
        # [m], [m s⁻¹]
        dh, exfilt = water_table_change(soil, net_flux, sy, i, dt)
        # [m] = [m] - [m]
        zi = zi_prev - dh
        # [m] > [m]
        if zi > d
            # [m³ s⁻¹] = ([m] * [m]) * [-] * ([m] - [m]) / [s]
            q_excess = (dw * dx) * sy * (zi - d) / dt
            # [m³ s⁻¹] = max([m³ s⁻¹] - [m³ s⁻¹], [m³ s⁻¹])
            q = max(q - q_excess, KIN_WAVE_MIN_FLOW)
        end
        # [m] = clamp([m], [m], [m])
        zi = clamp(zi, 0.0, d)
        # constrain water table depth change to 0.1 m per (sub) timestep based on first `zi`
        # computation
        max_delta_zi = 0.1
        its = Int(ceil(round(abs(zi - zi_prev) / max_delta_zi; sigdigits = 12)))
        if its > 1
            dt_s = dt / its
            # [m³ s⁻¹]
            q_sum = 0.0
            # [m³ s⁻¹]
            exfilt_sum = 0.0
            # [m³ s⁻¹]
            net_flux_sum = 0.0
            # [m]
            zi_start = zi_prev
            for _ in 1:its
                # [m s⁻¹]
                celerity = ssf_celerity(zi_prev, slope, sy, kh_profile, i)
                # [m²] = ([s] / [m]) * [m³ s⁻¹] + [m³ s⁻¹] / [m s⁻¹] + [m³ s⁻¹] ([s] / [m])
                constant_term =
                    (dt_s / dx) * q_in + q_prev / celerity + q_net_bnds * (dt_s / dx)
                # [m³ s⁻¹]
                q = kw_ssf_newton_raphson(q_prev, constant_term, celerity, dt_s, dx)
                # constrain maximum lateral subsurface flow rate q
                # [m³ s⁻¹] = min([m³ s⁻¹], [m² s⁻¹] * [m])
                q = min(q, (q_max * dw))
                # estimate water table depth zi, exfiltration rate and constrain zi and
                # lower boundary q
                # [m s⁻¹] = ([m³ s⁻¹] + [m³ s⁻¹] - [m³ s⁻¹]) / ([m] * [m])
                net_flux = (q_in + q_net_bnds - q) / (dw * dx)
                # [m], [m s⁻¹]
                dh, exfilt = water_table_change(soil, net_flux, sy, i, dt_s)
                # [m] = [m] - [m]
                zi = zi_prev - dh
                # [m] > [m]
                if zi > d
                    # [m³ s⁻¹] = ([m] * [m]) * [-] * ([m] - [m]) / [s]
                    q_excess = (dw * dx) * sy * (zi - d) / dt_s
                    # [m³ s⁻¹] = max([m³ s⁻¹] - [m³ s⁻¹], [m³ s⁻¹])
                    q = max(q - q_excess, KIN_WAVE_MIN_FLOW)
                end
                # [m] = clamp([m], [m], [m])
                zi = clamp(zi, 0.0, d)
                update_ustorelayerdepth!(soil, zi_prev, zi, i)
                # [m s⁻¹] += [m s⁻¹]
                exfilt_sum += exfilt
                # [m s⁻¹] += [m s⁻¹]
                net_flux_sum += net_flux
                # [m s⁻¹] += [m s⁻¹]
                q_sum += q
                # [m s⁻¹] = [m s⁻¹]
                q_prev = q
                # [m] = [m]
                zi_prev = zi
            end
            # [m³ s⁻¹] = [m³ s⁻¹] / [-]
            q = q_sum / its
            # [m s⁻¹] = [m s⁻¹] / [-]
            exfilt = exfilt_sum / its
            # [m s⁻¹] = [m s⁻¹] / [-]
            net_flux = net_flux_sum / its
            # [m] = [m] - [m]
            dh = zi_start - zi
        else
            update_ustorelayerdepth!(soil, zi_prev, zi, i)
        end
        return q, zi, exfilt, net_flux
    end
end

"""
    kinematic_wave_ssf(q_in, q_prev, zi_prev, q_net_bnds, slope, sy, d, dt, dx, dw, q_max, kh_profile, soil, i)

Kinematic wave for lateral subsurface flow for a single cell and timestep with a `KhLayered`
conductivity profile, using (average) hydraulic conductivity `kh`.

Return lateral subsurface flow `q`, water table depth `zi`, exfiltration rate `exfilt` and
net flux `net_flux`.
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
        # [m³ s⁻¹] = ([m³ s⁻¹] + [m³ s⁻¹]) / [-]
        q_ini = (q_prev + q_in) / 2.0
        # newton-raphson
        # [m s⁻¹]
        celerity = ssf_celerity(zi_prev, slope, sy, kh_profile, i)
        # [m²] = ([s] / [m]) * ([m³ s⁻¹] * [m³ s⁻¹]) + [m³ s⁻¹] / [m s⁻¹]
        constant_term = (dt / dx) * (q_in + q_net_bnds) + q_prev / celerity
        # [m³ s⁻¹]
        q = kw_ssf_newton_raphson(q_ini, constant_term, celerity, dt, dx)
        # constrain maximum lateral subsurface flow rate q
        # [m³ s⁻¹] = min([m³ s⁻¹], ([m² s⁻¹] * [m]))
        q = min(q, (q_max * dw))

        # estimate water table depth zi, exfiltration rate and constrain zi and lower
        # boundary q
        # [m s⁻¹] = ([m³ s⁻¹] + [m³ s⁻¹] - [m³ s⁻¹]) / ([m] * [m])
        net_flux = (q_in + q_net_bnds - q) / (dw * dx)
        # [m], [m s⁻¹]
        dh, exfilt = water_table_change(soil, net_flux, sy, i, dt)
        # [m] = [m] - [m]
        zi = zi_prev - dh
        # [-] = ([m s⁻¹] - [m s⁻¹]) * [s] / [m]
        sy_d = dh > 0.0 ? (net_flux - exfilt) * dt / dh : sy
        # [m] > [m]
        if zi > d
            # [m³ s⁻¹] = ([m] * [m]) * [-] * ([m] - [m]) / [s]
            q_excess = (dw * dx) * sy_d * (zi - d) / dt
            # [m³ s⁻¹] = max([m³ s⁻¹] - [m³ s⁻¹], [m³ s⁻¹])
            q = max(q - q_excess, KIN_WAVE_MIN_FLOW)
        end
        # [m] = clamp([m], [m], [m])
        zi = clamp(zi, 0.0, d)

        # update unsaturated zone
        update_ustorelayerdepth!(soil, zi_prev, zi, i)

        return q, zi, exfilt, net_flux
    end
end
