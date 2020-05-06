
function kinematic_wave_ssf(
    ssfin,
    ssfₜ₋₁,
    ziₜ₋₁,
    r,
    kh₀,
    β,
    θₑ,
    f,
    d,
    Δt,
    Δx,
    dw,
    ssfmax,
)

    ϵ = 1.0e-3
    max_iters = 3000

    if ssfin + ssfₜ₋₁ ≈ 0.0 && r <= 0.0
        return 0.0, d, 0.0
    else
        #initial estimate
        ssf = (ssfₜ₋₁ + ssfin) / 2.0
        count = 1

        # Estimate zi on the basis of the relation between subsurfacel flow and zi
        zi = log((f * ssf) / (dw * kh₀ * β) + exp(-f * d)) / -f
        # Reciprocal of derivative delta Q/ delta z_i, constrained w.r.t. neff on the basis of the continuity equation)
        Cn = (kh₀ * exp(-f * zi) * β) / θₑ
        # Term of the continuity equation for Newton-Raphson iteration for iteration 1
        # because celerity Cn is depending on zi, the increase or decrease of zi is moved to the recharge term of the continuity equation
        # then (1./Cn)*ssfₜ₋₁ can be replaced with (1./Cn)*ssf, and thus celerity and lateral flow rate ssf are then in line
        c =
            (Δt / Δx) * ssfin +
            (1.0 / Cn) * ssf +
            Δt * (r - (ziₜ₋₁ - zi) * θₑ * dw)

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
            c =
                (Δt / Δx) * ssfin +
                (1.0 / Cn) * ssf +
                Δt * (r - (ziₜ₋₁ - zi) * θₑ * dw)

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
        zi = max(0, zi)

        return ssf, zi, exfilt, count, fQ

    end
end
