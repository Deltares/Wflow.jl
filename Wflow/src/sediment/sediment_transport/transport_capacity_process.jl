"""
    mask_transport_capacity(
        transport_capacity,
        reservoirs,
        rivers,
    )

Mask transport capacity values for         reservoirs,
 and rivers.

# Arguments
- `transport_capacity` (total sediment transport capacity [t dt-1])
- `reservoirs` (reservoirs mask [-])
- `rivers` (rivers mask [-])

# Output
- `transport_capacity` (masked total sediment transport capacity [t dt-1])
"""
function mask_transport_capacity(transport_capacity, reservoirs, rivers)
    # Full deposition in rivers to switch to river concept
    if rivers
        tc = 0.0
        # Sediment flux in reservoirs will all reach the river
    elseif reservoirs
        tc = 1e9
    else
        tc = transport_capacity
    end

    return tc
end

"""
    limit_and_convert_transport_capacity(
        transport_capacity,
        q,
        waterlevel,
        width,
        length,
        dt,
    )

Limit to stremaflow and not debris flow and convert transport in ton/m3 to ton.

# Arguments
- `transport_capacity` (total sediment transport capacity [t m-3])
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `width` (drain width [m])
- `length` (drain length [m])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1])
"""
function limit_and_convert_transport_capacity(
    transport_capacity,
    q,
    waterlevel,
    width,
    length,
    dt,
)
    # 1285 g/L: boundary between streamflow and debris flow (Costa, 1988)
    transport_capacity = min(transport_capacity, 1.285)
    # Transport capacity [ton]
    transport_capacity *= waterlevel * width * length + q * dt

    return transport_capacity
end

"""
    transport_capacity_govers(
        q,
        waterlevel,
        c_govers,
        n_govers,
        density,
        slope,
        width,
        reservoirs,
        rivers,
        dt,
    )

Total sediment transport capacity based on Govers.

# Arguments
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `c_govers` (Govers transport capacity coefficient [-])
- `n_govers` (Govers transport capacity exponent [-])
- `density` (sediment density [kg m-3])
- `slope` (slope [-])
- `width` (drain width [m])
- `reservoirs` (reservoirs mask [-])
- `rivers` (rivers mask [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1])
"""
function transport_capacity_govers(
    q,
    waterlevel,
    c_govers,
    n_govers,
    density,
    slope,
    width,
    reservoirs,
    rivers,
    dt,
)
    # Transport capacity from govers 1990
    sinslope = sin(atan(slope)) #slope in radians
    # Unit stream power
    if waterlevel > 0.0
        velocity = q / (width * waterlevel) #m/s
    else
        velocity = 0.0
    end
    omega = 10 * sinslope * 100 * velocity #cm/s
    if omega > 0.4
        TCf = c_govers * (omega - 0.4)^(n_govers) * density #kg/m3
    else
        TCf = 0.0
    end
    transport_capacity = TCf * q * dt * 1e-3 #[ton/cell]

    # Mask transport capacity values for reservoirs and rivers
    transport_capacity = mask_transport_capacity(transport_capacity, reservoirs, rivers)

    return transport_capacity
end

"""
    transport_capacity_yalin(
        q,
        waterlevel,
        density,
        d50,
        slope,
        width,
        reservoirs,
        rivers,
        dt,
    )

Total sediment transport capacity based on Yalin.

# Arguments
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m-3])
- `d50` (median grain size [m])
- `slope` (slope [-])
- `width` (drain width [m])
- `reservoirs` (reservoirs mask [-])
- `rivers` (rivers mask [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1])
"""
function transport_capacity_yalin(
    q,
    waterlevel,
    density,
    d50,
    slope,
    width,
    reservoirs,
    rivers,
    dt,
)
    sinslope = sin(atan(slope)) #slope in radians
    # Transport capacity from Yalin without particle differentiation
    delta =
        max((waterlevel * sinslope / (d50 * 0.001 * (density / 1000 - 1)) / 0.06 - 1), 0.0)
    alphay = delta * 2.45 / (0.001 * density)^0.4 * 0.06^(0.5)
    if q > 0.0 && alphay != 0.0
        TC = (
            width / q *
            (density - 1000) *
            d50 *
            0.001 *
            (9.81 * waterlevel * sinslope) *
            0.635 *
            delta *
            (1 - log(1 + alphay) / (alphay))
        ) # [kg/m3]
        transport_capacity = TC * q * dt * 1e-3 #[ton]
    else
        transport_capacity = 0.0
    end

    # Mask transport capacity values for reservoirs and rivers
    transport_capacity = mask_transport_capacity(transport_capacity, reservoirs, rivers)

    return transport_capacity
end

"""
    transportability_yalin_differentiation(
        waterlevel,
        density,
        dm_clay,
        dm_silt,
        dm_sand,
        dm_sagg,
        dm_lagg,
        slope,
    )

Total flow transportability based on Yalin with particle differentiation.

# Arguments
- `waterlevel` (water level [m])
- `density` (sediment density [kg m-3])
- `dm_clay` (clay median grain size [m])
- `dm_silt` (silt median grain size [m])
- `dm_sand` (sand median grain size [m])
- `dm_sagg` (small aggregates median grain size [m])
- `dm_lagg` (large aggregates median grain size [m])
- `slope` (slope [-])

# Output
- `dtot` (total transportability of the flow [-])
"""
function transportability_yalin_differentiation(
    waterlevel,
    density,
    dm_clay,
    dm_silt,
    dm_sand,
    dm_sagg,
    dm_lagg,
    slope,
)
    sinslope = sin(atan(slope)) #slope in radians    
    # Delta parameter of Yalin for each particle class
    delta = waterlevel * sinslope / (1e-6 * (density / 1000 - 1)) / 0.06
    dclay = max(1 / dm_clay * delta - 1, 0.0)
    dsilt = max(1 / dm_silt * delta - 1, 0.0)
    dsand = max(1 / dm_sand * delta - 1, 0.0)
    dsagg = max(1 / dm_sagg * delta - 1, 0.0)
    dlagg = max(1 / dm_lagg * delta - 1, 0.0)
    # Total transportability
    dtot = dclay + dsilt + dsand + dsagg + dlagg

    return dtot
end

"""
    transport_capacity_yalin_differentiation(
        q,
        waterlevel,
        density,
        dm,
        slope,
        width,
        reservoirs,
        rivers,
        dtot,
        dt,
    )

Transport capacity for a specific grain size based on Yalin with particle differentiation.

# Arguments
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m-3])
- `dm` (median grain size [m])
- `slope` (slope [-])
- `width` (drain width [m])
- `reservoirs` (reservoirs mask [-])
- `rivers` (rivers mask [-])
- `dtot` (total flow transportability [t dt-1])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1])
"""
function transport_capacity_yalin_differentiation(
    q,
    waterlevel,
    density,
    dm,
    slope,
    width,
    reservoirs,
    rivers,
    dtot,
    dt,
)
    sinslope = sin(atan(slope)) #slope in radians
    # Transport capacity from Yalin with particle differentiation
    # Delta parameter of Yalin for the specific particle class
    delta = waterlevel * sinslope / (1e-6 * (density / 1000 - 1)) / 0.06
    d_part = max(1 / dm * delta - 1, 0.0)

    if q > 0.0
        TCa = width / q * (density - 1000) * 1e-6 * (9.81 * waterlevel * sinslope)
    else
        TCa = 0.0
    end

    TCb = 2.45 / (0.001 * density)^0.4 * 0.06^0.5

    if dtot != 0.0 && d_part != 0.0
        TC =
            TCa * dm * d_part / dtot *
            0.635 *
            d_part *
            (1 - log(1 + d_part * TCb) / d_part * TCb) # [kg/m3]
        transport_capacity = TC * q * dt * 1e-3 #[ton]
    else
        transport_capacity = 0.0
    end

    # Mask transport capacity values for reservoirs and rivers
    transport_capacity = mask_transport_capacity(transport_capacity, reservoirs, rivers)

    return transport_capacity
end

"""
    function trasnport_capacity_bagnold(
        q,
        waterlevel,
        c_bagnold,
        e_bagnold,
        width,
        length,
        dt,
    )

Total sediment transport capacity based on Bagnold.

# Arguments
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `c_bagnold` (Bagnold transport capacity coefficient [-])
- `e_bagnold` (Bagnold transport capacity exponent [-])
- `width` (drain width [m])
- `length` (drain length [m])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1])
"""
function transport_capacity_bagnold(q, waterlevel, c_bagnold, e_bagnold, width, length, dt)
    # Transport capacity from Bagnold
    if waterlevel > 0.0
        # Transport capacity [tons/m3]
        transport_capacity = c_bagnold * (q / (waterlevel * width))^e_bagnold
        transport_capacity = limit_and_convert_transport_capacity(
            transport_capacity,
            q,
            waterlevel,
            width,
            length,
            dt,
        )
    else
        transport_capacity = 0.0
    end

    return transport_capacity
end

"""
    function trasnport_capacity_engelund(
        q,
        waterlevel,
        density,
        d50,
        width,
        length,
        slope,
        dt,
    )

Total sediment transport capacity based on Engelund and Hansen.

# Arguments
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m-3])
- `d50` (median grain size [m])
- `width` (drain width [m])
- `length` (drain length [m])
- `slope` (slope [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1])
"""
function transport_capacity_engelund(q, waterlevel, density, d50, width, length, slope, dt)
    # Transport capacity from Engelund and Hansen
    if waterlevel > 0.0
        # Hydraulic radius of the river [m] (rectangular channel)
        hydrad = waterlevel * width / (width + 2 * waterlevel)
        vshear = sqrt(9.81 * hydrad * slope)

        # Flow velocity [m/s]
        velocity = (q / (waterlevel * width))

        # Concentration by weight
        cw =
            density / 1000 * 0.05 * velocity * vshear^3 /
            ((density / 1000 - 1)^2 * 9.81^2 * d50 * hydrad)
        cw = min(1.0, cw)

        # Transport capacity [tons/m3]
        transport_capacity = cw / (cw + (1 - cw) * density / 1000) * density / 1000
        transport_capacity = max(transport_capacity, 0.0)
        # Transport capacity [tons]
        transport_capacity = limit_and_convert_transport_capacity(
            transport_capacity,
            q,
            waterlevel,
            width,
            length,
            dt,
        )

    else
        transport_capacity = 0.0
    end

    return transport_capacity
end

"""
    function trasnport_capacity_kodatie(
        q,
        waterlevel,
        a_kodatie,
        b_kodatie,
        c_kodatie,
        d_kodatie,
        width,
        length,
        slope,
        dt,
    )

Total sediment transport capacity based on Kodatie.

# Arguments
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `a_kodatie` (Kodatie transport capacity coefficient [-])
- `b_kodatie` (Kodatie transport capacity coefficient [-])
- `c_kodatie` (Kodatie transport capacity coefficient [-])
- `d_kodatie` (Kodatie transport capacity coefficient [-])
- `width` (drain width [m])
- `slope` (slope [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1])
"""
function transport_capacity_kodatie(
    q,
    waterlevel,
    a_kodatie,
    b_kodatie,
    c_kodatie,
    d_kodatie,
    width,
    length,
    slope,
    dt,
)
    # Transport capacity from Kodatie
    if waterlevel > 0.0
        # Flow velocity [m/s]
        velocity = (q / (waterlevel * width))

        # Concentration
        transport_capacity =
            a_kodatie * velocity^b_kodatie * waterlevel^c_kodatie * slope^d_kodatie

        # Transport capacity [tons/m3]
        transport_capacity *= width / (q * dt)
        transport_capacity = limit_and_convert_transport_capacity(
            transport_capacity,
            q,
            waterlevel,
            width,
            length,
            dt,
        )

    else
        transport_capacity = 0.0
    end

    return transport_capacity
end

"""
    function trasnport_capacity_yang(
        q,
        waterlevel,
        density,
        d50,
        width,
        length,
        slope,
        dt,
    )

Total sediment transport capacity based on Yang sand and gravel equations.

# Arguments
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m-3])
- `d50` (median grain size [m])
- `width` (drain width [m])
- `length` (drain length [m])
- `slope` (slope [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1])
"""
function transport_capacity_yang(q, waterlevel, density, d50, width, length, slope, dt)
    # Transport capacity from Yang
    omegas = 411 * d50^2 / 3600
    # Hydraulic radius of the river [m] (rectangular channel)
    hydrad = waterlevel * width / (width + 2 * waterlevel)
    # Critical shear stress velocity
    vshear = sqrt(9.81 * hydrad * slope)
    var1 = vshear * d50 / 1000 / (1.16 * 1e-6)
    var2 = omegas * d50 / 1000 / (1.16 * 1e-6)
    vcr = ifelse(var1 >= 70.0, 2.05 * omegas, omegas * (2.5 / (log10(var1) - 0.06) + 0.66))
    vcr = min(vcr, 0.0)

    # Sand equation
    if (width * waterlevel) > vcr && d50 < 2.0
        logcppm = (
            5.435 - 0.286 * log10(var2) - 0.457 * log10(vshear / omegas) + 1.799 -
            0.409 * log10(var2) -
            0.314 *
            log10(vshear / omegas) *
            log10((q / (width * waterlevel) - vcr) * slope / omegas)
        )
        # Gravel equation
    elseif (width * waterlevel) > vcr && d50 >= 2.0
        logcppm = (
            6.681 - 0.633 * log10(var2) - 4.816 * log10(vshear / omegas) + 2.784 -
            0.305 * log10(var2) -
            0.282 *
            log10(vshear / omegas) *
            log10((q / (width * waterlevel) - vcr) * slope / omegas)
        )
    else
        logcppm = 0.0
    end

    # Sediment concentration by weight
    cw = 10^logcppm * 1e-6
    # Transport capacity [tons/m3]
    transport_capacity = cw / (cw + (1 - cw) * density / 1000) * density / 1000
    transport_capacity = max(transport_capacity, 0.0)
    # Transport capacity [tons]
    transport_capacity = limit_and_convert_transport_capacity(
        transport_capacity,
        q,
        waterlevel,
        width,
        length,
        dt,
    )

    return transport_capacity
end

"""
    function trasnport_capacity_molinas(
        q,
        waterlevel,
        density,
        d50,
        width,
        length,
        slope,
        dt,
    )

Total sediment transport capacity based on Molinas and Wu.

# Arguments
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m-3])
- `d50` (median grain size [m])
- `width` (drain width [m])
- `length` (drain length [m])
- `slope` (slope [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1])
"""
function transport_capacity_molinas(q, waterlevel, density, d50, width, length, slope, dt)
    # Transport capacity from Molinas and Wu
    if waterlevel > 0.0
        # Flow velocity [m/s]
        velocity = (q / (waterlevel * width))
        omegas = 411 * d50^2 / 3600

        # PSI parameter
        psi = (
            velocity^3 / (
                (density / 1000 - 1) *
                9.81 *
                waterlevel *
                omegas *
                log10(1000 * waterlevel / d50)^2
            )
        )
        # Concentration by weight
        cw = 1430 * (0.86 + psi^0.5) * psi^1.5 / (0.016 + psi) * 1e-6
        # Transport capacity [tons/m3]
        transport_capacity = cw / (cw + (1 - cw) * density / 1000) * density / 1000
        transport_capacity = max(transport_capacity, 0.0)
        # Transport capacity [tons]
        transport_capacity = limit_and_convert_transport_capacity(
            transport_capacity,
            q,
            waterlevel,
            width,
            length,
            dt,
        )

    else
        transport_capacity = 0.0
    end

    return transport_capacity
end
