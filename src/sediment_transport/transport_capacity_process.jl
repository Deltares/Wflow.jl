"""
    mask_transport_capacity(
        transport_capacity,
        waterbodies,
        rivers,
    )

Mask transport capacity values for waterbodies and rivers.

# Arguments
- `transport_capacity` (total sediment transport capacity [t dt-1])
- `waterbodies` (waterbodies mask [-])
- `rivers` (rivers mask [-])

# Output
- `transport_capacity` (masked total sediment transport capacity [t dt-1])
"""
function mask_transport_capacity(transport_capacity, waterbodies, rivers)
    # Full deposition in rivers to switch to river concept
    if rivers
        tc = 0.0
        # Sediment flux in waterbodies will all reach the river
    elseif waterbodies
        tc = 1e9
    else
        tc = transport_capacity
    end

    return tc
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
        waterbodies,
        rivers,
        ts,
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
- `waterbodies` (waterbodies mask [-])
- `rivers` (rivers mask [-])
- `ts` (time step [s])

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
    waterbodies,
    rivers,
    ts,
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
    transport_capacity = TCf * q * ts * 1e-3 #[ton/cell]

    # Mask transport capacity values for waterbodies and rivers
    transport_capacity = mask_transport_capacity(transport_capacity, waterbodies, rivers)

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
        waterbodies,
        rivers,
        ts,
    )

Total sediment transport capacity based on Yalin.

# Arguments
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m-3])
- `d50` (median grain size [m])
- `slope` (slope [-])
- `width` (drain width [m])
- `waterbodies` (waterbodies mask [-])
- `rivers` (rivers mask [-])
- `ts` (time step [s])

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
    waterbodies,
    rivers,
    ts,
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
        transport_capacity = TC * q * ts * 1e-3 #[ton]
    else
        transport_capacity = 0.0
    end

    # Mask transport capacity values for waterbodies and rivers
    transport_capacity = mask_transport_capacity(transport_capacity, waterbodies, rivers)

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
        waterbodies,
        rivers,
        dtot,
        ts,
    )

Transport capacity for a specific grain size based on Yalin with particle differentiation.

# Arguments
- `q` (discharge [m3 s-1])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m-3])
- `dm` (median grain size [m])
- `slope` (slope [-])
- `width` (drain width [m])
- `waterbodies` (waterbodies mask [-])
- `rivers` (rivers mask [-])
- `dtot` (total flow transportability [t dt-1])
- `ts` (time step [s])

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
    waterbodies,
    rivers,
    dtot,
    ts,
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
        transport_capacity = TC * q * ts * 1e-3 #[ton]
    else
        transport_capacity = 0.0
    end

    # Mask transport capacity values for waterbodies and rivers
    transport_capacity = mask_transport_capacity(transport_capacity, waterbodies, rivers)

    return transport_capacity
end