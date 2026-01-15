"""
    mask_transport_capacity(
        transport_capacity,
        reservoirs,
        rivers,
        dt
    )

Mask transport capacity values for reservoirs,
 and rivers.

# Arguments
- `transport_capacity` (total sediment transport capacity [t dt⁻¹ => kg s⁻¹])
- `reservoirs` (reservoirs mask [-])
- `rivers` (rivers mask [-])
- `dt` (timestep [s])

# Output
- `transport_capacity` (masked total sediment transport capacity [t dt-1 => kg s⁻¹])
"""
function mask_transport_capacity(transport_capacity, reservoirs, rivers, dt)
    # Full deposition in rivers to switch to river concept
    if rivers
        tc = 0.0
        # Sediment flux in reservoirs will all reach the river
    elseif reservoirs
        # 1e9 t dt⁻¹
        tc = to_SI(1e9, TON_PER_DT; dt_val = dt)
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
- `transport_capacity_density` (total sediment transport capacity [t m⁻³ => kg m⁻³])
- `q` (discharge [m³ s⁻¹])
- `waterlevel` (water level [m])
- `width` (drain width [m])
- `length` (drain length [m])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1 => kg s⁻¹])
"""
function limit_and_convert_transport_capacity(
    transport_capacity_density,
    q,
    waterlevel,
    width,
    length,
    dt,
)
    # 1285 g/L: boundary between streamflow and debris flow (Costa, 1988)
    # [kg m⁻³] = min([kg m⁻³], [kg m⁻³])
    transport_capacity_vol = min(transport_capacity_density, to_SI(1285, GRAM_PER_L))
    # Transport capacity
    # [kg s⁻¹] = [kg m⁻³]  * ([m] * [m] * [m] / [s] + [m³ s⁻¹])
    transport_capacity = transport_capacity_vol * (waterlevel * width * length / dt + q)

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
- `q` (discharge [m³ s⁻¹])
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
- `transport_capacity` (total sediment transport capacity [t dt⁻¹ => kg s⁻¹])
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
    sinslope = sin_slope(slope)
    if waterlevel > 0
        # [m s⁻¹] = [m³ s⁻¹] / ([m] * [m])
        velocity = q / (width * waterlevel)
        # [cm s⁻¹]
        omega = from_SI(10 * sinslope * velocity, CM_PER_S)
        if omega > 0.4
            # [kg m⁻³]
            TCf = c_govers * (omega - 0.4)^n_govers * density
            # [kg s⁻¹]
            transport_capacity = TCf * q
        else
            transport_capacity = 0.0
        end
    else
        transport_capacity = 0.0
    end
    # Mask transport capacity values for reservoirs and rivers
    transport_capacity = mask_transport_capacity(transport_capacity, reservoirs, rivers, dt)

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
- `q` (discharge [m³ s⁻¹])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m⁻³])
- `d50` (median grain size [mm => m])
- `slope` (slope [-])
- `width` (drain width [m])
- `reservoirs` (reservoirs mask [-])
- `rivers` (rivers mask [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt⁻¹ => kg s⁻¹])
"""
function transport_capacity_yalin(
    q::Float64,
    waterlevel::Float64,
    density::Float64,
    d50::Float64,
    slope::Float64,
    width::Float64,
    reservoirs::Bool,
    rivers::Bool,
    dt::Float64,
)
    sinslope = sin_slope(slope) #slope in radians
    # Transport capacity from Yalin without particle differentiation
    delta =
        max((waterlevel * sinslope / (d50 * (density / WATER_DENSITY - 1)) / 0.06 - 1), 0.0)
    alphay = delta * 2.45 / (1e-3 * density)^(2 // 5) * sqrt(0.06)
    if q > 0.0 && alphay != 0.0
        # [kg m⁻³]
        TC = (
            width / q *
            (density - 1000) *
            d50 *
            (GRAVITATIONAL_ACCELERATION * waterlevel * sinslope) *
            0.635 *
            delta *
            (1 - log(1 + alphay) / (alphay))
        )
        # [kg s⁻¹] = [kg m⁻³] * [m³ s⁻¹]
        transport_capacity = TC * q
    else
        transport_capacity = 0.0
    end

    # Mask transport capacity values for reservoirs and rivers
    transport_capacity = mask_transport_capacity(transport_capacity, reservoirs, rivers, dt)

    return transport_capacity
end

const WATER_DENSITY = 1e3 # [kg m⁻³]
const WATER_KINEMATIC_VISCOSITY = 1.16e-6 # [m² s⁻¹]

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
- `density` (sediment density [kg m⁻³])
- `dm_clay` (clay median grain size [μm => m])
- `dm_silt` (silt median grain size [μm => m])
- `dm_sand` (sand median grain size [μm => m])
- `dm_sagg` (small aggregates median grain size [μm => m])
- `dm_lagg` (large aggregates median grain size [μm => m])
- `slope` (slope [-])

# Output
- `dtot` (total transportability of the flow [-])
"""
function transportability_yalin_differentiation(
    waterlevel::Float64,
    density::Float64,
    dm_clay::Float64,
    dm_silt::Float64,
    dm_sand::Float64,
    dm_sagg::Float64,
    dm_lagg::Float64,
    slope::Float64,
)
    sinslope = sin_slope(slope) #slope in radians
    # Delta parameter of Yalin for each particle class
    delta = waterlevel * sinslope / (density / WATER_DENSITY - 1) / 0.06
    dclay = max(delta / dm_clay - 1, 0.0)
    dsilt = max(delta / dm_silt - 1, 0.0)
    dsand = max(delta / dm_sand - 1, 0.0)
    dsagg = max(delta / dm_sagg - 1, 0.0)
    dlagg = max(delta / dm_lagg - 1, 0.0)
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
- `q` (discharge [m³ s⁻¹])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m⁻³])
- `dm` (median grain size [m])
- `slope` (slope [-])
- `width` (drain width [m])
- `reservoirs` (reservoirs mask [-])
- `rivers` (rivers mask [-])
- `dtot` (total flow transportability [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt⁻¹])
"""
function transport_capacity_yalin_differentiation(
    q::Float64,
    waterlevel::Float64,
    density::Float64,
    dm::Float64,
    slope::Float64,
    width::Float64,
    reservoirs::Bool,
    rivers::Bool,
    dtot::Float64,
    dt::Float64,
)
    sinslope = sin_slope(slope) #slope in radians
    # Transport capacity from Yalin with particle differentiation
    # Delta parameter of Yalin for the specific particle class
    delta = waterlevel * sinslope / (1e-6 * (density / WATER_DENSITY - 1)) / 0.06
    d_part = max(delta / dm - 1, 0.0)

    if q > 0.0
        TCa =
            width / q *
            (density - 1000) *
            1e-6 *
            (GRAVITATIONAL_ACCELERATION * waterlevel * sinslope)
    else
        TCa = 0.0
    end

    TCb = 2.45 * sqrt(0.06) / density^(2 // 5)

    if dtot != 0.0 && d_part != 0.0
        TC =
            TCa * dm * d_part / dtot *
            0.635 *
            d_part *
            (1 - log(1 + d_part * TCb) / d_part * TCb) # [kg/m3]
        transport_capacity = TC * q * dt
    else
        transport_capacity = 0.0
    end

    # Mask transport capacity values for reservoirs and rivers
    transport_capacity = mask_transport_capacity(transport_capacity, reservoirs, rivers, dt)

    return transport_capacity
end

"""
    function transport_capacity_bagnold(
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
- `q` (discharge [m³ s⁻¹])
- `waterlevel` (water level [m])
- `c_bagnold` (Bagnold transport capacity coefficient [-])
- `e_bagnold` (Bagnold transport capacity exponent [-])
- `width` (drain width [m])
- `length` (drain length [m])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt⁻¹ => kg s⁻¹])
"""
function transport_capacity_bagnold(
    q::Float64,
    waterlevel::Float64,
    c_bagnold::Float64,
    e_bagnold::Float64,
    width::Float64,
    length::Float64,
    dt::Float64,
)
    # Transport capacity from Bagnold
    if waterlevel > 0.0
        # Transport capacity [kg s⁻¹]
        transport_capacity =
            to_SI(c_bagnold * (q / (waterlevel * width))^e_bagnold, TON_PER_M3)
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
    function transport_capacity_engelund(
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
- `q` (discharge [m³ s⁻¹])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m⁻³])
- `d50` (median grain size [m])
- `width` (drain width [m])
- `length` (drain length [m])
- `slope` (slope [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt⁻¹ => kg s⁻¹])
"""
function transport_capacity_engelund(
    q::Float64,
    waterlevel::Float64,
    density::Float64,
    d50::Float64,
    width::Float64,
    length::Float64,
    slope::Float64,
    dt::Float64,
)
    # Transport capacity from Engelund and Hansen
    if waterlevel > 0.0
        # Hydraulic radius of the river [m] (rectangular channel)
        # [m] = [m] * [m] / ([m] + [-] * [m])
        hydrad = waterlevel * width / (width + 2 * waterlevel)
        # [m s⁻¹] = sqrt([m s⁻²] * [m] * [-])
        vshear = sqrt(GRAVITATIONAL_ACCELERATION * hydrad * slope)

        # Flow velocity
        # [m s⁻¹] = [m³ s⁻¹] / ([m] * [m])
        velocity = q / (waterlevel * width)

        # Concentration by weight
        # [-] = [kg m⁻³] / [kg m⁻³]
        cw_ = density / WATER_DENSITY
        # [-] = [-] * [m s⁻¹] * [m s⁻¹]³ / (([-] - [-])^2 * [m s⁻²]^2 * [m] * [m])
        cw = 0.05 * velocity * vshear^3 / ((cw_ - 1)^2 * GRAVITATIONAL_ACCELERATION^2)
        cw = min(1.0, cw)

        # Transport capacity [kg m⁻³]
        # [kg m⁻³] = [-] / ([-] + ([-] - [-]) * [-]) * [kg m⁻³]
        transport_capacity_density = cw / (cw + (1 - cw) * cw_) * density
        transport_capacity_density = max(transport_capacity_density, 0.0)
        # Transport capacity [kg s⁻¹]
        transport_capacity = limit_and_convert_transport_capacity(
            transport_capacity_density,
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
- `q` (discharge [m³ s⁻¹])
- `waterlevel` (water level [m])
- `a_kodatie` (Kodatie transport capacity coefficient [-])
- `b_kodatie` (Kodatie transport capacity coefficient [-])
- `c_kodatie` (Kodatie transport capacity coefficient [-])
- `d_kodatie` (Kodatie transport capacity coefficient [-])
- `width` (drain width [m])
- `slope` (slope [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1 => kg s⁻¹])
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
        # [m s⁻¹] = [m³ s⁻¹] / ([m] * [m])
        velocity = q / (waterlevel * width)

        # Concentration
        # [kg m⁻¹]
        transport_capacity_concentration =
            a_kodatie * velocity^b_kodatie * waterlevel^c_kodatie * slope^d_kodatie

        # Transport capacity
        # [kg m⁻³] = [kg m⁻¹] * [m] / ([m³ s⁻¹] * [s])
        transport_capacity_density = transport_capacity_concentration * width / (q * dt)
        # [kg s⁻¹]
        transport_capacity = limit_and_convert_transport_capacity(
            transport_capacity_density,
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
- `q` (discharge [m³ s⁻¹])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m⁻³])
- `d50` (median grain size [mm => m])
- `width` (drain width [m])
- `length` (drain length [m])
- `slope` (slope [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt-1 => kg s⁻¹])
"""
function transport_capacity_yang(
    q::Float64,
    waterlevel::Float64,
    density::Float64,
    d50::Float64,
    width::Float64,
    length::Float64,
    slope::Float64,
    dt::Float64,
)
    # Transport capacity from Yang
    # [m s⁻¹]
    omegas = fall_velocity(d50)
    # Hydraulic radius of the river [m] (rectangular channel)
    # [m] = [m] * [m] / ([m] + [-] * [m])
    hydrad = waterlevel * width / (width + 2 * waterlevel)
    # Critical shear stress velocity
    # [m s⁻¹] = sqrt([m s⁻²] * [m] * [-])
    vshear = sqrt(g_gravity * hydrad * slope)
    # [-] = [m s⁻¹] * [m] / [m² s⁻¹]
    var1 = vshear * d50 / WATER_KINEMATIC_VISCOSITY
    # [-] = [m s⁻¹] * [m] / [m² s⁻¹]
    var2 = omegas * d50 / WATER_KINEMATIC_VISCOSITY
    vcr = ifelse(var1 >= 70.0, 2.05 * omegas, omegas * (2.5 / (log10(var1) - 0.06) + 0.66))
    vcr = min(vcr, 0.0)

    if width * waterlevel > vcr
        if d50 < 2.0
            # Sand equation
            logcppm = (
                5.435 - 0.286 * log10(var2) - 0.457 * log10(vshear / omegas) + 1.799 -
                0.409 * log10(var2) -
                0.314 *
                log10(vshear / omegas) *
                log10((q / (width * waterlevel) - vcr) * slope / omegas)
            )
        else
            # Gravel equation
            logcppm = (
                6.681 - 0.633 * log10(var2) - 4.816 * log10(vshear / omegas) + 2.784 -
                0.305 * log10(var2) -
                0.282 *
                log10(vshear / omegas) *
                log10((q / (width * waterlevel) - vcr) * slope / omegas)
            )
        end
    else
        logcppm = 0.0
    end

    # Sediment concentration by weight from parts per million
    # [-]
    cw = to_SI(10^logcppm, PPM)
    # Transport capacity
    # [kg m⁻³] = [-] / ([-] + ([-] - [-]) * [kg m⁻³] / [kg m⁻³]) * [kg m⁻³]
    transport_capacity_density = cw / (cw + (1 - cw) * density / WATER_DENSITY) * density
    transport_capacity_density = max(transport_capacity_density, 0.0)
    # [kg s⁻¹]
    transport_capacity = limit_and_convert_transport_capacity(
        transport_capacity_density,
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
- `q` (discharge [m³ s-1])
- `waterlevel` (water level [m])
- `density` (sediment density [kg m⁻³])
- `d50` (median grain size [m])
- `width` (drain width [m])
- `length` (drain length [m])
- `slope` (slope [-])
- `dt` (time step [s])

# Output
- `transport_capacity` (total sediment transport capacity [t dt⁻¹ => kg s⁻¹])
"""
function transport_capacity_molinas(q, waterlevel, density, d50, width, length, slope, dt)
    # Transport capacity from Molinas and Wu
    if waterlevel > 0.0
        # Flow velocity
        # [m s⁻¹] = ([m³ s⁻¹] / ([m] * [m]))
        velocity = (q / (waterlevel * width))
        # [m s⁻¹]
        omegas = fall_velocity(d50)

        # PSI parameter
        psi = (
            velocity^3 / (
                (density / WATER_DENSITY - 1) *
                GRAVITATIONAL_ACCELERATION *
                waterlevel *
                omegas *
                log10(waterlevel / d50)^2
            )
        )
        # Concentration by weight from parts per million
        # [-]
        cw = to_SI(1430 * (0.86 + sqrt(psi)) * psi^(3 // 2) / (0.016 + psi), PPM)
        # Transport capacity
        # [kg m⁻³] = [-] / ([-] + ([-] - [-]) * [kg m⁻³] / [kg m⁻³]) * [kg m⁻³]
        transport_capacity_density =
            cw / (cw + (1 - cw) * density / WATER_DENSITY) * density
        transport_capacity_density = max(transport_capacity_density, 0.0)
        # [kg s⁻¹]
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
