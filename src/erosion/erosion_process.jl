"""
    rainfall_erosion_eurosem(
        precip,
        interception,
        waterlevel,
        soil_detachability,
        eurosem_exponent,
        canopyheight,
        canopygapfraction,
        soilcover_fraction,
        area,
        dt,
    )

Rainfall erosion model based on EUROSEM.

# Arguments
- `precip` (precipitation [mm Δt⁻¹])
- `interception` (interception [mm Δt⁻¹])
- `waterlevel` (water level [m])
- `soil_detachability` (soil detachability [-])
- `eurosem_exponent` (EUROSEM exponent [-])
- `canopyheight` (canopy height [m])
- `canopygapfraction` (canopy gap fraction [-])
- `soilcover_fraction` (soil cover fraction [-])
- `area` (area [m2])
- `dt` (timestep [seconds])

# Output
- `rainfall_erosion` (soil loss [t Δt⁻¹])
"""
function rainfall_erosion_eurosem(
    precip,
    interception,
    waterlevel,
    soil_detachability,
    eurosem_exponent,
    canopyheight,
    canopygapfraction,
    soilcover_fraction,
    area,
    dt,
)
    # calculate rainfall intensity [mm/h]
    rintnsty = precip / (dt / 3600)
    # Kinetic energy of direct throughfall [J/m2/mm]
    # kedir = max(11.87 + 8.73 * log10(max(0.0001, rintnsty)),0.0) #basis used in USLE
    kedir = max(8.95 + 8.44 * log10(max(0.0001, rintnsty)), 0.0) #variant used in most distributed mdoels
    # Kinetic energy of leaf drainage [J/m2/mm]
    pheff = 0.5 * canopyheight
    keleaf = max((15.8 * pheff^0.5) - 5.87, 0.0)

    #Depths of rainfall (total, leaf drianage, direct) [mm]
    rdtot = precip
    rdleaf = rdtot * 0.1 * canopygapfraction #stemflow
    rddir = max(rdtot - rdleaf - interception, 0.0) #throughfall

    #Total kinetic energy by rainfall [J/m2]
    ketot = (rddir * kedir + rdleaf * keleaf) * 0.001
    # Rainfall / splash erosion [g/m2]
    sedspl = soil_detachability * ketot * exp(-eurosem_exponent * waterlevel)
    sedspl = sedspl * area * 1e-6 # ton/cell

    # Remove the impervious area
    sedspl = sedspl * (1.0 - soilcover_fraction)
    return rainfall_erosion
end

"""
    rainfall_erosion_answers(
        precip,
        usle_k,
        usle_c,
        area,
        dt,
    )

Rainfall erosion model based on ANSWERS.

# Arguments
- `precip` (precipitation [mm Δt⁻¹])
- `usle_k` (USLE soil erodibility [t ha-1 mm-1])
- `usle_c` (USLE cover and management factor [-])
- `soilcover_fraction` (soil cover fraction [-])
- `area` (area [m2])
- `dt` (timestep [seconds])

# Output
- `rainfall_erosion` (soil loss [t Δt⁻¹])
"""
function rainfall_erosion_answers(precip, usle_k, usle_c, area, dt)
    # calculate rainfall intensity [mm/min]
    rintnsty = precip / (ts / 60)
    # splash erosion [kg/min]
    sedspl = 0.108 * usle_c * usle_k * area * rintnsty^2
    # [ton/timestep]
    sedspl = sedspl * (dt / 60) * 1e-3
    return rainfall_erosion
end

"""
    overland_flow_erosion_answers(
        overland_flow,
        waterlevel,
        usle_k,
        usle_c,
        answers_k,
        slope,
        soilcover_fraction,
        area,
        dt,
    )

Overland flow erosion model based on ANSWERS.

# Arguments
- `overland_flow` (overland flow [m3 s-1])
- `waterlevel` (water level [m])
- `usle_k` (USLE soil erodibility [t ha-1 mm-1])
- `usle_c` (USLE cover and management factor [-])
- `answers_k` (ANSWERS overland flow factor [-])
- `slope` (slope [-])
- `area` (area [m2])
- `dt` (timestep [seconds])

# Output
- `overland_flow_erosion` (soil loss [t Δt⁻¹])
"""
function overland_flow_erosion_answers(
    overland_flow,
    usle_k,
    usle_c,
    answers_k,
    slope,
    area,
    dt,
)
    # Overland flow rate [m2/min]
    qr_land = overland_flow * 60 / ((area) / 2)
    # Sine of the slope
    sinslope = sin(atan(slope))

    # Overland flow erosion [kg/min]
    # For a wide range of slope, it is better to use the sine of slope rather than tangeant
    erosion = answers_k * usle_c * usle_k * area * sinslope * qr_land
    # [ton/timestep]
    erosion = erosion * (dt / 60) * 1e-3
    return erosion
end

"""
    total_soil_erosion(
        rainfall_erosion,
        overland_flow_erosion,
        clay_fraction,
        silt_fraction,
        sand_fraction,
        sagg_fraction,
        lagg_fraction,
    )

Calculate total soil erosion and particle differentiation.

# Arguments
- `rainfall_erosion` (soil loss from rainfall erosion [t Δt⁻¹])
- `overland_flow_erosion` (soil loss from overland flow erosion [t Δt⁻¹])
- `clay_fraction` (clay fraction [-])
- `silt_fraction` (silt fraction [-])
- `sand_fraction` (sand fraction [-])
- `sagg_fraction` (small aggregates fraction [-])
- `lagg_fraction` (large aggregates fraction [-])

# Output
- `soil_erosion` (total soil loss [t Δt⁻¹])
- `clay_erosion` (clay loss [t Δt⁻¹])
- `silt_erosion` (silt loss [t Δt⁻¹])
- `sand_erosion` (sand loss [t Δt⁻¹])
- `sagg_erosion` (small aggregates loss [t Δt⁻¹])
- `lagg_erosion` (large aggregates loss [t Δt⁻¹])
"""
function total_soil_erosion(
    rainfall_erosion,
    overland_flow_erosion,
    clay_fraction,
    silt_fraction,
    sand_fraction,
    sagg_fraction,
    lagg_fraction,
)
    # Total soil erosion
    soil_erosion = rainfall_erosion + overland_flow_erosion
    # Particle differentiation
    clay_erosion = soil_erosion * clay_fraction
    silt_erosion = soil_erosion * silt_fraction
    sand_erosion = soil_erosion * sand_fraction
    sagg_erosion = soil_erosion * sagg_fraction
    lagg_erosion = soil_erosion * lagg_fraction
    return soil_erosion,
    clay_erosion,
    silt_erosion,
    sand_erosion,
    sagg_erosion,
    lagg_erosion
end