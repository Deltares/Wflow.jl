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
- `precip` (precipitation [mm dt⁻¹ => m s⁻¹])
- `interception` (interception [mm dt⁻¹ => m s⁻¹])
- `waterlevel` (water level [m])
- `soil_detachability` (soil detachability [g J⁻¹ => kg J⁻¹])
- `eurosem_exponent` (EUROSEM exponent [-])
- `canopyheight` (canopy height [m])
- `canopygapfraction` (canopy gap fraction [-])
- `soilcover_fraction` (soil cover fraction [-])
- `area` (area [m²])
- `dt` (timestep [s])

# Output
- `rainfall_erosion` (soil loss [t dt⁻¹ => kg s⁻¹])
"""
function rainfall_erosion_eurosem(
    precip::Float64,
    interception::Float64,
    waterlevel::Float64,
    soil_detachability::Float64,
    eurosem_exponent::Float64,
    canopyheight::Float64,
    canopygapfraction::Float64,
    soilcover_fraction::Float64,
    area::Float64,
    dt::Float64,
)
    # Precipitation expressed in unit expected by model
    rainfall_intensity = from_SI(precip, MM_PER_HOUR)
    # Kinetic energy of direct throughfall
    # E_kin_direct = max(11.87 + 8.73 * log10(max(0.0001, rainfall_intensity)),0.0) #basis used in USLE
    # [J m⁻² mm⁻¹]
    E_kin_direct = max(8.95 + 8.44 * log10(max(0.0001, rainfall_intensity)), 0.0) #variant used in most distributed models
    # Kinetic energy of leaf drainage
    pheff = 0.5 * canopyheight
    # [J m⁻² mm⁻¹]
    E_kin_leaf = max((15.8 * sqrt(pheff)) - 5.87, 0.0)

    # Depths of rainfall (total, leaf drainage, direct)
    # [mm] = [mm h⁻¹] * [h]
    rainfall_depth_total = rainfall_intensity * from_SI(dt, HOUR)
    # [mm] = [mm] * [-] * [-]
    rainfall_depth_leaf = rainfall_depth_total * 0.1 * canopygapfraction # stemflow
    # [mm]
    intercepted = from_SI(interception * dt, MM)
    # [mm]
    rainfall_depth_direct =
        max(rainfall_depth_total - rainfall_depth_leaf - intercepted, 0.0) # throughfall

    # Total kinetic energy by rainfall
    # [J m⁻²] = [mm] * [J m⁻² mm⁻¹] + [mm] * [J m⁻² mm⁻¹]
    E_kin_tot = rainfall_depth_direct * E_kin_direct + rainfall_depth_leaf * E_kin_leaf
    # [kg s⁻¹] = [m²] * [kg J⁻¹] * [J m⁻²] * exp([m⁻¹] * [m]) / [s]
    rainfall_erosion =
        area * soil_detachability * E_kin_tot * exp(-eurosem_exponent * waterlevel) / dt

    # Remove the impervious area
    # [kg s⁻¹] = [kg s⁻¹] * [-]
    rainfall_erosion *= 1.0 - soilcover_fraction
    return rainfall_erosion
end

"""
    rainfall_erosion_answers(
        precip,
        usle_k,
        usle_c,
        answers_rainfall_factor,
        area,
        dt,
    )

Rainfall erosion model based on ANSWERS.

# Arguments
- `precip` (precipitation [mm dt⁻¹ => m s⁻¹])
- `usle_k` (USLE soil erodibility [-], treated as unitless but is actually [t ha-1 h-1 MJ-1 mm-1])
- `usle_c` (USLE cover and management factor [-])
- `answers_rainfall_factor` (ANSWERS rainfall erosion factor [-], treated as unitless but could have units)
- `area` (area [m²])
- `dt` (timestep [s])

# Output
- `rainfall_erosion` (soil loss [t dt⁻¹ => kg s⁻¹])
"""
function rainfall_erosion_answers(
    precip::Float64,
    usle_k::Float64,
    usle_c::Float64,
    answers_rainfall_factor::Float64,
    area::Float64,
)
    # calculate rainfall intensity
    rainfall_intensity = from_SI(precip, MM_PER_MIN)
    # splash erosion [kg min⁻¹]
    # The units here are hard to track because these equations are derived
    # empirically
    rainfall_erosion =
        answers_rainfall_factor * usle_c * usle_k * area * rainfall_intensity^2
    return to_SI(rainfall_erosion, KG_PER_MIN)
end

"""
    overland_flow_erosion_answers(
        overland_flow,
        waterlevel,
        usle_k,
        usle_c,
        answers_overland_flow_factor,
        slope,
        soilcover_fraction,
        area,
        dt,
    )

Overland flow erosion model based on ANSWERS.

# Arguments
- `overland_flow` (overland flow [m³ s⁻¹])
- `usle_k` (USLE soil erodibility [-], treated as unitless but is actually [t ha-1 h-1 MJ-1 mm-1])
- `usle_c` (USLE cover and management factor [-])
- `answers_overland_flow_factor` (ANSWERS overland flow factor [-], treated as unitless but could have units)
- `slope` (slope [-])
- `area` (area [m²])
- `dt` (timestep [s])

# Output
- `overland_flow_erosion` (soil loss [t dt⁻¹ => kg s⁻¹])
"""
function overland_flow_erosion_answers(
    overland_flow::Float64,
    usle_k::Float64,
    usle_c::Float64,
    answers_overland_flow_factor::Float64,
    slope::Float64,
    area::Float64,
)
    # Overland flow rate [m² min⁻¹]
    qr_land = from_SI(overland_flow, M3_PER_MIN) / sqrt(area)
    sinslope = sin_slope(slope)

    # Overland flow erosion [kg min⁻¹]
    # For a wide range of slope, it is better to use the sine of slope rather than tangeant
    erosion = answers_overland_flow_factor * usle_c * usle_k * area * sinslope * qr_land
    return to_SI(erosion, KG_PER_MIN)
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
- `rainfall_erosion` (soil loss from rainfall erosion [t dt⁻¹ => kg s⁻¹])
- `overland_flow_erosion` (soil loss from overland flow erosion [t dt⁻¹ => kg s⁻¹])
- `clay_fraction` (clay fraction [-])
- `silt_fraction` (silt fraction [-])
- `sand_fraction` (sand fraction [-])
- `sagg_fraction` (small aggregates fraction [-])
- `lagg_fraction` (large aggregates fraction [-])

# Output
- `soil_erosion` (total soil loss [t dt⁻¹ => kg s⁻¹])
- `clay_erosion` (clay loss [t dt⁻¹ => kg s⁻¹])
- `silt_erosion` (silt loss [t dt⁻¹ => kg s⁻¹])
- `sand_erosion` (sand loss [t dt⁻¹ => kg s⁻¹])
- `sagg_erosion` (small aggregates loss [t dt⁻¹ => kg s⁻¹])
- `lagg_erosion` (large aggregates loss [t dt⁻¹ => kg s⁻¹])
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
    # [kg s⁻¹] = [kg s⁻¹] + [kg s⁻¹]
    soil_erosion = rainfall_erosion + overland_flow_erosion
    # Particle differentiation
    # [kg s⁻¹] = # [kg s⁻¹] * [-]
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

"""
    function river_erosion_julian_torres(
        waterlevel,
        d50,
        width,
        length,
        slope,
        dt,
    )

River erosion model based on Julian Torres.
Repartition of the effective shear stress between the bank and the bed from Knight et al. 1984 [%]

# Arguments
- `waterlevel` (water level [m])
- `d50` (median grain size [m])
- `width` (width [m])
- `length` (length [m])
- `slope` (slope [-])
- `dt` (timestep [seconds])

# Output
- `bed` (potential river erosion [t dt⁻¹ => kg s⁻¹])
- `bank` (potential bank erosion [t dt⁻¹ => kg s⁻¹])
"""
function river_erosion_julian_torres(
    waterlevel::Float64,
    d50::Float64,
    width::Float64,
    length::Float64,
    slope::Float64,
    dt::Float64,
)
    if waterlevel > 0.0
        # Bed and Bank from Shields diagram, Da Silva & Yalin (2017)
        E_ = (2.65 - 1) * GRAVITATIONAL_ACCELERATION
        E = 10 * d50 * cbrt(E_)
        TCrbed =
            E_ *
            d50 *
            (0.13 * E^(-0.392) * exp(-0.015 * E^2) + 0.045 * (1 - exp(-0.068 * E)))
        TCrbank = TCrbed
        # kd from Hanson & Simon 2001
        kdbank = 0.2 * inv(sqrt(TCrbank)) * 1e-6
        kdbed = 0.2 * inv(sqrt(TCrbed)) * 1e-6

        # Hydraulic radius of the river [m] (rectangular channel)
        hydrad = waterlevel * width / (width + 2 * waterlevel)

        # Repartition of the effective shear stress between the bank and the Bed
        SFbank = exp(-3.23 * log10(width / waterlevel + 3) + 6.146)
        # Effective shear stress on river bed and banks [N/m2]
        TEffbank =
            1000 * GRAVITATIONAL_ACCELERATION * hydrad * slope * SFbank / 100 *
            (1 + width / (2 * waterlevel))
        TEffbed =
            1000 *
            GRAVITATIONAL_ACCELERATION *
            hydrad *
            slope *
            (1 - SFbank / 100) *
            (1 + 2 * waterlevel / width)

        # Potential erosion rates of the bed and bank [t/cell/timestep]
        #(assuming only one bank is eroding)
        Tex = max(TEffbank - TCrbank, 0.0)
        # 1.4 is bank default bulk density
        ERbank = kdbank * Tex * length * waterlevel * 1.4 * dt
        # 1.5 is bed default bulk density
        ERbed = kdbed * (TEffbed - TCrbed) * length * width * 1.5 * dt

        # Potential maximum bed/bank erosion
        bed = max(ERbed, 0.0)
        bank = max(ERbank, 0.0)

    else
        bed = 0.0
        bank = 0.0
    end

    return to_SI(bed, TON_PER_DT; dt_val = dt), to_SI(bank, TON_PER_DT; dt_val = dt)
end

"""
    function river_erosion_store(
        excess_sediment,
        store,
    )

River erosion of the previously deposited sediment.

# Arguments
- `excess_sediment` (excess sediment [tdt⁻¹ => kg s⁻¹])
- `store` (sediment store [t => kg])
- `dt` (timestep [s])

# Output
- `erosion` (river erosion [t dt⁻¹ => kg s⁻¹])
- `excess_sediment` (updated excess sediment [tdt⁻¹ => kg s⁻¹])
- `store` (updated sediment store [t => kg])
"""
function river_erosion_store(excess_sediment, store, dt)
    # River erosion of the previously deposited sediment

    # [kg s⁻¹] = min([kg] / [s], [kg s⁻¹])
    erosion = min(store / dt, excess_sediment)
    # Update the excess sediment and the sediment store
    # [kg s⁻¹] -= [kg s⁻¹]
    excess_sediment -= erosion
    # [kg] -= [kg s⁻¹] * [s]
    store -= erosion * dt
    return erosion, excess_sediment, store
end
