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
    rainfall_intensity = from_SI(precip, Unit(; mm = 1, h = -1))
    # Kinetic energy of direct throughfall [J m⁻² mm⁻¹]
    # E_kin_direct = max(11.87 + 8.73 * log10(max(0.0001, rainfall_intensity)),0.0) #basis used in USLE
    E_kin_direct = max(8.95 + 8.44 * log10(max(0.0001, rainfall_intensity)), 0.0) #variant used in most distributed models
    # Kinetic energy of leaf drainage [J m⁻² mm⁻¹]
    pheff = 0.5 * canopyheight
    E_kin_leaf = max((15.8 * sqrt(pheff)) - 5.87, 0.0)

    # Depths of rainfall (total, leaf drainage, direct) [mm]
    rainfall_depth_total = rainfall_intensity * from_SI(dt, Unit(; h = 1))
    rainfall_depth_leaf = rainfall_depth_total * 0.1 * canopygapfraction # stemflow
    intercepted = from_SI(interception * dt, Unit(; mm = 1))
    rainfall_depth_direct =
        max(rainfall_depth_total - rainfall_depth_leaf - intercepted, 0.0) # throughfall

    # Total kinetic energy by rainfall [J m⁻²]
    E_kin_tot = rainfall_depth_direct * E_kin_direct + rainfall_depth_leaf * E_kin_leaf
    # Rainfall / splash erosion [g m⁻²]
    # [kg] = [m²] * [kg J⁻¹] * [J m⁻²] * exp([m⁻¹] * [m])
    rainfall_erosion =
        area * soil_detachability * E_kin_tot * exp(-eurosem_exponent * waterlevel)

    # Remove the impervious area
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
- `precip` (precipitation [mm dt⁻¹])
- `usle_k` (USLE soil erodibility [t ha-1 mm-1])
- `usle_c` (USLE cover and management factor [-])
- `answers_rainfall_factor` (ANSWERS rainfall erosion factor [-])
- `area` (area [m2])
- `dt` (timestep [seconds])

# Output
- `rainfall_erosion` (soil loss [tdt⁻¹])
"""
function rainfall_erosion_answers(
    precip::Float64,
    usle_k::Float64,
    usle_c::Float64,
    answers_rainfall_factor::Float64,
    area::Float64,
    dt::Float64,
)
    # calculate rainfall intensity [mm/min]
    rintnsty = precip / (dt / 60)
    # splash erosion [kg/min]
    rainfall_erosion = answers_rainfall_factor * usle_c * usle_k * area * rintnsty^2
    # [ton/timestep]
    rainfall_erosion = rainfall_erosion * (dt / 60) * 1e-3
    return rainfall_erosion
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
- `overland_flow` (overland flow [m3 s-1])
- `waterlevel` (water level [m])
- `usle_k` (USLE soil erodibility [t ha-1 mm-1])
- `usle_c` (USLE cover and management factor [-])
- `answers_overland_flow_factor` (ANSWERS overland flow factor [-])
- `slope` (slope [-])
- `area` (area [m2])
- `dt` (timestep [seconds])

# Output
- `overland_flow_erosion` (soil loss [tdt⁻¹])
"""
function overland_flow_erosion_answers(
    overland_flow,
    usle_k,
    usle_c,
    answers_overland_flow_factor,
    slope,
    area,
    dt,
)
    # Overland flow rate [m2/min]
    qr_land = overland_flow * 60 / sqrt.(area)
    # Sine of the slope
    sinslope = slope / sqrt(1 + slope^2)

    # Overland flow erosion [kg/min]
    # For a wide range of slope, it is better to use the sine of slope rather than tangeant
    erosion = answers_overland_flow_factor * usle_c * usle_k * area * sinslope * qr_land
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
- `rainfall_erosion` (soil loss from rainfall erosion [tdt⁻¹])
- `overland_flow_erosion` (soil loss from overland flow erosion [tdt⁻¹])
- `clay_fraction` (clay fraction [-])
- `silt_fraction` (silt fraction [-])
- `sand_fraction` (sand fraction [-])
- `sagg_fraction` (small aggregates fraction [-])
- `lagg_fraction` (large aggregates fraction [-])

# Output
- `soil_erosion` (total soil loss [tdt⁻¹])
- `clay_erosion` (clay loss [tdt⁻¹])
- `silt_erosion` (silt loss [tdt⁻¹])
- `sand_erosion` (sand loss [tdt⁻¹])
- `sagg_erosion` (small aggregates loss [tdt⁻¹])
- `lagg_erosion` (large aggregates loss [tdt⁻¹])
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
- `bed` (potential river erosion [tdt⁻¹ => kg s⁻¹])
- `bank` (potential bank erosion [tdt⁻¹ => kg s⁻¹])
"""
function river_erosion_julian_torres(waterlevel, d50, width, length, slope, dt)
    if waterlevel > 0.0
        # Bed and Bank from Shields diagram, Da Silva & Yalin (2017)
        E_ = (2.65 - 1) * g_gravity
        E = (E_ * (d50 * 1e-3)^3 / 1e-12)^0.33
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
            g_gravity * hydrad * slope * SFbank / 100 * (1 + width / (2 * waterlevel))
        TEffbed =
            g_gravity * hydrad * slope * (1 - SFbank / 100) * (1 + 2 * waterlevel / width)

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

    return bed, bank
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

# Output
- `erosion` (river erosion [tdt⁻¹ => kg s⁻¹])
- `excess_sediment` (updated excess sediment [tdt⁻¹ => kg s⁻¹])
- `store` (updated sediment store [t => kg])
"""
function river_erosion_store(excess_sediment, store)
    # River erosion of the previously deposited sediment
    erosion = min(store, excess_sediment)
    # Update the excess sediment and the sediment store
    excess_sediment -= erosion
    store -= erosion
    return erosion, excess_sediment, store
end
