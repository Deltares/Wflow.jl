
"""
    rainfall_interception_gash(cmax, e_r, canopygapfraction, precipitation, canopystorage, maxevap)

Interception according to the Gash model (for daily timesteps). `cmax` is the maximum canopy
storage and `e_r` is the ratio of the average evaporation from the wet canopy and the
average precipitation intensity on a saturated canopy.
"""
function rainfall_interception_gash(
    cmax,
    e_r,
    canopygapfraction,
    precipitation,
    canopystorage,
    maxevap,
    dt,
)
    # TODO: add other rainfall interception method (lui)
    # TODO: include subdaily Gash model
    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopygapfraction)
    # [-]
    fraction_stemflow = min(0.1 * canopygapfraction, 1.0 - canopygapfraction)
    fraction_interception = 1.0 - canopygapfraction - fraction_stemflow

    # [m s⁻¹] (from cmax / dt)
    precipitation_saturation = if fraction_interception > e_r
        (-cmax / e_r) * log(1.0 - e_r / fraction_interception) / dt
    else
        0.0
    end

    # large storms P > P_sat
    largestorms = precipitation > precipitation_saturation

    # I_wet: Interception during wetting up of the canopy
    # I_sat: Evaporation during canopy saturation
    # I_dry: Evaporation after rainfall stops
    # [m s⁻¹]
    if largestorms
        I_wet = fraction_interception * precipitation_saturation - cmax / dt
        I_sat = e_r * (precipitation - precipitation_saturation)
        I_dry = cmax / dt
    else
        I_wet = precipitation * fraction_interception
        I_sat = 0.0
        I_dry = 0.0
    end

    I_trunc = 0.0

    # [m s⁻¹] = [-] * [m s⁻¹]
    stemflow = fraction_stemflow * precipitation

    throughfall = precipitation - I_wet - I_dry - I_sat - I_trunc - stemflow
    interception = I_wet + I_dry + I_sat + I_trunc

    # Now corect for area without any Interception (say open water Cmax -- zero)
    cmaxzero = (cmax <= 0.0)

    if cmaxzero
        throughfall = precipitation
        interception = 0.0
        stemflow = 0.0
    end

    # Now corect for maximum potential evap
    if interception > maxevap
        canopy_drainage = interception - maxevap
        interception = maxevap
    else
        canopy_drainage = 0
    end

    # Add surpluss to the throughfall
    throughfall += canopy_drainage

    return throughfall, interception, stemflow, canopystorage
end

"""
    rainfall_interception_modrut(precipitation, potential_evaporation, canopystorage, canopygapfraction, cmax)

Interception according to a modified Rutter model. The model is solved explicitly and there
is no drainage below `cmax`.
"""
function rainfall_interception_modrut(
    precipitation,
    potential_evaporation,
    canopy_storage,
    canopygapfraction,
    cmax,
    dt,
)

    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopygapfraction)
    fraction_stemflow = min(0.1 * canopygapfraction, 1.0 - canopygapfraction)

    # Amount of precipitation that falls on the canopy
    # [m s⁻¹] = [-] * [m s⁻¹]
    precipitation_canopy = (1.0 - canopygapfraction - fraction_stemflow) * precipitation

    # Calculate throughfall and stemflow
    # [m s⁻¹] = [-] * [m s⁻¹]
    throughfall = canopygapfraction * precipitation
    stemflow = fraction_stemflow * precipitation

    # Canopy storage cannot be larger than cmax, no gravity drainage below that. This check
    # is required because cmax can change over time
    if canopy_storage > cmax
        # [m]
        canopystorage_excess = canopy_storage - cmax
        canopy_storage = cmax

        # [m s⁻¹] = [m] / [s]
        throughfall += canopystorage_excess / dt
    end

    # Add the precipitation that falls on the canopy to the store
    # [m] += [m s⁻¹] * [s]
    canopy_storage += precipitation_canopy * dt

    # Evaporation, make sure the store does not get negative
    # [m s⁻¹] = min([m] / [s], [m s⁻¹])
    canopy_evaporation = min(canopy_storage / dt, potential_evaporation) # interception rate
    # [m] -= [m s⁻¹] * [s]
    canopy_storage -= canopy_evaporation * dt

    # Drain the canopystorage again if needed
    if canopy_storage > cmax
        # [m]
        canopystorage_excess = canopy_storage - cmax
        canopy_storage = cmax

        # [m s⁻¹] = [m] / [s]
        throughfall += canopy_storage_excess / dt
    end

    return throughfall, canopy_evaporation, stemflow, canopy_storage
end
