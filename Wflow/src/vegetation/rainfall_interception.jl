
"""
    rainfall_interception_gash(cmax, e_r, canopy_gap_fraction, precipitation, canopy_storage, max_evaporation)

Interception according to the Gash model (for daily timesteps). `cmax` is the maximum canopy
storage and `e_r` is the ratio of the average evaporation from the wet canopy and the
average precipitation intensity on a saturated canopy.
"""
function rainfall_interception_gash(
    cmax,
    e_r,
    canopy_gap_fraction,
    precipitation,
    canopy_storage,
    max_evaporation,
)
    # TODO: add other rainfall interception method (lui)
    # TODO: include subdaily Gash model
    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopy_gap_fraction)
    if cmax > 0.0
        if canopy_gap_fraction < inv(1.1)
            # [-]
            fraction_stemflow = 0.1 * canopy_gap_fraction
            fraction_interception = 1.0 - 1.1 * canopy_gap_fraction # > 0

            # [m s⁻¹] (from cmax / dt)
            precipitation_saturation = if e_r > fraction_interception
                0.0
            else
                -cmax / e_r * log(1.0 - e_r / fraction_interception)
            end
        else
            fraction_stemflow = 1.0 - canopy_gap_fraction
            fraction_interception = 0
            precipitation_saturation = 0
        end

        # large storms P > P_sat
        large_storms = precipitation > precipitation_saturation

        # [m s⁻¹]
        if large_storms
            iwet = fraction_interception * precipitation_saturation - cmax
            isat = e_r * (precipitation - precipitation_saturation)
            idry = cmax
            interception = iwet + isat + idry
        else
            iwet = fraction_interception * precipitation
            interception = iwet
        end

        # [m s⁻¹]
        stem_flow = fraction_stemflow * precipitation
        throughfall = precipitation - interception - stem_flow

        if interception > max_evaporation
            canopy_drainage = interception - max_evaporation
            interception = max_evaporation
            throughfall += canopy_drainage
        end
    else
        throughfall = precipitation
        interception = 0.0
        stem_flow = 0.0
    end
    return throughfall, interception, stem_flow, canopy_storage
end

"""
    rainfall_interception_modrut(precipitation, potential_evaporation, canopy_storage, canopy_gap_fraction, cmax)

Interception according to a modified Rutter model. The model is solved explicitly and there
is no drainage below `cmax`.
"""
function rainfall_interception_modrut(
    precipitation,
    potential_evaporation,
    canopy_storage,
    canopy_gap_fraction,
    cmax,
    dt,
)
    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopy_gap_fraction)
    if canopy_gap_fraction < inv(1.1)
        # [-] = [-] * [-]
        fraction_stemflow = 0.1 * canopy_gap_fraction
        # [m s⁻¹] = [-] * [m s⁻¹]
        precipitation_canopy =
            (1.0 - canopy_gap_fraction - fraction_stemflow) * precipitation
    else
        # [-] = [-] - [-]
        fraction_stemflow = 1.0 - canopy_gap_fraction
        precipitation_canopy = 0.0
    end

    # [m s⁻¹] = [-] * [m s⁻¹]
    stemflow = fraction_stemflow * precipitation
    throughfall = canopy_gap_fraction * precipitation

    # Canopystorage cannot be larger than cmax, no gravity drainage below that. This check
    # is required because cmax can change over time
    if canopy_storage > cmax
        # [m]
        canopy_drainage = canopy_storage - cmax
        canopy_storage = cmax

        # [m s⁻¹] += [m] / [s]
        throughfall += canopy_drainage / dt
    end

    # Add the precipitation that falls on the canopy to the store
    # [m] += [m s⁻¹] * [s]
    canopy_storage += precipitation_canopy * dt

    # Evaporation, make sure the store does not get negative
    max_evaporation = canopy_storage / dt
    if potential_evaporation > max_evaporation
        canopy_evaporation = max_evaporation
        canopy_storage = 0.0
    else
        canopy_evaporation = potential_evaporation
        # [m] -= [m s⁻¹] * [s]
        canopy_storage -= canopy_evaporation * dt
    end

    # Drain the canopy_storage again if needed
    if canopy_storage > cmax
        # [m s⁻¹] = [m] / [s]
        canopy_drainage = (canopy_storage - cmax) / dt
        canopy_storage = cmax
        throughfall += canopy_drainage
    end

    return throughfall, canopy_evaporation, stemflow, canopy_storage
end
