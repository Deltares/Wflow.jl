
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
            pt = 0.1 * canopy_gap_fraction
            pfrac = 1.0 - 1.1 * canopy_gap_fraction # > 0

            if e_r > pfrac
                p_sat = 0.0
            else
                p_sat = -cmax / e_r * log(1.0 - e_r / pfrac)
            end
        else
            pt = 1.0 - canopy_gap_fraction
            pfrac = 0
            p_sat = 0
        end

        # large storms P > P_sat
        large_storms = precipitation > p_sat

        if large_storms
            iwet = pfrac * p_sat - cmax
            isat = e_r * (precipitation - p_sat)
            idry = cmax
            interception = iwet + isat + idry
        else
            iwet = pfrac * precipitation
            interception = iwet
        end

        stem_flow = pt * precipitation
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
)
    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopy_gap_fraction)
    if canopy_gap_fraction < inv(1.1)
        pt = 0.1 * canopy_gap_fraction
        precip_canopy = (1.0 - canopy_gap_fraction - pt) * precipitation
    else
        pt = 1.0 - canopy_gap_fraction
        precip_canopy = 0.0
    end

    stemflow = precipitation * pt
    throughfall = canopy_gap_fraction * precipitation

    # Canopystorage cannot be larger than cmax, no gravity drainage below that. This check
    # is required because cmax can change over time
    if canopy_storage > cmax
        canopy_drainage = canopy_storage - cmax
        canopy_storage = cmax
        throughfall += canopy_drainage
    end

    # Add the precipitation that falls on the canopy to the store
    canopy_storage += precip_canopy

    # Evaporation, make sure the store does not get negative
    if potential_evaporation > canopy_storage
        canopy_evap = canopy_storage
        canopy_storage = 0.0
    else
        canopy_evap = potential_evaporation
        canopy_storage -= potential_evaporation
    end

    # Drain the canopy_storage again if needed
    if canopy_storage > cmax
        canopy_drainage = canopy_storage - cmax
        canopy_storage = cmax
        throughfall += canopy_drainage
    end

    return throughfall, canopy_evap, stemflow, canopy_storage
end
