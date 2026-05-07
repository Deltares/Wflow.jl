
"""
    rainfall_interception_gash(maximum_canopy_storage, evaporation_to_precipitation_ratio, canopy_gap_fraction, precipitation, canopy_storage, max_evaporation)

Interception according to the Gash model (for daily timesteps). `maximum_canopy_storage` is the maximum canopy
storage and `evaporation_to_precipitation_ratio` is the ratio of the average evaporation from the wet canopy and the
average precipitation intensity on a saturated canopy.
"""
function rainfall_interception_gash(
    maximum_canopy_storage,
    evaporation_to_precipitation_ratio,
    canopy_gap_fraction,
    precipitation,
    canopy_storage,
    max_evaporation,
)
    # TODO: add other rainfall interception method (lui)
    # TODO: include subdaily Gash model
    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopy_gap_fraction)
    if maximum_canopy_storage > 0.0
        if canopy_gap_fraction < inv(1.1)
            pt = 0.1 * canopy_gap_fraction
            pfrac = 1.0 - 1.1 * canopy_gap_fraction # > 0

            if evaporation_to_precipitation_ratio > pfrac
                p_sat = 0.0
            else
                p_sat =
                    -maximum_canopy_storage / evaporation_to_precipitation_ratio *
                    log(1.0 - evaporation_to_precipitation_ratio / pfrac)
            end
        else
            pt = 1.0 - canopy_gap_fraction
            pfrac = 0
            p_sat = 0
        end

        # large storms P > P_sat
        large_storms = precipitation > p_sat

        if large_storms
            iwet = pfrac * p_sat - maximum_canopy_storage
            isat = evaporation_to_precipitation_ratio * (precipitation - p_sat)
            idry = maximum_canopy_storage
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
    rainfall_interception_modrut(precipitation, potential_evaporation, canopy_storage, canopy_gap_fraction, maximum_canopy_storage)

Interception according to a modified Rutter model. The model is solved explicitly and there
is no drainage below `maximum_canopy_storage`.
"""
function rainfall_interception_modrut(
    precipitation,
    potential_evaporation,
    canopy_storage,
    canopy_gap_fraction,
    maximum_canopy_storage,
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

    # Canopystorage cannot be larger than maximum_canopy_storage, no gravity drainage below that. This check
    # is required because maximum_canopy_storage can change over time
    if canopy_storage > maximum_canopy_storage
        canopy_drainage = canopy_storage - maximum_canopy_storage
        canopy_storage = maximum_canopy_storage
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
    if canopy_storage > maximum_canopy_storage
        canopy_drainage = canopy_storage - maximum_canopy_storage
        canopy_storage = maximum_canopy_storage
        throughfall += canopy_drainage
    end

    return throughfall, canopy_evap, stemflow, canopy_storage
end
