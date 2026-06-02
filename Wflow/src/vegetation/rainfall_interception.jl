
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
    dt,
)
    # TODO: add other rainfall interception method (lui)
    # TODO: include subdaily Gash model
    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopy_gap_fraction)
    if maximum_canopy_storage > 0.0
        if canopy_gap_fraction < inv(1.1)
            fraction_stemflow = 0.1 * canopy_gap_fraction
            fraction_interception = 1.0 - 1.1 * canopy_gap_fraction # > 0

            precipitation_saturation = if evaporation_to_precipitation_ratio > fraction_interception
                0.0
            else
                -maximum_canopy_storage / (evaporation_to_precipitation_ratio * dt) * log(1.0 - evaporation_to_precipitation_ratio / fraction_interception)
            end
        else
            fraction_stemflow = 1.0 - canopy_gap_fraction
            fraction_interception = 0
            precipitation_saturation = 0
        end

        # large storms P > P_sat
        large_storms = precipitation > precipitation_saturation

        if large_storms
            iwet = fraction_interception * precipitation_saturation - maximum_canopy_storage / dt
            isat = evaporation_to_precipitation_ratio * (precipitation - precipitation_saturation)
            idry = maximum_canopy_storage / dt
            interception = iwet + isat + idry
        else
            iwet = fraction_interception * precipitation
            interception = iwet
        end

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
    dt,
)
    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopy_gap_fraction)
    if canopy_gap_fraction < inv(1.1)
        fraction_stemflow = 0.1 * canopy_gap_fraction
        precipitation_canopy =
            (1.0 - canopy_gap_fraction - fraction_stemflow) * precipitation
    else
        fraction_stemflow = 1.0 - canopy_gap_fraction
        precipitation_canopy = 0.0
    end

    stemflow = fraction_stemflow * precipitation
    throughfall = canopy_gap_fraction * precipitation

    # Canopystorage cannot be larger than maximum_canopy_storage, no gravity drainage below that. This check
    # is required because maximum_canopy_storage can change over time
    if canopy_storage > maximum_canopy_storage
        canopy_drainage = canopy_storage - maximum_canopy_storage
        canopy_storage = maximum_canopy_storage

        throughfall += canopy_drainage / dt
    end

    # Add the precipitation that falls on the canopy to the store
    canopy_storage += precipitation_canopy * dt

    # Evaporation, make sure the store does not get negative
    max_evaporation = canopy_storage / dt
    if potential_evaporation > max_evaporation
        canopy_evaporation = max_evaporation
        canopy_storage = 0.0
    else
        canopy_evaporation = potential_evaporation
        canopy_storage -= canopy_evaporation * dt
    end

    # Drain the canopy_storage again if needed
    if canopy_storage > maximum_canopy_storage
        canopy_drainage = canopy_storage - maximum_canopy_storage
        canopy_storage = maximum_canopy_storage

        throughfall += canopy_drainage / dt
    end

    return throughfall, canopy_evaporation, stemflow, canopy_storage
end
