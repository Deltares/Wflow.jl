
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
)
    # TODO: add other rainfall interception method (lui)
    # TODO: include subdaily Gash model
    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopygapfraction)
    pt = min(0.1 * canopygapfraction, 1.0 - canopygapfraction)
    pfrac = 1.0 - canopygapfraction - pt
    p_sat = (-cmax / e_r) * log(1.0 - min((e_r / pfrac), 1.0))
    p_sat = isinf(p_sat) ? 0.0 : p_sat

    # large storms P > P_sat
    largestorms = precipitation > p_sat

    iwet = largestorms ? (pfrac * p_sat) - cmax : precipitation * pfrac
    isat = largestorms ? (e_r) * (precipitation - p_sat) : 0.0
    idry = largestorms ? cmax : 0.0
    itrunc = 0.0

    stemflow = pt * precipitation

    throughfall = precipitation - iwet - idry - isat - itrunc - stemflow
    interception = iwet + idry + isat + itrunc

    # Now corect for area without any Interception (say open water Cmax -- zero)
    cmaxzero = cmax <= 0.0
    throughfall = cmaxzero ? precipitation : throughfall
    interception = cmaxzero ? 0.0 : interception
    stemflow = cmaxzero ? 0.0 : stemflow

    # Now corect for maximum potential evap
    canopy_drainage = interception > maxevap ? interception - maxevap : 0.0
    interception = min(interception, maxevap)

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
    canopystorage,
    canopygapfraction,
    cmax,
)

    # TODO: improve computation of stemflow partitioning coefficient pt (0.1 * canopygapfraction)
    pt = min(0.1 * canopygapfraction, 1.0 - canopygapfraction)

    # Amount of p that falls on the canopy
    precip_canopy = (1.0 - canopygapfraction - pt) * precipitation

    # Canopystorage cannot be larger than cmax, no gravity drainage below that. This check
    # is required because cmax can change over time
    canopy_drainage1 = canopystorage > cmax ? canopystorage - cmax : 0.0
    canopystorage -= canopy_drainage1

    # Add the precipitation that falls on the canopy to the store
    canopystorage += precip_canopy

    # Evaporation, make sure the store does not get negative
    canopy_evap = min(canopystorage, potential_evaporation) # interception rate
    canopystorage -= canopy_evap

    # Drain the canopystorage again if needed
    canopy_drainage2 = canopystorage > cmax ? canopystorage - cmax : 0.0
    canopystorage -= canopy_drainage2

    # Calculate throughfall and stemflow
    throughfall = canopy_drainage1 + canopy_drainage2 + canopygapfraction * precipitation
    stemflow = precipitation * pt

    return throughfall, canopy_evap, stemflow, canopystorage
end
