Base.@kwdef struct SimpleReservoir{T}
    maxvolume::T                                # maximum storage (above which water is spilled) [m³]
    area::T                                     # reservoir area [m²]
    maxrelease::T                               # maximum amount that can be released if below spillway [m³ s⁻¹]
    demand::T                                   # minimum (environmental) flow requirement downstream of the reservoir [m³ s⁻¹]
    targetminfrac::T                            # target minimum full fraction (of max storage) [-]
    targetfullfrac::T                           # target fraction full (of max storage) [-]
    volume::T = targetfullfrac * maxvolume      # reservoir volume [m³]
    inflow::T = mv                              # inflow into reservoir [m³]
    outflow::T = mv                             # outflow from reservoir [m³ s⁻¹]
    percfull::T = mv                            # fraction full (of max storage) [-]
    demandrelease::T = mv                       # minimum (environmental) flow released from reservoir [m³ s⁻¹]
    precipitation::T = mv                       # average precipitation for reservoir area [mm]
    evaporation::T = mv                         # average evaporation for reservoir area [mm]
end

statenames(::SimpleReservoir) = (:volume,)

function update(res::SimpleReservoir, inflow, p, pet, timestepsecs)

    vol = (
        res.volume + (inflow * timestepsecs) + (p / 1000.0) * res.area -
        (pet / 1000.0) * res.area
    )

    percfull = vol / res.maxvolume
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(percfull, a = res.targetminfrac, c = 30.0)
    demandrelease = min(fac * res.demand * timestepsecs, vol)
    vol = vol - demandrelease

    wantrel = max(0.0, vol - (res.maxvolume * res.targetfullfrac))
    # Assume extra maximum Q if spilling
    overflow_q = max((vol - res.maxvolume), 0.0)
    torelease = min(wantrel, overflow_q + res.maxrelease * timestepsecs - demandrelease)
    vol = vol - torelease
    outflow = torelease + demandrelease
    percfull = vol / res.maxvolume


    return setproperties(
        res,
        (
            outflow = outflow / timestepsecs,
            inflow = inflow,
            demandrelease = demandrelease / timestepsecs,
            percfull = percfull,
            precipitation = p * basetimestep.value / timestepsecs,
            evaporation = pet * basetimestep.value / timestepsecs,
            volume = vol,
        ),
    )


end
