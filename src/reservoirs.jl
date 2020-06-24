Base.@kwdef struct SimpleReservoir{T}
    maxvolume::T
    Δt::T
    area::T
    maxrelease::T
    demand::T
    targetminfrac::T
    targetfullfrac::T
    volume::T = targetfullfrac * maxvolume
    inflow::T = mv
    outflow::T = mv
    percfull::T = mv
    demandrelease::T = mv
    precipitation::T = mv
    evaporation::T = mv
end

function update(res::SimpleReservoir, inflow, p, pet)

    # dummy meteo variables pet and p
    pet = 4.0
    p = 3.0

    its = max(ceil(res.Δt / 21600.0), 1)
    _outflow = 0.0
    _demandrelease = 0.0
    percfull = 0.0
    vol = 0.0
    for n = 1:its
        vol = (
            res.volume + (inflow * res.Δt / its) + (p / its / 1000.0) * res.area - (pet / its / 1000.0) * res.area
        )

        percfull = vol / res.maxvolume
        # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
        fac = scurve(percfull, a = res.targetminfrac, c = 30.0)
        demandrelease = min(fac * res.demand * res.Δt / its, vol)
        vol = vol - demandrelease

        wantrel = max(0.0, vol - (res.maxvolume * res.targetfullfrac))
        # Assume extra maximum Q if spilling
        overflow_q = max((vol - res.maxvolume), 0.0)
        torelease =
            min(wantrel, overflow_q + res.maxrelease * res.Δt / its - demandrelease)
        vol = vol - torelease
        outflow = torelease + demandrelease
        percfull = vol / res.maxvolume

        _outflow = _outflow + outflow
        _demandrelease = _demandrelease + demandrelease

    end

    return setproperties(
        res,
        (
            outflow = _outflow / res.Δt,
            inflow = inflow,
            demandrelease = _demandrelease / res.Δt,
            percfull = percfull,
            precipitation = p,
            evaporation = pet,
            volume = vol,
        ),
    )


end
