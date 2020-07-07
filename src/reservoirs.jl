Base.@kwdef struct SimpleReservoir{T}
    maxvolume::T                                # maximum storage (above which water is spilled) [m³]
    area::T                                     # reservoir area [m²]
    maxrelease::T                               # maximum amount that can be released if below spillway [m³ s⁻¹]
    demand::T                                   # minimum (environmental) flow requirement downstream of the reservoir [m³ s⁻¹]
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

"""
    statenames(::Type{SimpleReservoir})

Returns Array{Symbol,1} for extracting model state fields.
"""
function statenames(::Type{SimpleReservoir})

    states = [:volume,]
    # TODO: (warm) states read from netcdf file or cold state (reinit=1, setting in ini file)

end

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
            precipitation = p,
            evaporation = pet,
            volume = vol,
        ),
    )


end
