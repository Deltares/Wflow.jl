Base.@kwdef struct SimpleReservoir{T}
    maxvolume::Vector{T}                                # maximum storage (above which water is spilled) [m³]
    area::Vector{T}                                     # reservoir area [m²]
    maxrelease::Vector{T}                               # maximum amount that can be released if below spillway [m³ s⁻¹]
    demand::Vector{T}                                   # minimum (environmental) flow requirement downstream of the reservoir [m³ s⁻¹]
    targetminfrac::Vector{T}                            # target minimum full fraction (of max storage) [-]
    targetfullfrac::Vector{T}                           # target fraction full (of max storage) [-]
    volume::Vector{T} = targetfullfrac .* maxvolume     # reservoir volume [m³]
    inflow::Vector{T} = fill(mv, length(area))          # inflow into reservoir [m³]
    outflow::Vector{T} = fill(mv, length(area))         # outflow from reservoir [m³ s⁻¹]
    percfull::Vector{T} = fill(mv, length(area))        # fraction full (of max storage) [-]
    demandrelease::Vector{T} = fill(mv, length(area))   # minimum (environmental) flow released from reservoir [m³ s⁻¹]
    precipitation::Vector{T} = fill(mv, length(area))   # average precipitation for reservoir area [mm]
    evaporation::Vector{T} = fill(mv, length(area))     # average evaporation for reservoir area [mm]

    function SimpleReservoir{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

"""
Update a single reservoir at position `i`.

This is called from within the kinematic wave loop, therefore updating only for a single
element rather than all at once.
"""
function update(res::SimpleReservoir, i, inflow, timestepsecs)

    vol = (
        res.volume[i] +
        (inflow * timestepsecs) +
        (res.precipitation[i] * (timestepsecs / basetimestep.value) / 1000.0) * res.area[i] -
        (res.evaporation[i] * (timestepsecs / basetimestep.value) / 1000.0) * res.area[i]
    )

    percfull = vol / res.maxvolume[i]
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(percfull, a = res.targetminfrac[i], c = 30.0)
    demandrelease = min(fac * res.demand[i] * timestepsecs, vol)
    vol = vol - demandrelease

    wantrel = max(0.0, vol - (res.maxvolume[i] * res.targetfullfrac[i]))
    # Assume extra maximum Q if spilling
    overflow_q = max((vol - res.maxvolume[i]), 0.0)
    torelease = min(wantrel, overflow_q + res.maxrelease[i] * timestepsecs - demandrelease)
    vol = vol - torelease
    outflow = torelease + demandrelease
    percfull = vol / res.maxvolume[i]

    # update values in place
    res.outflow[i] = outflow / timestepsecs
    res.inflow[i] = inflow
    res.demandrelease[i] = demandrelease / timestepsecs
    res.percfull[i] = percfull
    res.volume[i] = vol

    return res
end

Base.@kwdef struct NaturalLake{T}
    loc_id::Vector{Int}                     # location id of lake outlet
    lowerlake_ind::Vector{Int}              # Index of lower lake (linked lakes)
    area::Vector{T}                         # lake area [m²]
    threshold::Vector{T}                    # water level threshold H₀ [m] below that level outflow is zero
    storfunc::Vector{Int}                   # type of lake storage curve, 1: S = AH, 2: S = f(H) from lake data and interpolation
    outflowfunc::Vector{Int}                # type of lake rating curve, 1: Q = f(H) from lake data and interpolation, 2: General Q = b(H - H₀)ᵉ, 3: Case of Puls Approach Q = b(H - H₀)²
    b::Vector{T}                            # rating curve coefficient
    e::Vector{T}                            # rating curve exponent
    sh::Vector{DataFrame}                   # data for storage curve
    hq::Vector{DataFrame}                   # data for rating curve
    avg_waterlevel::Vector{T}               # average water level [m] (cold state)
    waterlevel::Vector{T} = copy(avg_waterlevel) # waterlevel H [m] of lake
    inflow::Vector{T} = fill(mv, length(area))   # inflow to the lake [m³ s⁻¹]
    storage::Vector{T} = initialize_storage(storfunc, area, waterlevel, sh) # storage lake [m³]
    outflow::Vector{T} = fill(mv, length(area))        # outflow lake [m³ s⁻¹]
    precipitation::Vector{T} = fill(mv, length(area))  # average precipitation for lake area [mm]
    evaporation::Vector{T} = fill(mv, length(area))    # average evaporation for lake area [mm]

    function NaturalLake{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

"Determine the initial storage depending on the storage function"
function initialize_storage(storfunc, area, waterlevel, sh)
    storage = similar(area)
    for i in eachindex(storage)
        if storfunc[i] == 1
            storage[i] = area[i] * waterlevel[i]
        else
            storage[i] = interpolate_lineare(waterlevel[i], sh[i].H, sh[i].S)
        end
    end
    return storage
end

statenames(::NaturalLake) = ("waterlevel_lake",)

"""
    statenames(::NaturalLake)
    statenames(::SimpleReservoir)
    statenames(::LateralSSF)
    statenames(::SurfaceFlow)

Get a tuple of symbols representing the fields that are model states.
"""
function statenames end

function interpolate_linear(x, xp, fp)
    if x <= min(xp)
        return min(xp)
    elseif x >= max(xp)
        return max(xp)
    else
        idx = findall(i -> i <= x, xp)
        i1 = last(idx)
        i2 = i1 + 1
        return fp[i1] * (1.0 - (x - xp[i1]) / (xp[i2] - xp[i1])) +
               fp[i2] * (x - xp[i1]) / (xp[i2] - xp[i1])
    end
end

"""
Update a single lake at position `i`.

This is called from within the kinematic wave loop, therefore updating only for a single
element rather than all at once.
"""
function update(lake::NaturalLake, i, inflow, doy, timestepsecs)

    lo = lake.lowerlake_ind[i]
    has_lowerlake = lo != 0

    col = max(doy + 1, 366)
    ### Modified Puls Approach (Burek et al., 2013, LISFLOOD) ###
    # outflowfunc = 3
    # Calculate lake factor and SI parameter
    if lake.outflowfunc[i] == 3
        lakefactor = lake.area[i] / (timestepsecs * pow(lake.b[i], 0.5))
        si_factor =
            (
                lake.storage[i] +
                (lake.precipitation[i] - lake.evaporation[i]) *
                (timestepsecs / basetimestep.value) *
                lake.area[i] / 1000.0
            ) / timestepsecs + inflow
        #Adjust SIFactor for ResThreshold != 0
        si_factor_adj = si_factor - lake.area[i] * lake.threshold[i] / timestepsecs
        #Calculate the new lake outflow/waterlevel/storage

        if si_factor_adj > 0.0
            outflow = pow(
                -lakefactor + pow((pow(lakefactor, 2.0) + 2.0 * si_factor_adj), 0.5),
                2.0,
            )
        else
            outflow = 0.0
        end
        storage = (si_factor - outflow) * timestepsecs
        waterlevel = storage / lake.area[i]
    end

    ### Linearisation for specific storage/rating curves ###
    if lake.outflowfunc[i] == 1 || lake.outflowfunc[i] == 2

        diff_wl = has_lowerlake ? lake.waterlevel[i] - lake.waterlevel[lo] : 0.0

        if lake.outflowfunc[i] == 1
            outflow =
                interpolate_linear(lake.waterlevel[i], lake.hq[i][!, 1], lake.hq[i][!, col])
        else
            if diff_wl >= 0.0
                outflow = max(
                    pow((lake.b[i] * (lake.waterlevel[i] - lake.threshold[i])), lake.e[i]),
                    0.0,
                )
            else
                outflow = min(
                    pow(
                        (-1.0 * lake.b[i] * (lake.waterlevel[lo] - lake.threshold[i])),
                        lake.e[i],
                    ),
                    0.0,
                )
            end
        end

        storage =
            lake.storage[i] +
            inflow * timestepsecs +
            (lake.precipitation[i] / 1000.0) *
            (timestepsecs / basetimestep.value) *
            lake.area[i] -
            (lake.evaporation[i] / 1000.0) *
            (timestepsecs / basetimestep.value) *
            lake.area[i] - outflow * timestepsecs

        waterlevel = if lake.storfunc[i] == 1
            lake.waterlevel[i] + (storage - lake.storage[i]) / lake.area[i]
        else
            interpolate_linear(storage, lake.sh[i].S, lake.sh[i].H)
        end

        # update lower lake (linked lakes) in case flow from lower lake to upper lake occurs
        if diff_wl < 0.0
            lowerlake_storage = lake.storage[lo] + outflow * timestepsecs

            lowerlake_waterlevel = if lake.storfunc[lo] == 1
                lake.waterlevel[lo] + (lowerlake_storage - lake.storage[lo]) / lake.area[lo]
            else
                interpolate_linear(lowerlake_storage, lake.sh[lo].S, lake.sh[lo].H)
            end

            # update values for the lower lake in place
            lake.outflow[lo] = -outflow
            lake.storage[lo] = lowerlake_storage
            lake.waterlevel[lo] = lowerlake_waterlevel
        end

    end

    # update values in place
    lake.outflow[i] = outflow
    lake.waterlevel[i] = waterlevel
    lake.inflow[i] = inflow
    lake.storage[i] = storage

    return lake
end
