Base.@kwdef struct NaturalLake{T}
    loc_id::Int                   # location id of lake outlet
    lowerlake_ind::Int            # Index of lower lake (linked lakes)
    area::T                         # lake area [m²]
    threshold::T                    # water level threshold H₀ [m] below that level outflow is zero
    storfunc::Int                 # type of lake storage curve, 1: S = AH, 2: S = f(H) from lake data and interpolation
    outflowfunc::Int              # type of lake rating curve, 1: Q = f(H) from lake data and interpolation, 2: General Q = b(H - H₀)ᵉ, 3: Case of Puls Approach Q = b(H - H₀)²
    b::T                            # rating curve coefficient
    e::T                            # rating curve exponent
    sh::DataFrame                   # data for storage curve
    hq::DataFrame                   # data for rating curve
    avg_waterlevel::T               # average water level [m] (cold state)
    waterlevel::T = avg_waterlevel  # waterlevel H [m] of lake
    inflow::T = mv                  # inflow to the lake [m³ s⁻¹]
    storage::T =
        storfunc == 1 ? area * waterlevel : interpolate_linear(waterlevel, sh.H, sh.S) # storage lake [m³]
    outflow::T = mv                 # outflow lake [m³ s⁻¹]
    precipitation::T = mv           # average precipitation for lake area [mm]
    evaporation::T = mv             # average evaporation for lake area [mm]
end

statenames(::NaturalLake) = (:waterlevel,)

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

function update(lake::NaturalLake, inflow, p, pet, doy, timestepsecs; lower_lake = nothing)

    col = max(doy + 1, 366)
    ### Modified Puls Approach (Burek et al., 2013, LISFLOOD) ###
    # outflowfunc = 3
    # Calculate lake factor and SI parameter
    if lake.outflowfunc == 3
        lakefactor = lake.area / (timestepsecs * pow(lake.b, 0.5))
        si_factor = (lake.storage + (p - pet) * lake.area / 1000.0) / timestepsecs + inflow
        #Adjust SIFactor for ResThreshold != 0
        si_factor_adj = si_factor - lake.area * lake.threshold / timestepsecs
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
        waterlevel = storage / lake.area
    end

    ### Linearisation for specific storage/rating curves ###
    if lake.outflowfunc == 1 || lake.outflowfunc == 2

        diff_wl = lower_lake != nothing ? lake.waterlevel - lower_lake.waterlevel : 0.0

        if lake.outflowfunc == 1
            outflow = interpolate_linear(lake.waterlevel, lake.hq[!, 1], lake.hq[!, col])
        else
            if diff_wl >= 0.0
                outflow =
                    max(pow((lake.b * (lake.waterlevel - lake.threshold)), lake.e), 0.0)
            else
                outflow = min(
                    pow((-1.0 * lake.b * (lower_lake.waterlevel - lake.threshold)), lake.e),
                    0.0,
                )
            end
        end

        storage =
            lake.storage + inflow * timestepsecs + (p / 1000.0) * lake.area -
            (pet / 1000.0) * lake.area - outflow * timestepsecs
        waterlevel =
            lake.storfunc == 1 ? lake.waterlevel + (storage - lake.storage) / lake.area :
            interpolate_linear(storage, lake.sh.S, lake.sh.H)

        # update lower lake (linked lakes) in case flow from lower lake to upper lake occurs
        if diff_wl < 0.0
            lowerlake_storage = lowerlake_storage + outflow * timestepsecs
            lowerlake_waterlevel = lower_lake.storfunc == 1 ?
                lower_lake.waterlevel +
            (lowerlake_storage - lower_lake.storage) / lower_lake.area :
                interpolate_linear(lowerlake_storage, lower_lake.sh.S, lower_lake.sh.H)

            setproperties(
                lower_lake,
                (
                    outflow = -outflow,
                    storage = lowerlake_storage,
                    waterlevel = lowerlake_waterlevel,
                ),
            )
        end

    end

    return setproperties(
        lake,
        (
            outflow = outflow,
            waterlevel = waterlevel,
            inflow = inflow,
            storage = storage,
            precipitation = p * basetimestep.value / timestepsecs,
            evaporation = pet * basetimestep.value / timestepsecs,
        ),
    )

end
