@get_units @with_kw struct SimpleReservoir{T}
    maxvolume::Vector{T} | "m3"                         # maximum storage (above which water is spilled) [m³]
    area::Vector{T} | "m2"                              # reservoir area [m²]
    maxrelease::Vector{T} | "m3 s-1"                    # maximum amount that can be released if below spillway [m³ s⁻¹]
    demand::Vector{T} | "m3 s-1"                        # minimum (environmental) flow requirement downstream of the reservoir [m³ s⁻¹]
    targetminfrac::Vector{T} | "-"                      # target minimum full fraction (of max storage) [-]
    targetfullfrac::Vector{T} | "-"                     # target fraction full (of max storage) [-]
    volume::Vector{T} | "m3"                            # reservoir volume [m³]
    inflow::Vector{T} | "m3"                        # inflow into reservoir [m³]
    outflow::Vector{T} | "m3 s-1"                       # outflow from reservoir [m³ s⁻¹]
    totaloutflow::Vector{T} | "m3"                      # total outflow from reservoir [m³]
    percfull::Vector{T} | "-"                           # fraction full (of max storage) [-]
    demandrelease::Vector{T} | "m3 s-1"                 # minimum (environmental) flow released from reservoir [m³ s⁻¹]
    precipitation::Vector{T}                            # average precipitation for reservoir area [mm]
    evaporation::Vector{T}                              # average evaporation for reservoir area [mm]

    function SimpleReservoir{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

statevars(::SimpleReservoir) = (:volume,)

function initialize_simple_reservoir(config, nc, inds_riv, nriv, pits)
    # read only reservoir data if reservoirs true
    # allow reservoirs only in river cells
    # note that these locations are only the reservoir outlet pixels
    reslocs = ncread(
        nc,
        param(config, "input.lateral.river.reservoir.locs");
        sel = inds_riv,
        type = Int,
        fill = 0,
    )

    # this holds the same ids as reslocs, but covers the entire reservoir
    rescoverage_2d = ncread(
        nc,
        param(config, "input.lateral.river.reservoir.areas");
        allow_missing = true,
    )
    # for each reservoir, a list of 2D indices, needed for getting the mean precipitation
    inds_res_cov = Vector{CartesianIndex{2}}[]

    rev_inds_reservoir = zeros(Int, size(rescoverage_2d))

    # construct a map from the rivers to the reservoirs and
    # a map of the reservoirs to the 2D model grid
    resindex = fill(0, nriv)
    inds_res = CartesianIndex{2}[]
    rescounter = 0
    for (i, ind) in enumerate(inds_riv)
        res_id = reslocs[i]
        if res_id > 0
            push!(inds_res, ind)
            rescounter += 1
            resindex[i] = rescounter
            rev_inds_reservoir[ind] = rescounter

            # get all indices related to this reservoir outlet
            # done in this loop to ensure that the order is equal to the order in the
            # SimpleReservoir struct
            res_cov = findall(isequal(res_id), rescoverage_2d)
            push!(inds_res_cov, res_cov)
        end
    end

    resdemand = ncread(
        nc,
        param(config, "input.lateral.river.reservoir.demand");
        sel = inds_res,
        type = Float64,
        fill = 0,
    )
    resmaxrelease = ncread(
        nc,
        param(config, "input.lateral.river.reservoir.maxrelease");
        sel = inds_res,
        type = Float64,
        fill = 0,
    )
    resmaxvolume = ncread(
        nc,
        param(config, "input.lateral.river.reservoir.maxvolume");
        sel = inds_res,
        type = Float64,
        fill = 0,
    )
    resarea = ncread(
        nc,
        param(config, "input.lateral.river.reservoir.area");
        sel = inds_res,
        type = Float64,
        fill = 0,
    )
    res_targetfullfrac = ncread(
        nc,
        param(config, "input.lateral.river.reservoir.targetfullfrac");
        sel = inds_res,
        type = Float64,
        fill = 0,
    )
    res_targetminfrac = ncread(
        nc,
        param(config, "input.lateral.river.reservoir.targetminfrac");
        sel = inds_res,
        type = Float64,
        fill = 0,
    )

    # for surface water routing reservoir locations are considered pits in the flow network
    # all upstream flow goes to the river and flows into the reservoir
    pits[inds_res] .= true

    n = length(resarea)
    reservoirs = SimpleReservoir{Float64}(
        demand = resdemand,
        maxrelease = resmaxrelease,
        maxvolume = resmaxvolume,
        area = resarea,
        targetfullfrac = res_targetfullfrac,
        targetminfrac = res_targetminfrac,
        volume = res_targetfullfrac .* resmaxvolume,
        inflow = fill(mv, n),
        outflow = fill(mv, n),
        totaloutflow = fill(mv, n),
        percfull = fill(mv, n),
        demandrelease = fill(mv, n),
        precipitation = fill(mv, n),
        evaporation = fill(mv, n),
    )

    return reservoirs,
    resindex,
    (
        indices_outlet = inds_res,
        indices_coverage = inds_res_cov,
        reverse_indices = rev_inds_reservoir,
    ),
    pits
end

"""
Update a single reservoir at position `i`.

This is called from within the kinematic wave loop, therefore updating only for a single
element rather than all at once.
"""
function update(res::SimpleReservoir, i, inflow, timestepsecs)

    vol = max(
        0.0,
        (
            res.volume[i] +
            (inflow * timestepsecs) +
            (res.precipitation[i] * (timestepsecs / tosecond(basetimestep)) / 1000.0) *
            res.area[i] -
            (res.evaporation[i] * (timestepsecs / tosecond(basetimestep)) / 1000.0) *
            res.area[i]
        ),
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
    res.inflow[i] += inflow * timestepsecs
    res.totaloutflow[i] += outflow
    res.demandrelease[i] = demandrelease / timestepsecs
    res.percfull[i] = percfull
    res.volume[i] = vol

    return res
end

@get_units @with_kw struct NaturalLake{T}
    lowerlake_ind::Vector{Int} | "-"        # Index of lower lake (linked lakes)
    area::Vector{T} | "m2"                  # lake area [m²]
    threshold::Vector{T} | "m"              # water level threshold H₀ [m] below that level outflow is zero
    storfunc::Vector{Int} | "-"             # type of lake storage curve, 1: S = AH, 2: S = f(H) from lake data and interpolation
    outflowfunc::Vector{Int} | "-"          # type of lake rating curve, 1: Q = f(H) from lake data and interpolation, 2: General Q = b(H - H₀)ᵉ, 3: Case of Puls Approach Q = b(H - H₀)²
    b::Vector{T} | "m3/2 s-1 (if e=3/2)"    # rating curve coefficient
    e::Vector{T} | "-"                      # rating curve exponent
    sh::Vector{Union{SH, Missing}}          # data for storage curve
    hq::Vector{Union{HQ, Missing}}          # data for rating curve
    waterlevel::Vector{T} | "m"             # waterlevel H [m] of lake
    inflow::Vector{T} | "m3"                # inflow to the lake [m³]
    storage::Vector{T} | "m3"               # storage lake [m³]
    outflow::Vector{T} | "m3 s-1"           # outflow lake [m³ s⁻¹]
    totaloutflow::Vector{T} | "m3"          # total outflow lake [m³] 
    precipitation::Vector{T}                # average precipitation for lake area [mm]
    evaporation::Vector{T}                  # average evaporation for lake area [mm]

    function NaturalLake{T}(args...) where {T}
        equal_size_vectors(args)
        return new(args...)
    end
end

function initialize_natural_lake(config, static_path, nc, inds_riv, nriv, pits)
    # read only lake data if lakes true
    # allow lakes only in river cells
    # note that these locations are only the lake outlet pixels
    lakelocs = ncread(
        nc,
        param(config, "input.lateral.river.lake.locs");
        sel = inds_riv,
        type = Int,
        fill = 0,
    )

    # this holds the same ids as lakelocs, but covers the entire lake
    lakecoverage_2d =
        ncread(nc, param(config, "input.lateral.river.lake.areas"); allow_missing = true)
    # for each lake, a list of 2D indices, needed for getting the mean precipitation
    inds_lake_cov = Vector{CartesianIndex{2}}[]

    rev_inds_lake = zeros(Int, size(lakecoverage_2d))

    # construct a map from the rivers to the lakes and
    # a map of the lakes to the 2D model grid
    lakeindex = fill(0, nriv)
    inds_lake = CartesianIndex{2}[]
    lakecounter = 0
    for (i, ind) in enumerate(inds_riv)
        lake_id = lakelocs[i]
        if lake_id > 0
            push!(inds_lake, ind)
            lakecounter += 1
            lakeindex[i] = lakecounter
            rev_inds_lake[ind] = lakecounter

            # get all indices related to this lake outlet
            # done in this loop to ensure that the order is equal to the order in the
            # NaturalLake struct
            lake_cov = findall(isequal(lake_id), lakecoverage_2d)
            push!(inds_lake_cov, lake_cov)
        end
    end

    lakearea = ncread(
        nc,
        param(config, "input.lateral.river.lake.area");
        sel = inds_lake,
        type = Float64,
        fill = 0,
    )
    lake_b = ncread(
        nc,
        param(config, "input.lateral.river.lake.b");
        sel = inds_lake,
        type = Float64,
        fill = 0,
    )
    lake_e = ncread(
        nc,
        param(config, "input.lateral.river.lake.e");
        sel = inds_lake,
        type = Float64,
        fill = 0,
    )
    lake_threshold = ncread(
        nc,
        param(config, "input.lateral.river.lake.threshold");
        sel = inds_lake,
        type = Float64,
        fill = 0,
    )
    linked_lakelocs = ncread(
        nc,
        param(config, "input.lateral.river.lake.linkedlakelocs");
        sel = inds_lake,
        type = Int,
        fill = 0,
    )
    lake_storfunc = ncread(
        nc,
        param(config, "input.lateral.river.lake.storfunc");
        sel = inds_lake,
        type = Int,
        fill = 0,
    )
    lake_outflowfunc = ncread(
        nc,
        param(config, "input.lateral.river.lake.outflowfunc");
        sel = inds_lake,
        type = Int,
        fill = 0,
    )
    lake_waterlevel = ncread(
        nc,
        param(config, "input.lateral.river.lake.waterlevel");
        sel = inds_lake,
        type = Float64,
        fill = 0,
    )

    # for surface water routing lake locations are considered pits in the flow network
    # all upstream flow goes to the river and flows into the lake
    pits[inds_lake] .= true

    # This is currently the same length as all river cells, but will be the
    # length of all lake cells. To do that we need to introduce a mapping.
    n_lakes = length(inds_lake)

    sh = Vector{Union{SH, Missing}}(missing, n_lakes)
    hq = Vector{Union{HQ, Missing}}(missing, n_lakes)
    for i = 1:n_lakes
        lakeloc = lakelocs[i]

        if linked_lakelocs[i] > 0
            linked_lakelocs[i] = i
        else
            linked_lakelocs[i] = 0
        end

        if lake_storfunc[i] == 2
            csv_path = joinpath(static_path, "lake_sh_$lakeloc.csv")
            sh[i] = read_sh_csv(csv_path)
        end

        if lake_outflowfunc[i] == 1
            csv_path = joinpath(static_path, "lake_hq_$lakeloc.csv")
            hq[i] = read_hq_csv(csv_path)
        end

        if lake_outflowfunc[i] == 3 && lake_storfunc[i] != 1
            @warn("For the modified pulse approach (LakeOutflowFunc = 3) the LakeStorFunc should be 1")
        end
    end
    n = length(lakearea)
    lakes = NaturalLake{Float64}(
        lowerlake_ind = linked_lakelocs,
        area = lakearea,
        threshold = lake_threshold,
        storfunc = lake_storfunc,
        outflowfunc = lake_outflowfunc,
        b = lake_b,
        e = lake_e,
        waterlevel = lake_waterlevel,
        sh = sh,
        hq = hq,
        inflow = fill(mv, n),
        storage = initialize_storage(lake_storfunc, lakearea, lake_waterlevel, sh),
        outflow = fill(mv, n),
        totaloutflow = fill(mv, n),
        precipitation = fill(mv, n),
        evaporation = fill(mv, n),
    )

    return lakes,
    lakeindex,
    (
        indices_outlet = inds_lake,
        indices_coverage = inds_lake_cov,
        reverse_indices = rev_inds_lake,
    ),
    pits
end


"Determine the initial storage depending on the storage function"
function initialize_storage(storfunc, area, waterlevel, sh)
    storage = similar(area)
    for i in eachindex(storage)
        if storfunc[i] == 1
            storage[i] = area[i] * waterlevel[i]
        else
            storage[i] = interpolate_linear(waterlevel[i], sh[i].H, sh[i].S)
        end
    end
    return storage
end

statevars(::NaturalLake) = (:waterlevel,)

"""
    statevars(::SBM; snow=false)
    statevars(::LateralSSF)
    statevars(::SurfaceFlow)
    statevars(::NaturalLake)
    statevars(::SimpleReservoir)

Get a tuple of symbols representing the fields that are model states.
"""
function statevars end

function interpolate_linear(x, xp, fp)
    if x <= minimum(xp)
        return minimum(xp)
    elseif x >= maximum(xp)
        return maximum(xp)
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

    col = max(doy, 365)
    ### Modified Puls Approach (Burek et al., 2013, LISFLOOD) ###
    # outflowfunc = 3
    # Calculate lake factor and SI parameter
    if lake.outflowfunc[i] == 3
        lakefactor = lake.area[i] / (timestepsecs * pow(lake.b[i], 0.5))
        si_factor =
            (
                lake.storage[i] +
                (lake.precipitation[i] - lake.evaporation[i]) *
                (timestepsecs / tosecond(basetimestep)) *
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
                interpolate_linear(lake.waterlevel[i], lake.hq[i].H, lake.hq[i].Q[:, col])
        else
            if diff_wl >= 0.0
                outflow = max(
                    lake.b[i] * pow((lake.waterlevel[i] - lake.threshold[i]), lake.e[i]),
                    0.0,
                )
            else
                outflow = min(
                    -1.0 *
                    lake.b[i] *
                    pow((lake.waterlevel[lo] - lake.threshold[i]), lake.e[i]),
                    0.0,
                )
            end
        end

        storage =
            lake.storage[i] +
            inflow * timestepsecs +
            (lake.precipitation[i] / 1000.0) *
            (timestepsecs / tosecond(basetimestep)) *
            lake.area[i] -
            (lake.evaporation[i] / 1000.0) *
            (timestepsecs / tosecond(basetimestep)) *
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
    lake.outflow[i] = max(outflow, 0.0) # for a linked lake flow can be negative
    lake.waterlevel[i] = waterlevel
    lake.inflow[i] += inflow * timestepsecs
    lake.totaloutflow[i] += outflow * timestepsecs
    lake.storage[i] = storage

    return lake
end
