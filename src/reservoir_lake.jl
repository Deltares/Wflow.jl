
@get_units @grid_loc @with_kw struct ReservoirParameters{T}
    dt::T                                   # Model time step [s]
    maxvolume::Vector{T} | "m3"             # maximum storage (above which water is spilled) [m³]
    area::Vector{T} | "m2"                  # reservoir area [m²]
    maxrelease::Vector{T} | "m3 s-1"        # maximum amount that can be released if below spillway [m³ s⁻¹]
    demand::Vector{T} | "m3 s-1"            # minimum (environmental) flow requirement downstream of the reservoir [m³ s⁻¹]
    targetminfrac::Vector{T} | "-"          # target minimum full fraction (of max storage) [-]
    targetfullfrac::Vector{T} | "-"         # target fraction full (of max storage
end

function ReservoirParameters(config, nc, inds_riv, nriv, pits, dt)
    # read only reservoir data if reservoirs true
    # allow reservoirs only in river cells
    # note that these locations are only the reservoir outlet pixels
    reslocs = ncread(
        nc,
        config,
        "lateral.river.reservoir.locs";
        optional = false,
        sel = inds_riv,
        type = Int,
        fill = 0,
    )

    # this holds the same ids as reslocs, but covers the entire reservoir
    rescoverage_2d = ncread(
        nc,
        config,
        "lateral.river.reservoir.areas";
        optional = false,
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
        config,
        "lateral.river.reservoir.demand";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    resmaxrelease = ncread(
        nc,
        config,
        "lateral.river.reservoir.maxrelease";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    resmaxvolume = ncread(
        nc,
        config,
        "lateral.river.reservoir.maxvolume";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    resarea = ncread(
        nc,
        config,
        "lateral.river.reservoir.area";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    res_targetfullfrac = ncread(
        nc,
        config,
        "lateral.river.reservoir.targetfullfrac";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    res_targetminfrac = ncread(
        nc,
        config,
        "lateral.river.reservoir.targetminfrac";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )

    # for surface water routing reservoir locations are considered pits in the flow network
    # all upstream flow goes to the river and flows into the reservoir
    pits[inds_res] .= true

    reservoir_indices = (
        indices_outlet = inds_res,
        indices_coverage = inds_res_cov,
        reverse_indices = rev_inds_reservoir,
    )

    parameters = ReservoirParameters{Float}(;
        dt = dt,
        demand = resdemand,
        maxrelease = resmaxrelease,
        maxvolume = resmaxvolume,
        area = resarea,
        targetfullfrac = res_targetfullfrac,
        targetminfrac = res_targetminfrac,
    )

    return parameters, reservoir_indices, resindex, pits
end

@get_units @grid_loc @with_kw struct ReservoirVariables{T}
    volume::Vector{T} | "m3"                # reservoir volume [m³]
    outflow::Vector{T} | "m3 s-1"           # outflow from reservoir [m³ s⁻¹]
    totaloutflow::Vector{T} | "m3"          # total outflow from reservoir [m³]
    percfull::Vector{T} | "-"               # fraction full (of max storage) [-]
    demandrelease::Vector{T} | "m3 s-1"     # minimum (environmental) flow released from reservoir [m³ s⁻¹]
    actevap::Vector{T}                      # average actual evaporation for reservoir area [mm Δt⁻¹]
end

function ReservoirVariables(n, parameters)
    (; targetfullfrac, maxvolume) = parameters
    variables = ReservoirVariables{Float}(;
        volume = targetfullfrac .* maxvolume,
        outflow = fill(mv, n),
        totaloutflow = fill(mv, n),
        percfull = fill(mv, n),
        demandrelease = fill(mv, n),
        actevap = fill(mv, n),
    )
    return variables
end

@get_units @grid_loc @with_kw struct ReservoirBC{T}
    inflow::Vector{T} | "m3"        # total inflow into reservoir [m³]
    precipitation::Vector{T}        # average precipitation for reservoir area [mm Δt⁻¹]
    evaporation::Vector{T}          # average potential evaporation for reservoir area [mm Δt⁻¹]
end

function ReservoirBC(n)
    bc = ReservoirBC{Float}(;
        inflow = fill(mv, n),
        precipitation = fill(mv, n),
        evaporation = fill(mv, n),
    )
    return bc
end

@with_kw struct SimpleReservoir{T}
    boundary_conditions::ReservoirBC{T}
    parameters::ReservoirParameters{T}
    variables::ReservoirVariables{T}
end

function SimpleReservoir(config, nc, inds_riv, nriv, pits, dt)
    parameters, reservoir_indices, resindex, pits =
        ReservoirParameters(config, nc, inds_riv, nriv, pits, dt)

    n = length(parameters.area)
    @info "Read `$n` reservoir locations."

    variables = ReservoirVariables(n, parameters)
    boundary_conditions = ReservoirBC(n)
    reservoirs = SimpleReservoir{Float}(; boundary_conditions, parameters, variables)

    return reservoirs, resindex, reservoir_indices, pits
end

"""
Update a single reservoir at position `i`.

This is called from within the kinematic wave loop, therefore updating only for a single
element rather than all at once.
"""
function update!(res::SimpleReservoir, i, inflow, timestepsecs)
    res_bc = res.boundary_conditions
    res_p = res.parameters
    res_v = res.variables

    # limit lake evaporation based on total available volume [m³]
    precipitation =
        0.001 * res_bc.precipitation[i] * (timestepsecs / res_p.dt) * res_p.area[i]
    available_volume = res_v.volume[i] + inflow * timestepsecs + precipitation
    evap = 0.001 * res_bc.evaporation[i] * (timestepsecs / res_p.dt) * res_p.area[i]
    actevap = min(available_volume, evap) # [m³/timestepsecs]

    vol = res_v.volume[i] + (inflow * timestepsecs) + precipitation - actevap
    vol = max(vol, 0.0)

    percfull = vol / res_p.maxvolume[i]
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(percfull, res_p.targetminfrac[i], Float(1.0), Float(30.0))
    demandrelease = min(fac * res_p.demand[i] * timestepsecs, vol)
    vol = vol - demandrelease

    wantrel = max(0.0, vol - (res_p.maxvolume[i] * res_p.targetfullfrac[i]))
    # Assume extra maximum Q if spilling
    overflow_q = max((vol - res_p.maxvolume[i]), 0.0)
    torelease =
        min(wantrel, overflow_q + res_p.maxrelease[i] * timestepsecs - demandrelease)
    vol = vol - torelease
    outflow = torelease + demandrelease
    percfull = vol / res_p.maxvolume[i]

    # update values in place
    res_v.outflow[i] = outflow / timestepsecs
    res_bc.inflow[i] += inflow * timestepsecs
    res_v.totaloutflow[i] += outflow
    res_v.demandrelease[i] = demandrelease / timestepsecs
    res_v.percfull[i] = percfull
    res_v.volume[i] = vol
    res_v.actevap[i] += 1000.0 * (actevap / res_p.area[i])

    return nothing
end

@get_units @grid_loc @with_kw struct LakeParameters{T}
    dt::T                                           # Model time step [s]
    lowerlake_ind::Vector{Int} | "-"                # Index of lower lake (linked lakes)
    area::Vector{T} | "m2"                          # lake area [m²]
    maxstorage::Vector{Union{T, Missing}} | "m3"    # lake maximum storage from rating curve 1 [m³]
    threshold::Vector{T} | "m"                      # water level threshold H₀ [m] below that level outflow is zero
    storfunc::Vector{Int} | "-"                     # type of lake storage curve, 1: S = AH, 2: S = f(H) from lake data and interpolation
    outflowfunc::Vector{Int} | "-"                  # type of lake rating curve, 1: Q = f(H) from lake data and interpolation, 2: General Q = b(H - H₀)ᵉ, 3: Case of Puls Approach Q = b(H - H₀)²
    b::Vector{T} | "m3/2 s-1 (if e=3/2)"            # rating curve coefficient
    e::Vector{T} | "-"                              # rating curve exponent
    sh::Vector{Union{SH, Missing}}                  # data for storage curve
    hq::Vector{Union{HQ, Missing}}                  # data for rating curve    
end

function LakeParameters(config, nc, inds_riv, nriv, pits, dt)
    # read only lake data if lakes true
    # allow lakes only in river cells
    # note that these locations are only the lake outlet pixels
    lakelocs_2d = ncread(
        nc,
        config,
        "lateral.river.lake.locs";
        optional = false,
        type = Int,
        fill = 0,
    )
    lakelocs = lakelocs_2d[inds_riv]

    # this holds the same ids as lakelocs, but covers the entire lake
    lakecoverage_2d = ncread(
        nc,
        config,
        "lateral.river.lake.areas";
        optional = false,
        allow_missing = true,
    )
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
        config,
        "lateral.river.lake.area";
        optional = false,
        sel = inds_lake,
        type = Float,
        fill = 0,
    )
    lake_b = ncread(
        nc,
        config,
        "lateral.river.lake.b";
        optional = false,
        sel = inds_lake,
        type = Float,
        fill = 0,
    )
    lake_e = ncread(
        nc,
        config,
        "lateral.river.lake.e";
        optional = false,
        sel = inds_lake,
        type = Float,
        fill = 0,
    )
    lake_threshold = ncread(
        nc,
        config,
        "lateral.river.lake.threshold";
        optional = false,
        sel = inds_lake,
        type = Float,
        fill = 0,
    )
    linked_lakelocs = ncread(
        nc,
        config,
        "lateral.river.lake.linkedlakelocs";
        sel = inds_lake,
        defaults = 0,
        type = Int,
        fill = 0,
    )
    lake_storfunc = ncread(
        nc,
        config,
        "lateral.river.lake.storfunc";
        optional = false,
        sel = inds_lake,
        type = Int,
        fill = 0,
    )
    lake_outflowfunc = ncread(
        nc,
        config,
        "lateral.river.lake.outflowfunc";
        optional = false,
        sel = inds_lake,
        type = Int,
        fill = 0,
    )
    lake_waterlevel = ncread(
        nc,
        config,
        "lateral.river.lake.waterlevel";
        optional = false,
        sel = inds_lake,
        type = Float,
        fill = 0,
    )

    # for surface water routing lake locations are considered pits in the flow network
    # all upstream flow goes to the river and flows into the lake
    pits[inds_lake] .= true

    # This is currently the same length as all river cells, but will be the
    # length of all lake cells. To do that we need to introduce a mapping.
    n_lakes = length(inds_lake)
    lakelocs = lakelocs_2d[inds_lake]

    @info "Read `$n_lakes` lake locations."

    sh = Vector{Union{SH, Missing}}(missing, n_lakes)
    hq = Vector{Union{HQ, Missing}}(missing, n_lakes)
    lowerlake_ind = fill(0, n_lakes)
    # lake CSV parameter files are expected in the same directory as path_static
    path = dirname(input_path(config, config.input.path_static))

    for i in 1:n_lakes
        lakeloc = lakelocs[i]
        if linked_lakelocs[i] > 0
            lowerlake_ind[i] = only(findall(x -> x == linked_lakelocs[i], lakelocs))
        end

        if lake_storfunc[i] == 2
            csv_path = joinpath(path, "lake_sh_$lakeloc.csv")
            @info(
                "Read a storage curve from CSV file `$csv_path`, for lake location `$lakeloc`"
            )
            sh[i] = read_sh_csv(csv_path)
        end

        if lake_outflowfunc[i] == 1
            csv_path = joinpath(path, "lake_hq_$lakeloc.csv")
            @info(
                "Read a rating curve from CSV file `$csv_path`, for lake location `$lakeloc`"
            )
            hq[i] = read_hq_csv(csv_path)
        end

        if lake_outflowfunc[i] == 3 && lake_storfunc[i] != 1
            @warn(
                "For the modified pulse approach (LakeOutflowFunc = 3) the LakeStorFunc should be 1"
            )
        end
    end
    parameters = LakeParameters{Float}(;
        dt = dt,
        lowerlake_ind = lowerlake_ind,
        area = lakearea,
        maxstorage = maximum_storage(lake_storfunc, lake_outflowfunc, lakearea, sh, hq),
        threshold = lake_threshold,
        storfunc = lake_storfunc,
        outflowfunc = lake_outflowfunc,
        b = lake_b,
        e = lake_e,
        sh = sh,
        hq = hq,
    )
    lake_indices = (
        indices_outlet = inds_lake,
        indices_coverage = inds_lake_cov,
        reverse_indices = rev_inds_lake,
    )
    return parameters, lake_indices, lakeindex, lake_Waterlevel, pits
end

@get_units @grid_loc @with_kw struct LakeVariables{T}
    waterlevel::Vector{T} | "m"                 # waterlevel H [m] of lake
    storage::Vector{T} | "m3"                   # storage lake [m³]
    outflow::Vector{T} | "m3 s-1"               # outflow lake [m³ s⁻¹]
    totaloutflow::Vector{T} | "m3"              # total outflow lake [m³]
    actevap::Vector{T}                          # average actual evapotranspiration for lake area [mm Δt⁻¹] 
end

function LakeVariables(n, lake_waterlevel)
    variables = LakeVariables{Float}(;
        waterlevel = lake_waterlevel,
        inflow = fill(mv, n),
        storage = initialize_storage(lake_storfunc, lakearea, lake_waterlevel, sh),
        outflow = fill(mv, n),
        totaloutflow = fill(mv, n),
        actevap = fill(mv, n),
    )
    return variables
end

@get_units @grid_loc @with_kw struct LakeBC{T}
    inflow::Vector{T} | "m3"                    # inflow to the lake [m³]
    precipitation::Vector{T}                    # average precipitation for lake area [mm Δt⁻¹]
    evaporation::Vector{T}                      # average potential evaporation for lake area [mm Δt⁻¹]
end

function LakeBC(n)
    bc = LakeBC{Float}(;
        inflow = fill(mv, n),
        precipitation = fill(mv, n),
        evaporation = fill(mv, n),
    )
    return bc
end

@with_kw struct Lake{T}
    boundary_conditions::LakeBC{T}
    parameters::LakeParameters{T}
    variables::LakeVariables{T}
end

function Lake(config, nc, inds_riv, nriv, pits, dt)
    parameters, lake_indices, lakeindex, lake_waterlevel, pits =
        LakeParameters(config, nc, inds_riv, nriv, pits, dt)

    n = length(parameters.area)
    variables = LakeVariables(n, lake_waterlevel)
    boundary_conditions = LakeBC(n)

    lakes = Lake{Float}(; boundary_conditions, parameters, variables)

    return lakes, lakeindex, lake_indices, pits
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

"Determine the water level depending on the storage function"
function waterlevel(storfunc, area, storage, sh)
    waterlevel = similar(area)
    for i in eachindex(storage)
        if storfunc[i] == 1
            waterlevel[i] = storage[i] / area[i]
        else
            waterlevel[i] = interpolate_linear(storage[i], sh[i].S, sh[i].H)
        end
    end
    return waterlevel
end

"Determine the maximum storage for lakes with a rating curve of type 1"
function maximum_storage(storfunc, outflowfunc, area, sh, hq)
    maxstorage = Vector{Union{Float, Missing}}(missing, length(area))
    # maximum storage is based on the maximum water level (H) value in the H-Q table
    for i in eachindex(maxstorage)
        if outflowfunc[i] == 1
            if storfunc[i] == 2
                maxstorage[i] = interpolate_linear(maximum(hq[i].H), sh[i].H, sh[i].S)
            else
                maxstorage[i] = area[i] * maximum(hq[i].H)
            end
        end
    end
    return maxstorage
end

function interpolate_linear(x, xp, fp)
    if x <= minimum(xp)
        return minimum(fp)
    elseif x >= maximum(xp)
        return maximum(fp)
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
function update!(lake::Lake, i, inflow, doy, timestepsecs)
    lake_bc = lake.boundary_conditions
    lake_p = lake.parameters
    lake_v = lake.variables

    lo = lake_p.lowerlake_ind[i]
    has_lowerlake = lo != 0

    # limit lake evaporation based on total available volume [m³]
    precipitation =
        0.001 * lake_bc.precipitation[i] * (timestepsecs / lake_p.dt) * lake_p.area[i]
    available_volume = lake_v.storage[i] + inflow * timestepsecs + precipitation
    evap = 0.001 * lake_bc.evaporation[i] * (timestepsecs / lake_p.dt) * lake_p.area[i]
    actevap = min(available_volume, evap) # [m³/timestepsecs]

    ### Modified Puls Approach (Burek et al., 2013, LISFLOOD) ###
    # outflowfunc = 3
    # Calculate lake factor and SI parameter
    if lake_p.outflowfunc[i] == 3
        lakefactor = lake_p.area[i] / (timestepsecs * pow(lake_p.b[i], 0.5))
        si_factor = (lake_v.storage[i] + precipitation - actevap) / timestepsecs + inflow
        # Adjust SIFactor for ResThreshold != 0
        si_factor_adj = si_factor - lake_p.area[i] * lake_p.threshold[i] / timestepsecs
        # Calculate the new lake outflow/waterlevel/storage
        if si_factor_adj > 0.0
            quadratic_sol_term =
                -lakefactor + pow((pow(lakefactor, 2.0) + 4.0 * si_factor_adj), 0.5)
            if quadratic_sol_term > 0.0
                outflow = pow(0.5 * quadratic_sol_term, 2.0)
            else
                outflow = 0.0
            end
        else
            outflow = 0.0
        end
        outflow = min(outflow, si_factor)
        storage = (si_factor - outflow) * timestepsecs
        waterlevel = storage / lake_p.area[i]
    end

    ### Linearisation for specific storage/rating curves ###
    if lake_p.outflowfunc[i] == 1 || lake_p.outflowfunc[i] == 2
        diff_wl = has_lowerlake ? lake_v.waterlevel[i] - lake_v.waterlevel[lo] : 0.0

        storage_input =
            (lake_v.storage[i] + precipitation - actevap) / timestepsecs + inflow
        if lake_p.outflowfunc[i] == 1
            outflow = interpolate_linear(
                lake_v.waterlevel[i],
                lake_p.hq[i].H,
                lake_p.hq[i].Q[:, doy],
            )
            outflow = min(outflow, storage_input)
        else
            if diff_wl >= 0.0
                if lake_v.waterlevel[i] > lake_p.threshold[i]
                    dh = lake_v.waterlevel[i] - lake_p.threshold[i]
                    outflow = lake_p.b[i] * pow(dh, lake_p.e[i])
                    maxflow = (dh * lake_p.area[i]) / timestepsecs
                    outflow = min(outflow, maxflow)
                else
                    outflow = Float(0)
                end
            else
                if lake_v.waterlevel[lo] > lake_p.threshold[i]
                    dh = lake_v.waterlevel[lo] - lake_p.threshold[i]
                    outflow = -1.0 * lake_p.b[i] * pow(dh, lake_p.e[i])
                    maxflow = (dh * lake_p.area[lo]) / timestepsecs
                    outflow = max(outflow, -maxflow)
                else
                    outflow = Float(0)
                end
            end
        end
        storage = (storage_input - outflow) * timestepsecs

        # update storage and outflow for lake with rating curve of type 1.
        if lake_p.outflowfunc[i] == 1
            overflow = max(0.0, (storage - lake_p.maxstorage[i]) / timestepsecs)
            storage = min(storage, lake_p.maxstorage[i])
            outflow = outflow + overflow
        end

        waterlevel = if lake_p.storfunc[i] == 1
            lake_v.waterlevel[i] + (storage - lake_v.storage[i]) / lake_p.area[i]
        else
            interpolate_linear(storage, lake_p.sh[i].S, lake_p.sh[i].H)
        end

        # update lower lake (linked lakes) in case flow from lower lake to upper lake occurs
        if diff_wl < 0.0
            lowerlake_storage = lake_v.storage[lo] + outflow * timestepsecs

            lowerlake_waterlevel = if lake_p.storfunc[lo] == 1
                lake_v.waterlevel[lo] +
                (lowerlake_storage - lake_v.storage[lo]) / lake_p.area[lo]
            else
                interpolate_linear(lowerlake_storage, lake_p.sh[lo].S, lake_p.sh[lo].H)
            end

            # update values for the lower lake in place
            lake_v.outflow[lo] = -outflow
            lake_v.totaloutflow[lo] += -outflow * timestepsecs
            lake_v.storage[lo] = lowerlake_storage
            lake_v.waterlevel[lo] = lowerlake_waterlevel
        end
    end

    # update values in place
    lake_v.outflow[i] = max(outflow, 0.0) # for a linked lake flow can be negative
    lake_v.waterlevel[i] = waterlevel
    lake_bc.inflow[i] += inflow * timestepsecs
    lake_v.totaloutflow[i] += outflow * timestepsecs
    lake_v.storage[i] = storage
    lake_v.actevap[i] += 1000.0 * (actevap / lake_p.area[i])

    return nothing
end
