"Struct for storing lake model parameters"
@get_units @grid_loc @with_kw struct LakeParameters{T}
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

"Initialize lake model parameters"
function LakeParameters(config, dataset, inds_riv, nriv, pits)
    # read only lake data if lakes true
    # allow lakes only in river cells
    # note that these locations are only the lake outlet pixels
    lakelocs_2d = ncread(
        dataset,
        config,
        "routing.river_flow.lake.locs";
        optional = false,
        type = Int,
        fill = 0,
    )
    lakelocs = lakelocs_2d[inds_riv]

    # this holds the same ids as lakelocs, but covers the entire lake
    lakecoverage_2d = ncread(
        dataset,
        config,
        "routing.river_flow.lake.areas";
        optional = false,
        allow_missing = true,
    )
    # for each lake, a list of 2D indices, needed for getting the mean precipitation
    inds_lake_cov = Vector{CartesianIndex{2}}[]

    rev_inds_lake = zeros(Int, size(lakecoverage_2d))

    # construct a map from the rivers to the lakes and
    # a map of the lakes to the 2D model grid
    inds_lake_map2river = fill(0, nriv)
    inds_lake = CartesianIndex{2}[]
    lakecounter = 0
    for (i, ind) in enumerate(inds_riv)
        lake_id = lakelocs[i]
        if lake_id > 0
            push!(inds_lake, ind)
            lakecounter += 1
            inds_lake_map2river[i] = lakecounter
            rev_inds_lake[ind] = lakecounter

            # get all indices related to this lake outlet
            # done in this loop to ensure that the order is equal to the order in the
            # NaturalLake struct
            lake_cov = findall(isequal(lake_id), lakecoverage_2d)
            push!(inds_lake_cov, lake_cov)
        end
    end

    lakearea = ncread(
        dataset,
        config,
        "routing.river_flow.lake.area";
        optional = false,
        sel = inds_lake,
        type = FLOAT,
        fill = 0,
    )
    lake_b = ncread(
        dataset,
        config,
        "routing.river_flow.lake.b";
        optional = false,
        sel = inds_lake,
        type = FLOAT,
        fill = 0,
    )
    lake_e = ncread(
        dataset,
        config,
        "routing.river_flow.lake.e";
        optional = false,
        sel = inds_lake,
        type = FLOAT,
        fill = 0,
    )
    lake_threshold = ncread(
        dataset,
        config,
        "routing.river_flow.lake.threshold";
        optional = false,
        sel = inds_lake,
        type = FLOAT,
        fill = 0,
    )
    linked_lakelocs = ncread(
        dataset,
        config,
        "routing.river_flow.lake.linkedlakelocs";
        sel = inds_lake,
        defaults = 0,
        type = Int,
        fill = 0,
    )
    lake_storfunc = ncread(
        dataset,
        config,
        "routing.river_flow.lake.storfunc";
        optional = false,
        sel = inds_lake,
        type = Int,
        fill = 0,
    )
    lake_outflowfunc = ncread(
        dataset,
        config,
        "routing.river_flow.lake.outflowfunc";
        optional = false,
        sel = inds_lake,
        type = Int,
        fill = 0,
    )
    lake_waterlevel = ncread(
        dataset,
        config,
        "routing.river_flow.lake.waterlevel";
        optional = false,
        sel = inds_lake,
        type = FLOAT,
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
    parameters = LakeParameters{FLOAT}(;
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
    lake_network = (
        indices_outlet = inds_lake,
        indices_coverage = inds_lake_cov,
        reverse_indices = rev_inds_lake,
        river_indices = findall(x -> x ≠ 0, inds_lake_map2river),
    )
    return parameters, lake_network, inds_lake_map2river, lake_waterlevel, pits
end

"Struct for storing Lake model parameters"
@get_units @grid_loc @with_kw struct LakeVariables{T}
    waterlevel::Vector{T} | "m"                 # waterlevel H [m] of lake
    storage::Vector{T} | "m3"                   # storage lake [m³]
    outflow::Vector{T} | "m3 s-1"               # outflow of lake outlet [m³ s⁻¹]
    outflow_av::Vector{T} | "m3 s-1"            # average outflow lake [m³ s⁻¹] for model timestep Δt (including flow from lower to upper lake)
    actevap::Vector{T}                          # average actual evapotranspiration for lake area [mm Δt⁻¹] 
end

"Initialize lake model variables"
function LakeVariables(n, lake_waterlevel)
    variables = LakeVariables{FLOAT}(;
        waterlevel = lake_waterlevel,
        inflow = fill(MISSING_VALUE, n),
        storage = initialize_storage(lake_storfunc, lakearea, lake_waterlevel, sh),
        outflow = fill(MISSING_VALUE, n),
        outflow_av = fill(MISSING_VALUE, n),
        actevap = fill(MISSING_VALUE, n),
    )
    return variables
end

"Struct for storing lake model boundary conditions"
@get_units @grid_loc @with_kw struct LakeBC{T}
    inflow::Vector{T} | "m3 s-1"                # inflow to the lake [m³ s⁻¹] for model timestep Δt
    precipitation::Vector{T}                    # average precipitation for lake area [mm Δt⁻¹]
    evaporation::Vector{T}                      # average potential evaporation for lake area [mm Δt⁻¹]
end

"Initialize lake model boundary conditions"
function LakeBC(n)
    bc = LakeBC{FLOAT}(;
        inflow = fill(MISSING_VALUE, n),
        precipitation = fill(MISSING_VALUE, n),
        evaporation = fill(MISSING_VALUE, n),
    )
    return bc
end

"Lake model"
@with_kw struct Lake{T}
    boundary_conditions::LakeBC{T}
    parameters::LakeParameters{T}
    variables::LakeVariables{T}
end

"Initialize lake model"
function Lake(dataset, config, indices_river, n_river_cells, pits)
    parameters, lake_network, inds_lake_map2river, lake_waterlevel, pits =
        LakeParameters(dataset, config, indices_river, n_river_cells, pits)

    n_lakes = length(parameters.area)
    variables = LakeVariables(n_lakes, lake_waterlevel)
    boundary_conditions = LakeBC(n_lakes)

    lake = Lake{FLOAT}(; boundary_conditions, parameters, variables)

    return lake, lake_network, inds_lake_map2river, pits
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
    maxstorage = Vector{Union{FLOAT, Missing}}(missing, length(area))
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
function update!(model::Lake, i, inflow, doy, dt, dt_forcing)
    lake_bc = model.boundary_conditions
    lake_p = model.parameters
    lake_v = model.variables

    lo = lake_p.lowerlake_ind[i]
    has_lowerlake = lo != 0

    # limit lake evaporation based on total available volume [m³]
    precipitation = 0.001 * lake_bc.precipitation[i] * (dt / dt_forcing) * lake_p.area[i]
    available_volume = lake_v.storage[i] + inflow * dt + precipitation
    evap = 0.001 * lake_bc.evaporation[i] * (dt / dt_forcing) * lake_p.area[i]
    actevap = min(available_volume, evap) # [m³/dt]

    ### Modified Puls Approach (Burek et al., 2013, LISFLOOD) ###
    # outflowfunc = 3
    # Calculate lake factor and SI parameter
    if lake_p.outflowfunc[i] == 3
        lakefactor = lake_p.area[i] / (dt * pow(lake_p.b[i], 0.5))
        si_factor = (lake_v.storage[i] + precipitation - actevap) / dt + inflow
        # Adjust SIFactor for ResThreshold != 0
        si_factor_adj = si_factor - lake_p.area[i] * lake_p.threshold[i] / dt
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
        storage = (si_factor - outflow) * dt
        waterlevel = storage / lake_p.area[i]
    end

    ### Linearisation for specific storage/rating curves ###
    if lake_p.outflowfunc[i] == 1 || lake_p.outflowfunc[i] == 2
        diff_wl = has_lowerlake ? lake_v.waterlevel[i] - lake_v.waterlevel[lo] : 0.0

        storage_input = (lake_v.storage[i] + precipitation - actevap) / dt + inflow
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
                    maxflow = (dh * lake_p.area[i]) / dt
                    outflow = min(outflow, maxflow)
                else
                    outflow = FLOAT(0)
                end
            else
                if lake_v.waterlevel[lo] > lake_p.threshold[i]
                    dh = lake_v.waterlevel[lo] - lake_p.threshold[i]
                    outflow = -1.0 * lake_p.b[i] * pow(dh, lake_p.e[i])
                    maxflow = (dh * lake_p.area[lo]) / dt
                    outflow = max(outflow, -maxflow)
                else
                    outflow = FLOAT(0)
                end
            end
        end
        storage = (storage_input - outflow) * dt

        # update storage and outflow for lake with rating curve of type 1.
        if lake_p.outflowfunc[i] == 1
            overflow = max(0.0, (storage - lake_p.maxstorage[i]) / dt)
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
            lowerlake_storage = lake_v.storage[lo] + outflow * dt

            lowerlake_waterlevel = if lake_p.storfunc[lo] == 1
                lake_v.waterlevel[lo] +
                (lowerlake_storage - lake_v.storage[lo]) / lake_p.area[lo]
            else
                interpolate_linear(lowerlake_storage, lake_p.sh[lo].S, lake_p.sh[lo].H)
            end

            # update values for the lower lake in place
            lake_v.outflow[lo] = -outflow
            lake_v.outflow_av[lo] += -outflow * dt
            lake_v.storage[lo] = lowerlake_storage
            lake_v.waterlevel[lo] = lowerlake_waterlevel
        end
    end

    # update values in place
    lake_v.outflow[i] = outflow
    lake_v.waterlevel[i] = waterlevel
    lake_bc.inflow[i] += inflow * dt
    lake_v.outflow_av[i] += outflow * dt
    lake_v.storage[i] = storage
    lake_v.actevap[i] += 1000.0 * (actevap / lake_p.area[i])

    return nothing
end
