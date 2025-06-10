"Struct for storing lake model parameters"
@with_kw struct LakeParameters
    lowerlake_ind::Vector{Int}                  # Index of lower lake (linked lakes) [-]
    area::Vector{Float64}                       # lake area [m²]
    maxstorage::Vector{Union{Float64, Missing}} # lake maximum storage from rating curve 1 [m³]
    threshold::Vector{Float64}                  # water level threshold H₀ [m] below that level outflow is zero
    storfunc::Vector{Int}                       # type of lake storage curve, 1: S = AH, 2: S = f(H) from lake data and interpolation
    outflowfunc::Vector{Int}                    # type of lake rating curve, 1: Q = f(H) from lake data and interpolation, 2: General Q = b(H - H₀)ᵉ, 3: Case of Puls Approach Q = b(H - H₀)²
    b::Vector{Float64}                          # rating curve coefficient [m3/2 s-1] (if e=3/2)
    e::Vector{Float64}                          # rating curve exponent [-]
    sh::Vector{Union{SH, Missing}}              # data for storage curve
    hq::Vector{Union{HQ, Missing}}              # data for rating curve    
end

"Initialize lake model parameters"
function LakeParameters(dataset::NCDataset, config::Config, network::NetworkWaterBody)
    (; indices_outlet) = network

    lens = lens_input_parameter(config, "lake_surface__area"; optional = false)
    lakearea = ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
    lens = lens_input_parameter(
        config,
        "lake_water__rating_curve_coefficient";
        optional = false,
    )
    lake_b = ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
    lens =
        lens_input_parameter(config, "lake_water__rating_curve_exponent"; optional = false)
    lake_e = ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
    lens = lens_input_parameter(
        config,
        "lake_water_flow_threshold-level__elevation";
        optional = false,
    )
    lake_threshold =
        ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
    lens = lens_input(config, "lake~lower_location__count")
    linked_lakelocs = ncread(
        dataset,
        config,
        lens;
        sel = indices_outlet,
        defaults = 0,
        type = Int,
        fill = 0,
    )
    lens = lens_input_parameter(
        config,
        "lake_water__storage_curve_type_count";
        optional = false,
    )
    lake_storfunc =
        ncread(dataset, config, lens; sel = indices_outlet, type = Int, fill = 0)
    lens = lens_input_parameter(
        config,
        "lake_water__rating_curve_type_count";
        optional = false,
    )
    lake_outflowfunc =
        ncread(dataset, config, lens; sel = indices_outlet, type = Int, fill = 0)
    lens = lens_input_parameter(
        config,
        "lake_water_surface__initial_elevation";
        optional = false,
    )
    lake_waterlevel =
        ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)

    n_lakes = length(indices_outlet)
    lakelocs = get_waterbody_locs(dataset, config, indices_outlet, "lake")

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
    parameters = LakeParameters(;
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

    return parameters, lake_waterlevel
end

"Struct for storing Lake model parameters"
@with_kw struct LakeVariables
    waterlevel::Vector{Float64}       # waterlevel H [m] of lake
    waterlevel_av::Vector{Float64}    # average waterlevel H [m] of lake for model timestep Δt
    storage::Vector{Float64}          # storage lake [m³]
    storage_av::Vector{Float64}       # average storage lake for model timestep Δt [m³]
    outflow::Vector{Float64}          # outflow of lake outlet [m³ s⁻¹]
    outflow_av::Vector{Float64}       # average outflow lake [m³ s⁻¹] for model timestep Δt (including flow from lower to upper lake)
    actevap::Vector{Float64}          # average actual evapotranspiration for lake area [mm Δt⁻¹] 
end

"Initialize lake model variables"
function LakeVariables(n::Int, parameters::LakeParameters, lake_waterlevel::Vector{Float64})
    (; storfunc, area, sh) = parameters
    variables = LakeVariables(;
        waterlevel = lake_waterlevel,
        waterlevel_av = fill(MISSING_VALUE, n),
        storage = initialize_storage(storfunc, area, lake_waterlevel, sh),
        storage_av = fill(MISSING_VALUE, n),
        outflow = fill(MISSING_VALUE, n),
        outflow_av = fill(MISSING_VALUE, n),
        actevap = fill(MISSING_VALUE, n),
    )
    return variables
end

"Struct for storing lake model boundary conditions"
@with_kw struct LakeBC
    inflow::Vector{Float64}           # inflow to the lake [m³ s⁻¹] for model timestep Δt
    precipitation::Vector{Float64}    # average precipitation for lake area [mm Δt⁻¹]
    evaporation::Vector{Float64}      # average potential evaporation for lake area [mm Δt⁻¹]
end

"Initialize lake model boundary conditions"
function LakeBC(n::Int)
    bc = LakeBC(;
        inflow = fill(MISSING_VALUE, n),
        precipitation = fill(MISSING_VALUE, n),
        evaporation = fill(MISSING_VALUE, n),
    )
    return bc
end

"Lake model"
@with_kw struct Lake
    boundary_conditions::LakeBC
    parameters::LakeParameters
    variables::LakeVariables
end

"Initialize lake model"
function Lake(dataset::NCDataset, config::Config, network::NetworkWaterBody)
    parameters, lake_waterlevel = LakeParameters(dataset, config, network)

    n_lakes = length(parameters.area)
    variables = LakeVariables(n_lakes, parameters, lake_waterlevel)
    boundary_conditions = LakeBC(n_lakes)

    lake = Lake(; boundary_conditions, parameters, variables)

    return lake
end

"Determine the initial storage depending on the storage function"
function initialize_storage(
    storfunc::Vector{Int},
    area::Vector{Float64},
    waterlevel::Vector{Float64},
    sh::Vector{Union{SH, Missing}},
)
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
function waterlevel(storfunc::Int, area::Float64, storage::Float64, sh::Union{SH, Missing})
    if storfunc == 1
        waterlevel = storage / area
    else
        waterlevel = interpolate_linear(storage, sh.S, sh.H)
    end
    return waterlevel
end

"Determine the maximum storage for lakes with a rating curve of type 1"
function maximum_storage(
    storfunc::Vector{Int},
    outflowfunc::Vector{Int},
    area::Vector{Float64},
    sh::Vector{Union{SH, Missing}},
    hq::Vector{Union{HQ, Missing}},
)
    maxstorage = Vector{Union{Float64, Missing}}(missing, length(area))
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
function update!(
    model::Lake,
    i::Int,
    inflow::Float64,
    doy::Int,
    dt::Float64,
    dt_forcing::Float64,
)
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
                    outflow = Float64(0)
                end
            else
                if lake_v.waterlevel[lo] > lake_p.threshold[i]
                    dh = lake_v.waterlevel[lo] - lake_p.threshold[i]
                    outflow = -1.0 * lake_p.b[i] * pow(dh, lake_p.e[i])
                    maxflow = (dh * lake_p.area[lo]) / dt
                    outflow = max(outflow, -maxflow)
                else
                    outflow = Float64(0)
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
    # instantaneous variables
    lake_v.outflow[i] = outflow
    lake_v.waterlevel[i] = waterlevel
    lake_v.storage[i] = storage
    # average variables (here accumulated for model timestep Δt)
    lake_bc.inflow[i] += inflow * dt
    lake_v.outflow_av[i] += outflow * dt
    lake_v.storage_av[i] += storage * dt
    lake_v.waterlevel_av[i] += waterlevel * dt
    lake_v.actevap[i] += 1000.0 * (actevap / lake_p.area[i])

    return nothing
end
