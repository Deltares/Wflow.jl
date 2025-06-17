"Struct for storing reservoir model parameters"
@with_kw struct ReservoirParameters
    # type of reservoir storage curve, 1: S = AH, 2: S = f(H) from reservoir data and
    # interpolation
    storfunc::Vector{Int}
    # type of reservoir rating curve, 1: Q = f(H) from reservoir data and interpolation, 2:
    # General Q = b(H - H₀)ᵉ, 3: Case of Puls Approach Q = b(H - H₀)², 4: Simple reservoir
    outflowfunc::Vector{Int}
    # reservoir area [m²]
    area::Vector{Float64}
    # index of lower reservoir (linked reservoirs) [-]
    lower_reservoir_ind::Vector{Int} = fill(0, length(area))
    # reservoir maximum storage for rating curve types 1 and 4 [m³]
    maxstorage::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # water level threshold H₀ [m] below that level outflow is zero
    threshold::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # rating curve coefficient [m3/2 s-1] (if e=3/2)
    b::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # rating curve exponent [-]
    e::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # data for storage curve
    sh::Vector{Union{SH, Missing}} = Vector{Union{SH, Missing}}(missing, length(area))
    # data for rating curve
    hq::Vector{Union{HQ, Missing}} = Vector{Union{HQ, Missing}}(missing, length(area))
    # column index of rating curve data hq
    col_index_hq::Vector{Int} = [1]
    # maximum amount that can be released if below spillway for rating curve type 4 [m³ s⁻¹]
    maxrelease::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # minimum (environmental) flow requirement downstream of the reservoir for rating curve
    # type 4 [m³ s⁻¹]
    demand::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # target minimum full fraction (of max storage) for rating curve type 4 [-]
    targetminfrac::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # target fraction full (of max storage) for rating curve type 4 [-]
    targetfullfrac::Vector{Float64} = fill(MISSING_VALUE, length(area))
end

"Initialize reservoir model parameters"
function ReservoirParameters(dataset::NCDataset, config::Config, network::NetworkReservoir)
    (; indices_outlet) = network

    lens = lens_input_parameter(config, "reservoir_surface__area"; optional = false)
    area = ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
    lens = lens_input_parameter(
        config,
        "reservoir_water_surface__initial_elevation";
        optional = false,
    )
    waterlevel =
        ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
    lens = lens_input_parameter(
        config,
        "reservoir_water__storage_curve_type_count";
        optional = false,
    )
    storfunc = ncread(dataset, config, lens; sel = indices_outlet, type = Int, fill = 0)
    lens = lens_input_parameter(
        config,
        "reservoir_water__rating_curve_type_count";
        optional = false,
    )
    outflowfunc = ncread(dataset, config, lens; sel = indices_outlet, type = Int, fill = 0)

    lens = lens_input(config, "reservoir~lower_location__count")
    linked_reslocs = ncread(
        dataset,
        config,
        lens;
        sel = indices_outlet,
        defaults = 0,
        type = Int,
        fill = 0,
    )

    n_reservoirs = length(area)
    lens = lens_input(config, "reservoir_location__count"; optional = false)
    reslocs = ncread(dataset, config, lens; sel = indices_outlet, type = Int, fill = 0)
    @info "Read `$n_reservoirs` reservoir locations."

    parameters = ReservoirParameters(; area, outflowfunc, storfunc)

    if 2 in outflowfunc || 3 in outflowfunc
        lens = lens_input_parameter(
            config,
            "reservoir_water_flow_threshold-level__elevation";
            optional = false,
        )
        threshold =
            ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
        lens = lens_input_parameter(
            config,
            "reservoir_water__rating_curve_coefficient";
            optional = false,
        )
        b = ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
        lens = lens_input_parameter(
            config,
            "reservoir_water__rating_curve_exponent";
            optional = false,
        )
        e = ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
    end
    if 4 in outflowfunc
        lens = lens_input_parameter(
            config,
            "reservoir_water_demand~required~downstream__volume_flow_rate";
            optional = false,
        )
        demand =
            ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
        lens = lens_input_parameter(
            config,
            "reservoir_water_release-below-spillway__max_volume_flow_rate";
            optional = false,
        )
        maxrelease =
            ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
        lens = lens_input_parameter(config, "reservoir_water__max_volume"; optional = false)
        maxstorage =
            ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
        lens = lens_input_parameter(
            config,
            "reservoir_water~full-target__volume_fraction";
            optional = false,
        )
        targetfullfrac =
            ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
        lens = lens_input_parameter(
            config,
            "reservoir_water~min-target__volume_fraction";
            optional = false,
        )
        targetminfrac =
            ncread(dataset, config, lens; sel = indices_outlet, type = Float64, fill = 0)
    end

    # reservoir CSV parameter files are expected in the same directory as path_static
    path = dirname(input_path(config, config.input.path_static))
    for i in 1:n_reservoirs
        resloc = reslocs[i]
        if linked_reslocs[i] > 0
            @reset parameters.lower_reservoir_ind[i] =
                only(findall(x -> x == linked_reslocs[i], reslocs))
        end

        if storfunc[i] == 2
            csv_path = joinpath(path, "reservoir_sh_$resloc.csv")
            @info(
                "Read a storage curve from CSV file `$csv_path`, for reservoir location `$resloc`"
            )
            @reset parameters.sh[i] = read_sh_csv(csv_path)
        end

        if outflowfunc[i] == 1
            csv_path = joinpath(path, "reservoir_hq_$resloc.csv")
            @info(
                "Read a rating curve from CSV file `$csv_path`, for reservoir location `$resloc`"
            )
            @reset parameters.hq[i] = read_hq_csv(csv_path)
            @reset parameters.maxstorage[i] = maximum_storage(parameters, i)
        elseif outflowfunc[i] == 2 || outflowfunc[i] == 3
            @reset parameters.threshold[i] = threshold[i]
            @reset parameters.b[i] = b[i]
            @reset parameters.e[i] = e[i]
        elseif outflowfunc[i] == 4
            @reset parameters.demand[i] = demand[i]
            @reset parameters.maxrelease[i] = maxrelease[i]
            @reset parameters.maxstorage[i] = maxstorage[i]
            @reset parameters.targetfullfrac[i] = targetfullfrac[i]
            @reset parameters.targetminfrac[i] = targetminfrac[i]
        end

        if outflowfunc[i] == 3 && storfunc[i] != 1
            @warn(
                "For the modified puls approach (outflowfunc = 3) the storfunc should be 1"
            )
        end
    end

    return parameters, waterlevel
end

"Struct for storing reservoir model variables"
@with_kw struct ReservoirVariables
    # waterlevel H [m] of reservoir
    waterlevel::Vector{Float64}
    # average waterlevel H [m] of reservoir for model timestep Δt
    waterlevel_av::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # reservoir storage [m³]
    storage::Vector{Float64}
    # average reservoir storage [m³] for model timestep Δt
    storage_av::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # outflow from reservoir [m³ s⁻¹]
    outflow::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # average outflow from reservoir [m³ s⁻¹] for model timestep Δt
    outflow_av::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # average actual evaporation for reservoir area [mm Δt⁻¹]
    actevap::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # fraction full (of max storage) for rating curve type 4 [-]
    percfull::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # minimum (environmental) flow released from reservoir for rating curve type 4 [m³ s⁻¹]
    demandrelease::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
end

"Initialize reservoir model variables"
function ReservoirVariables(parameters::ReservoirParameters, waterlevel::Vector{Float64})
    (; storfunc, area, sh) = parameters
    variables = ReservoirVariables(;
        waterlevel,
        storage = initialize_storage(storfunc, area, waterlevel, sh),
    )
    return variables
end

"Struct for storing reservoir model boundary conditions"
@with_kw struct ReservoirBC
    inflow::Vector{Float64}               # inflow into reservoir [m³ s⁻¹] for model timestep Δt
    precipitation::Vector{Float64}        # average precipitation for reservoir area [mm Δt⁻¹]
    evaporation::Vector{Float64}          # average potential evaporation for reservoir area [mm Δt⁻¹]
end

"Initialize reservoir model boundary conditions"
function ReservoirBC(n::Int)
    bc = ReservoirBC(;
        inflow = fill(MISSING_VALUE, n),
        precipitation = fill(MISSING_VALUE, n),
        evaporation = fill(MISSING_VALUE, n),
    )
    return bc
end

"Reservoir model"
@with_kw struct Reservoir
    boundary_conditions::ReservoirBC
    parameters::ReservoirParameters
    variables::ReservoirVariables
end

"Initialize reservoir model `SimpleReservoir`"
function Reservoir(dataset::NCDataset, config::Config, network::NetworkReservoir)
    parameters, waterlevel = ReservoirParameters(dataset, config, network)
    variables = ReservoirVariables(parameters, waterlevel)
    boundary_conditions = ReservoirBC(length(parameters.area))
    reservoir = Reservoir(; boundary_conditions, parameters, variables)

    return reservoir
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

"Determine the maximum storage for reservoirs with a rating curve of type 1"
function maximum_storage(parameters::ReservoirParameters, i::Int)
    (; storfunc, hq, sh, area) = parameters

    # maximum storage is based on the maximum water level (H) value in the H-Q table
    if storfunc[i] == 2
        maxstorage = interpolate_linear(maximum(hq[i].H), sh[i].H, sh[i].S)
    else
        maxstorage = area[i] * maximum(hq[i].H)
    end

    return maxstorage
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

"Update the column index of reservoir rating curve HQ data"
function update_index_hq!(reservoir::Reservoir, clock::Clock)
    (; outflowfunc, col_index_hq) = reservoir.parameters
    if 1 in outflowfunc
        col_index_hq[1] = julian_day(clock.time - clock.dt)
    end
    return nothing
end
update_index_hq!(reservoir, clock::Clock) = nothing

"Update reservoir with rating curve type (`ouflowfunc`) 4 for a single timestep"
function update_reservoir_simple(
    model::Reservoir,
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = model.parameters
    res_v = model.variables
    (; precipitation, actevap, inflow) = boundary_vars

    storage = res_v.storage[i] + (inflow * dt) + precipitation - actevap
    storage = max(storage, 0.0)

    percfull = storage / res_p.maxstorage[i]
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(percfull, res_p.targetminfrac[i], 1.0, 30.0)
    demandrelease = min(fac * res_p.demand[i] * dt, storage)
    storage = storage - demandrelease

    wantrel = max(0.0, storage - (res_p.maxstorage[i] * res_p.targetfullfrac[i]))
    # Assume extra maximum Q if spilling
    overflow_q = max((storage - res_p.maxstorage[i]), 0.0)
    torelease = min(wantrel, overflow_q + res_p.maxrelease[i] * dt - demandrelease)
    storage = storage - torelease
    outflow = torelease + demandrelease
    percfull = storage / res_p.maxstorage[i]

    outflow /= dt
    demandrelease /= dt

    return outflow, storage, demandrelease, percfull
end

"""
Update reservoir with rating curve type (`ouflowfunc`) 3 (Modified Puls approach) for a
single timestep.
"""
function update_reservoir_modified_puls(
    model::Reservoir,
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = model.parameters
    res_v = model.variables
    (; precipitation, actevap, inflow) = boundary_vars

    res_factor = res_p.area[i] / (dt * pow(res_p.b[i], 0.5))
    si_factor = (res_v.storage[i] + precipitation - actevap) / dt + inflow
    # Adjust si_factor for reservoir threshold != 0
    si_factor_adj = si_factor - res_p.area[i] * res_p.threshold[i] / dt
    # Calculate the new reservoir outflow/waterlevel/storage
    if si_factor_adj > 0.0
        quadratic_sol_term =
            -res_factor + pow((pow(res_factor, 2.0) + 4.0 * si_factor_adj), 0.5)
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
    return outflow, storage
end

"Update reservoir with rating curve type (`ouflowfunc`) 1 (HQ data) for a single timestep."
function update_reservoir_hq(
    model::Reservoir,
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = model.parameters
    res_v = model.variables
    (; precipitation, actevap, inflow) = boundary_vars

    storage_input = (res_v.storage[i] + precipitation - actevap) / dt + inflow
    outflow = interpolate_linear(
        res_v.waterlevel[i],
        res_p.hq[i].H,
        res_p.hq[i].Q[:, res_p.col_index_hq[1]],
    )
    outflow = min(outflow, storage_input)

    storage = (storage_input - outflow) * dt

    overflow = max(0.0, (storage - res_p.maxstorage[i]) / dt)
    storage = min(storage, res_p.maxstorage[i])
    outflow = outflow + overflow

    return outflow, storage
end

"Update reservoir with rating curve type (`ouflowfunc`) 2 (free weir) for a single timestep."
function update_reservoir_free_weir(
    model::Reservoir,
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = model.parameters
    res_v = model.variables
    (; precipitation, actevap, inflow) = boundary_vars

    lo = res_p.lower_reservoir_ind[i]
    has_lower_res = lo != 0
    diff_wl = has_lower_res ? res_v.waterlevel[i] - res_v.waterlevel[lo] : 0.0

    storage_input = (res_v.storage[i] + precipitation - actevap) / dt + inflow

    if diff_wl >= 0.0
        if res_v.waterlevel[i] > res_p.threshold[i]
            dh = res_v.waterlevel[i] - res_p.threshold[i]
            outflow = res_p.b[i] * pow(dh, res_p.e[i])
            maxflow = (dh * res_p.area[i]) / dt
            outflow = min(outflow, maxflow)
        else
            outflow = Float64(0)
        end
    else
        if res_v.waterlevel[lo] > res_p.threshold[i]
            dh = res_v.waterlevel[lo] - res_p.threshold[i]
            outflow = -1.0 * res_p.b[i] * pow(dh, res_p.e[i])
            maxflow = (dh * res_p.area[lo]) / dt
            outflow = max(outflow, -maxflow)
        else
            outflow = Float64(0)
        end
    end
    storage = (storage_input - outflow) * dt

    # update lower reservoir (linked reservoirs) in case flow from lower reservoir to upper reservoir occurs
    if diff_wl < 0.0
        lower_res_storage = res_v.storage[lo] + outflow * dt

        lower_res_waterlevel = if res_p.storfunc[lo] == 1
            res_v.waterlevel[lo] + (lower_res_storage - res_v.storage[lo]) / res_p.area[lo]
        else
            interpolate_linear(lower_res_storage, res_p.sh[lo].S, res_p.sh[lo].H)
        end

        # update values for the lower reservoir in place
        res_v.outflow[lo] = -outflow
        res_v.outflow_av[lo] += -outflow * dt
        res_v.storage[lo] = lower_res_storage
        res_v.waterlevel[lo] = lower_res_waterlevel
    end
    return outflow, storage
end

"""
Update a single reservoir at position `i`.

This is called from within the river routing scheme, therefore updating only for a single
element rather than all at once.
"""
function update!(
    model::Reservoir,
    i::Int,
    inflow::Float64,
    dt::Float64,
    dt_forcing::Float64,
)
    res_bc = model.boundary_conditions
    res_p = model.parameters
    res_v = model.variables

    # limit reservoir evaporation based on total available volume [m³]
    precipitation = 0.001 * res_bc.precipitation[i] * (dt / dt_forcing) * res_p.area[i]
    available_storage = res_v.storage[i] + inflow * dt + precipitation
    evap = 0.001 * res_bc.evaporation[i] * (dt / dt_forcing) * res_p.area[i]
    actevap = min(available_storage, evap) # [m³/dt]

    boundary_vars = (; precipitation, actevap, inflow)

    if res_p.outflowfunc[i] == 1
        outflow, storage = update_reservoir_hq(model, i, boundary_vars, dt)
    elseif res_p.outflowfunc[i] == 2
        outflow, storage = update_reservoir_free_weir(model, i, boundary_vars, dt)
    elseif res_p.outflowfunc[i] == 3
        outflow, storage = update_reservoir_modified_puls(model, i, boundary_vars, dt)
    elseif res_p.outflowfunc[i] == 4
        outflow, storage, demandrelease, percfull =
            update_reservoir_simple(model, i, boundary_vars, dt)
    end

    waterlevel = if res_p.storfunc[i] == 1
        res_v.waterlevel[i] + (storage - res_v.storage[i]) / res_p.area[i]
    else
        interpolate_linear(storage, res_p.sh[i].S, res_p.sh[i].H)
    end

    # update values in place
    # instantaneous variables
    res_v.storage[i] = storage
    res_v.waterlevel[i] = waterlevel
    res_v.outflow[i] = outflow
    if res_p.outflowfunc[i] == 4
        res_v.demandrelease[i] = demandrelease
        res_v.percfull[i] = percfull
    end
    # average variables (here accumulated for model timestep Δt)
    res_bc.inflow[i] += inflow * dt
    res_v.outflow_av[i] += outflow * dt
    res_v.storage_av[i] += storage * dt
    res_v.waterlevel_av[i] += waterlevel * dt
    res_v.actevap[i] += 1000.0 * (actevap / res_p.area[i])

    return nothing
end