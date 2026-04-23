# Explicit integers added because the enumerators are 0-based by default
@enumx ReservoirProfileType linear = 1 interpolation = 2
@enumx ReservoirOutflowType rating_curve = 1 free_weir = 2 modified_puls = 3 simple = 4

"Struct for storing reservoir model parameters"
@with_kw struct ReservoirParameters
    # reservoir location id
    id::Vector{Int}
    # type of reservoir storage curve, 1: S = AH, 2: S = f(H) from reservoir data and
    # interpolation
    storfunc::Vector{ReservoirProfileType.T}
    # type of reservoir rating curve, 1: Q = f(H) from reservoir data and interpolation, 2:
    # General Q = b(H - H₀)ᵉ, 3: Case of Puls Approach Q = b(H - H₀)², 4: Simple reservoir
    outflowfunc::Vector{ReservoirOutflowType.T}
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
    sh::Vector{Union{SH,Missing}} = Vector{Union{SH,Missing}}(missing, length(area))
    # data for rating curve
    hq::Vector{Union{HQ,Missing}} = Vector{Union{HQ,Missing}}(missing, length(area))
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
    (; outlet_indices_2d) = network

    area =
        ncread(dataset, config, "reservoir_surface__area", Routing; sel=outlet_indices_2d)
    waterlevel = ncread(
        dataset,
        config,
        "reservoir_water_surface__initial_elevation",
        Routing;
        sel=outlet_indices_2d,
    )
    storfunc = ncread(
        dataset,
        config,
        "reservoir_water__storage_curve_type_count",
        Routing;
        sel=outlet_indices_2d,
    )
    storfunc = to_enumx.(ReservoirProfileType.T, storfunc)
    outflowfunc = ncread(
        dataset,
        config,
        "reservoir_water__rating_curve_type_count",
        Routing;
        sel=outlet_indices_2d,
    )
    outflowfunc = to_enumx.(ReservoirOutflowType.T, outflowfunc)
    linked_reslocs = ncread(
        dataset,
        config,
        "reservoir_lower_location__count",
        Routing;
        sel=outlet_indices_2d,
    )

    n_reservoirs = length(area)
    reslocs = ncread(
        dataset,
        config,
        "reservoir_location__count",
        Routing;
        sel=outlet_indices_2d,
    )
    @info "Read `$n_reservoirs` reservoir locations."

    parameters = ReservoirParameters(; id=reslocs, area, outflowfunc, storfunc)

    if ReservoirOutflowType.free_weir in outflowfunc ||
       ReservoirOutflowType.modified_puls in outflowfunc
        threshold = ncread(
            dataset,
            config,
            "reservoir_water_flow_threshold_level__elevation",
            Routing;
            sel=outlet_indices_2d,
        )
        b = ncread(
            dataset,
            config,
            "reservoir_water__rating_curve_coefficient",
            Routing;
            sel=outlet_indices_2d,
        )
        e = ncread(
            dataset,
            config,
            "reservoir_water__rating_curve_exponent",
            Routing;
            sel=outlet_indices_2d,
        )
    end
    if ReservoirOutflowType.simple in outflowfunc
        demand = ncread(
            dataset,
            config,
            "reservoir_water_demand__required_downstream_volume_flow_rate",
            Routing;
            sel=outlet_indices_2d,
        )
        maxrelease = ncread(
            dataset,
            config,
            "reservoir_water_release_below_spillway__max_volume_flow_rate",
            Routing;
            sel=outlet_indices_2d,
        )
        maxstorage = ncread(
            dataset,
            config,
            "reservoir_water__max_volume",
            Routing;
            sel=outlet_indices_2d,
        )
        targetfullfrac = ncread(
            dataset,
            config,
            "reservoir_water__target_full_volume_fraction",
            Routing;
            sel=outlet_indices_2d,
        )
        targetminfrac = ncread(
            dataset,
            config,
            "reservoir_water__target_min_volume_fraction",
            Routing;
            sel=outlet_indices_2d,
        )
    end

    # reservoir CSV parameter files are expected in the same directory as path_static
    path = dirname(input_path(config, config.input.path_static))
    for reservoir_idx in 1:n_reservoirs
        resloc = reslocs[reservoir_idx]
        if linked_reslocs[reservoir_idx] > 0
            @reset parameters.lower_reservoir_ind[reservoir_idx] =
                only(findall(x -> x == linked_reslocs[reservoir_idx], reslocs))
        end

        if storfunc[reservoir_idx] == ReservoirProfileType.interpolation
            csv_path = joinpath(path, "reservoir_sh_$resloc.csv")
            @info(
                "Read a storage curve from CSV file `$csv_path`, for reservoir location `$resloc`"
            )
            @reset parameters.sh[reservoir_idx] = read_sh_csv(csv_path)
        end

        if outflowfunc[reservoir_idx] == ReservoirOutflowType.rating_curve
            csv_path = joinpath(path, "reservoir_hq_$resloc.csv")
            @info(
                "Read a rating curve from CSV file `$csv_path`, for reservoir location `$resloc`"
            )
            parameters.hq[reservoir_idx] = read_hq_csv(csv_path)
            parameters.maxstorage[reservoir_idx] =
                maximum_storage(parameters, reservoir_idx)
        elseif outflowfunc[reservoir_idx] == ReservoirOutflowType.free_weir ||
               outflowfunc[reservoir_idx] == ReservoirOutflowType.modified_puls
            parameters.threshold[reservoir_idx] = threshold[reservoir_idx]
            parameters.b[reservoir_idx] = b[reservoir_idx]
            parameters.e[reservoir_idx] = e[reservoir_idx]
        elseif outflowfunc[reservoir_idx] == ReservoirOutflowType.simple
            parameters.demand[reservoir_idx] = demand[reservoir_idx]
            parameters.maxrelease[reservoir_idx] = maxrelease[reservoir_idx]
            parameters.maxstorage[reservoir_idx] = maxstorage[reservoir_idx]
            parameters.targetfullfrac[reservoir_idx] = targetfullfrac[reservoir_idx]
            parameters.targetminfrac[reservoir_idx] = targetminfrac[reservoir_idx]
        end

        if outflowfunc[reservoir_idx] == ReservoirOutflowType.modified_puls &&
           storfunc[reservoir_idx] != ReservoirProfileType.linear
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
    # reservoir storage [m³]
    storage::Vector{Float64}
    # outflow from reservoir [m³ s⁻¹]
    outflow::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # average outflow from reservoir [m³ s⁻¹] for model timestep Δt
    outflow_av::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # observed outflow from reservoir [m³ s⁻¹]
    outflow_obs::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # average actual evaporation for reservoir area [mm Δt⁻¹]
    actevap::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
end

"Initialize reservoir model variables"
function ReservoirVariables(
    dataset::NCDataset,
    config::Config,
    network::NetworkReservoir,
    parameters::ReservoirParameters,
    waterlevel::Vector{Float64},
)
    (; storfunc, area, sh) = parameters
    (; outlet_indices_2d) = network
    outflow_obs = ncread(
        dataset,
        config,
        "reservoir_water__outgoing_observed_volume_flow_rate",
        Routing;
        sel=outlet_indices_2d,
    )
    variables = ReservoirVariables(;
        waterlevel,
        storage=initialize_storage(storfunc, area, waterlevel, sh),
        outflow_obs,
    )
    return variables
end

"Struct for storing reservoir model boundary conditions"
@with_kw struct ReservoirBC
    n_reservoirs::Int
    inflow_subsurface::Vector{Float64} = fill(MISSING_VALUE, n_reservoirs)    # inflow from subsurface flow into reservoir [m³ s⁻¹]
    inflow_overland::Vector{Float64} = fill(MISSING_VALUE, n_reservoirs)      # inflow from overland flow into reservoir [m³ s⁻¹]
    inflow::Vector{Float64} = fill(MISSING_VALUE, n_reservoirs)               # total inflow into reservoir [m³ s⁻¹] for model timestep Δt
    external_inflow::Vector{Float64}                               # external inflow (abstraction/supply/demand) [m³ s⁻¹]
    actual_external_abstraction_av::Vector{Float64} = zeros(n_reservoirs)  # actual abstraction from external negative inflow [m³ s⁻¹]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n_reservoirs)        # average precipitation for reservoir area [mm Δt⁻¹]
    evaporation::Vector{Float64} = fill(MISSING_VALUE, n_reservoirs)          # average potential evaporation for reservoir area [mm Δt⁻¹]
end

"Initialize reservoir model boundary conditions"
function ReservoirBC(dataset::NCDataset, config::Config, network::NetworkReservoir)
    (; outlet_indices_2d) = network
    external_inflow = ncread(
        dataset,
        config,
        "reservoir_water__external_inflow_volume_flow_rate",
        Routing;
        sel=outlet_indices_2d,
    )
    n_reservoirs = length(outlet_indices_2d)
    bc = ReservoirBC(; n_reservoirs, external_inflow)
    return bc
end

"Reservoir model"
@with_kw struct ReservoirModel
    boundary_conditions::ReservoirBC
    parameters::ReservoirParameters
    variables::ReservoirVariables
end

"Initialize reservoir model `SimpleReservoir`"
function ReservoirModel(dataset::NCDataset, config::Config, network::NetworkReservoir)
    parameters, waterlevel = ReservoirParameters(dataset, config, network)
    variables = ReservoirVariables(dataset, config, network, parameters, waterlevel)
    boundary_conditions = ReservoirBC(dataset, config, network)
    reservoir = ReservoirModel(; boundary_conditions, parameters, variables)

    return reservoir
end

"Determine the water level depending on the storage function"
function waterlevel(
    storfunc::ReservoirProfileType.T,
    area::Float64,
    storage::Float64,
    sh::Union{SH,Missing},
)
    if storfunc == ReservoirProfileType.linear
        waterlevel = storage / area
    else # storfunc == ReservoirProfileType.interpolation
        waterlevel = interpolate_linear(storage, sh.S, sh.H)
    end
    return waterlevel
end

"Determine the maximum storage for reservoirs with a rating curve of type 1"
function maximum_storage(parameters::ReservoirParameters, reservoir_idx::Int)
    (; storfunc, hq, sh, area) = parameters

    # maximum storage is based on the maximum water level (H) value in the H-Q table
    if storfunc[reservoir_idx] == ReservoirProfileType.interpolation
        maxstorage = interpolate_linear(
            maximum(hq[reservoir_idx].H),
            sh[reservoir_idx].H,
            sh[reservoir_idx].S,
        )
    else # storfunc[reservoir_idx] == ReservoirProfileType.linear
        maxstorage = area[reservoir_idx] * maximum(hq[reservoir_idx].H)
    end

    return maxstorage
end

"Determine the initial storage depending on the storage function"
function initialize_storage(
    storfunc::Vector{ReservoirProfileType.T},
    area::Vector{Float64},
    waterlevel::Vector{Float64},
    sh::Vector{Union{SH,Missing}},
)
    storage = similar(area)
    for reservoir_idx in eachindex(storage)
        if storfunc[reservoir_idx] == ReservoirProfileType.linear
            storage[reservoir_idx] = area[reservoir_idx] * waterlevel[reservoir_idx]
        else # storfunc[reservoir_idx] == ReservoirProfileType.interpolation
            storage[reservoir_idx] = interpolate_linear(
                waterlevel[reservoir_idx],
                sh[reservoir_idx].H,
                sh[reservoir_idx].S,
            )
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
function update_index_hq!(reservoir_model::ReservoirModel, clock::Clock)
    (; outflowfunc, col_index_hq) = reservoir_model.parameters
    if ReservoirOutflowType.rating_curve in outflowfunc
        col_index_hq[1] = julian_day(clock.time - clock.dt)
    end
    return nothing
end
update_index_hq!(reservoir_model::Any, clock::Clock) = nothing

"Update reservoir with rating curve type (`ouflowfunc`) 4 for a single timestep"
function update_reservoir_simple(
    reservoir_model::ReservoirModel,
    reservoir_idx::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables
    (; precipitation, actevap, inflow) = boundary_vars

    storage = res_v.storage[reservoir_idx] + (inflow * dt) + precipitation - actevap
    storage = max(storage, 0.0)

    percfull = storage / res_p.maxstorage[reservoir_idx]
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(percfull, res_p.targetminfrac[reservoir_idx], 1.0, 30.0)
    demandrelease = min(fac * res_p.demand[reservoir_idx] * dt, storage)
    storage -= demandrelease

    wantrel = max(
        0.0,
        storage - (res_p.maxstorage[reservoir_idx] * res_p.targetfullfrac[reservoir_idx]),
    )
    # Assume extra maximum Q if spilling
    overflow_q = max((storage - res_p.maxstorage[reservoir_idx]), 0.0)
    torelease =
        min(wantrel, overflow_q + res_p.maxrelease[reservoir_idx] * dt - demandrelease)
    storage -= torelease
    outflow = torelease + demandrelease
    outflow /= dt

    return outflow, storage
end

"""
Update reservoir with rating curve type (`ouflowfunc`) 3 (Modified Puls approach) for a
single timestep.
"""
function update_reservoir_modified_puls(
    reservoir_model::ReservoirModel,
    reservoir_idx::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables
    (; precipitation, actevap, inflow) = boundary_vars

    res_factor = res_p.area[reservoir_idx] / (dt * pow(res_p.b[reservoir_idx], 0.5))
    si_factor =
        max((res_v.storage[reservoir_idx] + precipitation - actevap) / dt + inflow, 0.0) # prevent negative values
    # Adjust si_factor for reservoir threshold != 0
    si_factor_adj =
        si_factor - res_p.area[reservoir_idx] * res_p.threshold[reservoir_idx] / dt
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
    reservoir_model::ReservoirModel,
    reservoir_idx::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables
    (; precipitation, actevap, inflow) = boundary_vars

    storage_input =
        max((res_v.storage[reservoir_idx] + precipitation - actevap) / dt + inflow, 0.0) # prevent negative values

    outflow = interpolate_linear(
        res_v.waterlevel[reservoir_idx],
        res_p.hq[reservoir_idx].H,
        res_p.hq[reservoir_idx].Q[:, res_p.col_index_hq[1]],
    )
    outflow = min(outflow, storage_input)

    storage = (storage_input - outflow) * dt

    overflow = max(0.0, (storage - res_p.maxstorage[reservoir_idx]) / dt)
    storage = min(storage, res_p.maxstorage[reservoir_idx])
    outflow += overflow

    return outflow, storage
end

"Update reservoir with rating curve type (`ouflowfunc`) 2 (free weir) for a single timestep."
function update_reservoir_free_weir(
    reservoir_model::ReservoirModel,
    reservoir_idx::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables
    (; precipitation, actevap, inflow) = boundary_vars

    lo = res_p.lower_reservoir_ind[reservoir_idx]
    has_lower_res = lo != 0
    diff_wl = has_lower_res ? res_v.waterlevel[reservoir_idx] - res_v.waterlevel[lo] : 0.0

    storage_input =
        max((res_v.storage[reservoir_idx] + precipitation - actevap) / dt + inflow, 0.0) # prevent negative values

    if diff_wl >= 0.0
        if res_v.waterlevel[reservoir_idx] > res_p.threshold[reservoir_idx]
            dh = res_v.waterlevel[reservoir_idx] - res_p.threshold[reservoir_idx]
            outflow = res_p.b[reservoir_idx] * pow(dh, res_p.e[reservoir_idx])
            maxflow = (dh * res_p.area[reservoir_idx]) / dt
            outflow = min(outflow, maxflow)
        else
            outflow = Float64(0)
        end
    else
        if res_v.waterlevel[lo] > res_p.threshold[reservoir_idx]
            dh = res_v.waterlevel[lo] - res_p.threshold[reservoir_idx]
            outflow = -1.0 * res_p.b[reservoir_idx] * pow(dh, res_p.e[reservoir_idx])
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

        lower_res_waterlevel = if res_p.storfunc[lo] == ReservoirProfileType.linear
            res_v.waterlevel[lo] + (lower_res_storage - res_v.storage[lo]) / res_p.area[lo]
        else # res_p.storfunc[lo] == ReservoirProfileType.interpolation
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

"Update reservoir using observed outflow for a single timestep."
function update_reservoir_outflow_obs(
    reservoir_model::ReservoirModel,
    reservoir_idx::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_v = reservoir_model.variables
    (; precipitation, actevap, inflow) = boundary_vars

    storage_input =
        max((res_v.storage[reservoir_idx] + precipitation - actevap) / dt + inflow, 0.0) # prevent negative values
    outflow = min(res_v.outflow_obs[reservoir_idx], storage_input)
    storage = (storage_input - outflow) * dt
    return outflow, storage
end

"""
Update a single reservoir at position `reservoir_idx`.

This is called from within the river routing scheme, therefore updating only for a single
element rather than all at once.
"""
function update_reservoir_model!(
    reservoir_model::ReservoirModel,
    reservoir_idx::Int,
    inflow::Float64,
    dt::Float64,
    dt_forcing::Float64,
)
    res_bc = reservoir_model.boundary_conditions
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables

    # limit reservoir evaporation based on total available volume [m³]
    precipitation =
        0.001 *
        res_bc.precipitation[reservoir_idx] *
        (dt / dt_forcing) *
        res_p.area[reservoir_idx]
    available_storage = res_v.storage[reservoir_idx] + inflow * dt + precipitation
    evap =
        0.001 *
        res_bc.evaporation[reservoir_idx] *
        (dt / dt_forcing) *
        res_p.area[reservoir_idx]
    actevap = min(available_storage, evap) # [m³/dt]

    boundary_vars = (; precipitation, actevap, inflow)
    update_reservoir_args = (reservoir_model, reservoir_idx, boundary_vars, dt)

    if !isnan(res_v.outflow_obs[reservoir_idx])
        outflow, storage = update_reservoir_outflow_obs(update_reservoir_args...)
    elseif res_p.outflowfunc[reservoir_idx] == ReservoirOutflowType.rating_curve
        outflow, storage = update_reservoir_hq(update_reservoir_args...)
    elseif res_p.outflowfunc[reservoir_idx] == ReservoirOutflowType.free_weir
        outflow, storage = update_reservoir_free_weir(update_reservoir_args...)
    elseif res_p.outflowfunc[reservoir_idx] == ReservoirOutflowType.modified_puls
        outflow, storage = update_reservoir_modified_puls(update_reservoir_args...)
    elseif res_p.outflowfunc[reservoir_idx] == ReservoirOutflowType.simple
        outflow, storage = update_reservoir_simple(update_reservoir_args...)
    end

    waterlevel = if res_p.storfunc[reservoir_idx] == ReservoirProfileType.linear
        res_v.waterlevel[reservoir_idx] +
        (storage - res_v.storage[reservoir_idx]) / res_p.area[reservoir_idx]
    else # res_p.storfunc[reservoir_idx] == ReservoirProfileType.interpolation
        interpolate_linear(storage, res_p.sh[reservoir_idx].S, res_p.sh[reservoir_idx].H)
    end

    # update values in place
    # instantaneous variables
    res_v.storage[reservoir_idx] = storage
    res_v.waterlevel[reservoir_idx] = waterlevel
    res_v.outflow[reservoir_idx] = outflow

    # average variables (here accumulated for model timestep Δt)
    res_bc.inflow[reservoir_idx] += inflow * dt
    res_v.outflow_av[reservoir_idx] += outflow * dt
    res_v.actevap[reservoir_idx] += 1000.0 * (actevap / res_p.area[reservoir_idx])

    return nothing
end

"Generate log message for using observed outflow at reservoir locations"
function log_message_observed_outflow(reservoir_model::ReservoirModel)
    not_nan = findall(x -> !isnan(x), reservoir_model.variables.outflow_obs)
    if isempty(not_nan)
        msg = "Observed outflow is not used for any reservoir location"
    else
        ids = reservoir_model.parameters.id[not_nan]
        msg = "Observed outflow is used for reservoir location ids $ids"
    end
    return msg
end

"Check if observed outflow is used for reservoirs"
function using_observed_outflow(
    reservoir_model::Union{ReservoirModel,Nothing},
    config::Config,
)
    par = "reservoir_water__outgoing_observed_volume_flow_rate"
    check =
        !isnothing(reservoir_model) &&
        (haskey(config.input.forcing, par) || haskey(config.input.cyclic, par))
    return check
end
