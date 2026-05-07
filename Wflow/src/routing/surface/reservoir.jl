# Explicit integers added because the enumerators are 0-based by default
@enumx ReservoirProfileType linear = 1 interpolation = 2
@enumx ReservoirOutflowType rating_curve = 1 free_weir = 2 modified_puls = 3 simple = 4

"Struct for storing reservoir model parameters"
@with_kw struct ReservoirParameters
    # reservoir location id
    id::Vector{Int}
    # type of reservoir storage curve, 1: S = AH, 2: S = f(H) from reservoir data and
    # interpolation
    storage_curve_type::Vector{ReservoirProfileType.T}
    # type of reservoir rating curve, 1: Q = f(H) from reservoir data and interpolation, 2:
    # General Q = b(H - H₀)ᵉ, 3: Case of Puls Approach Q = b(H - H₀)², 4: Simple reservoir
    outflow_curve_type::Vector{ReservoirOutflowType.T}
    # reservoir area [m²]
    area::Vector{Float64}
    # index of lower reservoir (linked reservoirs) [-]
    lower_reservoir_ind::Vector{Int} = zeros(Int, length(area))
    # reservoir maximum storage for rating curve types 1 and 4 [m³]
    maximum_storage::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # water level threshold H₀ [m] below that level outflow is zero
    threshold::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # rating curve coefficient [m3/2 s-1] (if rating_curve_exponent=3/2)
    rating_curve_coefficient::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # rating curve exponent [-]
    rating_curve_exponent::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # data for storage curve
    storage_height_curve::Vector{Union{SH, Missing}} =
        Vector{Union{SH, Missing}}(missing, length(area))
    # data for rating curve
    height_discharge_curve::Vector{Union{HQ, Missing}} =
        Vector{Union{HQ, Missing}}(missing, length(area))
    # column index of rating curve data height_discharge_curve
    col_index_hq::Vector{Int} = [1]
    # maximum amount that can be released if below spillway for rating curve type 4 [m³ s⁻¹]
    maximum_release::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # minimum (environmental) flow requirement downstream of the reservoir for rating curve
    # type 4 [m³ s⁻¹]
    demand::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # target minimum full fraction (of max storage) for rating curve type 4 [-]
    target_minimum_fraction::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # target fraction full (of max storage) for rating curve type 4 [-]
    target_full_fraction::Vector{Float64} = fill(MISSING_VALUE, length(area))
end

"Initialize reservoir model parameters"
function ReservoirParameters(dataset::NCDataset, config::Config, network::NetworkReservoir)
    (; indices_outlet) = network

    area = ncread(dataset, config, "reservoir_surface__area", Routing; sel = indices_outlet)
    waterlevel = ncread(
        dataset,
        config,
        "reservoir_water_surface__initial_elevation",
        Routing;
        sel = indices_outlet,
    )
    storage_curve_type = ncread(
        dataset,
        config,
        "reservoir_water__storage_curve_type_count",
        Routing;
        sel = indices_outlet,
    )
    storage_curve_type = to_enumx.(ReservoirProfileType.T, storage_curve_type)
    outflow_curve_type = ncread(
        dataset,
        config,
        "reservoir_water__rating_curve_type_count",
        Routing;
        sel = indices_outlet,
    )
    outflow_curve_type = to_enumx.(ReservoirOutflowType.T, outflow_curve_type)
    linked_reslocs = ncread(
        dataset,
        config,
        "reservoir_lower_location__count",
        Routing;
        sel = indices_outlet,
    )

    n_reservoirs = length(area)
    reslocs =
        ncread(dataset, config, "reservoir_location__count", Routing; sel = indices_outlet)
    @info "Read `$n_reservoirs` reservoir locations."

    parameters =
        ReservoirParameters(; id = reslocs, area, outflow_curve_type, storage_curve_type)

    if ReservoirOutflowType.free_weir in outflow_curve_type ||
       ReservoirOutflowType.modified_puls in outflow_curve_type
        threshold = ncread(
            dataset,
            config,
            "reservoir_water_flow_threshold_level__elevation",
            Routing;
            sel = indices_outlet,
        )
        rating_curve_coefficient = ncread(
            dataset,
            config,
            "reservoir_water__rating_curve_coefficient",
            Routing;
            sel = indices_outlet,
        )
        rating_curve_exponent = ncread(
            dataset,
            config,
            "reservoir_water__rating_curve_exponent",
            Routing;
            sel = indices_outlet,
        )
    end
    if ReservoirOutflowType.simple in outflow_curve_type
        demand = ncread(
            dataset,
            config,
            "reservoir_water_demand__required_downstream_volume_flow_rate",
            Routing;
            sel = indices_outlet,
        )
        maximum_release = ncread(
            dataset,
            config,
            "reservoir_water_release_below_spillway__max_volume_flow_rate",
            Routing;
            sel = indices_outlet,
        )
        maximum_storage = ncread(
            dataset,
            config,
            "reservoir_water__max_volume",
            Routing;
            sel = indices_outlet,
        )
        target_full_fraction = ncread(
            dataset,
            config,
            "reservoir_water__target_full_volume_fraction",
            Routing;
            sel = indices_outlet,
        )
        target_minimum_fraction = ncread(
            dataset,
            config,
            "reservoir_water__target_min_volume_fraction",
            Routing;
            sel = indices_outlet,
        )
    end

    # reservoir CSV parameter files are expected in the same directory as path_static
    path = dirname(input_path(config, config.input.path_static))
    for i in 1:n_reservoirs
        resloc = reslocs[i]
        if linked_reslocs[i] > 0
            @reset parameters.lower_reservoir_ind[i] =
                only(findall(x -> x == linked_reslocs[i], reslocs))
        end

        if storage_curve_type[i] == ReservoirProfileType.interpolation
            csv_path = joinpath(path, "reservoir_sh_$resloc.csv")
            @info(
                "Read a storage curve from CSV file `$csv_path`, for reservoir location `$resloc`"
            )
            @reset parameters.storage_height_curve[i] = read_sh_csv(csv_path)
        end

        if outflow_curve_type[i] == ReservoirOutflowType.rating_curve
            csv_path = joinpath(path, "reservoir_hq_$resloc.csv")
            @info(
                "Read a rating curve from CSV file `$csv_path`, for reservoir location `$resloc`"
            )
            parameters.height_discharge_curve[i] = read_hq_csv(csv_path)
            parameters.maximum_storage[i] = maximum_storage(parameters, i)
        elseif outflow_curve_type[i] == ReservoirOutflowType.free_weir ||
               outflow_curve_type[i] == ReservoirOutflowType.modified_puls
            parameters.threshold[i] = threshold[i]
            parameters.rating_curve_coefficient[i] = rating_curve_coefficient[i]
            parameters.rating_curve_exponent[i] = rating_curve_exponent[i]
        elseif outflow_curve_type[i] == ReservoirOutflowType.simple
            parameters.demand[i] = demand[i]
            parameters.maximum_release[i] = maximum_release[i]
            parameters.maximum_storage[i] = maximum_storage[i]
            parameters.target_full_fraction[i] = target_full_fraction[i]
            parameters.target_minimum_fraction[i] = target_minimum_fraction[i]
        end

        if outflow_curve_type[i] == ReservoirOutflowType.modified_puls &&
           storage_curve_type[i] != ReservoirProfileType.linear
            @warn(
                "For the modified puls approach (outflow_curve_type = 3) the storage_curve_type should be 1"
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
    actual_evaporation::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
end

"Initialize reservoir model variables"
function ReservoirVariables(
    dataset::NCDataset,
    config::Config,
    network::NetworkReservoir,
    parameters::ReservoirParameters,
    waterlevel::Vector{Float64},
)
    (; storage_curve_type, area, storage_height_curve) = parameters
    (; indices_outlet) = network
    outflow_obs = ncread(
        dataset,
        config,
        "reservoir_water__outgoing_observed_volume_flow_rate",
        Routing;
        sel = indices_outlet,
    )
    variables = ReservoirVariables(;
        waterlevel,
        storage = initialize_storage(
            storage_curve_type,
            area,
            waterlevel,
            storage_height_curve,
        ),
        outflow_obs,
    )
    return variables
end

"Struct for storing reservoir model boundary conditions"
@with_kw struct ReservoirBC
    n::Int
    inflow_subsurface::Vector{Float64} = fill(MISSING_VALUE, n)    # inflow from subsurface flow into reservoir [m³ s⁻¹]
    inflow_overland::Vector{Float64} = fill(MISSING_VALUE, n)      # inflow from overland flow into reservoir [m³ s⁻¹]
    inflow::Vector{Float64} = fill(MISSING_VALUE, n)               # total inflow into reservoir [m³ s⁻¹] for model timestep Δt
    external_inflow::Vector{Float64}                               # external inflow (abstraction/supply/demand) [m³ s⁻¹]
    actual_external_abstraction_av::Vector{Float64} = zeros(n)  # actual abstraction from external negative inflow [m³ s⁻¹]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)        # average precipitation for reservoir area [mm Δt⁻¹]
    evaporation::Vector{Float64} = fill(MISSING_VALUE, n)          # average potential evaporation for reservoir area [mm Δt⁻¹]
end

"Initialize reservoir model boundary conditions"
function ReservoirBC(dataset::NCDataset, config::Config, network::NetworkReservoir)
    (; indices_outlet) = network
    external_inflow = ncread(
        dataset,
        config,
        "reservoir_water__external_inflow_volume_flow_rate",
        Routing;
        sel = indices_outlet,
    )
    n = length(indices_outlet)
    bc = ReservoirBC(; n, external_inflow)
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
    storage_curve_type::ReservoirProfileType.T,
    area::Float64,
    storage::Float64,
    storage_height_curve::Union{SH, Missing},
)
    if storage_curve_type == ReservoirProfileType.linear
        waterlevel = storage / area
    else # storage_curve_type == ReservoirProfileType.interpolation
        waterlevel =
            interpolate_linear(storage, storage_height_curve.S, storage_height_curve.H)
    end
    return waterlevel
end

"Determine the maximum storage for reservoirs with a rating curve of type 1"
function maximum_storage(parameters::ReservoirParameters, i::Int)
    (; storage_curve_type, height_discharge_curve, storage_height_curve, area) = parameters

    # maximum storage is based on the maximum water level (H) value in the H-Q table
    if storage_curve_type[i] == ReservoirProfileType.interpolation
        max_storage = interpolate_linear(
            maximum(height_discharge_curve[i].H),
            storage_height_curve[i].H,
            storage_height_curve[i].S,
        )
    else # storage_curve_type[i] == ReservoirProfileType.linear
        max_storage = area[i] * maximum(height_discharge_curve[i].H)
    end

    return max_storage
end

"Determine the initial storage depending on the storage function"
function initialize_storage(
    storage_curve_type::Vector{ReservoirProfileType.T},
    area::Vector{Float64},
    waterlevel::Vector{Float64},
    storage_height_curve::Vector{Union{SH, Missing}},
)
    storage = similar(area)
    for i in eachindex(storage)
        if storage_curve_type[i] == ReservoirProfileType.linear
            storage[i] = area[i] * waterlevel[i]
        else # storage_curve_type[i] == ReservoirProfileType.interpolation
            storage[i] = interpolate_linear(
                waterlevel[i],
                storage_height_curve[i].H,
                storage_height_curve[i].S,
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
    (; outflow_curve_type, col_index_hq) = reservoir_model.parameters
    if ReservoirOutflowType.rating_curve in outflow_curve_type
        col_index_hq[1] = julian_day(clock.time - clock.dt)
    end
    return nothing
end
update_index_hq!(reservoir_model::Any, clock::Clock) = nothing

"Update reservoir with rating curve type (`ouflowfunc`) 4 for a single timestep"
function update_reservoir_simple(
    reservoir_model::ReservoirModel,
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables
    (; precipitation, actual_evaporation, inflow) = boundary_vars

    storage = res_v.storage[i] + (inflow * dt) + precipitation - actual_evaporation
    storage = max(storage, 0.0)

    percfull = storage / res_p.maximum_storage[i]
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(percfull, res_p.target_minimum_fraction[i], 1.0, 30.0)
    demandrelease = min(fac * res_p.demand[i] * dt, storage)
    storage -= demandrelease

    wantrel = max(0.0, storage - (res_p.maximum_storage[i] * res_p.target_full_fraction[i]))
    # Assume extra maximum Q if spilling
    overflow_q = max((storage - res_p.maximum_storage[i]), 0.0)
    torelease = min(wantrel, overflow_q + res_p.maximum_release[i] * dt - demandrelease)
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
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables
    (; precipitation, actual_evaporation, inflow) = boundary_vars

    res_factor = res_p.area[i] / (dt * pow(res_p.rating_curve_coefficient[i], 0.5))
    si_factor =
        max((res_v.storage[i] + precipitation - actual_evaporation) / dt + inflow, 0.0) # prevent negative values
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
    reservoir_model::ReservoirModel,
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables
    (; precipitation, actual_evaporation, inflow) = boundary_vars

    storage_input =
        max((res_v.storage[i] + precipitation - actual_evaporation) / dt + inflow, 0.0) # prevent negative values

    outflow = interpolate_linear(
        res_v.waterlevel[i],
        res_p.height_discharge_curve[i].H,
        res_p.height_discharge_curve[i].Q[:, res_p.col_index_hq[1]],
    )
    outflow = min(outflow, storage_input)

    storage = (storage_input - outflow) * dt

    overflow = max(0.0, (storage - res_p.maximum_storage[i]) / dt)
    storage = min(storage, res_p.maximum_storage[i])
    outflow += overflow

    return outflow, storage
end

"Update reservoir with rating curve type (`ouflowfunc`) 2 (free weir) for a single timestep."
function update_reservoir_free_weir(
    reservoir_model::ReservoirModel,
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables
    (; precipitation, actual_evaporation, inflow) = boundary_vars

    lo = res_p.lower_reservoir_ind[i]
    has_lower_res = lo != 0
    diff_wl = has_lower_res ? res_v.waterlevel[i] - res_v.waterlevel[lo] : 0.0

    storage_input =
        max((res_v.storage[i] + precipitation - actual_evaporation) / dt + inflow, 0.0) # prevent negative values

    if diff_wl >= 0.0
        if res_v.waterlevel[i] > res_p.threshold[i]
            dh = res_v.waterlevel[i] - res_p.threshold[i]
            outflow =
                res_p.rating_curve_coefficient[i] * pow(dh, res_p.rating_curve_exponent[i])
            maxflow = (dh * res_p.area[i]) / dt
            outflow = min(outflow, maxflow)
        else
            outflow = Float64(0)
        end
    else
        if res_v.waterlevel[lo] > res_p.threshold[i]
            dh = res_v.waterlevel[lo] - res_p.threshold[i]
            outflow =
                -1.0 *
                res_p.rating_curve_coefficient[i] *
                pow(dh, res_p.rating_curve_exponent[i])
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

        lower_res_waterlevel =
            if res_p.storage_curve_type[lo] == ReservoirProfileType.linear
                res_v.waterlevel[lo] +
                (lower_res_storage - res_v.storage[lo]) / res_p.area[lo]
            else # res_p.storage_curve_type[lo] == ReservoirProfileType.interpolation
                interpolate_linear(
                    lower_res_storage,
                    res_p.storage_height_curve[lo].S,
                    res_p.storage_height_curve[lo].H,
                )
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
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    res_v = reservoir_model.variables
    (; precipitation, actual_evaporation, inflow) = boundary_vars

    storage_input =
        max((res_v.storage[i] + precipitation - actual_evaporation) / dt + inflow, 0.0) # prevent negative values
    outflow = min(res_v.outflow_obs[i], storage_input)
    storage = (storage_input - outflow) * dt
    return outflow, storage
end

"""
Update a single reservoir at position `i`.

This is called from within the river routing scheme, therefore updating only for a single
element rather than all at once.
"""
function update_reservoir_model!(
    reservoir_model::ReservoirModel,
    i::Int,
    inflow::Float64,
    dt::Float64,
    dt_forcing::Float64,
)
    res_bc = reservoir_model.boundary_conditions
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables

    # limit reservoir evaporation based on total available volume [m³]
    precipitation = 0.001 * res_bc.precipitation[i] * (dt / dt_forcing) * res_p.area[i]
    available_storage = res_v.storage[i] + inflow * dt + precipitation
    evap = 0.001 * res_bc.evaporation[i] * (dt / dt_forcing) * res_p.area[i]
    actual_evaporation = min(available_storage, evap) # [m³/dt]

    boundary_vars = (; precipitation, actual_evaporation, inflow)
    update_reservoir_args = (reservoir_model, i, boundary_vars, dt)

    if !isnan(res_v.outflow_obs[i])
        outflow, storage = update_reservoir_outflow_obs(update_reservoir_args...)
    elseif res_p.outflow_curve_type[i] == ReservoirOutflowType.rating_curve
        outflow, storage = update_reservoir_hq(update_reservoir_args...)
    elseif res_p.outflow_curve_type[i] == ReservoirOutflowType.free_weir
        outflow, storage = update_reservoir_free_weir(update_reservoir_args...)
    elseif res_p.outflow_curve_type[i] == ReservoirOutflowType.modified_puls
        outflow, storage = update_reservoir_modified_puls(update_reservoir_args...)
    elseif res_p.outflow_curve_type[i] == ReservoirOutflowType.simple
        outflow, storage = update_reservoir_simple(update_reservoir_args...)
    end

    waterlevel = if res_p.storage_curve_type[i] == ReservoirProfileType.linear
        res_v.waterlevel[i] + (storage - res_v.storage[i]) / res_p.area[i]
    else # res_p.storage_curve_type[i] == ReservoirProfileType.interpolation
        interpolate_linear(
            storage,
            res_p.storage_height_curve[i].S,
            res_p.storage_height_curve[i].H,
        )
    end

    # update values in place
    # instantaneous variables
    res_v.storage[i] = storage
    res_v.waterlevel[i] = waterlevel
    res_v.outflow[i] = outflow

    # average variables (here accumulated for model timestep Δt)
    res_bc.inflow[i] += inflow * dt
    res_v.outflow_av[i] += outflow * dt
    res_v.actual_evaporation[i] += 1000.0 * (actual_evaporation / res_p.area[i])

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
    reservoir_model::Union{ReservoirModel, Nothing},
    config::Config,
)
    par = "reservoir_water__outgoing_observed_volume_flow_rate"
    check =
        !isnothing(reservoir_model) &&
        (haskey(config.input.forcing, par) || haskey(config.input.cyclic, par))
    return check
end
