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
    # rating curve coefficient [m³⁻ᵉ s⁻¹]
    rating_curve_coefficient::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # rating curve exponent [-]
    rating_curve_exponent::Vector{Float64} = fill(MISSING_VALUE, length(area))
    # data for storage curve
    storage_waterlevel_curve::Vector{Union{SH, Missing}} =
        Vector{Union{SH, Missing}}(missing, length(area))
    # data for rating curve
    waterlevel_discharge_curve::Vector{Union{HQ, Missing}} =
        Vector{Union{HQ, Missing}}(missing, length(area))
    # column index of rating curve data hq
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
    for reservoir_idx in 1:n_reservoirs
        resloc = reslocs[reservoir_idx]
        if linked_reslocs[reservoir_idx] > 0
            parameters.lower_reservoir_ind[reservoir_idx] =
                only(findall(x -> x == linked_reslocs[reservoir_idx], reslocs))
        end

        if storage_curve_type[reservoir_idx] == ReservoirProfileType.interpolation
            csv_path = joinpath(path, "reservoir_sh_$resloc.csv")
            @info(
                "Read a storage curve from CSV file `$csv_path`, for reservoir location `$resloc`"
            )
            parameters.storage_waterlevel_curve[reservoir_idx] = read_sh_csv(csv_path)
        end

        if outflow_curve_type[reservoir_idx] == ReservoirOutflowType.rating_curve
            csv_path = joinpath(path, "reservoir_hq_$resloc.csv")
            @info(
                "Read a rating curve from CSV file `$csv_path`, for reservoir location `$resloc`"
            )
            parameters.waterlevel_discharge_curve[reservoir_idx] = read_hq_csv(csv_path)
            parameters.maximum_storage[reservoir_idx] =
                maximum_storage(parameters, reservoir_idx)
        elseif outflow_curve_type[reservoir_idx] == ReservoirOutflowType.free_weir ||
               outflow_curve_type[reservoir_idx] == ReservoirOutflowType.modified_puls
            parameters.threshold[reservoir_idx] = threshold[reservoir_idx]
            parameters.rating_curve_coefficient[reservoir_idx] =
                rating_curve_coefficient[reservoir_idx]
            parameters.rating_curve_exponent[reservoir_idx] =
                rating_curve_exponent[reservoir_idx]
        elseif outflow_curve_type[reservoir_idx] == ReservoirOutflowType.simple
            parameters.demand[reservoir_idx] = demand[reservoir_idx]
            parameters.maximum_release[reservoir_idx] = maximum_release[reservoir_idx]
            parameters.maximum_storage[reservoir_idx] = maximum_storage[reservoir_idx]
            parameters.target_full_fraction[reservoir_idx] =
                target_full_fraction[reservoir_idx]
            parameters.target_minimum_fraction[reservoir_idx] =
                target_minimum_fraction[reservoir_idx]
        end

        if outflow_curve_type[reservoir_idx] == ReservoirOutflowType.modified_puls &&
           storage_curve_type[reservoir_idx] != ReservoirProfileType.linear
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
    # Cumulative outflow from reservoir [m³] for model timestep dt
    outflow_cumulative::Vector{Float64} = zeros(length(waterlevel))
    # average outflow from reservoir [m³ s⁻¹] for model timestep dt
    outflow_average::Vector{Float64} = zeros(length(waterlevel))
    # observed outflow from reservoir [m³ s⁻¹]
    outflow_obs::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # cumulative actual evaporation for reservoir area [m]
    actevap_cumulative::Vector{Float64} = zeros(length(waterlevel))
    # average actual evaporation for reservoir area [m s⁻¹]
    actevap_average::Vector{Float64} = zeros(length(waterlevel))
end

"Initialize reservoir model variables"
function ReservoirVariables(
    dataset::NCDataset,
    config::Config,
    network::NetworkReservoir,
    parameters::ReservoirParameters,
    waterlevel::Vector{Float64},
)
    (; storage_curve_type, area, storage_waterlevel_curve) = parameters
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
            storage_waterlevel_curve,
        ),
        outflow_obs,
    )
    return variables
end

"Struct for storing reservoir model boundary conditions"
@with_kw struct ReservoirBC
    n::Int
    # inflow from subsurface flow into reservoir [m³ s⁻¹]
    inflow_subsurface::Vector{Float64} = fill(MISSING_VALUE, n)
    # inflow from overland flow into reservoir [m³ s⁻¹]
    inflow_overland::Vector{Float64} = fill(MISSING_VALUE, n)
    # cumulative inflow reservoir [m³] for model timestep dt
    inflow_cumulative::Vector{Float64} = zeros(n)
    # average inflow into reservoir [m³ s⁻¹] for model timestep dt
    inflow_average::Vector{Float64} = zeros(n)
    # external inflow (abstraction/supply/demand) [m³ s⁻¹]
    external_inflow::Vector{Float64} = zeros(n)
    # cumulative actual abstractoin from external negative flow [m³]
    actual_external_abstraction_cumulative::Vector{Float64} = zeros(n)
    # average actual abstraction from external negative inflow [m³ s⁻¹]
    actual_external_abstraction_average::Vector{Float64} = zeros(n)
    # average precipitation for reservoir area [m s⁻¹]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)
    # average potential evaporation for reservoir area [m s⁻¹]
    evaporation::Vector{Float64} = fill(MISSING_VALUE, n)
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
    storage_waterlevel_curve::Union{SH, Missing},
)
    if storage_curve_type == ReservoirProfileType.linear
        waterlevel = storage / area
    else # storage_curve_type == ReservoirProfileType.interpolation
        waterlevel = interpolate_linear(
            storage,
            storage_waterlevel_curve.S,
            storage_waterlevel_curve.H,
        )
    end
    return waterlevel
end

"Determine the maximum storage for reservoirs with a rating curve of type 1"
function maximum_storage(parameters::ReservoirParameters, i::Int)
    (; storage_curve_type, waterlevel_discharge_curve, storage_waterlevel_curve, area) =
        parameters

    # maximum storage is based on the maximum water level (H) value in the H-Q table
    if storage_curve_type[i] == ReservoirProfileType.interpolation
        maximum_storage = interpolate_linear(
            maximum(waterlevel_discharge_curve[i].H),
            storage_waterlevel_curve[i].H,
            storage_waterlevel_curve[i].S,
        )
    else # storage_curve_type[i] == ReservoirProfileType.linear
        maximum_storage = area[i] * maximum(waterlevel_discharge_curve[i].H)
    end

    return maximum_storage
end

"Determine the initial storage depending on the storage function"
function initialize_storage(
    storage_curve_type::Vector{ReservoirProfileType.T},
    area::Vector{Float64},
    waterlevel::Vector{Float64},
    storage_waterlevel_curve::Vector{Union{SH, Missing}},
)
    storage = similar(area)
    for reservoir_idx in eachindex(storage)
        if storage_curve_type[reservoir_idx] == ReservoirProfileType.linear
            storage[reservoir_idx] = area[reservoir_idx] * waterlevel[reservoir_idx]
        else # storage_curve_type[reservoir_idx] == ReservoirProfileType.interpolation
            storage[reservoir_idx] = interpolate_linear(
                waterlevel[reservoir_idx],
                storage_waterlevel_curve[reservoir_idx].H,
                storage_waterlevel_curve[reservoir_idx].S,
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
    (;
        maximum_storage,
        target_minimum_fraction,
        target_full_fraction,
        demand,
        maximum_release,
    ) = reservoir_model.parameters
    res_v = reservoir_model.variables
    (; precipitation, evaporation, inflow) = boundary_vars

    storage = res_v.storage[i] + (inflow + precipitation - evaporation) * dt
    storage = max(storage, 0.0)

    fill_fraction = storage / maximum_storage[i]
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(fill_fraction, target_minimum_fraction[i], 1.0, 30.0)
    demand_release = min(fac * demand[i], storage / dt)
    storage -= demand_release * dt
    release_wanted = max(0.0, (storage - maximum_storage[i] * target_full_fraction[i]) / dt)
    # Assume extra maximum Q if spilling
    overflow_q = max(0.0, (storage - maximum_storage[i]) / dt)
    release_realized = min(release_wanted, overflow_q + maximum_release[i] - demand_release)
    storage -= release_realized * dt
    outflow = release_realized + demand_release

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
    (; area, threshold, rating_curve_coefficient) = reservoir_model.parameters
    (; storage) = reservoir_model.variables
    (; precipitation, evaporation, inflow) = boundary_vars

    res_factor = area[i] / (dt * sqrt(rating_curve_coefficient[i]))
    si_factor = storage[i] / dt + precipitation - evaporation + inflow
    # Adjust si_factor for reservoir threshold != 0
    si_factor_adj = si_factor - area[i] * threshold[i] / dt
    # Calculate the new reservoir outflow/waterlevel/storage
    if si_factor_adj > 0.0
        quadratic_sol_term = -res_factor + sqrt((res_factor^2 + 4 * si_factor_adj))
        if quadratic_sol_term > 0.0
            outflow = 0.25 * quadratic_sol_term^2
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
    (; waterlevel_discharge_curve, col_index_hq, maximum_storage) =
        reservoir_model.parameters
    (; storage, waterlevel) = reservoir_model.variables
    (; precipitation, evaporation, inflow) = boundary_vars

    storage_input = storage[i] / dt + precipitation - evaporation + inflow
    outflow = interpolate_linear(
        waterlevel[i],
        waterlevel_discharge_curve[i].H,
        waterlevel_discharge_curve[i].Q[:, col_index_hq[1]],
    )
    outflow = min(outflow, storage_input)
    storage = (storage_input - outflow) * dt

    overflow = max(0.0, (storage - maximum_storage[i]) / dt)
    storage = min(storage, maximum_storage[i])
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
    (;
        threshold,
        rating_curve_coefficient,
        rating_curve_exponent,
        area,
        storage_curve_type,
        storage_waterlevel_curve,
        lower_reservoir_ind,
    ) = reservoir_model.parameters
    res_v = reservoir_model.variables
    (; waterlevel) = res_v
    (; precipitation, evaporation, inflow) = boundary_vars

    lo = lower_reservoir_ind[i]
    has_lower_res = (lo != 0)
    diff_wl = has_lower_res ? waterlevel[i] - waterlevel[lo] : 0.0

    storage_input = max(res_v.storage[i] / dt + precipitation - evaporation + inflow, 0.0)

    if diff_wl >= 0.0
        if res_v.waterlevel[i] > threshold[i]
            dh = waterlevel[i] - threshold[i]
            outflow = rating_curve_coefficient[i] * pow(dh, rating_curve_exponent[i])
            maxflow = dh * area[i] / dt
            outflow = min(outflow, maxflow)
        else
            outflow = 0.0
        end
    else
        if waterlevel[lo] > threshold[i]
            dh = waterlevel[lo] - threshold[i]
            outflow = -rating_curve_coefficient[i] * pow(dh, rating_curve_exponent[i])
            maxflow = dh * area[lo] / dt
            outflow = max(outflow, -maxflow)
        else
            outflow = 0.0
        end
    end
    storage = (storage_input - outflow) * dt

    # update lower reservoir (linked reservoirs) in case flow from lower reservoir to upper reservoir occurs
    if diff_wl < 0.0
        lower_res_storage = res_v.storage[lo] + outflow * dt
        lower_res_waterlevel = if storage_curve_type[lo] == ReservoirProfileType.linear
            waterlevel[lo] + (lower_res_storage - storage[lo]) / area[lo]
        else # ReservoirProfileType.interpolation
            interpolate_linear(
                lower_res_storage,
                storage_waterlevel_curve[lo].S,
                storage_waterlevel_curve[lo].H,
            )
        end

        # update values for the lower reservoir in place
        res_v.outflow[lo] = -outflow
        res_v.outflow_cumulative[lo] -= outflow * dt
        res_v.storage[lo] = lower_res_storage
        waterlevel[lo] = lower_res_waterlevel
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
    (; storage, outflow_obs) = reservoir_model.variables
    (; precipitation, evaporation, inflow) = boundary_vars
    storage_input = max(storage[i] / dt + precipitation - evaporation + inflow, 0.0)
    outflow = min(outflow_obs[i], storage_input)
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
)
    res_bc = reservoir_model.boundary_conditions
    res_p = reservoir_model.parameters
    res_v = reservoir_model.variables

    # limit reservoir evaporation based on total available volume [m³]
    precipitation = res_bc.precipitation[i] * res_p.area[i]
    available_storage = res_v.storage[i] + (inflow + precipitation) * dt
    potential_evaporation = res_bc.evaporation[i] * res_p.area[i]
    evaporation = min(available_storage / dt, potential_evaporation)

    boundary_vars = (; precipitation, evaporation, inflow)
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
            res_p.storage_waterlevel_curve[i].S,
            res_p.storage_waterlevel_curve[i].H,
        )
    end

    # update values in place
    # instantaneous variables
    res_v.storage[i] = storage
    res_v.waterlevel[i] = waterlevel
    res_v.outflow[i] = outflow

    # average variables (here accumulated for model timestep dt)
    res_bc.inflow_cumulative[i] += inflow * dt
    res_v.outflow_cumulative[i] += outflow * dt
    res_v.actevap_cumulative[i] += evaporation / res_p.area[i] * dt
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
