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
    lower_reservoir_ind::Vector{Int} = zeros(Int, length(area))
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

    area = ncread(
        dataset,
        config,
        "reservoir_surface__area",
        Routing;
        optional = false,
        sel = indices_outlet,
        type = Float64,
        fill = 0,
    )
    waterlevel = ncread(
        dataset,
        config,
        "reservoir_water_surface__initial_elevation",
        Routing;
        optional = false,
        sel = indices_outlet,
        type = Float64,
        fill = 0,
    )
    storfunc = ncread(
        dataset,
        config,
        "reservoir_water__storage_curve_type_count",
        Routing;
        optional = false,
        sel = indices_outlet,
        type = Int,
        fill = 0,
    )
    storfunc = to_enumx.(ReservoirProfileType.T, storfunc)
    outflowfunc = ncread(
        dataset,
        config,
        "reservoir_water__rating_curve_type_count",
        Routing;
        optional = false,
        sel = indices_outlet,
        type = Int,
        fill = 0,
    )
    outflowfunc = to_enumx.(ReservoirOutflowType.T, outflowfunc)
    linked_reslocs = ncread(
        dataset,
        config,
        "reservoir_lower_location__count",
        Routing;
        sel = indices_outlet,
        defaults = 0,
        type = Int,
        fill = 0,
    )

    n_reservoirs = length(area)
    reslocs = ncread(
        dataset,
        config,
        "reservoir_location__count",
        Routing;
        optional = false,
        sel = indices_outlet,
        type = Int,
        fill = 0,
    )
    @info "Read `$n_reservoirs` reservoir locations."

    parameters = ReservoirParameters(; id = reslocs, area, outflowfunc, storfunc)

    if ReservoirOutflowType.free_weir in outflowfunc ||
       ReservoirOutflowType.modified_puls in outflowfunc
        threshold = ncread(
            dataset,
            config,
            "reservoir_water_flow_threshold_level__elevation",
            Routing;
            optional = false,
            sel = indices_outlet,
            type = Float64,
            fill = 0,
        )
        b = ncread(
            dataset,
            config,
            "reservoir_water__rating_curve_coefficient",
            Routing;
            optional = false,
            sel = indices_outlet,
            type = Float64,
            fill = 0,
        )
        e = ncread(
            dataset,
            config,
            "reservoir_water__rating_curve_exponent",
            Routing;
            optional = false,
            sel = indices_outlet,
            type = Float64,
            fill = 0,
        )
    end
    if ReservoirOutflowType.simple in outflowfunc
        demand = ncread(
            dataset,
            config,
            "reservoir_water_demand__required_downstream_volume_flow_rate",
            Routing;
            optional = false,
            sel = indices_outlet,
            type = Float64,
            fill = 0,
        )
        maxrelease = ncread(
            dataset,
            config,
            "reservoir_water_release_below_spillway__max_volume_flow_rate",
            Routing;
            optional = false,
            sel = indices_outlet,
            type = Float64,
            fill = 0,
        )
        maxstorage = ncread(
            dataset,
            config,
            "reservoir_water__max_volume",
            Routing;
            optional = false,
            sel = indices_outlet,
            type = Float64,
            fill = 0,
        )
        targetfullfrac = ncread(
            dataset,
            config,
            "reservoir_water__target_full_volume_fraction",
            Routing;
            optional = false,
            sel = indices_outlet,
            type = Float64,
            fill = 0,
        )
        targetminfrac = ncread(
            dataset,
            config,
            "reservoir_water__target_min_volume_fraction",
            Routing;
            optional = false,
            sel = indices_outlet,
            type = Float64,
            fill = 0,
        )
    end

    # reservoir CSV parameter files are expected in the same directory as path_static
    path = dirname(input_path(config, config.input.path_static))
    for i in 1:n_reservoirs
        resloc = reslocs[i]
        if linked_reslocs[i] > 0
            parameters.lower_reservoir_ind[i] =
                only(findall(x -> x == linked_reslocs[i], reslocs))
        end

        if storfunc[i] == ReservoirProfileType.interpolation
            csv_path = joinpath(path, "reservoir_sh_$resloc.csv")
            @info(
                "Read a storage curve from CSV file `$csv_path`, for reservoir location `$resloc`"
            )
            parameters.sh[i] = read_sh_csv(csv_path)
        end

        if outflowfunc[i] == ReservoirOutflowType.rating_curve
            csv_path = joinpath(path, "reservoir_hq_$resloc.csv")
            @info(
                "Read a rating curve from CSV file `$csv_path`, for reservoir location `$resloc`"
            )
            parameters.hq[i] = read_hq_csv(csv_path)
            parameters.maxstorage[i] = maximum_storage(parameters, i)
        elseif outflowfunc[i] == ReservoirOutflowType.free_weir ||
               outflowfunc[i] == ReservoirOutflowType.modified_puls
            parameters.threshold[i] = threshold[i]
            parameters.b[i] = b[i]
            parameters.e[i] = e[i]
        elseif outflowfunc[i] == ReservoirOutflowType.simple
            parameters.demand[i] = demand[i]
            parameters.maxrelease[i] = maxrelease[i]
            parameters.maxstorage[i] = maxstorage[i]
            parameters.targetfullfrac[i] = targetfullfrac[i]
            parameters.targetminfrac[i] = targetminfrac[i]
        end

        if outflowfunc[i] == ReservoirOutflowType.modified_puls &&
           storfunc[i] != ReservoirProfileType.linear
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
    # average outflow from reservoir [m³ s⁻¹] for model timestep dt
    outflow_av::AverageVector = AverageVector(; n = length(waterlevel))
    # observed outflow from reservoir [m³ s⁻¹]
    outflow_obs::Vector{Float64} = fill(MISSING_VALUE, length(waterlevel))
    # average actual evaporation for reservoir area [mm dt⁻¹ => m s⁻¹]
    actevap::AverageVector = AverageVector(; n = length(waterlevel))
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
    (; indices_outlet) = network
    outflow_obs = ncread(
        dataset,
        config,
        "reservoir_water__outgoing_observed_volume_flow_rate",
        Routing;
        sel = indices_outlet,
        defaults = MISSING_VALUE,
        type = Float64,
        fill = MISSING_VALUE,
    )
    variables = ReservoirVariables(;
        waterlevel,
        storage = initialize_storage(storfunc, area, waterlevel, sh),
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
    # total inflow into reservoir [m³ s⁻¹] for model timestep dt
    inflow::AverageVector = AverageVector(; n)
    # external inflow (abstraction/supply/demand) [m³ s⁻¹]
    external_inflow::Vector{Float64} = zeros(n)
    # cumulative actual abstraction from external negative inflow [m³]
    actual_external_abstraction_av::AverageVector = AverageVector(; n)
    # average precipitation for reservoir area [mm dt⁻¹ => m s⁻¹]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)
    # average potential evaporation for reservoir area [mm dt⁻¹ => m s⁻¹]
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
        defaults = 0.0,
        type = Float64,
    )
    n = length(indices_outlet)
    bc = ReservoirBC(; n, external_inflow)
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
    variables = ReservoirVariables(dataset, config, network, parameters, waterlevel)
    boundary_conditions = ReservoirBC(dataset, config, network)
    reservoir = Reservoir(; boundary_conditions, parameters, variables)

    return reservoir
end

"Determine the water level depending on the storage function"
function waterlevel(
    storfunc::ReservoirProfileType.T,
    area::Float64,
    storage::Float64,
    sh::Union{SH, Missing},
)
    if storfunc == ReservoirProfileType.linear
        waterlevel = storage / area
    else # storfunc == ReservoirProfileType.interpolation
        waterlevel = interpolate_linear(storage, sh.S, sh.H)
    end
    return waterlevel
end

"Determine the maximum storage for reservoirs with a rating curve of type 1"
function maximum_storage(parameters::ReservoirParameters, i::Int)
    (; storfunc, hq, sh, area) = parameters

    # maximum storage is based on the maximum water level (H) value in the H-Q table
    if storfunc[i] == ReservoirProfileType.interpolation
        maxstorage = interpolate_linear(maximum(hq[i].H), sh[i].H, sh[i].S)
    else # storfunc[i] == ReservoirProfileType. linear
        maxstorage = area[i] * maximum(hq[i].H)
    end

    return maxstorage
end

"Determine the initial storage depending on the storage function"
function initialize_storage(
    storfunc::Vector{ReservoirProfileType.T},
    area::Vector{Float64},
    waterlevel::Vector{Float64},
    sh::Vector{Union{SH, Missing}},
)
    storage = similar(area)
    for i in eachindex(storage)
        if storfunc[i] == ReservoirProfileType.linear
            storage[i] = area[i] * waterlevel[i]
        else # storfunc[i] == ReservoirProfileType.interpolation
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
    if ReservoirOutflowType.rating_curve in outflowfunc
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
    (; maxstorage, targetminfrac, targetfullfrac, demand, maxrelease) = model.parameters
    (; storage) = model.variables
    (; precipitation, evaporation, inflow) = boundary_vars

    # [m³] = [m³] + ([m³ s⁻¹] + [m³ s⁻¹] + [m³ s⁻¹]) * [s]
    storage = storage[i] + (inflow + precipitation - evaporation) * dt
    storage = max(storage, 0.0)

    # [-] = [m³] / [m³]
    fill_fraction = storage / maxstorage[i]
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(fill_fraction, targetminfrac[i], 1.0, 30.0)
    # [m³ s⁻¹] = min([-] * [m³ s⁻¹], [m³] / [s])
    demand_release = min(fac * demand[i], storage / dt)
    # [m³] -= [m³ s⁻¹] * [s]
    storage -= demand_release * dt
    # [m³ s⁻¹] = max([m³ s⁻¹], ([m³] - [m³] * [-]) / [s])
    release_wanted = max(0.0, (storage - maxstorage[i] * targetfullfrac[i]) / dt)
    # Assume extra maximum Q if spilling
    # [m³ s⁻¹] = max([m³ s⁻¹], ([m³] - [m³]) / [s])
    overflow_q = max(0.0, (storage - maxstorage[i]) / dt)
    # [m³ s⁻¹] = min([m³ s⁻¹], [m³ s⁻¹] + [m³ s⁻¹] - [m³ s⁻¹])
    release_realized = min(release_wanted, overflow_q + maxrelease[i] - demand_release)
    # [m³] -= [m³ s⁻¹] * [s]
    storage -= release_realized * dt
    # [m³ s⁻¹] = [m³ s⁻¹] + [m³ s⁻¹]
    outflow = release_realized + demand_release

    return outflow, storage
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
    (; area, threshold, b) = model.parameters
    (; storage) = model.variables
    (; precipitation, evaporation, inflow) = boundary_vars

    # [m³ᐟ² s⁻¹ᐟ²] = [m²] / ([s] * sqrt([m s⁻¹]))
    res_factor = area[i] / (dt * sqrt(b[i]))
    # [m³ s⁻¹] = [m³] / [s] + [m³ s⁻¹] - [m³ s⁻¹] + [m³ s⁻¹]
    si_factor = storage[i] / dt + precipitation - evaporation + inflow
    # Adjust si_factor for reservoir threshold != 0
    # [m³ s⁻¹] = [m³ s⁻¹] - [m²] * [m] / [s]
    si_factor_adj = si_factor - area[i] * threshold[i] / dt
    # Calculate the new reservoir outflow/waterlevel/storage
    if si_factor_adj > 0.0
        # [m³ᐟ² s⁻¹ᐟ²] = -[m³ᐟ² s⁻¹ᐟ²] + sqrt([m³ᐟ² s⁻¹ᐟ²]^2 + [-] * [m³ s⁻¹])
        quadratic_sol_term = -res_factor + sqrt((res_factor^2 + 4 * si_factor_adj))
        if quadratic_sol_term > 0.0
            # [m³ s⁻¹] = [-] * [m³ᐟ² s⁻¹ᐟ²]^2
            outflow = 0.25 * quadratic_sol_term^2
        else
            outflow = 0.0
        end
    else
        outflow = 0.0
    end
    # [m³ s⁻¹] = min([m³ s⁻¹], [m³ s⁻¹])
    outflow = min(outflow, si_factor)
    # [m³] = ([m³ s⁻¹] - [m³ s⁻¹]) * dt
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
    (; hq, col_index_hq, maxstorage) = model.parameters
    (; storage, waterlevel) = model.variables
    (; precipitation, evaporation, inflow) = boundary_vars

    # [m³ s⁻¹] = [m³] / [s] + [m³ s⁻¹] - [m³ s⁻¹] + [m³ s⁻¹]
    storage_input = storage[i] / dt + precipitation - evaporation + inflow
    # [m³ s⁻¹]
    outflow = interpolate_linear(waterlevel[i], hq[i].H, hq[i].Q[:, col_index_hq[1]])
    # [m³ s⁻¹] = min([m³ s⁻¹], [m³ s⁻¹])
    outflow = min(outflow, storage_input)
    # [m³] = ([m³ s⁻¹] - [m³ s⁻¹]) * [s]
    storage = (storage_input - outflow) * dt

    # [m³ s⁻¹] = max([m³ s⁻¹], ([m³] - [m³]) / [s])
    overflow = max(0.0, (storage - maxstorage[i]) / dt)
    storage = min(storage, maxstorage[i])
    outflow += overflow

    return outflow, storage
end

"Update reservoir with rating curve type (`ouflowfunc`) 2 (free weir) for a single timestep."
function update_reservoir_free_weir(
    model::Reservoir,
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    (; threshold, b, e, area, storfunc, sh, lower_reservoir_ind) = model.parameters
    res_v = model.variables
    (; waterlevel) = res_v
    (; precipitation, evaporation, inflow) = boundary_vars

    # [-]
    lo = lower_reservoir_ind[i]
    has_lower_res = (lo != 0)
    # [m]
    diff_wl = has_lower_res ? waterlevel[i] - waterlevel[lo] : 0.0

    # [m³ s⁻¹] = [m³] / [s] + [m³ s⁻¹] - [m³ s⁻¹] + [m³ s⁻¹]
    storage_input = res_v.storage[i] / dt + precipitation - evaporation + inflow

    if diff_wl >= 0.0
        if res_v.waterlevel[i] > threshold[i]
            # [m]
            dh = waterlevel[i] - threshold[i]
            # [m³ s⁻¹] = [m³⁻ᵉ s⁻¹] * [m]ᵉ
            outflow = b[i] * pow(dh, e[i])
            # [m³ s⁻¹] = [m] * [m²] / [s]
            maxflow = dh * area[i] / dt
            # [m³ s⁻¹] = min([m³ s⁻¹], [m³ s⁻¹])
            outflow = min(outflow, maxflow)
        else
            outflow = 0.0
        end
    else
        if waterlevel[lo] > threshold[i]
            # [m]
            dh = waterlevel[lo] - threshold[i]
            # [m³ s⁻¹] = -[m³⁻ᵉ s⁻¹] * [m]ᵉ
            outflow = -b[i] * pow(dh, e[i])
            # [m³ s⁻¹] = [m] * [m²] / [s]
            maxflow = dh * area[lo] / dt
            # [m³ s⁻¹] = max([m³ s⁻¹], [m³ s⁻¹])
            outflow = max(outflow, -maxflow)
        else
            outflow = 0.0
        end
    end
    # [m³] = ([m³ s⁻¹] - [m³ s⁻¹]) * [s]
    storage = (storage_input - outflow) * dt

    # update lower reservoir (linked reservoirs) in case flow from lower reservoir to upper reservoir occurs
    if diff_wl < 0.0
        # [m³] = [m³] + [m³ s⁻¹] * [s]
        lower_res_storage = res_v.storage[lo] + outflow * dt
        # [m]
        lower_res_waterlevel = if storfunc[lo] == ReservoirProfileType.linear
            # [m] + ([m³] - [m³]) / [m²]
            waterlevel[lo] + (lower_res_storage - storage[lo]) / area[lo]
        else # ReservoirProfileType.interpolation
            interpolate_linear(lower_res_storage, sh[lo].S, sh[lo].H)
        end

        # update values for the lower reservoir in place
        res_v.outflow[lo] = -outflow
        add_to_cumulative!(res_v.outflow_av, lo, -outflow * dt)
        res_v.storage[lo] = lower_res_storage
        waterlevel[lo] = lower_res_waterlevel
    end
    return outflow, storage
end

"Update reservoir using observed outflow for a single timestep."
function update_reservoir_outflow_obs(
    model::Reservoir,
    i::Int,
    boundary_vars::NamedTuple,
    dt::Float64,
)
    (; storage, outflow_obs) = model.variables
    (; precipitation, evaporation, inflow) = boundary_vars
    # [m³ s⁻¹] = [m³] / [s] + [m³ s⁻¹] - [m³ s⁻¹] + [m³ s⁻¹]
    storage_input = storage[i] / dt + precipitation - evaporation + inflow
    # [m³ s⁻¹] = min([m³ s⁻¹], [m³ s⁻¹])
    outflow = min(outflow_obs[i], storage_input)
    # [m³] = ([m³ s⁻¹] - [m³ s⁻¹]) * [s]
    storage = (storage_input - outflow) * dt
    return outflow, storage
end

"""
Update a single reservoir at position `i`.

This is called from within the river routing scheme, therefore updating only for a single
element rather than all at once.
"""
function update!(model::Reservoir, i::Int, inflow::Float64, dt::Float64)
    res_bc = model.boundary_conditions
    res_p = model.parameters
    res_v = model.variables

    # limit reservoir evaporation based on total available volume [m³]
    # [m³ s⁻¹] = [m s⁻¹] * [m²]
    precipitation = res_bc.precipitation[i] * res_p.area[i]
    # [m³] = [m³] + ([m³ s⁻¹] + [m³ s⁻¹]) * [s]
    available_storage = res_v.storage[i] + (inflow + precipitation) * dt
    # [m³ s⁻¹] = [m s⁻¹] * [m²]
    potential_evaporation = res_bc.evaporation[i] * res_p.area[i]
    # [m³ s⁻¹] = min([m³] / [s], [m³ s⁻¹])
    evaporation = min(available_storage / dt, potential_evaporation)

    boundary_vars = (; precipitation, evaporation, inflow)
    update_reservoir_args = (model, i, boundary_vars, dt)
    if !isnan(res_v.outflow_obs[i])
        outflow, storage = update_reservoir_outflow_obs(update_reservoir_args...)
    elseif res_p.outflowfunc[i] == ReservoirOutflowType.rating_curve
        outflow, storage = update_reservoir_hq(update_reservoir_args...)
    elseif res_p.outflowfunc[i] == ReservoirOutflowType.free_weir
        outflow, storage = update_reservoir_free_weir(update_reservoir_args...)
    elseif res_p.outflowfunc[i] == ReservoirOutflowType.modified_puls
        outflow, storage = update_reservoir_modified_puls(update_reservoir_args...)
    elseif res_p.outflowfunc[i] == ReservoirOutflowType.simple
        outflow, storage = update_reservoir_simple(update_reservoir_args...)
    end

    waterlevel = if res_p.storfunc[i] == ReservoirProfileType.linear
        res_v.waterlevel[i] + (storage - res_v.storage[i]) / res_p.area[i]
    else # res_p.storfunc[i] == ReservoirProfileType.interpolation
        interpolate_linear(storage, res_p.sh[i].S, res_p.sh[i].H)
    end

    # update values in place
    # instantaneous variables
    res_v.storage[i] = storage
    res_v.waterlevel[i] = waterlevel
    res_v.outflow[i] = outflow

    # average variables (here accumulated for model timestep dt)
    add_to_cumulative!(res_bc.inflow, i, inflow * dt)
    add_to_cumulative!(res_v.outflow_av, i, outflow * dt)
    add_to_cumulative!(res_v.actevap, i, evaporation * dt / res_p.area[i])
    return nothing
end

"Generate log message for using observed outflow at reservoir locations"
function log_message_observed_outflow(reservoir::Reservoir)
    not_nan = findall(x -> !isnan(x), reservoir.variables.outflow_obs)
    if isempty(not_nan)
        msg = "Observed outflow is not used for any reservoir location"
    else
        ids = reservoir.parameters.id[not_nan]
        msg = "Observed outflow is used for reservoir location ids $ids"
    end
    return msg
end

"Check if observed outflow is used for reservoirs"
function using_observed_outflow(reservoir::Union{Reservoir, Nothing}, config::Config)
    par = "reservoir_water__outgoing_observed_volume_flow_rate"
    check =
        !isnothing(reservoir) &&
        (haskey(config.input.forcing, par) || haskey(config.input.cyclic, par))
    return check
end
