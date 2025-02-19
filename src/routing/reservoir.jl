"Struct for storing reservoir model parameters"
@with_kw struct ReservoirParameters{T}
    maxstorage::Vector{T}       # maximum storage (above which water is spilled) [m³]
    area::Vector{T}             # reservoir area [m²]
    maxrelease::Vector{T}       # maximum amount that can be released if below spillway [m³ s⁻¹]
    demand::Vector{T}           # minimum (environmental) flow requirement downstream of the reservoir [m³ s⁻¹]
    targetminfrac::Vector{T}    # target minimum full fraction (of max storage) [-]
    targetfullfrac::Vector{T}   # target fraction full (of max storage
end

"Initialize reservoir model parameters"
function ReservoirParameters(dataset, config, indices_river, n_river_cells, pits)
    # read only reservoir data if reservoirs true
    # allow reservoirs only in river cells
    # note that these locations are only the reservoir outlet pixels
    lens = lens_input(config, "reservoir_location__count"; optional = false)
    reslocs = ncread(dataset, config, lens; sel = indices_river, type = Int, fill = 0)

    # this holds the same ids as reslocs, but covers the entire reservoir
    lens = lens_input(config, "reservoir_area__count"; optional = false)
    rescoverage_2d = ncread(dataset, config, lens; allow_missing = true)
    # for each reservoir, a list of 2D indices, needed for getting the mean precipitation
    inds_res_cov = Vector{CartesianIndex{2}}[]

    rev_inds_reservoir = zeros(Int, size(rescoverage_2d))

    # construct a map from the rivers to the reservoirs and
    # a map of the reservoirs to the 2D model grid
    inds_reservoir_map2river = fill(0, n_river_cells)
    inds_res = CartesianIndex{2}[]
    rescounter = 0
    for (i, ind) in enumerate(indices_river)
        res_id = reslocs[i]
        if res_id > 0
            push!(inds_res, ind)
            rescounter += 1
            inds_reservoir_map2river[i] = rescounter
            rev_inds_reservoir[ind] = rescounter

            # get all indices related to this reservoir outlet
            # done in this loop to ensure that the order is equal to the order in the
            # SimpleReservoir struct
            res_cov = findall(isequal(res_id), rescoverage_2d)
            push!(inds_res_cov, res_cov)
        end
    end
    lens = lens_input_parameter(
        config,
        "reservoir_water_demand~required~downstream__volume_flow_rate";
        optional = false,
    )
    resdemand = ncread(dataset, config, lens; sel = inds_res, type = Float, fill = 0)
    lens = lens_input_parameter(
        config,
        "reservoir_water_release-below-spillway__max_volume_flow_rate";
        optional = false,
    )
    resmaxrelease = ncread(dataset, config, lens; sel = inds_res, type = Float, fill = 0)
    lens = lens_input_parameter(config, "reservoir_water__max_volume"; optional = false)
    resmaxstorage = ncread(dataset, config, lens; sel = inds_res, type = Float, fill = 0)
    lens = lens_input_parameter(config, "reservoir_surface__area"; optional = false)
    resarea = ncread(dataset, config, lens; sel = inds_res, type = Float, fill = 0)
    lens = lens_input_parameter(
        config,
        "reservoir_water~full-target__volume_fraction";
        optional = false,
    )
    res_targetfullfrac =
        ncread(dataset, config, lens; sel = inds_res, type = Float, fill = 0)
    lens = lens_input_parameter(
        config,
        "reservoir_water~min-target__volume_fraction";
        optional = false,
    )
    res_targetminfrac =
        ncread(dataset, config, lens; sel = inds_res, type = Float, fill = 0)

    # for surface water routing reservoir locations are considered pits in the flow network
    # all upstream flow goes to the river and flows into the reservoir
    pits[inds_res] .= true

    reservoir_network = (
        indices_outlet = inds_res,
        indices_coverage = inds_res_cov,
        reverse_indices = rev_inds_reservoir,
        river_indices = findall(x -> x ≠ 0, inds_reservoir_map2river),
    )

    parameters = ReservoirParameters{Float}(;
        demand = resdemand,
        maxrelease = resmaxrelease,
        maxstorage = resmaxstorage,
        area = resarea,
        targetfullfrac = res_targetfullfrac,
        targetminfrac = res_targetminfrac,
    )

    return parameters, reservoir_network, inds_reservoir_map2river, pits
end

"Struct for storing reservoir model variables"
@with_kw struct ReservoirVariables{T}
    storage::Vector{T}          # reservoir storage [m³]
    storage_av::Vector{T}       # average reservoir storage [m³] for model timestep Δt
    outflow::Vector{T}          # outflow from reservoir [m³ s⁻¹]
    outflow_av::Vector{T}       # average outflow from reservoir [m³ s⁻¹] for model timestep Δt
    percfull::Vector{T}         # fraction full (of max storage) [-]
    demandrelease::Vector{T}    # minimum (environmental) flow released from reservoir [m³ s⁻¹]
    actevap::Vector{T}          # average actual evaporation for reservoir area [mm Δt⁻¹]
end

"Initialize reservoir model variables"
function ReservoirVariables(n, parameters)
    (; targetfullfrac, maxstorage) = parameters
    variables = ReservoirVariables{Float}(;
        storage = targetfullfrac .* maxstorage,
        storage_av = fill(mv, n),
        outflow = fill(mv, n),
        outflow_av = fill(mv, n),
        percfull = fill(mv, n),
        demandrelease = fill(mv, n),
        actevap = fill(mv, n),
    )
    return variables
end

"Struct for storing reservoir model boundary conditions"
@with_kw struct ReservoirBC{T}
    inflow::Vector{T}               # inflow into reservoir [m³ s⁻¹] for model timestep Δt
    precipitation::Vector{T}        # average precipitation for reservoir area [mm Δt⁻¹]
    evaporation::Vector{T}          # average potential evaporation for reservoir area [mm Δt⁻¹]
end

"Initialize reservoir model boundary conditions"
function ReservoirBC(n)
    bc = ReservoirBC{Float}(;
        inflow = fill(mv, n),
        precipitation = fill(mv, n),
        evaporation = fill(mv, n),
    )
    return bc
end

"Reservoir model `SimpleReservoir`"
@with_kw struct SimpleReservoir{T}
    boundary_conditions::ReservoirBC{T}
    parameters::ReservoirParameters{T}
    variables::ReservoirVariables{T}
end

"Initialize reservoir model `SimpleReservoir`"
function SimpleReservoir(dataset, config, indices_river, n_river_cells, pits)
    parameters, reservoir_network, inds_reservoir_map2river, pits =
        ReservoirParameters(dataset, config, indices_river, n_river_cells, pits)

    n_reservoirs = length(parameters.area)
    @info "Read `$n_reservoirs` reservoir locations."

    variables = ReservoirVariables(n_reservoirs, parameters)
    boundary_conditions = ReservoirBC(n_reservoirs)
    reservoir = SimpleReservoir{Float}(; boundary_conditions, parameters, variables)

    return reservoir, reservoir_network, inds_reservoir_map2river, pits
end

"""
Update a single reservoir at position `i`.

This is called from within the kinematic wave loop, therefore updating only for a single
element rather than all at once.
"""
function update!(model::SimpleReservoir, i, inflow, dt, dt_forcing)
    res_bc = model.boundary_conditions
    res_p = model.parameters
    res_v = model.variables

    # limit lake evaporation based on total available volume [m³]
    precipitation = 0.001 * res_bc.precipitation[i] * (dt / dt_forcing) * res_p.area[i]
    available_storage = res_v.storage[i] + inflow * dt + precipitation
    evap = 0.001 * res_bc.evaporation[i] * (dt / dt_forcing) * res_p.area[i]
    actevap = min(available_storage, evap) # [m³/dt]

    storage = res_v.storage[i] + (inflow * dt) + precipitation - actevap
    storage = max(storage, 0.0)

    percfull = storage / res_p.maxstorage[i]
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(percfull, res_p.targetminfrac[i], Float(1.0), Float(30.0))
    demandrelease = min(fac * res_p.demand[i] * dt, storage)
    storage = storage - demandrelease

    wantrel = max(0.0, storage - (res_p.maxstorage[i] * res_p.targetfullfrac[i]))
    # Assume extra maximum Q if spilling
    overflow_q = max((storage - res_p.maxstorage[i]), 0.0)
    torelease = min(wantrel, overflow_q + res_p.maxrelease[i] * dt - demandrelease)
    storage = storage - torelease
    outflow = torelease + demandrelease
    percfull = storage / res_p.maxstorage[i]

    # update values in place
    # instantaneous variables
    res_v.demandrelease[i] = demandrelease / dt
    res_v.percfull[i] = percfull
    res_v.storage[i] = storage
    res_v.outflow[i] = outflow / dt
    # average variables (here accumulated for model timestep Δt)
    res_bc.inflow[i] += inflow * dt
    res_v.outflow_av[i] += outflow
    res_v.storage_av[i] += storage * dt
    res_v.actevap[i] += 1000.0 * (actevap / res_p.area[i])

    return nothing
end