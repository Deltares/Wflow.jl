"Struct for storing reservoir model parameters"
@get_units @grid_loc @with_kw struct ReservoirParameters{T}
    maxvolume::Vector{T} | "m3"             # maximum storage (above which water is spilled) [m³]
    area::Vector{T} | "m2"                  # reservoir area [m²]
    maxrelease::Vector{T} | "m3 s-1"        # maximum amount that can be released if below spillway [m³ s⁻¹]
    demand::Vector{T} | "m3 s-1"            # minimum (environmental) flow requirement downstream of the reservoir [m³ s⁻¹]
    targetminfrac::Vector{T} | "-"          # target minimum full fraction (of max storage) [-]
    targetfullfrac::Vector{T} | "-"         # target fraction full (of max storage
end

"Initialize reservoir model parameters"
function ReservoirParameters(dataset, config, indices_river, n_river_cells, pits)
    # read only reservoir data if reservoirs true
    # allow reservoirs only in river cells
    # note that these locations are only the reservoir outlet pixels
    reslocs = ncread(
        dataset,
        config,
        "lateral.river.reservoir.locs";
        optional = false,
        sel = indices_river,
        type = Int,
        fill = 0,
    )

    # this holds the same ids as reslocs, but covers the entire reservoir
    rescoverage_2d = ncread(
        dataset,
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

    resdemand = ncread(
        dataset,
        config,
        "lateral.river.reservoir.demand";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    resmaxrelease = ncread(
        dataset,
        config,
        "lateral.river.reservoir.maxrelease";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    resmaxvolume = ncread(
        dataset,
        config,
        "lateral.river.reservoir.maxvolume";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    resarea = ncread(
        dataset,
        config,
        "lateral.river.reservoir.area";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    res_targetfullfrac = ncread(
        dataset,
        config,
        "lateral.river.reservoir.targetfullfrac";
        optional = false,
        sel = inds_res,
        type = Float,
        fill = 0,
    )
    res_targetminfrac = ncread(
        dataset,
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

    reservoir_network = (
        indices_outlet = inds_res,
        indices_coverage = inds_res_cov,
        reverse_indices = rev_inds_reservoir,
        river_indices = findall(x -> x ≠ 0, inds_reservoir_map2river),
    )

    parameters = ReservoirParameters{Float}(;
        demand = resdemand,
        maxrelease = resmaxrelease,
        maxvolume = resmaxvolume,
        area = resarea,
        targetfullfrac = res_targetfullfrac,
        targetminfrac = res_targetminfrac,
    )

    return parameters, reservoir_network, inds_reservoir_map2river, pits
end

"Struct for storing reservoir model variables"
@get_units @grid_loc @with_kw struct ReservoirVariables{T}
    volume::Vector{T} | "m3"                # reservoir volume [m³]
    outflow::Vector{T} | "m3 s-1"           # outflow from reservoir [m³ s⁻¹]
    outflow_av::Vector{T} | "m3 s-1"        # average outflow from reservoir [m³ s⁻¹] for model timestep Δt
    percfull::Vector{T} | "-"               # fraction full (of max storage) [-]
    demandrelease::Vector{T} | "m3 s-1"     # minimum (environmental) flow released from reservoir [m³ s⁻¹]
    actevap::Vector{T}                      # average actual evaporation for reservoir area [mm Δt⁻¹]
end

"Initialize reservoir model variables"
function ReservoirVariables(n, parameters)
    (; targetfullfrac, maxvolume) = parameters
    variables = ReservoirVariables{Float}(;
        volume = targetfullfrac .* maxvolume,
        outflow = fill(mv, n),
        outflow_av = fill(mv, n),
        percfull = fill(mv, n),
        demandrelease = fill(mv, n),
        actevap = fill(mv, n),
    )
    return variables
end

"Struct for storing reservoir model boundary conditions"
@get_units @grid_loc @with_kw struct ReservoirBC{T}
    inflow::Vector{T} | "m3 s-1"    # inflow into reservoir [m³ s⁻¹] for model timestep Δt
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
    available_volume = res_v.volume[i] + inflow * dt + precipitation
    evap = 0.001 * res_bc.evaporation[i] * (dt / dt_forcing) * res_p.area[i]
    actevap = min(available_volume, evap) # [m³/dt]

    vol = res_v.volume[i] + (inflow * dt) + precipitation - actevap
    vol = max(vol, 0.0)

    percfull = vol / res_p.maxvolume[i]
    # first determine minimum (environmental) flow using a simple sigmoid curve to scale for target level
    fac = scurve(percfull, res_p.targetminfrac[i], Float(1.0), Float(30.0))
    demandrelease = min(fac * res_p.demand[i] * dt, vol)
    vol = vol - demandrelease

    wantrel = max(0.0, vol - (res_p.maxvolume[i] * res_p.targetfullfrac[i]))
    # Assume extra maximum Q if spilling
    overflow_q = max((vol - res_p.maxvolume[i]), 0.0)
    torelease = min(wantrel, overflow_q + res_p.maxrelease[i] * dt - demandrelease)
    vol = vol - torelease
    outflow = torelease + demandrelease
    percfull = vol / res_p.maxvolume[i]

    # update values in place
    res_v.outflow[i] = outflow / dt
    res_bc.inflow[i] += inflow * dt
    res_v.outflow_av[i] += outflow
    res_v.demandrelease[i] = demandrelease / dt
    res_v.percfull[i] = percfull
    res_v.volume[i] = vol
    res_v.actevap[i] += 1000.0 * (actevap / res_p.area[i])

    return nothing
end