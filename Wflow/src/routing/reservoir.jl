"Struct for storing reservoir model parameters"
@with_kw struct ReservoirParameters{T <: DenseArray{Float}}
    maxstorage::T                     # maximum storage (above which water is spilled) [m³]
    area::T                           # reservoir area [m²]
    maxrelease::T                     # maximum amount that can be released if below spillway [m³ s⁻¹]
    demand::T                         # minimum (environmental) flow requirement downstream of the reservoir [m³ s⁻¹]
    targetminfrac::T                  # target minimum full fraction (of max storage) [-]
    targetfullfrac::T                 # target fraction full (of max storage
end

function Adapt.adapt_structure(to, from::ReservoirParameters)
    return ReservoirParameters(
        adapt(to, from.maxstorage),
        adapt(to, from.area),
        adapt(to, from.maxrelease),
        adapt(to, from.demand),
        adapt(to, from.targetminfrac),
        adapt(to, from.targetfullfrac),
    )
end

"Initialize reservoir model parameters"
function ReservoirParameters(dataset::NCDataset, config::Config, network::NetworkWaterBody)
    (; indices_outlet) = network

    lens = lens_input_parameter(
        config,
        "reservoir_water_demand~required~downstream__volume_flow_rate";
        optional = false,
    )
    resdemand = ncread(dataset, config, lens; sel = indices_outlet, type = Float, fill = 0)
    lens = lens_input_parameter(
        config,
        "reservoir_water_release-below-spillway__max_volume_flow_rate";
        optional = false,
    )
    resmaxrelease =
        ncread(dataset, config, lens; sel = indices_outlet, type = Float, fill = 0)
    lens = lens_input_parameter(config, "reservoir_water__max_volume"; optional = false)
    resmaxstorage =
        ncread(dataset, config, lens; sel = indices_outlet, type = Float, fill = 0)
    lens = lens_input_parameter(config, "reservoir_surface__area"; optional = false)
    resarea = ncread(dataset, config, lens; sel = indices_outlet, type = Float, fill = 0)
    lens = lens_input_parameter(
        config,
        "reservoir_water~full-target__volume_fraction";
        optional = false,
    )
    res_targetfullfrac =
        ncread(dataset, config, lens; sel = indices_outlet, type = Float, fill = 0)
    lens = lens_input_parameter(
        config,
        "reservoir_water~min-target__volume_fraction";
        optional = false,
    )
    res_targetminfrac =
        ncread(dataset, config, lens; sel = indices_outlet, type = Float, fill = 0)

    parameters = ReservoirParameters(;
        demand = resdemand,
        maxrelease = resmaxrelease,
        maxstorage = resmaxstorage,
        area = resarea,
        targetfullfrac = res_targetfullfrac,
        targetminfrac = res_targetminfrac,
    )

    return parameters
end

"Struct for storing reservoir model variables"
@with_kw struct ReservoirVariables{T <: DenseArray{Float}}
    storage::T                        # reservoir storage [m³]
    storage_av::T                     # average reservoir storage [m³] for model timestep Δt
    outflow::T                        # outflow from reservoir [m³ s⁻¹]
    outflow_av::T                     # average outflow from reservoir [m³ s⁻¹] for model timestep Δt
    percfull::T                       # fraction full (of max storage) [-]
    demandrelease::T                  # minimum (environmental) flow released from reservoir [m³ s⁻¹]
    actevap::T                        # average actual evaporation for reservoir area [mm Δt⁻¹]
end

function Adapt.adapt_structure(to, from::ReservoirVariables)
    return ReservoirVariables(
        adapt(to, from.storage),
        adapt(to, from.storage_av),
        adapt(to, from.outflow),
        adapt(to, from.outflow_av),
        adapt(to, from.percfull),
        adapt(to, from.demandrelease),
        adapt(to, from.actevap),
    )
end

"Initialize reservoir model variables"
function ReservoirVariables(n::Int, parameters::ReservoirParameters)
    (; targetfullfrac, maxstorage) = parameters
    variables = ReservoirVariables(;
        storage = targetfullfrac .* maxstorage,
        storage_av = fill(MISSING_VALUE, n),
        outflow = fill(MISSING_VALUE, n),
        outflow_av = fill(MISSING_VALUE, n),
        percfull = fill(MISSING_VALUE, n),
        demandrelease = fill(MISSING_VALUE, n),
        actevap = fill(MISSING_VALUE, n),
    )
    return variables
end

"Struct for storing reservoir model boundary conditions"
@with_kw struct ReservoirBC{T <: DenseArray{Float}}
    inflow::T                             # inflow into reservoir [m³ s⁻¹] for model timestep Δt
    precipitation::T                      # average precipitation for reservoir area [mm Δt⁻¹]
    evaporation::T                        # average potential evaporation for reservoir area [mm Δt⁻¹]
end

function Adapt.adapt_structure(to, from::ReservoirBC)
    return ReservoirBC(
        adapt(to, from.inflow),
        adapt(to, from.precipitation),
        adapt(to, from.evaporation),
    )
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

"Reservoir model `SimpleReservoir`"
@with_kw struct SimpleReservoir{T <: DenseArray{Float}}
    boundary_conditions::ReservoirBC{T}
    parameters::ReservoirParameters{T}
    variables::ReservoirVariables{T}
end

function Adapt.adapt_structure(to, from::SimpleReservoir)
    return SimpleReservoir(
        adapt(to, from.boundary_conditions),
        adapt(to, from.parameters),
        adapt(to, from.variables),
    )
end

"Initialize reservoir model `SimpleReservoir`"
function SimpleReservoir(dataset::NCDataset, config::Config, network::NetworkWaterBody)
    parameters = ReservoirParameters(dataset, config, network)

    n_reservoirs = Int(length(parameters.area))
    @info "Read `$n_reservoirs` reservoir locations."

    variables = ReservoirVariables(n_reservoirs, parameters)
    boundary_conditions = ReservoirBC(n_reservoirs)
    reservoir = SimpleReservoir(; boundary_conditions, parameters, variables)

    return reservoir
end

"""
Update a single reservoir at position `i`.

This is called from within the kinematic wave loop, therefore updating only for a single
element rather than all at once.
"""
function update!(
    model::SimpleReservoir,
    i::Int,
    inflow::Float,
    dt::Float,
    dt_forcing::Float,
)
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