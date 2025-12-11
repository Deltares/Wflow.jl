
struct NoMassBalance <: AbstractMassBalance end

"""
Store mass balance error results (balance error and relative error), and storage at previous
time step `storage_prev` required for mass balance computation.
"""
@with_kw struct MassBalance <: AbstractMassBalance
    n::Int # number of cells/nodes
    storage_prev::Vector{Float64} = zeros(n)
    error::Vector{Float64} = zeros(n)
    relative_error::Vector{Float64} = zeros(n)
end

"""
Water mass balance error results (balance error and relative error) for river, overland,
subsurface and reservoir flow routing.
"""
@with_kw struct FlowRoutingMassBalance{
    R <: AbstractMassBalance,
    RT <: AbstractMassBalance,
} <: AbstractMassBalance
    river_water_balance::R
    reservoir_water_balance::RT
    overland_water_balance::MassBalance
    subsurface_water_balance::MassBalance
end

"""
Store water mass balance error results (balance error and relative error) computed for each
model time step Δt for a hydrological model.

- `land_water_balance`: Water mass balance results for a land hydrology model. This is
    largely based on vertical fluxes and lateral fluxes that control for example total soil
    water storage and snow storage (lateral snow transport).
- `routing`: Water mass balance results for river, overland, subsurface and reservoir flow
    routing.
"""
@with_kw struct HydrologicalMassBalance <: AbstractMassBalance
    land_water_balance::MassBalance
    routing::FlowRoutingMassBalance
end

"Initialize `HydrologicalMassBalance` for storing water mass balance error results."
function HydrologicalMassBalance(domain::Domain, config::Config)
    (; river_routing, land_routing, water_mass_balance__flag) = config.model
    if water_mass_balance__flag
        n_land = length(domain.land.network.indices)
        n_river = length(domain.river.network.indices)
        if config.model.reservoir__flag
            n_reservoir = length(domain.reservoir.network.indices_outlet)
            reservoir_water_balance = MassBalance(; n = n_reservoir)
        else
            reservoir_water_balance = NoMassBalance()
        end
        if land_routing == RoutingType.local_inertial &&
           river_routing == RoutingType.local_inertial
            river_water_balance = NoMassBalance()
        else
            river_water_balance = MassBalance(; n = n_river)
        end
        routing = FlowRoutingMassBalance(;
            river_water_balance,
            reservoir_water_balance,
            overland_water_balance = MassBalance(; n = n_land),
            subsurface_water_balance = MassBalance(; n = n_land),
        )
        mass_balance = HydrologicalMassBalance(;
            land_water_balance = MassBalance(; n = n_land),
            routing,
        )
    else
        mass_balance = NoMassBalance()
    end
    return mass_balance
end

"Compute mass balance error and relative error."
function compute_mass_balance_error(
    total_in::Float64,
    total_out::Float64,
    storage_rate::Float64,
)
    error = total_in - total_out - storage_rate
    average_flow_rate = (total_in + total_out) / 2.0
    relative_error = iszero(average_flow_rate) ? 0.0 : error / average_flow_rate
    return error, relative_error
end

"Compute total storage of land hydrolology model `LandHydrologySBM` at index `i`."
function compute_total_storage(model::LandHydrologySBM, i::Int)
    (; total_soilwater_storage) = model.soil.variables
    (; canopy_storage) = model.interception.variables
    (; snow, glacier, demand) = model

    snow_storage = get_snow_storage(snow)[i] + get_snow_water(snow)[i]
    glacier_storage = get_glacier_store(glacier)[i] * get_glacier_fraction(glacier)[i]
    paddy_storage = get_water_depth(demand.paddy)[i]

    total_storage =
        total_soilwater_storage[i] +
        canopy_storage[i] +
        snow_storage +
        glacier_storage +
        paddy_storage

    return total_storage
end

"Compute total storage of land hydrology model `LandHydrologySBM`."
function compute_total_storage!(model::LandHydrologySBM, water_balance::MassBalance)
    (; storage_prev) = water_balance
    for i in eachindex(storage_prev)
        storage_prev[i] = compute_total_storage(model, i)
    end
    return nothing
end

"""
    get_storage(model::LocalInertialRiverFlow, i)
    get_storage(model::KinWaveRiverFlow, i)

Return storage of a river flow model at index `i`. For `LocalInertialRiverFlow` floodplain
storage is added to river storage if an optional floodplain is included.
"""
function get_storage(model::LocalInertialRiverFlow, i)
    (; storage) = model.variables
    if isnothing(model.floodplain)
        return storage[i]
    else
        total_storage = storage[i] + model.floodplain.variables.storage[i]
        return total_storage
    end
end
get_storage(model::KinWaveRiverFlow, i) = model.variables.storage[i]

"""
Save river (+ floodplain) storage at previous time step as `storage_prev` of river
`water_balance`.
"""
function storage_prev!(model::AbstractRiverFlowModel, water_balance::MassBalance)
    (; storage_prev) = water_balance
    for i in eachindex(storage_prev)
        storage_prev[i] = get_storage(model, i)
    end
    return nothing
end

"""
Save reservoir storage at previous time step as `storage_prev` of reservoir `water_balance`.
"""
function storage_prev!(reservoir::Reservoir, water_balance::MassBalance)
    water_balance.storage_prev .= reservoir.variables.storage
end

"""
    storage_prev!(model, ::HydrologicalMassBalance)
    storage_prev!(model, ::NoMassBalance)

Save storage at previous time step as `storage_prev` for river, reservoir, overland and
subsurface flow routing and `land` hydrology model as part of water mass balance
`HydrologicalMassBalance`. For `NoMassBalance` storage at previous time step is not
required.
"""
function storage_prev!(model, ::HydrologicalMassBalance)
    (;
        river_water_balance,
        reservoir_water_balance,
        overland_water_balance,
        subsurface_water_balance,
    ) = model.mass_balance.routing
    (; land_water_balance) = model.mass_balance
    (; land) = model
    (; overland_flow, river_flow, subsurface_flow) = model.routing
    (; reservoir) = model.routing.river_flow.boundary_conditions

    compute_total_storage!(land, land_water_balance)
    storage_prev!(river_flow, river_water_balance)
    overland_water_balance.storage_prev .= overland_flow.variables.storage
    subsurface_water_balance.storage_prev .= get_storage(subsurface_flow)
    storage_prev!(reservoir, reservoir_water_balance)
    return nothing
end

function storage_prev!(model, ::NoMassBalance)
    return nothing
end

"Compute total incoming vertical flux of land hydrology `SBM` at index `i`."
function vertical_in(model::LandHydrologySBM, i::Int)
    (; precipitation) = model.atmospheric_forcing
    (; allocation) = model
    total_in = precipitation[i] + get_irrigation_allocated(allocation)[i]
    return total_in
end

"Compute total outgoing vertical flux of land hydrology `SBM` at index `i`."
function vertical_out(model::LandHydrologySBM, i::Int)
    (; allocation) = model
    (; net_runoff, actevap, actleakage) = model.soil.variables
    (; net_runoff_river) = model.runoff.variables
    total_out =
        net_runoff[i] +
        actevap[i] +
        net_runoff_river[i] +
        actleakage[i] +
        get_groundwater_abstraction_flux(allocation)[i]
    return total_out
end

"""
Compute water mass balance error and relative error for `land` hydrology `SBM` of model type
`SbmModel`.
"""
function compute_land_hydrology_balance!(model::AbstractModel{<:SbmModel})
    (; storage_prev, error, relative_error) = model.mass_balance.land_water_balance
    (; snow) = model.land
    (; area) = model.domain.land.parameters
    (; subsurface_flow) = model.routing

    for i in eachindex(storage_prev)
        f_conv = (model.clock.dt / BASETIMESTEP) / (area[i] * 0.001)
        subsurface_flux_in = get_inflow(subsurface_flow)[i] * f_conv
        total_in = subsurface_flux_in + vertical_in(model.land, i) + get_snow_in(snow)[i]

        subsurface_flux_out = get_outflow(subsurface_flow)[i] * f_conv
        vertical_flux_out = vertical_out(model.land, i)
        total_out = subsurface_flux_out + vertical_flux_out + get_snow_out(snow)[i]
        storage = compute_total_storage(model.land, i)
        storage_rate = storage - storage_prev[i]
        error[i], relative_error[i] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    return nothing
end

"""
Compute water mass balance error and relative error for `land` hydrology `SBM` of model type
`SbmGwfModel`. Errors for subsurface flow constant head boundaries are set at zero.
"""
function compute_land_hydrology_balance!(model::AbstractModel{<:SbmGwfModel})
    (; storage_prev, error, relative_error) = model.mass_balance.land_water_balance
    (; snow) = model.land
    (; area) = model.domain.land.parameters
    (; subsurface_flow) = model.routing

    # exclude recharge from computing total incoming and outgoing boundary fluxes for
    # groundwaterflow, other boundaries are required for the total soil water balance.
    boundaries_flow_in, boundaries_flow_out =
        sum_boundary_fluxes(subsurface_flow; exclude = Recharge)

    for i in eachindex(storage_prev)
        f_conv = (model.clock.dt / BASETIMESTEP) / (area[i] * 0.001)
        subsurface_flux_in = get_inflow(subsurface_flow)[i] * f_conv
        total_in =
            subsurface_flux_in +
            vertical_in(model.land, i) +
            get_snow_in(snow)[i] +
            boundaries_flow_in[i] * f_conv
        subsurface_flux_out = get_outflow(subsurface_flow)[i] * f_conv
        vertical_flux_out = vertical_out(model.land, i)
        total_out =
            subsurface_flux_out +
            vertical_flux_out +
            get_snow_out(snow)[i] +
            boundaries_flow_out[i] * f_conv
        storage = compute_total_storage(model.land, i)
        storage_rate = storage - storage_prev[i]
        error[i], relative_error[i] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    k = subsurface_flow.constanthead.index
    error[k] .= 0.0
    relative_error[k] .= 0.0
    return nothing
end

"""
    compute_flow_balance!(reservoir::Reservoir, water_balance::MassBalance, dt::Float64)
    compute_flow_balance!(reservoir::Nothing, water_balance::NoMassBalance, dt::Float64)

Compute reservoir water mass balance error and relative error if reservoirs are included.
"""
function compute_flow_balance!(
    reservoir::Reservoir,
    water_balance::MassBalance,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; storage, outflow_av, actevap) = reservoir.variables
    (; precipitation, inflow) = reservoir.boundary_conditions
    (; area) = reservoir.parameters

    for i in eachindex(storage_prev)
        total_in = inflow[i] + (precipitation[i] * 0.001 * area[i]) / dt
        total_out = outflow_av[i] + (actevap[i] * 0.001 * area[i]) / dt
        storage_rate = (storage[i] - storage_prev[i]) / dt
        error[i], relative_error[i] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    return nothing
end

compute_flow_balance!(reservoir::Nothing, water_balance::NoMassBalance, dt::Float64) =
    nothing

"Compute water mass balance error and relative error for river kinematic wave routing."
function compute_flow_balance!(
    river_flow::KinWaveRiverFlow,
    water_balance::MassBalance,
    network::NetworkRiver,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; inwater, external_inflow, actual_external_abstraction_av, abstraction) =
        river_flow.boundary_conditions
    (; qin_av, q_av, storage) = river_flow.variables

    for i in eachindex(storage_prev)
        total_in = inwater[i] + qin_av[i] + max(0.0, external_inflow[i])
        total_out = q_av[i] + actual_external_abstraction_av[i] + abstraction[i]
        storage_rate = (storage[i] - storage_prev[i]) / dt
        error[i], relative_error[i] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    return nothing
end

"""
Compute water mass balance error and relative error for river (and floodplain) local
inertial routing.
"""
function compute_flow_balance!(
    river_flow::LocalInertialRiverFlow,
    water_balance::MassBalance,
    network::NetworkRiver,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; inwater, external_inflow, actual_external_abstraction_av, abstraction) =
        river_flow.boundary_conditions
    (; edges_at_node) = network

    for i in river_flow.parameters.active_n
        total_in = 0.0
        total_out = 0.0
        q_src = sum_at(river_flow.variables.q_av, edges_at_node.src[i])
        total_in, total_out = add_inflow(total_in, total_out, [q_src, inwater[i]])
        total_in += max(0.0, external_inflow[i])
        q_dst = sum_at(river_flow.variables.q_av, edges_at_node.dst[i])
        total_in, total_out = add_outflow(total_in, total_out, q_dst)
        total_out += actual_external_abstraction_av[i] + abstraction[i]
        storage = river_flow.variables.storage[i]
        if !isnothing(river_flow.floodplain)
            storage += river_flow.floodplain.variables.storage[i]
        end
        storage_rate = (storage - storage_prev[i]) / dt
        error[i], relative_error[i] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    return nothing
end

"Helper function to add `inflow` to incoming flow `qin` or outgoing flow `qout`"
function add_inflow(qin, qout, inflow)
    for q in inflow
        if q > 0.0
            qin += q
        else
            qout -= q
        end
    end
    return qin, qout
end

"Helper function to add `outflow` to incoming flow `qin` or outgoing flow `qout`"
function add_outflow(qin, qout, outflow)
    for q in outflow
        if q > 0.0
            qout += q
        else
            qin -= q
        end
    end
    return qin, qout
end

"""
Compute water mass balance error and relative error for overland flow kinematic wave
routing.
"""
function compute_flow_balance!(
    overland_flow::KinWaveOverlandFlow,
    water_balance::MassBalance,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; inwater) = overland_flow.boundary_conditions
    (; qin_av, q_av, storage) = overland_flow.variables

    for i in eachindex(storage_prev)
        total_in = inwater[i] + qin_av[i]
        total_out = q_av[i]
        storage_rate = (storage[i] - storage_prev[i]) / dt
        error[i], relative_error[i] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    return nothing
end

"""
Compute water mass balance error and relative error for 1D river local inertial routing
(subgrid channel) and 2D overland flow local inertial routing. Water balance errors are
computed for each land cell (total storage) considering both river and overland flow.
"""
function compute_flow_balance!(
    river_flow::LocalInertialRiverFlow,
    overland_flow::LocalInertialOverlandFlow,
    water_balance::MassBalance,
    domain::Domain,
    dt::Float64,
)
    indices = domain.land.network.edge_indices
    inds_river = domain.land.network.river_indices
    (; storage_prev, error, relative_error) = water_balance
    (; river_location, reservoir_outlet) = domain.land.parameters
    (; edges_at_node) = domain.river.network
    (; qx_av, qy_av, storage) = overland_flow.variables
    (; runoff) = overland_flow.boundary_conditions
    (; external_inflow, actual_external_abstraction_av, abstraction) =
        river_flow.boundary_conditions

    for i in 1:(overland_flow.parameters.n)
        yd = indices.yd[i]
        xd = indices.xd[i]
        total_in = 0.0
        total_out = 0.0
        total_in, total_out =
            add_inflow(total_in, total_out, [qx_av[xd], qy_av[yd], runoff[i]])
        total_in, total_out = add_outflow(total_in, total_out, [qx_av[i], qy_av[i]])
        if river_location[i]
            reservoir_outlet[i] && continue
            k = inds_river[i]
            q_src = sum_at(river_flow.variables.q_av, edges_at_node.src[k])
            total_in, total_out = add_inflow(total_in, total_out, q_src)
            total_in += max(0.0, external_inflow[k])
            q_dst = sum_at(river_flow.variables.q_av, edges_at_node.dst[k])
            total_in, total_out = add_outflow(total_in, total_out, q_dst)
            total_out += actual_external_abstraction_av[k] + abstraction[k]
        end
        storage_rate = (storage[i] - storage_prev[i]) / dt
        error[i], relative_error[i] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
end

"""
Compute water mass balance error and relative error for subsurface flow kinematic wave
routing.
"""
function compute_flow_balance!(
    subsurface_flow::LateralSSF,
    water_balance::MassBalance,
    parameters::LandParameters,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; storage, ssfin, ssf, exfiltwater) = subsurface_flow.variables
    (; recharge) = subsurface_flow.boundary_conditions
    (; flow_length, area) = parameters

    f_conv = dt / tosecond(BASETIMESTEP)
    for i in eachindex(storage_prev)
        total_in = ssfin[i] * f_conv
        total_out = ssf[i] * f_conv + exfiltwater[i] * area[i]
        total_in, total_out = add_inflow(total_in, total_out, recharge[i] * flow_length[i])
        storage_rate = (storage[i] - storage_prev[i])
        error[i], relative_error[i] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
end

"""
Compute water mass balance error and relative error for single layer subsurface flow
`GroundwaterFlow` in an unconfined aquifer. Errors for subsurface flow constant head
boundaries are set at zero.
"""
function compute_flow_balance!(
    subsurface_flow::GroundwaterFlow{A},
    water_balance::MassBalance,
    parameters::LandParameters,
    dt::Float64,
) where {A <: UnconfinedAquifer}
    (; storage_prev, error, relative_error) = water_balance
    (; storage, q_in_av, q_out_av, exfiltwater) = subsurface_flow.aquifer.variables
    (; area) = subsurface_flow.aquifer.parameters

    n = length(storage_prev)
    flux_in = zeros(n)
    flux_out = zeros(n)
    flux_in, flux_out = sum_boundary_fluxes(subsurface_flow)

    f_conv = dt / tosecond(BASETIMESTEP)
    for i in eachindex(storage_prev)
        total_in = (q_in_av[i] + flux_in[i]) * f_conv
        total_out = f_conv * (q_out_av[i] + flux_out[i]) + exfiltwater[i] * area[i]
        storage_rate = (storage[i] - storage_prev[i])
        error[i], relative_error[i] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    index_const_head = subsurface_flow.constanthead.index
    error[index_const_head] .= 0.0
    relative_error[index_const_head] .= 0.0
    return nothing
end

"""
Compute water mass balance error and relative error for river, overland, reservoir and
subsurface flow routing.
"""
function compute_flow_routing_balance!(model)
    (; river_flow, overland_flow, subsurface_flow) = model.routing
    (; reservoir) = river_flow.boundary_conditions
    (;
        river_water_balance,
        reservoir_water_balance,
        overland_water_balance,
        subsurface_water_balance,
    ) = model.mass_balance.routing
    dt = tosecond(model.clock.dt)

    compute_flow_balance!(river_flow, river_water_balance, model.domain.river.network, dt)
    compute_flow_balance!(reservoir, reservoir_water_balance, dt)
    compute_flow_balance!(overland_flow, overland_water_balance, dt)
    compute_flow_balance!(
        subsurface_flow,
        subsurface_water_balance,
        model.domain.land.parameters,
        dt,
    )
end

"""
Compute water mass balance error and relative error for combined river and overland,
reservoir and subsurface flow routing. Combined river and overland flow consists of local
inertial 1D river flow routing (subgrid channel) and local inertial 2D overland flow
routing.
"""
function compute_flow_routing_balance!(
    model::Model{R},
) where {R <: Routing{<:LocalInertialOverlandFlow, <:LocalInertialRiverFlow}}
    (; river_flow, overland_flow, subsurface_flow) = model.routing
    (; reservoir) = river_flow.boundary_conditions
    (; overland_water_balance, reservoir_water_balance, subsurface_water_balance) =
        model.mass_balance.routing
    dt = tosecond(model.clock.dt)

    compute_flow_balance!(
        river_flow,
        overland_flow,
        overland_water_balance,
        model.domain,
        dt,
    )
    compute_flow_balance!(reservoir, reservoir_water_balance, dt)
    compute_flow_balance!(
        subsurface_flow,
        subsurface_water_balance,
        model.domain.land.parameters,
        dt,
    )
end

"""
    compute_mass_balance!(model, ::HydrologicalMassBalance)
    compute_mass_balance!(model, ::NoMassBalance)
Compute water mass balance error and relative error if `model` contains a mass balance
`HydrologicalMassBalance`, skip computations with mass balance `NoMassBalance`.
"""
function compute_mass_balance!(model, ::HydrologicalMassBalance)
    compute_land_hydrology_balance!(model)
    compute_flow_routing_balance!(model)
    return nothing
end

function compute_mass_balance!(model, ::NoMassBalance)
    return nothing
end
