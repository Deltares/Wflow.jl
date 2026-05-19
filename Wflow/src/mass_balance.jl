
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

function MassBalance(::AbstractSubsurfaceFlowModel, n::Int)
    mass_balance = MassBalance(; n, storage_prev = [])
    return mass_balance
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
function HydrologicalMassBalance(domain::Domain, subsurface_flow, config::Config)
    (; river_routing, land_routing, water_mass_balance__flag) = config.model
    if water_mass_balance__flag
        n_land_cells = length(domain.land.network.land_indices_2d)
        n_river_cells = length(domain.river.network.river_indices_2d)
        if config.model.reservoir__flag
            n_reservoir = length(domain.reservoir.network.outlet_indices_2d)
            reservoir_water_balance = MassBalance(; n = n_reservoir)
        else
            reservoir_water_balance = NoMassBalance()
        end
        if land_routing == RoutingType.local_inertial &&
           river_routing == RoutingType.local_inertial
            river_water_balance = NoMassBalance()
        else
            river_water_balance = MassBalance(; n = n_river_cells)
        end
        routing = FlowRoutingMassBalance(;
            river_water_balance,
            reservoir_water_balance,
            overland_water_balance = MassBalance(; n = n_land_cells),
            subsurface_water_balance = MassBalance(subsurface_flow, n_land_cells),
        )
        mass_balance = HydrologicalMassBalance(;
            land_water_balance = MassBalance(; n = n_land_cells),
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

"Compute total storage of land hydrolology model `LandHydrologySBM` at `cell_idx`."
function compute_total_storage(land_hydrology_model::LandHydrologySBM, cell_idx::Int)
    (; soil, interception, snow, glacier, demand) = land_hydrology_model
    (; total_soilwater_storage) = soil.variables
    (; canopy_storage) = interception.variables

    snow_storage = get_snow_storage(snow)[cell_idx] + get_snow_water(snow)[cell_idx]
    glacier_storage =
        get_glacier_store(glacier)[cell_idx] * get_glacier_fraction(glacier)[cell_idx]
    paddy_storage = get_water_depth(demand.paddy)[cell_idx]

    total_storage =
        total_soilwater_storage[cell_idx] +
        canopy_storage[cell_idx] +
        snow_storage +
        glacier_storage +
        paddy_storage

    return total_storage
end

"Compute total storage of land hydrology model `LandHydrologySBM`."
function compute_total_storage!(
    land_hydrology_model::LandHydrologySBM,
    water_balance::MassBalance,
)
    (; storage_prev) = water_balance
    for cell_idx in eachindex(storage_prev)
        storage_prev[cell_idx] = compute_total_storage(land_hydrology_model, cell_idx)
    end
    return nothing
end

"""
    get_storage(river_flow_model::LocalInertialRiverFlowModel, cell_idx)
    get_storage(river_flow_model::KinWaveRiverFlowModel, cell_idx)

Return storage of a river flow model at `cell_idx`. For `LocalInertialRiverFlowModel` floodplain
storage is added to river storage if an optional floodplain is included.
"""
function get_storage(river_flow_model::LocalInertialRiverFlowModel, cell_idx::Int)
    (; storage) = river_flow_model.variables
    if isnothing(river_flow_model.floodplain)
        return storage[cell_idx]
    else
        total_storage =
            storage[cell_idx] + river_flow_model.floodplain.variables.storage[cell_idx]
        return total_storage
    end
end
get_storage(river_flow_model::KinWaveRiverFlowModel, cell_idx) =
    river_flow_model.variables.storage[cell_idx]

"""
Save river (+ floodplain) storage at previous time step as `storage_prev` of river
`water_balance`.
"""
function storage_prev!(river_flow_model::AbstractRiverFlowModel, water_balance::MassBalance)
    (; storage_prev) = water_balance
    for cell_idx in eachindex(storage_prev)
        storage_prev[cell_idx] = get_storage(river_flow_model, cell_idx)
    end
    return nothing
end

"""
Save reservoir storage at previous time step as `storage_prev` of reservoir `water_balance`.
"""
function storage_prev!(reservoir_model::ReservoirModel, water_balance::MassBalance)
    water_balance.storage_prev .= reservoir_model.variables.storage
end

"""
    storage_prev!(model, ::HydrologicalMassBalance)
    storage_prev!(model, ::NoMassBalance)

Save storage at previous time step as `storage_prev` for river, reservoir and overland flow
routing and `land` hydrology model as part of water mass balance `HydrologicalMassBalance`.
For `NoMassBalance` storage at previous time step is not required.
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
    storage_prev!(reservoir, reservoir_water_balance)
    return nothing
end

function storage_prev!(model, ::NoMassBalance)
    return nothing
end

"Compute total incoming vertical flux of land hydrology `SBM` at `cell_idx`."
function vertical_in(land_hydrology_model::LandHydrologySBM, cell_idx::Int)
    (; atmospheric_forcing, allocation) = land_hydrology_model
    (; precipitation) = atmospheric_forcing
    total_in = precipitation[cell_idx] + get_irrigation_allocated(allocation)[cell_idx]
    return total_in
end

"Compute total outgoing vertical flux of land hydrology `SBM` at `cell_idx`."
function vertical_out(land_hydrology_model::LandHydrologySBM, cell_idx::Int)
    (; allocation, soil, runoff) = land_hydrology_model
    (; net_runoff, actevap, actleakage) = soil.variables
    (; net_runoff_river) = runoff.variables
    total_out =
        net_runoff[cell_idx] +
        actevap[cell_idx] +
        net_runoff_river[cell_idx] +
        actleakage[cell_idx] +
        get_groundwater_abstraction_flux(allocation)[cell_idx]
    return total_out
end

"""
Compute water mass balance error and relative error for `land` hydrology `SBM`. Errors for
subsurface flow constant head boundaries are set at zero.
"""
function compute_land_hydrology_balance!(
    model::AbstractModel{<:Union{SbmModel, SbmGwfModel}},
)
    (; land_water_balance) = model.mass_balance
    (; storage_prev, error, relative_error) = land_water_balance
    (; snow) = model.land
    (; area) = model.domain.land.parameters
    (; subsurface_flow) = model.routing

    # exclude recharge from computing total incoming and outgoing boundary fluxes for
    # groundwaterflow, other boundaries are required for the total soil water balance.
    boundaries_flow_in, boundaries_flow_out =
        sum_boundary_fluxes(subsurface_flow, model.domain; exclude = RechargeModel)

    for cell_idx in eachindex(storage_prev)
        f_conv = (model.clock.dt / BASETIMESTEP) / (area[cell_idx] * 0.001)
        subsurface_flux_in = subsurface_flow.variables.q_in_av[cell_idx] * f_conv
        total_in =
            subsurface_flux_in +
            vertical_in(model.land, cell_idx) +
            get_snow_in(snow)[cell_idx] +
            boundaries_flow_in[cell_idx] * f_conv
        subsurface_flux_out = subsurface_flow.variables.q_av[cell_idx] * f_conv
        vertical_flux_out = vertical_out(model.land, cell_idx)
        total_out =
            subsurface_flux_out +
            vertical_flux_out +
            get_snow_out(snow)[cell_idx] +
            boundaries_flow_out[cell_idx] * f_conv
        storage = compute_total_storage(model.land, cell_idx)
        storage_rate = storage - storage_prev[cell_idx]
        error[cell_idx], relative_error[cell_idx] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    constant_head_boundary_error!(subsurface_flow, land_water_balance)
    return nothing
end

"""
    compute_flow_balance!(reservoir_model::ReservoirModel, water_balance::MassBalance, dt::Float64)
    compute_flow_balance!(reservoir_model::Nothing, water_balance::NoMassBalance, dt::Float64)

Compute reservoir water mass balance error and relative error if reservoirs are included.
"""
function compute_flow_balance!(
    reservoir_model::ReservoirModel,
    water_balance::MassBalance,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; storage, outflow_av, actevap) = reservoir_model.variables
    (; precipitation, inflow) = reservoir_model.boundary_conditions
    (; area) = reservoir_model.parameters
    n_reservoirs = length(storage_prev)

    for reservoir_idx in 1:n_reservoirs
        total_in =
            inflow[reservoir_idx] +
            (precipitation[reservoir_idx] * 0.001 * area[reservoir_idx]) / dt
        total_out =
            outflow_av[reservoir_idx] +
            (actevap[reservoir_idx] * 0.001 * area[reservoir_idx]) / dt
        storage_rate = (storage[reservoir_idx] - storage_prev[reservoir_idx]) / dt
        error[reservoir_idx], relative_error[reservoir_idx] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    return nothing
end

compute_flow_balance!(reservoir_model::Nothing, water_balance::NoMassBalance, dt::Float64) =
    nothing

"Compute water mass balance error and relative error for river kinematic wave routing."
function compute_flow_balance!(
    river_flow_model::KinWaveRiverFlowModel,
    water_balance::MassBalance,
    network::NetworkRiver,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; inwater, external_inflow, actual_external_abstraction_av, abstraction) =
        river_flow_model.boundary_conditions
    (; qin_av, q_av, storage, n_cells) = river_flow_model.variables

    for cell_idx in 1:n_cells
        total_in =
            inwater[cell_idx] + qin_av[cell_idx] + max(0.0, external_inflow[cell_idx])
        total_out =
            q_av[cell_idx] +
            actual_external_abstraction_av[cell_idx] +
            abstraction[cell_idx]
        storage_rate = (storage[cell_idx] - storage_prev[cell_idx]) / dt
        error[cell_idx], relative_error[cell_idx] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    return nothing
end

"""
Compute water mass balance error and relative error for river (and floodplain) local
inertial routing.
"""
function compute_flow_balance!(
    river_flow_model::LocalInertialRiverFlowModel,
    water_balance::MassBalance,
    network::NetworkRiver,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; inwater, external_inflow, actual_external_abstraction_av, abstraction) =
        river_flow_model.boundary_conditions
    (; edges_at_node) = network

    for cell_idx in river_flow_model.parameters.active_n
        total_in = 0.0
        total_out = 0.0
        q_src = sum_at(river_flow_model.variables.q_av, edges_at_node.src[cell_idx])
        total_in, total_out = add_inflow(total_in, total_out, [q_src, inwater[cell_idx]])
        total_in += max(0.0, external_inflow[cell_idx])
        q_dst = sum_at(river_flow_model.variables.q_av, edges_at_node.dst[cell_idx])
        total_in, total_out = add_outflow(total_in, total_out, q_dst)
        total_out += actual_external_abstraction_av[cell_idx] + abstraction[cell_idx]
        storage = river_flow_model.variables.storage[cell_idx]
        if !isnothing(river_flow_model.floodplain)
            storage += river_flow_model.floodplain.variables.storage[cell_idx]
        end
        storage_rate = (storage - storage_prev[cell_idx]) / dt
        error[cell_idx], relative_error[cell_idx] =
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
    overland_flow_model::KinWaveOverlandFlowModel,
    water_balance::MassBalance,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; inwater) = overland_flow_model.boundary_conditions
    (; qin_av, q_av, storage, n_cells) = overland_flow_model.variables

    for cell_idx in 1:n_cells
        total_in = inwater[cell_idx] + qin_av[cell_idx]
        total_out = q_av[cell_idx]
        storage_rate = (storage[cell_idx] - storage_prev[cell_idx]) / dt
        error[cell_idx], relative_error[cell_idx] =
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
    river_flow_model::LocalInertialRiverFlowModel,
    overland_flow_model::LocalInertialOverlandFlowModel,
    water_balance::MassBalance,
    domain::Domain,
    dt::Float64,
)
    indices = domain.land.network.edge_indices
    inds_river = domain.land.network.river_cell_indices
    (; storage_prev, error, relative_error) = water_balance
    (; river_location, reservoir_outlet) = domain.land.parameters
    (; edges_at_node) = domain.river.network
    (; qx_av, qy_av, storage) = overland_flow_model.variables
    (; runoff) = overland_flow_model.boundary_conditions
    (; external_inflow, actual_external_abstraction_av, abstraction) =
        river_flow_model.boundary_conditions

    for cell_idx in 1:(overland_flow_model.parameters.n_cells)
        yd = indices.yd[cell_idx]
        xd = indices.xd[cell_idx]
        total_in = 0.0
        total_out = 0.0
        total_in, total_out =
            add_inflow(total_in, total_out, [qx_av[xd], qy_av[yd], runoff[cell_idx]])
        total_in, total_out =
            add_outflow(total_in, total_out, [qx_av[cell_idx], qy_av[cell_idx]])
        if river_location[cell_idx]
            reservoir_outlet[cell_idx] && continue
            river_cell_idx = inds_river[cell_idx]
            q_src =
                sum_at(river_flow_model.variables.q_av, edges_at_node.src[river_cell_idx])
            total_in, total_out = add_inflow(total_in, total_out, q_src)
            total_in += max(0.0, external_inflow[river_cell_idx])
            q_dst =
                sum_at(river_flow_model.variables.q_av, edges_at_node.dst[river_cell_idx])
            total_in, total_out = add_outflow(total_in, total_out, q_dst)
            total_out +=
                actual_external_abstraction_av[river_cell_idx] + abstraction[river_cell_idx]
        end
        storage_rate = (storage[cell_idx] - storage_prev[cell_idx]) / dt
        error[cell_idx], relative_error[cell_idx] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
end

function constant_head_boundary_error!(
    model::GroundwaterFlowModel,
    water_balance::MassBalance,
)
    (; error, relative_error) = water_balance

    index_const_head = model.constanthead.index
    error[index_const_head] .= 0.0
    relative_error[index_const_head] .= 0.0
    return nothing
end

constant_head_boundary_error!(model::LateralSSFModel, water_balance::MassBalance) = nothing

"Compute water mass balance error and relative error for a subsurface flow model."
function compute_flow_balance!(
    subsurface_flow_model::AbstractSubsurfaceFlowModel,
    water_balance::MassBalance,
    domain::Domain,
    dt::Float64,
)
    (; error, relative_error) = water_balance
    (; q_in_av, q_av, q_net_av, n_cells) = subsurface_flow_model.variables

    flux_in, flux_out = sum_boundary_fluxes(subsurface_flow_model, domain)

    f_conv = dt / tosecond(BASETIMESTEP)
    for cell_idx in 1:n_cells
        total_in = (q_in_av[cell_idx] + flux_in[cell_idx]) * f_conv
        total_out = f_conv * (q_av[cell_idx] + flux_out[cell_idx])
        storage_rate = q_net_av[cell_idx] * f_conv
        error[cell_idx], relative_error[cell_idx] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    constant_head_boundary_error!(subsurface_flow_model, water_balance)
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
    compute_flow_balance!(subsurface_flow, subsurface_water_balance, model.domain, dt)
end

"""
Compute water mass balance error and relative error for combined river and overland,
reservoir and subsurface flow routing. Combined river and overland flow consists of local
inertial 1D river flow routing (subgrid channel) and local inertial 2D overland flow
routing.
"""
function compute_flow_routing_balance!(
    model::Model{R},
) where {R <: Routing{<:LocalInertialOverlandFlowModel, <:LocalInertialRiverFlowModel}}
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
    compute_flow_balance!(subsurface_flow, subsurface_water_balance, model.domain, dt)
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
