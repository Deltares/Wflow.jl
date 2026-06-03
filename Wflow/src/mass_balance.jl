
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
model time step dt for a hydrological model.

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
            subsurface_water_balance = MassBalance(subsurface_flow, n_land),
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
function compute_total_storage(land_hydrology_model::LandHydrologySBM, i::Int)
    (; soil, interception, snow, glacier, demand) = land_hydrology_model
    (; total_soil_water_storage) = soil.variables
    (; canopy_storage) = interception.variables

    snow_storage = get_snow_storage(snow)[i] + get_snow_water(snow)[i]
    glacier_storage = get_glacier_store(glacier)[i] * get_glacier_fraction(glacier)[i]
    paddy_storage = get_water_depth(demand.paddy)[i]

    total_storage =
        total_soil_water_storage[i] +
        canopy_storage[i] +
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
    for i in eachindex(storage_prev)
        storage_prev[i] = compute_total_storage(land_hydrology_model, i)
    end
    return nothing
end

"""
    get_storage(river_flow_model::RiverFlowModel, i)

Return storage of a river flow model at index `i`, floodplain storage is added to river
storage if an optional floodplain is included.
"""
function get_storage(river_flow_model::AbstractRiverFlowModel, i)
    (; storage) = river_flow_model.variables
    if isnothing(river_flow_model.floodplain)
        return storage[i]
    else
        total_storage = storage[i] + river_flow_model.floodplain.variables.storage[i]
        return total_storage
    end
end

"""
Save river (+ floodplain) storage at previous time step as `storage_prev` of river
`water_balance`.
"""
function storage_prev!(river_flow_model::AbstractRiverFlowModel, water_balance::MassBalance)
    (; storage_prev) = water_balance
    for i in eachindex(storage_prev)
        storage_prev[i] = get_storage(river_flow_model, i)
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
    (; river_water_balance, reservoir_water_balance, overland_water_balance) =
        model.mass_balance.routing
    (; land_water_balance) = model.mass_balance
    (; land) = model
    (; overland_flow, river_flow) = model.routing
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

"Compute total incoming vertical flux of land hydrology `SBM` at index `i`."
function vertical_in(land_hydrology_model::LandHydrologySBM, i::Int)
    (; atmospheric_forcing, allocation) = land_hydrology_model
    (; precipitation) = atmospheric_forcing
    total_in = precipitation[i] + get_irrigation_allocated(allocation)[i]
    return total_in
end

"Compute total outgoing vertical flux of land hydrology `SBM` at index `i`."
function vertical_out(land_hydrology_model::LandHydrologySBM, i::Int)
    (; allocation, soil, runoff) = land_hydrology_model
    (; net_runoff, actual_evapotranspiration, actual_leakage) = soil.variables
    (; net_runoff_river) = runoff.variables
    total_out =
        net_runoff[i] +
        actual_evapotranspiration[i] +
        net_runoff_river[i] +
        actual_leakage[i] +
        get_groundwater_abstraction_flux(allocation)[i]
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
    dt = tosecond(model.clock.dt)

    # exclude recharge from computing total incoming and outgoing boundary fluxes for
    # groundwaterflow, other boundaries are required for the total soil water balance.
    boundaries_flow_in, boundaries_flow_out =
        sum_boundary_fluxes(subsurface_flow, model.domain; exclude = RechargeModel)

    for i in eachindex(storage_prev)
        area_i = area[i]
        subsurface_flux_in = subsurface_flow.variables.q_in_average[i] / area_i
        total_in =
            subsurface_flux_in +
            vertical_in(model.land, i) +
            get_snow_in(snow)[i] +
            boundaries_flow_in[i] / area[i]

        subsurface_flux_out = subsurface_flow.variables.q_average[i] / area_i
        vertical_flux_out = vertical_out(model.land, i)
        total_out =
            subsurface_flux_out +
            vertical_flux_out +
            get_snow_out(snow)[i] +
            boundaries_flow_out[i] / area[i]
        storage = compute_total_storage(model.land, i)
        storage_rate = (storage - storage_prev[i]) / dt
        error[i], relative_error[i] =
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
    (; storage, outflow_average, actevap_average, actevap_cumulative) =
        reservoir_model.variables
    (; precipitation, inflow_average, inflow_cumulative) =
        reservoir_model.boundary_conditions
    (; area) = reservoir_model.parameters

    @. inflow_average = inflow_cumulative / dt
    @. actevap_average = actevap_cumulative / dt

    for i in eachindex(storage_prev)
        total_in = inflow_average[i] + precipitation[i] * area[i]
        total_out = outflow_average[i] + actevap_average[i] * area[i]
        storage_rate = (storage[i] - storage_prev[i]) / dt
        error[i], relative_error[i] =
            compute_mass_balance_error(total_in, total_out, storage_rate)
    end
    return nothing
end

compute_flow_balance!(reservoir_model::Nothing, water_balance::NoMassBalance, dt::Float64) =
    nothing

"Compute water mass balance error and relative error for river kinematic wave routing."
function compute_flow_balance!(
    river_flow_model::RiverFlowModel{<:KinematicWave},
    water_balance::MassBalance,
    network::NetworkRiver,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; inwater, external_inflow, actual_external_abstraction_average, abstraction) =
        river_flow_model.boundary_conditions
    (; qin_average, q_average) = river_flow_model.variables

    for i in eachindex(storage_prev)
        total_in = inwater[i] + qin_average[i] + max(0.0, external_inflow[i])
        total_out = q_average[i] + actual_external_abstraction_average[i] + abstraction[i]
        storage = river_flow_model.variables.storage[i]
        if !isnothing(river_flow_model.floodplain)
            storage += river_flow_model.floodplain.variables.storage[i]
        end
        storage_rate = (storage - storage_prev[i]) / dt
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
    river_flow_model::RiverFlowModel{<:AbstractStaggeredRoutingMethod},
    water_balance::MassBalance,
    network::NetworkRiver,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; inwater, external_inflow, actual_external_abstraction_average, abstraction) =
        river_flow_model.boundary_conditions
    (; edges_at_node) = network

    (; q_average) = river_flow_model.variables

    for i in river_flow_model.parameters.active_n
        total_in = 0.0
        total_out = 0.0
        q_src = sum_at(q_average, edges_at_node.src[i])
        total_in, total_out = add_inflow(total_in, total_out, [q_src, inwater[i]])
        total_in += max(0.0, external_inflow[i])
        q_dst = sum_at(q_average, edges_at_node.dst[i])
        total_in, total_out = add_outflow(total_in, total_out, q_dst)
        total_out += actual_external_abstraction_average[i] + abstraction[i]
        storage = river_flow_model.variables.storage[i]

        if !isnothing(river_flow_model.floodplain)
            storage += river_flow_model.floodplain.variables.storage[i]
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
    overland_flow_model::OverlandFlowModel{<:KinematicWave},
    water_balance::MassBalance,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; inwater) = overland_flow_model.boundary_conditions
    (; qin_average, q_average, storage) = overland_flow_model.variables

    for i in eachindex(storage_prev)
        total_in = inwater[i] + qin_average[i]
        total_out = q_average[i]
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
    river_flow_model::RiverFlowModel{<:LocalInertial},
    overland_flow_model::OverlandFlowModel{<:LocalInertial},
    water_balance::MassBalance,
    domain::Domain,
    dt::Float64,
)
    indices = domain.land.network.edge_indices
    inds_river = domain.land.network.river_indices
    (; storage_prev, error, relative_error) = water_balance
    (; river_location, reservoir_outlet) = domain.land.parameters
    (; edges_at_node) = domain.river.network
    (; qx_average, qy_average, storage) = overland_flow_model.variables
    (; runoff) = overland_flow_model.boundary_conditions
    (; external_inflow, actual_external_abstraction_average, abstraction) =
        river_flow_model.boundary_conditions

    q_av_average = river_flow_model.variables.q_average
    actual_external_abstraction_av_average = actual_external_abstraction_average
    qx_av_average = qx_average
    qy_av_average = qy_average

    for i in 1:(overland_flow_model.parameters.n)
        idx_down = indices.idx_down[i]
        idx_left = indices.idx_left[i]
        total_in = 0.0
        total_out = 0.0
        total_in, total_out = add_inflow(
            total_in,
            total_out,
            [qx_av_average[idx_left], qy_av_average[idx_down], runoff[i]],
        )
        total_in, total_out =
            add_outflow(total_in, total_out, [qx_av_average[i], qy_av_average[i]])
        if river_location[i]
            reservoir_outlet[i] && continue
            k = inds_river[i]
            q_src = sum_at(q_av_average, edges_at_node.src[k])
            total_in, total_out = add_inflow(total_in, total_out, q_src)
            total_in += max(0.0, external_inflow[k])
            q_dst = sum_at(q_av_average, edges_at_node.dst[k])
            total_in, total_out = add_outflow(total_in, total_out, q_dst)
            total_out += actual_external_abstraction_av_average[k] + abstraction[k]
        end
        storage_rate = (storage[i] - storage_prev[i]) / dt
        error[i], relative_error[i] =
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
    (; q_in_average, q_average, q_net_average) = subsurface_flow_model.variables

    flux_in, flux_out = sum_boundary_fluxes(subsurface_flow_model, domain)

    for i in eachindex(q_net_average)
        total_in = q_in_average[i] + flux_in[i]
        total_out = q_average[i] + flux_out[i]
        storage_rate = q_net_average[i]
        error[i], relative_error[i] =
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

function compute_flow_routing_balance!(
    model::Model{R},
) where {
    R <: Routing{<:OverlandFlowModel{<:LocalInertial}, <:RiverFlowModel{<:LocalInertial}},
}
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
function compute_mass_balance!(model::AbstractModel, ::HydrologicalMassBalance)
    compute_land_hydrology_balance!(model)
    compute_flow_routing_balance!(model)
    return nothing
end

function compute_mass_balance!(model, ::NoMassBalance)
    return nothing
end
