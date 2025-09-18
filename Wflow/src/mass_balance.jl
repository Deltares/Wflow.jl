
struct NoMassBalance <: AbstractMassBalance end

@with_kw struct MassBalance
    n::Int
    storage_prev::Vector{Float64} = fill(MISSING_VALUE, n)
    error::Vector{Float64} = fill(MISSING_VALUE, n)
    relative_error::Vector{Float64} = fill(MISSING_VALUE, n)
end

@with_kw struct FlowRoutingMassBalance
    river::MassBalance
    reservoir::Union{MassBalance, NoMassBalance}
    land::MassBalance
    subsurface::MassBalance
end

@with_kw struct HydrologicalMassBalance <: AbstractMassBalance
    land::MassBalance
    routing::FlowRoutingMassBalance
end

function HydrologicalMassBalance(domain::Domain, modelsettings::NamedTuple)
    do_water_mass_balance = modelsettings.water_mass_balance
    if do_water_mass_balance
        do_reservoirs = modelsettings.reservoirs
        n_land = length(domain.land.network.indices)
        n_river = length(domain.river.network.indices)
        if do_reservoirs
            n_reservoir = length(domain.reservoir.network.indices_outlet)
            reservoir = MassBalance(; n = n_reservoir)
        else
            reservoir = NoMassBalance()
        end
        routing = FlowRoutingMassBalance(;
            river = MassBalance(; n = n_river),
            reservoir,
            land = MassBalance(; n = n_land),
            subsurface = MassBalance(; n = n_land),
        )
        mass_balance = HydrologicalMassBalance(; land = MassBalance(; n = n_land), routing)
    else
        mass_balance = NoMassBalance()
    end
    return mass_balance
end

function compute_mass_balance_error(total_input, total_output, storage_rate)
    error = total_input - total_output - storage_rate
    average_flow_rate = (total_input + total_output) / 2.0
    relative_error = iszero(average_flow_rate) ? 0.0 : error / average_flow_rate
    return error, relative_error
end

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

function compute_total_storage!(model::LandHydrologySBM, water_balance::MassBalance)
    (; storage_prev) = water_balance
    for i in eachindex(storage_prev)
        storage_prev[i] = compute_total_storage(model, i)
    end
    return nothing
end

function storage_prev!(model, reservoir::MassBalance)
    (; storage) = model.routing.river_flow.boundary_conditions.reservoir.variables
    reservoir.storage_prev .= storage
end

function storage_prev!(model, ::HydrologicalMassBalance)
    (; river, reservoir, land, subsurface) = model.mass_balance.routing

    compute_total_storage!(model.land, model.mass_balance.land)
    river.storage_prev .= model.routing.river_flow.variables.storage[1:(river.n)]
    land.storage_prev .= model.routing.overland_flow.variables.storage
    subsurface.storage_prev .= get_storage(model.routing.subsurface_flow)
    storage_prev!(model, reservoir)
    return nothing
end

function storage_prev!(model, ::NoMassBalance)
    return nothing
end

function compute_land_hydrology_balance!(
    model::AbstractModel{<:Union{SbmModel, SbmGwfModel}},
)
    (; storage_prev, error, relative_error) = model.mass_balance.land
    (; allocation, snow) = model.land
    (; net_runoff, actevap, actleakage) = model.land.soil.variables
    (; net_runoff_river) = model.land.runoff.variables
    (; precipitation) = model.land.atmospheric_forcing
    (; area) = model.domain.land.parameters
    (; subsurface_flow) = model.routing

    for i in eachindex(storage_prev)
        subsurface_flow_in = get_inflow(subsurface_flow)[i]
        subsurface_flux_in =
            1000.0 * subsurface_flow_in / area[i] * (model.clock.dt / BASETIMESTEP)
        total_input =
            subsurface_flux_in +
            precipitation[i] +
            get_irrigation_allocated(allocation)[i] +
            get_snow_in(snow)[i]

        subsurface_flow_out = get_outflow(subsurface_flow)[i]
        subsurface_flux_out =
            1000.0 * subsurface_flow_out / area[i] * (model.clock.dt / BASETIMESTEP)
        vertical_flux_out =
            net_runoff[i] +
            actevap[i] +
            net_runoff_river[i] +
            actleakage[i] +
            get_groundwater_abstraction_flux(allocation)[i]
        total_output = subsurface_flux_out + vertical_flux_out + get_snow_out(snow)[i]
        storage = compute_total_storage(model.land, i)

        storage_rate = storage - storage_prev[i]
        _error, _relative_error =
            compute_mass_balance_error(total_input, total_output, storage_rate)
        error[i] = _error
        relative_error[i] = _relative_error
    end
    return nothing
end

function compute_flow_balance!(
    reservoir::Reservoir,
    water_balance::MassBalance,
    i::Int,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; storage, outflow_av, actevap) = reservoir.variables
    (; precipitation, inflow) = reservoir.boundary_conditions
    (; area) = reservoir.parameters

    total_input = inflow[i] + (precipitation[i] * 0.001 * area[i]) / dt
    total_output = outflow_av[i] + (actevap[i] * 0.001 * area[i]) / dt
    storage_rate = (storage[i] - storage_prev[i]) / dt
    _error, _relative_error =
        compute_mass_balance_error(total_input, total_output, storage_rate)
    error[i] = _error
    relative_error[i] = _relative_error
end

function compute_flow_balance!(
    river_flow::KinWaveRiverFlow,
    routing_balance::FlowRoutingMassBalance,
    network::NetworkRiver,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = routing_balance.river
    (; reservoir) = river_flow.boundary_conditions
    (; inwater, external_inflow, actual_external_abstraction_av, abstraction) =
        river_flow.boundary_conditions
    (; qin_av, q_av, storage) = river_flow.variables

    for i in eachindex(storage_prev)
        res_index = network.reservoir_indices[i]
        if !isnothing(reservoir) && res_index != 0
            compute_flow_balance!(reservoir, routing_balance.reservoir, res_index, dt)
        end
        total_input = inwater[i] + qin_av[i] + max(0.0, external_inflow[i])
        total_output = q_av[i] + actual_external_abstraction_av[i] + abstraction[i]
        storage_rate = (storage[i] - storage_prev[i]) / dt
        _error, _relative_error =
            compute_mass_balance_error(total_input, total_output, storage_rate)
        error[i] = _error
        relative_error[i] = _relative_error
    end
    return nothing
end

function compute_flow_balance!(
    overland_flow::KinWaveOverlandFlow,
    water_balance::MassBalance,
    dt::Float64,
)
    (; storage_prev, error, relative_error) = water_balance
    (; inwater) = overland_flow.boundary_conditions
    (; qin_av, q_av, storage) = overland_flow.variables

    for i in eachindex(storage_prev)
        total_input = inwater[i] + qin_av[i]
        total_output = q_av[i]
        storage_rate = (storage[i] - storage_prev[i]) / dt
        _error, _relative_error =
            compute_mass_balance_error(total_input, total_output, storage_rate)
        error[i] = _error
        relative_error[i] = _relative_error
    end
    return nothing
end

function compute_flow_routing_balance!(model)
    (; river_flow, overland_flow) = model.routing
    (; routing) = model.mass_balance
    (; land) = routing
    (; network) = model.domain.river
    dt = tosecond(model.clock.dt)
    compute_flow_balance!(river_flow, routing, network, dt)
    compute_flow_balance!(overland_flow, land, dt)
end

function compute_mass_balance!(model, ::HydrologicalMassBalance)
    compute_land_hydrology_balance!(model)
    compute_flow_routing_balance!(model)
    return nothing
end

function compute_mass_balance!(model, ::NoMassBalance)
    return nothing
end
