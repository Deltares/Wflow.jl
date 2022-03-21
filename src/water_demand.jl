abstract type NonIrrigationWaterDemand end

@get_units @with_kw struct Industry{T} <: NonIrrigationWaterDemand
    demand_gross::Vector{T}   # gross industry water demand [mm Δt⁻¹]
    demand_net::Vector{T}     # net industry water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} # return flow fraction [-]
end

@get_units @with_kw struct Domestic{T} <: NonIrrigationWaterDemand
    demand_gross::Vector{T}   # gross domestic water demand [mm Δt⁻¹]
    demand_net::Vector{T}     # net domestic water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} # return flow fraction [-]
end

@get_units @with_kw struct Livestock{T} <: NonIrrigationWaterDemand
    demand_gross::Vector{T}  # gross livestock water demand [mm Δt⁻¹]
    demand_net::Vector{T}    # net livestock water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} # return flow fraction [-]
end

function update(non_irrigation::N) where {N<:NonIrrigationWaterDemand}
    for i = 1:length(non_irrigation.demand_gross)
        fraction = non_irrigation.demand_net[i] / non_irrigation.demand_gross[i]
        non_irrigation.returnflow_fraction[i] = 1.0 - fraction == Inf ? 0.0 : fraction
    end
    return nothing
end

function initialize_domestic_demand(nc, config, inds, Δt)

    demand_gross = ncread(
        nc,
        config,
        "vertical.domestic.demand_gross";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)
    demand_net = ncread(
        nc,
        config,
        "vertical.domestic.demand_net";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)

    domestic = Domestic{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = fill(mv, length(inds)),
    )
    
    return domestic
end

function initialize_industry_demand(nc, config, inds, Δt)

    demand_gross = ncread(
        nc,
        config,
        "vertical.industry.demand_gross";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)
    demand_net = ncread(
        nc,
        config,
        "vertical.industry.demand_net";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)

    industry = Industry{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = fill(mv, length(inds)),
    )
    
    return industry
end

function initialize_livestock_demand(nc, config, inds, Δt)

    demand_gross = ncread(
        nc,
        config,
        "vertical.livestock.demand_gross";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)
    demand_net = ncread(
        nc,
        config,
        "vertical.livestock.demand_net";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)

    livestock = Livestock{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = fill(mv, length(inds)),
    )
    
    return livestock
end