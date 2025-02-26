struct WaterBalanceVertical{T <: AbstractFloat} # Terms in [mm]
    storage_total::Vector{T}
    storage_total_prev::Vector{T}
    inflow_total::Vector{T}
    outflow_total::Vector{T}
    balance_error::Vector{T}
    relative_error::Vector{T}
end

function WaterBalanceVertical(n::Integer, float_type::Type{T}) where {T <: AbstractFloat}
    WaterBalanceVertical(ntuple(_ -> fill(float_type(Nan), n), 6)...)
end

function compute_water_balance!(model, water_balance::WaterBalanceVertical)::Nothing
    (;
        storage_total,
        storage_total_prev,
        inflow_total,
        outflow_total,
        balance_error,
        relative_error,
    ) = water_balance
    (; vertical, lateral) = model

    (; soil, runoff) = vertical
    (; satwaterdepth, ustoredepth, actevap) = soil.variables
    (; runoff_land, net_runoff_river) = runoff.variables

    (; subsurface) = lateral
    (; ssf, ssfin) = subsurface.variables
    (; flow_length, flow_width) = subsurface.parameters

    storage_total .= satwaterdepth
    storage_total .+= ustoredepth

    inflow_total .= vertical.atmospheric_forcing.precipitation
    inflow_total .+= runoff.variables.net_runoff_river
    inflow_total .+= 1000.0 * ssfin ./ (flow_length .* flow_width)

    outflow_total .= runoff_land
    outflow_total .+= net_runoff_river
    outflow_total .+= actevap
    outflow_total .+= 1000.0 * ssf ./ (flow_length .* flow_width)

    @. balance_error = (inflow_total - outflow_total - (storage_total - storage_total_prev))
    @. relative_error = balance_error / ((inflow_total + outflow_total) / 2)

    copyto!(storage_total_prev, storage_total)

    #TODO : What to do with results?

    return nothing
end

struct WaterBalanceOverLand{T <: AbstractFloat} # Terms in [?]
    volume_t_prev::Vector{T}
    volume_rate::Vector{T}
    water_in::Vector{T}
    water_out::Vector{T}
    balance_error::Vector{T}
    relative_error::Vector{T}
end

function WaterBalanceOverLand(n::Integer, float_type::Type{T}) where {T <: AbstractFloat}
    WaterBalanceOverLand(ntuple(_ -> fill(float_type(Nan), n), 6)...)
end

# function compute_water_balance!(model, water_balance::WaterBalanceOverLand)::Nothing
#     (; vertical, lateral) = model
#     (; volume_t_prev, volume_rate, water_in, water_out, balance_error, relative_error) = water_balance_overland

#     (; land, river) = lateral
#     (; volume) = land.variables
#     volume_t = volume

#     water_in .= lateral.land.qin_av
#     @. water_in += lateral.land.qlat_av * lateral.land.dl

#     water_out .= lateral.land.q_av

#     @. volume_rate = (volume_t - volume_t_previous) / vertical.dt

#     @. balance_error = (water_in - water_out) - volume_rate

#     #TODO : What to do with results?

#     return nothing
# end

struct WaterBalance{T}
    water_balance_vertical::WaterBalanceVertical{T}
    water_balance_over_land::WaterBalanceOverLand{T}
end

function WaterBalance(n::Integer, float_type::Type{T}) where {T <: AbstractFloat}
    WaterBalance(WaterBalanceVertical(n, float_type), WaterBalanceOverLand(n, float_type))
end