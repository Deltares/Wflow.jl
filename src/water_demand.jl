@get_units @with_kw struct Industry{T}
    demand_gross::Vector{T}   # gross industry water demand [mm Δt⁻¹]
    demand_net::Vector{T}     # net industry water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} # return flow fraction [-]
end

@get_units @with_kw struct Domestic{T}
    demand_gross::Vector{T}   # gross domestic water demand [mm Δt⁻¹]
    demand_net::Vector{T}     # net domestic water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} # return flow fraction [-]
end

@get_units @with_kw struct Livestock{T}
    demand_gross::Vector{T}  # gross livestock water demand [mm Δt⁻¹]
    demand_net::Vector{T}    # net livestock water demand [mm Δt⁻¹]
    returnflow_fraction::Vector{T} # return flow fraction [-]
end

@get_units @with_kw struct NonPaddy{T}
    demand_gross::Vector{T}
    irrigation_efficiency::Vector{T}
end

@get_units @with_kw struct Paddy{T}
    demand_gross::Vector{T}
    irrigation_efficiency::Vector{T}
    h_min::Vector{T}                    # minimum required water depth in the irrigated rice field [mm]
    h_max::Vector{T}                    # optimal water depth in the irrigated rice fields [mm]
    h_p::Vector{T}                      # water depth when rice field starts spilling water (overflow) [mm]
    h::Vector{T}                        # actual water depth in rice field
end

function set_returnflow_fraction(returnflow_fraction, demand_gross, demand_net)
    for i = 1:length(demand_gross)
        fraction = demand_net[i] / demand_gross[i]
        returnflow_fraction[i] = 1.0 - fraction
    end
    return returnflow_fraction
end

function initialize_domestic_demand(nc, config, inds, Δt)

    demand_gross =
        ncread(
            nc,
            config,
            "vertical.domestic.demand_gross";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)
    demand_net =
        ncread(
            nc,
            config,
            "vertical.domestic.demand_net";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)

    returnflow_fraction =
        set_returnflow_fraction(fill(mv, length(inds)), demand_gross, demand_net)

    domestic = Domestic{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_fraction,
    )

    return domestic
end

function initialize_industry_demand(nc, config, inds, Δt)

    demand_gross =
        ncread(
            nc,
            config,
            "vertical.industry.demand_gross";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)
    demand_net =
        ncread(
            nc,
            config,
            "vertical.industry.demand_net";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)

    returnflow_fraction =
        set_returnflow_fraction(fill(mv, length(inds)), demand_gross, demand_net)

    industry = Industry{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_fraction,
    )

    return industry
end

function initialize_livestock_demand(nc, config, inds, Δt)

    demand_gross =
        ncread(
            nc,
            config,
            "vertical.livestock.demand_gross";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)
    demand_net =
        ncread(
            nc,
            config,
            "vertical.livestock.demand_net";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (Δt / basetimestep)

    returnflow_fraction =
        set_returnflow_fraction(fill(mv, length(inds)), demand_gross, demand_net)

    livestock = Livestock{Float}(
        demand_gross = demand_gross,
        demand_net = demand_net,
        returnflow_fraction = returnflow_fraction,
    )

    return livestock
end

function initialize_paddy(nc, config, inds)

    h_min =
        ncread(nc, config, "vertical.paddy.h_min"; sel = inds, defaults = 0.0, type = Float)
    h_max =
        ncread(nc, config, "vertical.paddy.h_max"; sel = inds, defaults = 0.0, type = Float)
    h_p = ncread(nc, config, "vertical.paddy.h_p"; sel = inds, defaults = 0.0, type = Float)
    efficiency = ncread(
        nc,
        config,
        "vertical.paddy.irrigation_efficiency";
        sel = inds,
        defaults = 1.0,
        type = Float,
    )

    paddy = Paddy{Float}(
        demand_gross = fill(mv, length(inds)),
        irrigation_efficiency = efficiency,
        h_min = h_min,
        h_max = h_max,
        h_p = h_p,
        h = fill(0.0, length(inds)),
    )

    return paddy
end

function initialize_nonpaddy(nc, config, inds)

    efficiency = ncread(
        nc,
        config,
        "vertical.nonpaddy.irrigation_efficiency";
        sel = inds,
        defaults = 1.0,
        type = Float,
    )

    nonpaddy = NonPaddy{Float}(
        demand_gross = fill(mv, length(inds)),
        irrigation_efficiency = efficiency,
    )

    return nonpaddy
end
