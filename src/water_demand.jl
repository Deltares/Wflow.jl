@get_units @with_kw struct WaterDemand{T}
    domestic_gross::Vector{T}   # gross domestic water demand [mm Δt⁻¹]
    domestic_net::Vector{T}     # net domestic water demand [mm Δt⁻¹]
    industry_gross::Vector{T}   # gross industry water demand [mm Δt⁻¹]
    industry_net::Vector{T}     # net industry water demand [mm Δt⁻¹]
    livestock_gross::Vector{T}  # gross livestock water demand [mm Δt⁻¹]
    livestock_net::Vector{T}    # net livestock water demand [mm Δt⁻¹]
end

function initialize_water_demand(nc, config, inds, Δt)

    domestic_gross = ncread(
        nc,
        config.input,
        "vertical.water_demand.domestic_gross";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)
    domestic_net = ncread(
        nc,
        config.input,
        "vertical.water_demand.domestic_net";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)
    industry_gross = ncread(
        nc,
        config.input,
        "vertical.water_demand.industry_gross";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)
    industry_net = ncread(
        nc,
        config.input,
        "vertical.water_demand.industry_net";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)
    livestock_gross = ncread(
        nc,
        config.input,
        "vertical.water_demand.livestock_gross";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)
    livestock_net = ncread(
        nc,
        config.input,
        "vertical.water_demand.livestock_net";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)

    water_demand = WaterDemand{Float}(
        domestic_gross = domestic_gross,
        domestic_net = domestic_net,
        industry_gross = industry_gross,
        industry_net = industry_net,
        livestock_gross = livestock_gross,
        livestock_net = livestock_net,
    )

    return water_demand

end