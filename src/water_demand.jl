@get_units @with_kw struct WaterDemand{T}
    domestic_gross::Vector{T}   # gross domestic water demand [mm Δt⁻¹]
    domestic_netto::Vector{T}   # net domestic water demand [mm Δt⁻¹]
    industry_gross::Vector{T}   # gross industry water demand [mm Δt⁻¹]
    industry_netto::Vector{T}   # net industry water demand [mm Δt⁻¹]
    livestock_gross::Vector{T}  # gross livestock water demand [mm Δt⁻¹]
    livestock_netto::Vector{T}  # net livestock water demand [mm Δt⁻¹]
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
    domestic_netto = ncread(
        nc,
        config.input,
        "vertical.water_demand.domestic_netto";
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
    industry_netto = ncread(
        nc,
        config.input,
        "vertical.water_demand.industry_netto";
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
    livestock_netto = ncread(
        nc,
        config.input,
        "vertical.water_demand.livestock_netto";
        sel = inds,
        defaults = 0.0,
        type = Float,
    ).* (Δt / basetimestep)

    water_demand = WaterDemand{Float}(
        domestic_gross = domestic_gross,
        domestic_netto = domestic_netto,
        industry_gross = industry_gross,
        industry_netto = industry_netto,
        livestock_gross = livestock_gross,
        livestock_netto = livestock_netto,
    )

    return water_demand

end