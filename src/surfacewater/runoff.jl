abstract type AbstractRunoffModel{T} end

@get_units @grid_loc @with_kw struct OpenWaterRunoffVariables{T}
    # Runoff from river based on riverfrac [mm Δt⁻¹]
    runoff_river::Vector{T}
    # Net runoff from river [mm Δt⁻¹]
    net_runoff_river::Vector{T}
    # Runoff from land based on waterfrac [mm Δt⁻¹]
    runoff_land::Vector{T}
    # Actual evaporation from open water (land) [mm Δt⁻¹]
    ae_openw_l::Vector{T}
    # Actual evaporation from river [mm Δt⁻¹]
    ae_openw_r::Vector{T}
end

function OpenWaterRunoffVariables(
    n;
    runoff_river::Vector{T} = fill(mv, n),
    runoff_land::Vector{T} = fill(mv, n),
    ae_openw_l::Vector{T} = fill(mv, n),
    ae_openw_r::Vector{T} = fill(mv, n),
    net_runoff_river::Vector{T} = fill(mv, n),
) where {T}
    return OpenWaterRunoffVariables{T}(;
        runoff_river = runoff_river,
        runoff_land = runoff_land,
        ae_openw_l = ae_openw_l,
        ae_openw_r = ae_openw_r,
        net_runoff_river = net_runoff_river,
    )
end

@get_units @grid_loc @with_kw struct OpenWaterRunoffParameters{T}
    # Fraction of river [-]
    riverfrac::Vector{T} | "-"
    # Fraction of open water (excluding rivers) [-]
    waterfrac::Vector{T} | "-"
end

@get_units @grid_loc @with_kw struct OpenWaterRunoffBC{T}
    surface_water_flux::Vector{T}
    waterlevel_land::Vector{T} | "mm"
    waterlevel_river::Vector{T} | "mm"
end

function OpenWaterRunoffBC(
    n;
    surface_water_flux::Vector{T} = fill(mv, n),
    waterlevel_land::Vector{T} = fill(mv, n),
    waterlevel_river::Vector{T} = zeros(Float, n), # set to zero to account for cells outside river domain
) where {T}
    return OpenWaterRunoffBC{T}(;
        surface_water_flux = surface_water_flux,
        waterlevel_land = waterlevel_land,
        waterlevel_river = waterlevel_river,
    )
end

function OpenWaterRunoffParameters(nc, config, inds, riverfrac)
    # fraction open water
    waterfrac = ncread(
        nc,
        config,
        "vertical.runoff.parameters.waterfrac";
        sel = inds,
        defaults = 0.0,
        type = Float,
    )
    waterfrac = max.(waterfrac .- riverfrac, Float(0.0))
    params = OpenWaterRunoffParameters(; waterfrac = waterfrac, riverfrac = riverfrac)
    return params
end

@get_units @with_kw struct OpenWaterRunoff{T} <: AbstractRunoffModel{T}
    boundary_conditions::OpenWaterRunoffBC{T} | "-"
    parameters::OpenWaterRunoffParameters{T} | "-"
    variables::OpenWaterRunoffVariables{T} | "-"
end

function OpenWaterRunoff(nc, config, inds, riverfrac)
    n = length(riverfrac)
    vars = OpenWaterRunoffVariables(n)
    bc = OpenWaterRunoffBC(n)
    params = OpenWaterRunoffParameters(nc, config, inds, riverfrac)
    model =
        OpenWaterRunoff(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

get_surface_water_flux(snow::NoSnowModel, glacier, effective_precip) = effective_precip
function get_surface_water_flux(snow::AbstractSnowModel, glacier, effective_precip)
    return get_runoff(snow) .+ get_glacier_melt(glacier)
end

function update_boundary_conditions!(
    model::OpenWaterRunoff,
    external_models::NamedTuple,
    lateral,
    network,
)
    (; surface_water_flux, waterlevel_river, waterlevel_land) = model.boundary_conditions
    inds_riv = network.index_river
    (; snow, glacier, interception) = external_models
    (; throughfall, stemflow) = interception.variables

    surface_water_flux .= get_surface_water_flux(snow, glacier, throughfall .+ stemflow)

    # extract water levels h_av [m] from the land and river domains this is used to limit
    # open water evaporation
    waterlevel_land .= lateral.land.h_av .* 1000.0
    waterlevel_river[inds_riv] .= lateral.river.h_av .* 1000.0
end

function update!(model::OpenWaterRunoff, atmospheric_forcing::AtmosphericForcing)
    (; potential_evaporation) = atmospheric_forcing
    (; runoff_river, net_runoff_river, runoff_land, ae_openw_r, ae_openw_l) =
        model.variables
    (; riverfrac, waterfrac) = model.parameters
    (; surface_water_flux, waterlevel_river, waterlevel_land) = model.boundary_conditions

    @. runoff_river = min(1.0, riverfrac) * surface_water_flux
    @. runoff_land = min(1.0, waterfrac) * surface_water_flux
    @. ae_openw_r = min(waterlevel_river * riverfrac, riverfrac * potential_evaporation)
    @. ae_openw_l = min(waterlevel_land * waterfrac, waterfrac * potential_evaporation)
    @. net_runoff_river = runoff_river - ae_openw_r
end