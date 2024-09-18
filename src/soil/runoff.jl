
@get_units @grid_loc @with_kw struct SurfaceRunoffModelVars{T}
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

@get_units @grid_loc @with_kw struct SurfaceRunoffModelParameters{T}
    # Fraction of river [-]
    riverfrac::Vector{T} | "-"
    # Fraction of open water (excluding rivers) [-]
    waterfrac::Vector{T} | "-"
end

@get_units @grid_loc @with_kw struct SurfaceRunoffBC{T}
    surface_water_flux::Vector{T}
    waterlevel_land::Vector{T} | "mm"
    waterlevel_river::Vector{T} | "mm"
end

function surface_runoff_model_bc(n)
    bc = SurfaceRunoffBC(;
        surface_water_flux = fill(mv, n),
        waterlevel_land = fill(mv, n),
        waterlevel_river = zeros(Float, n), # set to zero to account for cells outside river domain
    )
    return bc
end

function initialize_surface_runoff_model_params(nc, config, inds, riverfrac)
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
    params = SurfaceRunoffModelParameters(; waterfrac = waterfrac, riverfrac = riverfrac)
    return params
end

@get_units @with_kw struct SurfaceRunoff{T}
    boundary_conditions::SurfaceRunoffBC{T} | "-"
    parameters::SurfaceRunoffModelParameters{T} | "-"
    variables::SurfaceRunoffModelVars{T} | "-"
end

function surface_runoff_model_vars(n)
    vars = SurfaceRunoffModelVars(;
        runoff_river = fill(mv, n),
        runoff_land = fill(mv, n),
        ae_openw_l = fill(mv, n),
        ae_openw_r = fill(mv, n),
        net_runoff_river = fill(mv, n),
    )
    return vars
end

function initialize_surface_runoff_model(nc, config, inds, riverfrac)
    n = length(riverfrac)
    vars = surface_runoff_model_vars(n)
    bc = surface_runoff_model_bc(n)
    params = initialize_surface_runoff_model_params(nc, config, inds, riverfrac)
    model = SurfaceRunoff(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

get_surface_water_flux(snow::NoSnowModel, glacier, effective_precip) = effective_precip
function get_surface_water_flux(snow::AbstractSnowModel, glacier, effective_precip)
    return get_runoff(snow) .+ get_glacier_melt(glacier)
end

function update_boundary_conditions!(
    model::SurfaceRunoff,
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

function update!(model::SurfaceRunoff, atmospheric_forcing::AtmosphericForcing)
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