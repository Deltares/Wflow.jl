abstract type AbstractRunoffModel{T} end

"Struct for storing open water runoff variables"
@with_kw struct OpenWaterRunoffVariables{T}
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

"Initialize open water runoff model variables"
function OpenWaterRunoffVariables(T::Type{<:AbstractFloat}, n::Int)
    return OpenWaterRunoffVariables{T}(;
        runoff_river = fill(MISSING_VALUE, n),
        runoff_land = fill(MISSING_VALUE, n),
        ae_openw_l = fill(MISSING_VALUE, n),
        ae_openw_r = fill(MISSING_VALUE, n),
        net_runoff_river = fill(MISSING_VALUE, n),
    )
end

"Struct for storing open water runoff parameters"
@with_kw struct OpenWaterRunoffParameters{T}
    # Fraction of river [-]
    riverfrac::Vector{T}
    # Fraction of open water (excluding rivers) [-]
    waterfrac::Vector{T}
end

"Initialize open water runoff parameters"
function OpenWaterRunoffParameters(dataset, config, indices, riverfrac)
    # fraction open water
    lens = lens_input_parameter(config, "land~water-covered__area_fraction")
    waterfrac = ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float64)
    waterfrac = max.(waterfrac .- riverfrac, 0.0)
    params = OpenWaterRunoffParameters(; waterfrac = waterfrac, riverfrac = riverfrac)
    return params
end

"Struct for storing open water runoff boundary conditions"
@with_kw struct OpenWaterRunoffBC{T}
    water_flux_surface::Vector{T} # [mm dt-1]
    waterdepth_land::Vector{T} # [mm]
    waterdepth_river::Vector{T} # [mm]
end

"Initialize open water runoff boundary conditions"
function OpenWaterRunoffBC(T::Type{<:AbstractFloat}, n::Int)
    return OpenWaterRunoffBC{T}(;
        water_flux_surface = fill(MISSING_VALUE, n),
        waterdepth_land = fill(MISSING_VALUE, n),
        waterdepth_river = zeros(T, n),
    )
end

"Open water runoff model"
@with_kw struct OpenWaterRunoff{T} <: AbstractRunoffModel{T}
    boundary_conditions::OpenWaterRunoffBC{T}
    parameters::OpenWaterRunoffParameters{T}
    variables::OpenWaterRunoffVariables{T}
end

"Initialize open water runoff model"
function OpenWaterRunoff(dataset, config, indices, riverfrac)
    n = length(riverfrac)
    vars = OpenWaterRunoffVariables(Float64, n)
    bc = OpenWaterRunoffBC(Float64, n)
    params = OpenWaterRunoffParameters(dataset, config, indices, riverfrac)
    model =
        OpenWaterRunoff(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

"Return the water flux at the surface (boundary condition) when snow is not modelled"
function get_water_flux_surface!(
    water_flux_surface,
    snow::NoSnowModel,
    glacier,
    interception,
)
    (; throughfall, stemflow) = interception.variables
    @. water_flux_surface = throughfall + stemflow
    return nothing
end

"Return the water flux at the surface (boundary condition) when snow is modelled"
function get_water_flux_surface!(
    water_flux_surface,
    snow::AbstractSnowModel,
    glacier,
    interception,
)
    water_flux_surface .=
        get_runoff(snow) .+ get_glacier_melt(glacier) .* get_glacier_fraction(glacier)
    return nothing
end

"Update boundary conditions of the open water runoff model for a single timestep"
function update_boundary_conditions!(
    model::OpenWaterRunoff,
    external_models::NamedTuple,
    routing,
    network,
)
    (; water_flux_surface, waterdepth_river, waterdepth_land) = model.boundary_conditions
    (; land_indices) = network.river
    (; snow, glacier, interception) = external_models

    get_water_flux_surface!(water_flux_surface, snow, glacier, interception)

    # extract water depth h [m] from the land and river routing, used to limit open water
    # evaporation
    waterdepth_land .= routing.overland_flow.variables.h .* 1000.0
    for (i, land_index) in enumerate(land_indices)
        waterdepth_river[land_index] = routing.river_flow.variables.h[i] * 1000.0
    end
    return nothing
end

"Update the open water runoff model for a single timestep"
function update!(model::OpenWaterRunoff, atmospheric_forcing::AtmosphericForcing)
    (; potential_evaporation) = atmospheric_forcing
    (; runoff_river, net_runoff_river, runoff_land, ae_openw_r, ae_openw_l) =
        model.variables
    (; riverfrac, waterfrac) = model.parameters
    (; water_flux_surface, waterdepth_river, waterdepth_land) = model.boundary_conditions

    @. runoff_river = min(1.0, riverfrac) * water_flux_surface
    @. runoff_land = min(1.0, waterfrac) * water_flux_surface
    @. ae_openw_r = min(waterdepth_river * riverfrac, riverfrac * potential_evaporation)
    @. ae_openw_l = min(waterdepth_land * waterfrac, waterfrac * potential_evaporation)
    @. net_runoff_river = runoff_river - ae_openw_r

    return nothing
end