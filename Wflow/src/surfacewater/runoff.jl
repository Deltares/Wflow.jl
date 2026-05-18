abstract type AbstractRunoffModel end

"Struct for storing open water runoff variables"
@with_kw struct OpenWaterRunoffVariables
    n_cells::Int
    # Runoff from river based on riverfrac [mm Δt⁻¹]
    runoff_river::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Net runoff from river [mm Δt⁻¹]
    net_runoff_river::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Runoff from land based on waterfrac [mm Δt⁻¹]
    runoff_land::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Actual evaporation from open water (land) [mm Δt⁻¹]
    ae_openw_l::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Actual evaporation from river [mm Δt⁻¹]
    ae_openw_r::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct for storing open water runoff boundary conditions"
@with_kw struct OpenWaterRunoffBC
    n_cells::Int
    water_flux_surface::Vector{Float64} = fill(MISSING_VALUE, n_cells) # [mm dt-1]
    waterdepth_land::Vector{Float64} = fill(MISSING_VALUE, n_cells) # [mm]
    waterdepth_river::Vector{Float64} = zeros(n_cells) # [mm]
end

"Open water runoff model"
@with_kw struct OpenWaterRunoff <: AbstractRunoffModel
    n_cells::Int
    boundary_conditions::OpenWaterRunoffBC = OpenWaterRunoffBC(; n_cells)
    variables::OpenWaterRunoffVariables = OpenWaterRunoffVariables(; n_cells)
end

"Return the water flux at the surface (boundary condition) when snow is not modelled"
function get_water_flux_surface!(
    water_flux_surface::Vector{Float64},
    snow::NoSnowModel,
    glacier::AbstractGlacierModel,
    interception::AbstractInterceptionModel,
)
    (; throughfall, stemflow) = interception.variables
    @. water_flux_surface = throughfall + stemflow
    return nothing
end

"Return the water flux at the surface (boundary condition) when snow is modelled"
function get_water_flux_surface!(
    water_flux_surface::Vector{Float64},
    snow::AbstractSnowModel,
    glacier::AbstractGlacierModel,
    interception::AbstractInterceptionModel,
)
    water_flux_surface .=
        get_runoff(snow) .+ get_glacier_melt(glacier) .* get_glacier_fraction(glacier)
    return nothing
end

"Update boundary conditions of the open water runoff model for a single timestep"
function update_bc_open_water_runoff_model!(
    open_water_runoff_model::OpenWaterRunoff,
    external_models::NamedTuple,
    routing::Routing,
    network::NetworkRiver,
)
    (; water_flux_surface, waterdepth_river, waterdepth_land) =
        open_water_runoff_model.boundary_conditions
    (; land_cell_indices_containing_river) = network
    (; snow, glacier, interception) = external_models

    get_water_flux_surface!(water_flux_surface, snow, glacier, interception)

    # extract water depth h [m] from the land and river routing, used to limit open water
    # evaporation
    waterdepth_land .= routing.overland_flow.variables.h .* 1000.0
    for (river_cell_idx, land_cell_idx) in enumerate(land_cell_indices_containing_river)
        waterdepth_river[land_cell_idx] =
            routing.river_flow.variables.h[river_cell_idx] * 1000.0
    end
    return nothing
end

"Update the open water runoff model for a single timestep"
function update_open_water_runoff_model!(
    open_water_runoff_model::OpenWaterRunoff,
    atmospheric_forcing::AtmosphericForcing,
    parameters::LandParameters,
)
    (; boundary_conditions, variables) = open_water_runoff_model
    (; water_flux_surface, waterdepth_river, waterdepth_land) = boundary_conditions
    (; runoff_river, net_runoff_river, runoff_land, ae_openw_r, ae_openw_l) = variables
    (; potential_evaporation) = atmospheric_forcing
    (; river_fraction, water_fraction) = parameters

    @. runoff_river = min(1.0, river_fraction) * water_flux_surface
    @. runoff_land = min(1.0, water_fraction) * water_flux_surface
    @. ae_openw_r =
        min(waterdepth_river * river_fraction, river_fraction * potential_evaporation)
    @. ae_openw_l =
        min(waterdepth_land * water_fraction, water_fraction * potential_evaporation)
    @. net_runoff_river = runoff_river - ae_openw_r

    return nothing
end
