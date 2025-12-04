abstract type AbstractRunoffModel end

"Struct for storing open water runoff variables"
@with_kw struct OpenWaterRunoffVariables
    n::Int
    # Runoff from river based on riverfrac [mm dt⁻¹ => m s⁻¹]
    runoff_river::Vector{Float64} = fill(MISSING_VALUE, n)
    # Net runoff from river [mm dt⁻¹ => m s⁻¹]
    net_runoff_river::Vector{Float64} = fill(MISSING_VALUE, n)
    # Runoff from land based on waterfrac [mm dt⁻¹ => m s⁻¹]
    runoff_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual evaporation from open water (land) [mm dt⁻¹ => m s⁻¹]
    ae_openw_l::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual evaporation from river [mm dt⁻¹ => m s⁻¹]
    ae_openw_r::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing open water runoff boundary conditions"
@with_kw struct OpenWaterRunoffBC
    n::Int
    # [mm dt-1 => m s⁻¹]
    water_flux_surface::Vector{Float64} = fill(MISSING_VALUE, n)
    # [mm => m]
    waterdepth_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # [mm => m]
    waterdepth_river::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Open water runoff model"
@with_kw struct OpenWaterRunoff <: AbstractRunoffModel
    n::Int
    boundary_conditions::OpenWaterRunoffBC = OpenWaterRunoffBC(; n)
    variables::OpenWaterRunoffVariables = OpenWaterRunoffVariables(; n)
end

"Return the water flux at the surface (boundary condition) when snow is not modelled"
function get_water_flux_surface!(
    water_flux_surface::Vector{Float64},
    snow::NoSnowModel,
    glacier::AbstractGlacierModel,
    interception::AbstractInterceptionModel,
)
    (; throughfall, stemflow) = interception.variables
    # [m s⁻¹] = [m s⁻¹] + [m s⁻¹]
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
    # [m s⁻¹] = [m s⁻¹] + [m s⁻¹] * [-]
    water_flux_surface .=
        get_runoff(snow) .+ get_glacier_melt(glacier) .* get_glacier_fraction(glacier)
    return nothing
end

"Update boundary conditions of the open water runoff model for a single timestep"
function update_boundary_conditions!(
    model::OpenWaterRunoff,
    external_models::NamedTuple,
    routing::Routing,
    network::NetworkRiver,
)
    (; water_flux_surface, waterdepth_river, waterdepth_land) = model.boundary_conditions
    (; land_indices) = network
    (; snow, glacier, interception) = external_models

    get_water_flux_surface!(water_flux_surface, snow, glacier, interception)

    # extract water depth h from the land and river routing, used to limit open water
    # evaporation
    waterdepth_land .= routing.overland_flow.variables.h
    for (i, land_index) in enumerate(land_indices)
        # [m] = [m]
        waterdepth_river[land_index] = routing.river_flow.variables.h[i]
    end
    return nothing
end

"Update the open water runoff model for a single timestep"
function update!(
    model::OpenWaterRunoff,
    atmospheric_forcing::AtmosphericForcing,
    parameters::LandParameters,
    dt::Number,
)
    (; potential_evaporation) = atmospheric_forcing
    (; runoff_river, net_runoff_river, runoff_land, ae_openw_r, ae_openw_l) =
        model.variables
    (; river_fraction, water_fraction) = parameters
    (; water_flux_surface, waterdepth_river, waterdepth_land) = model.boundary_conditions

    @. runoff_river = min(1.0, river_fraction) * water_flux_surface
    @. runoff_land = min(1.0, water_fraction) * water_flux_surface
    @. ae_openw_r =
        min(waterdepth_river * river_fraction, river_fraction * potential_evaporation)
    @. ae_openw_l =
        min(waterdepth_land * water_fraction, water_fraction * potential_evaporation)
    @. net_runoff_river = runoff_river - ae_openw_r

    return nothing
end
