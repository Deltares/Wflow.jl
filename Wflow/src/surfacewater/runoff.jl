abstract type AbstractRunoffModel end

"Struct for storing open water runoff variables"
@with_kw struct OpenWaterRunoffVariables
    # Runoff from river based on riverfrac [mm Δt⁻¹]
    runoff_river::Vector{Float}
    # Net runoff from river [mm Δt⁻¹]
    net_runoff_river::Vector{Float}
    # Runoff from land based on waterfrac [mm Δt⁻¹]
    runoff_land::Vector{Float}
    # Actual evaporation from open water (land) [mm Δt⁻¹]
    ae_openw_l::Vector{Float}
    # Actual evaporation from river [mm Δt⁻¹]
    ae_openw_r::Vector{Float}
end

"Initialize open water runoff model variables"
function OpenWaterRunoffVariables(n::Int)
    return OpenWaterRunoffVariables(;
        runoff_river = fill(MISSING_VALUE, n),
        runoff_land = fill(MISSING_VALUE, n),
        ae_openw_l = fill(MISSING_VALUE, n),
        ae_openw_r = fill(MISSING_VALUE, n),
        net_runoff_river = fill(MISSING_VALUE, n),
    )
end

"Struct for storing open water runoff boundary conditions"
@with_kw struct OpenWaterRunoffBC
    water_flux_surface::Vector{Float} # [mm dt-1]
    waterdepth_land::Vector{Float} # [mm]
    waterdepth_river::Vector{Float} # [mm]
end

"Initialize open water runoff boundary conditions"
function OpenWaterRunoffBC(n::Int)
    return OpenWaterRunoffBC(;
        water_flux_surface = fill(MISSING_VALUE, n),
        waterdepth_land = fill(MISSING_VALUE, n),
        waterdepth_river = zeros(Float, n),
    )
end

"Open water runoff model"
@with_kw struct OpenWaterRunoff <: AbstractRunoffModel
    boundary_conditions::OpenWaterRunoffBC
    variables::OpenWaterRunoffVariables
end

"Initialize open water runoff model"
function OpenWaterRunoff(n::Int)
    variables = OpenWaterRunoffVariables(n)
    boundary_conditions = OpenWaterRunoffBC(n)
    model = OpenWaterRunoff(; boundary_conditions, variables)
    return model
end

"Return the water flux at the surface (boundary condition) when snow is not modelled"
function get_water_flux_surface!(
    water_flux_surface::Vector{Float},
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
    water_flux_surface::Vector{Float},
    snow::AbstractSnowModel,
    glacier::AbstractGlacierModel,
    interception::AbstractInterceptionModel,
)
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

    # extract water depth h [m] from the land and river routing, used to limit open water
    # evaporation
    waterdepth_land .= routing.overland_flow.variables.h .* 1000.0
    for (i, land_index) in enumerate(land_indices)
        waterdepth_river[land_index] = routing.river_flow.variables.h[i] * 1000.0
    end
    return nothing
end

"Update the open water runoff model for a single timestep"
function update!(
    model::OpenWaterRunoff,
    atmospheric_forcing::AtmosphericForcing,
    parameters::LandParameters,
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