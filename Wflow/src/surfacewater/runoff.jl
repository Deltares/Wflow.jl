abstract type AbstractRunoffModel end

"Struct for storing open water runoff variables"
@with_kw struct OpenWaterRunoffVariables
    n::Int
    # Runoff from river based on riverfrac [m s⁻¹]
    runoff_river::Vector{Float64} = zeros(n)
    # Net runoff from river [m s⁻¹]
    net_runoff_river::Vector{Float64} = zeros(n)
    # Runoff from land based on waterfrac [m s⁻¹]
    runoff_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual evaporation from open water (land) [m s⁻¹]
    actual_open_water_evaporation_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual evaporation from river [m s⁻¹]
    actual_open_water_evaporation_river::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing open water runoff boundary conditions"
@with_kw struct OpenWaterRunoffBC
    n::Int
    # [m s⁻¹]
    water_flux_surface::Vector{Float64} = fill(MISSING_VALUE, n)
    # [m]
    waterdepth_land::Vector{Float64} = fill(MISSING_VALUE, n)
    # [m]
    waterdepth_river::Vector{Float64} = zeros(n)
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

function update_river_water_depth!(
        waterdepth_river::Vector{Float64},
        land_indices::Vector{Int},
        river_flow_model::AbstractRiverFlowModel,
    )
    # extract water depth h from the river (used to limit open water evaporation)
    for (i, land_index) in enumerate(land_indices)
        waterdepth_river[land_index] = river_flow_model.variables.h[i]
    end
    return nothing
end

function update_river_water_depth!(
        waterdepth_river::Vector{Float64},
        land_indices::Vector{Int},
        river_flow_model::NoRiverFlowModel,
    )
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
    (; land_indices) = network
    (; river_flow) = routing
    (; snow, glacier, interception) = external_models

    get_water_flux_surface!(water_flux_surface, snow, glacier, interception)
    update_river_water_depth!(waterdepth_river, land_indices, river_flow)

    # extract water depth h from the land (used to limit open water evaporation)
    waterdepth_land .= routing.overland_flow.variables.h

    return nothing
end

function update_river_runoff_model!(
        open_water_runoff_model::OpenWaterRunoff,
        river_flow_model::AbstractRiverFlowModel,
        potential_evaporation::Vector{Float64},
        river_fraction::Vector{Float64},
        dt::Float64,
    )

    (; boundary_conditions, variables) = open_water_runoff_model
    (; water_flux_surface, waterdepth_river) = boundary_conditions
    (;
        runoff_river,
        net_runoff_river,
        actual_open_water_evaporation_river,
    ) = variables

    @. runoff_river = min(1.0, river_fraction) * water_flux_surface
    @. actual_open_water_evaporation_river =
        river_fraction * min(waterdepth_river / dt, potential_evaporation)
    @. net_runoff_river = runoff_river - actual_open_water_evaporation_river
    return nothing
end

function update_river_runoff_model!(
        open_water_runoff_model::OpenWaterRunoff,
        river_flow_model::NoRiverFlowModel,
        potential_evaporation::Vector{Float64},
        river_fraction::Vector{Float64},
        dt::Float64,
    )
    return nothing
end

function update_land_runoff_model!(
        open_water_runoff_model::OpenWaterRunoff,
        potential_evaporation::Vector{Float64},
        water_fraction::Vector{Float64},
        dt::Float64,
    )
    (; boundary_conditions, variables) = open_water_runoff_model
    (; water_flux_surface, waterdepth_land) = boundary_conditions
    (;
        runoff_land,
        actual_open_water_evaporation_land,
    ) = variables

    @. runoff_land = min(1.0, water_fraction) * water_flux_surface
    @. actual_open_water_evaporation_land =
        water_fraction * min(waterdepth_land / dt, potential_evaporation)
    return nothing
end

"Update the open water runoff model for a single timestep"
function update_open_water_runoff_model!(
        open_water_runoff_model::OpenWaterRunoff,
        atmospheric_forcing::AtmosphericForcing,
        routing::Routing,
        parameters::LandParameters,
        dt::Float64,
    )
    (; potential_evaporation) = atmospheric_forcing
    (; river_fraction, water_fraction) = parameters
    (; river_flow) = routing

    update_river_runoff_model!(open_water_runoff_model, river_flow, potential_evaporation, river_fraction, dt)
    update_land_runoff_model!(open_water_runoff_model, potential_evaporation, water_fraction, dt)
    return nothing
end
