"Total soil erosion model"
@with_kw struct SoilLoss{RE, OFE, SE, T}
    atmospheric_forcing::AtmosphericForcing{T}
    hydrological_forcing::HydrologicalForcing{T}
    geometry::LandParameters{T}
    rainfall_erosion::RE
    overland_flow_erosion::OFE
    soil_erosion::SE
end

function SoilLoss(nc, config, inds)
    n = length(inds)

    atmospheric_forcing = AtmosphericForcing(n)
    hydrological_forcing = HydrologicalForcing(n)
    geometry = LandParameters(nc, config, inds)

    # Rainfall erosion
    rainfallerosionmodel = get(config.model, "rainfall_erosion", "answers")::String
    if rainfallerosionmodel == "answers"
        rainfall_erosion_model = RainfallErosionAnswersModel(nc, config, inds)
    elseif rainfallerosionmodel == "eurosem"
        rainfall_erosion_model = RainfallErosionEurosemModel(nc, config, inds)
    else
        error("Unknown rainfall erosion model: $rainfallerosionmodel")
    end

    # Overland flow erosion
    overlandflowerosionmodel = get(config.model, "overland_flow_erosion", "answers")::String
    if overlandflowerosionmodel == "answers"
        overland_flow_erosion_model = OverlandFlowErosionAnswersModel(nc, config, inds)
    else
        error("Unknown overland flow erosion model: $overlandflowerosionmodel")
        # overland_flow_erosion_model = NoOverlandFlowErosionModel()
    end

    # Total soil erosion and particle differentiation
    soil_erosion_model = SoilErosionModel(nc, config, inds)

    soil_loss = SoilLoss{
        typeof(rainfall_erosion_model),
        typeof(overland_flow_erosion_model),
        typeof(soil_erosion_model),
        Float,
    }(;
        atmospheric_forcing = atmospheric_forcing,
        hydrological_forcing = hydrological_forcing,
        geometry = geometry,
        rainfall_erosion = rainfall_erosion_model,
        overland_flow_erosion = overland_flow_erosion_model,
        soil_erosion = soil_erosion_model,
    )
    return soil_loss
end

function update!(model::SoilLoss, dt)
    (;
        atmospheric_forcing,
        hydrological_forcing,
        geometry,
        rainfall_erosion,
        overland_flow_erosion,
        soil_erosion,
    ) = model
    # Convert dt to integer
    ts = tosecond(dt)
    #TODO add interception/canopygapfraction calculation here for eurosem
    #need SBM refactor

    # Rainfall erosion
    update_boundary_conditions!(rainfall_erosion, atmospheric_forcing, hydrological_forcing)
    update!(rainfall_erosion, geometry, ts)
    # Overland flow erosion
    update_boundary_conditions!(overland_flow_erosion, hydrological_forcing)
    update!(overland_flow_erosion, geometry, ts)
    # Total soil erosion and particle differentiation
    update_boundary_conditions!(soil_erosion, rainfall_erosion, overland_flow_erosion)
    update!(soil_erosion)
end
