@get_units @with_kw struct SoilLoss{RE, OLE, SE, T}
    hydrometeo_forcing::HydrometeoForcing | "-"
    geometry::LandGeometry | "-"
    rainfall_erosion::RE | "-"
    overland_flow_erosion::OLE | "-"
    soil_erosion::SE | "-"
end

function initialize_soil_loss(nc, config, inds)
    n = length(inds)
    hydrometeo_forcing = initialize_hydrometeo_forcing(n)
    geometry = initialize_land_geometry(nc, config, inds)
    # Rainfall erosion
    rainfallerosionmodel = get(config.model, "rainfall_erosion", "answers")::String
    if rainfallerosionmodel == "answers"
        rainfall_erosion_model = initialize_answers_rainfall_erosion_model(nc, config, inds)
    elseif rainfallerosionmodel == "eurosem"
        rainfall_erosion_model = initialize_eurosem_rainfall_erosion_model(nc, config, inds)
    else
        error("Unknown rainfall erosion model: $rainfallerosionmodel")
    end

    # Overland flow erosion
    overlandflowerosionmodel = get(config.model, "overland_flow_erosion", "answers")::String
    if overlandflowerosionmodel == "answers"
        overland_flow_erosion_model =
            initialize_answers_overland_flow_erosion_model(nc, config, inds)
    else
        error("Unknown overland flow erosion model: $overlandflowerosionmodel")
        # overland_flow_erosion_model = NoOverlandFlowErosionModel()
    end

    # Total soil erosion and particle differentiation
    soil_erosion_model = initialize_soil_erosion_model(nc, config, inds)

    soil_loss = SoilLoss{
        typeof(rainfall_erosion_model),
        typeof(overland_flow_erosion_model),
        typeof(soil_erosion_model),
        Float,
    }(;
        hydrometeo_forcing = hydrometeo_forcing,
        geometry = geometry,
        rainfall_erosion = rainfall_erosion_model,
        overland_flow_erosion = overland_flow_erosion_model,
        soil_erosion = soil_erosion_model,
    )
    return soil_loss
end

function update!(model::SoilLoss, dt)
    # Convert dt to integer
    ts = tosecond(dt)
    #TODO add interception/canopygapfraction calculation here for eurosem
    #need SBM refactor
    # Rainfall erosion
    update!(model.rainfall_erosion, model.hydrometeo_forcing, model.area, ts)
    # Overland flow erosion
    update!(model.overland_flow_erosion, model.hydrometeo_forcing, model.geometry, ts)
    # Total soil erosion and particle differentiation
    re = get_rainfall_erosion(model.rainfall_erosion)
    ole = get_overland_flow_erosion(model.overland_flow_erosion)
    (; rainfall_erosion, overland_flow_erosion) = model.soil_erosion.boundary_conditions
    @. rainfall_erosion = re
    @. overland_flow_erosion = ole
    update!(model.soil_erosion)
end
