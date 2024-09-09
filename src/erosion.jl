@get_units @with_kw struct SoilErosion{RE, OLE, TE}
    rainfall_erosion::RE | "-"
    overland_flow_erosion::OLE | "-"
    total_erosion::TE | "-"
end

function initialize_soil_erosion(nc, config, inds)

    # Rainfall erosion
    rainfallerosionmodel = get(config.model, "rainfall_erosion", "answers")::String
    if rainfallerosionmodel == "answers"
        rainfall_erosion_model = initialize_rainfall_erosion_answers(nc, config, inds)
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
    total_erosion_model = initialize_soil_erosion_model(nc, config, inds)

    soil_erosion = SoilErosion{
        typeof(rainfall_erosion_model),
        typeof(overland_flow_erosion_model),
        typeof(total_erosion_model),
    }(;
        rainfall_erosion = rainfall_erosion_model,
        overland_flow_erosion = overland_flow_erosion_model,
        total_erosion = total_erosion_model,
    )
    return soil_erosion
end

function update!(model::SoilErosion, area, dt)
    #TODO add interception/canopygapfraction calculation here for eurosem
    #need SBM refactor
    # Rainfall erosion
    update!(model.rainfall_erosion, area, dt)
    # Overland flow erosion
    update!(model.overland_flow_erosion, area, dt)
    # Total soil erosion and particle differentiation
    (; rainfall_erosion, overland_flow_erosion) = model.total_erosion.boundary_conditions
    @. rainfall_erosion = get_rainfall_erosion(model.rainfall_erosion)
    @. overland_flow_erosion = get_overland_flow_erosion(model.overland_flow_erosion)
    update!(model.total_erosion)
end
