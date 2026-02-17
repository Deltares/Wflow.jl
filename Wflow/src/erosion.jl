"Soil loss model"
@with_kw struct SoilLoss{RE <: AbstractRainfallErosionModel} <: AbstractLandModel
    n::Int
    atmospheric_forcing::AtmosphericForcing = AtmosphericForcing(; n)
    hydrological_forcing::HydrologicalForcing = HydrologicalForcing(; n)
    rainfall_erosion::RE
    overland_flow_erosion::OverlandFlowErosionAnswersModel
    soil_erosion::SoilErosionModel
end

"Initialize soil loss model"
function SoilLoss(dataset::NCDataset, config::Config, indices::Vector{CartesianIndex{2}})
    (; rainfall_erosion, overland_flow_erosion) = config.model
    n = length(indices)

    # Rainfall erosion
    if rainfall_erosion == RainfallErosionType.answers
        rainfall_erosion = RainfallErosionAnswersModel(dataset, config, indices)
    elseif rainfall_erosion == RainfallErosionType.eurosem
        rainfall_erosion = RainfallErosionEurosemModel(dataset, config, indices)
    end

    # Overland flow erosion
    overland_flow_erosion = OverlandFlowErosionAnswersModel(dataset, config, indices)

    # Total soil erosion and particle differentiation
    soil_erosion = SoilErosionModel(dataset, config, indices)

    soil_loss = SoilLoss(; n, rainfall_erosion, overland_flow_erosion, soil_erosion)
    return soil_loss
end

"Update soil loss model for a single timestep"
function update!(model::SoilLoss, parameters::LandParameters, dt::Float64)
    (;
        atmospheric_forcing,
        hydrological_forcing,
        rainfall_erosion,
        overland_flow_erosion,
        soil_erosion,
    ) = model

    #TODO add interception/canopygapfraction calculation here for eurosem
    #need SBM refactor

    # Rainfall erosion
    update_boundary_conditions!(rainfall_erosion, atmospheric_forcing, hydrological_forcing)
    update!(rainfall_erosion, parameters, dt)
    # Overland flow erosion
    update_boundary_conditions!(overland_flow_erosion, hydrological_forcing)
    update!(overland_flow_erosion, parameters, dt)
    # Total soil erosion and particle differentiation
    update_boundary_conditions!(soil_erosion, rainfall_erosion, overland_flow_erosion)
    update!(soil_erosion)
end
