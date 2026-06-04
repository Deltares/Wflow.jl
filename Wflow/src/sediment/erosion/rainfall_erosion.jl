abstract type AbstractRainfallErosionModel end

"Struct for storing rainfall erosion model variables"
@with_kw struct RainfallErosionModelVariables
    n_cells::Int
    # Total soil erosion rate [kg s⁻¹] from rainfall (splash)
    soil_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct for storing EUROSEM rainfall erosion model boundary conditions"
@with_kw struct RainfallErosionEurosemBC
    n_cells::Int
    # precipitation [m s⁻¹]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Interception [m s⁻¹]
    interception::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Waterlevel on land [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct for storing EUROSEM rainfall erosion model parameters"
@with_kw struct RainfallErosionEurosemParameters
    # Soil detachability factor [kg (kg m² s⁻²)⁻¹]
    soil_detachability::Vector{Float64}
    # Exponent EUROSEM [m⁻¹]
    eurosem_exponent::Vector{Float64}
    # Canopy height [m]
    canopy_height::Vector{Float64}
    # Canopy gap fraction [-]
    canopy_gap_fraction::Vector{Float64}
    # Fraction of the soil that is covered (eg paved, snow, etc) [-]
    soilcover_fraction::Vector{Float64}
end

"Initialize EUROSEM rainfall erosion model parameters"
function RainfallErosionEurosemParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    soil_detachability = ncread(
        dataset,
        config,
        "soil_erosion__rainfall_soil_detachability_factor",
        SoilLossModel;
        sel = indices,
    )
    eurosem_exponent = ncread(
        dataset,
        config,
        "soil_erosion__eurosem_exponent",
        SoilLossModel;
        sel = indices,
    )
    canopy_height =
        ncread(dataset, config, "vegetation_canopy__height", SoilLossModel; sel = indices)
    canopy_gap_fraction = ncread(
        dataset,
        config,
        "vegetation_canopy__gap_fraction",
        SoilLossModel;
        sel = indices,
    )
    soilcover_fraction = ncread(
        dataset,
        config,
        "compacted_soil__area_fraction",
        SoilLossModel;
        sel = indices,
    )

    eurosem_parameters = RainfallErosionEurosemParameters(;
        soil_detachability,
        eurosem_exponent,
        canopy_height,
        canopy_gap_fraction,
        soilcover_fraction,
    )
    return eurosem_parameters
end

"EUROSEM rainfall erosion model"
@with_kw struct RainfallErosionEurosemModel <: AbstractRainfallErosionModel
    n_cells::Int
    boundary_conditions::RainfallErosionEurosemBC = RainfallErosionEurosemBC(; n_cells)
    parameters::RainfallErosionEurosemParameters
    variables::RainfallErosionModelVariables = RainfallErosionModelVariables(; n_cells)
end

"Initialize EUROSEM rainfall erosion model"
function RainfallErosionEurosemModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n_cells = length(indices)
    parameters = RainfallErosionEurosemParameters(dataset, config, indices)
    rainfall_erosion_model = RainfallErosionEurosemModel(; n_cells, parameters)
    return rainfall_erosion_model
end

"Update EUROSEM rainfall erosion model boundary conditions for a single timestep"
function update_bc_rainfall_erosion_model!(
    rainfall_erosion_model::RainfallErosionEurosemModel,
    atmospheric_forcing::AtmosphericForcing,
    hydrological_forcing::HydrologicalForcing,
)
    (; precipitation, interception, waterlevel) = rainfall_erosion_model.boundary_conditions
    @. precipitation = atmospheric_forcing.precipitation
    @. waterlevel = hydrological_forcing.waterlevel_land
    @. interception = hydrological_forcing.interception
end

"Update EUROSEM rainfall erosion model for a single timestep"
function update_rainfall_erosion_model!(
    rainfall_erosion_model::RainfallErosionEurosemModel,
    parameters::LandParameters,
    dt::Float64,
)
    (; precipitation, interception, waterlevel) = rainfall_erosion_model.boundary_conditions
    (;
        soil_detachability,
        eurosem_exponent,
        canopy_height,
        canopy_gap_fraction,
        soilcover_fraction,
    ) = rainfall_erosion_model.parameters
    (; soil_erosion_rate) = rainfall_erosion_model.variables

    n_cells = length(precipitation)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        soil_erosion_rate[cell_idx] = rainfall_erosion_eurosem(
            precipitation[cell_idx],
            interception[cell_idx],
            waterlevel[cell_idx],
            soil_detachability[cell_idx],
            eurosem_exponent[cell_idx],
            canopy_height[cell_idx],
            canopy_gap_fraction[cell_idx],
            soilcover_fraction[cell_idx],
            parameters.area[cell_idx],
            dt,
        )
    end
end

"Struct for storing ANSWERS rainfall erosion model boundary conditions"
@with_kw struct RainfallErosionAnswersBC
    n_cells::Int
    # precipitation [m s⁻¹]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct for storing ANSWERS rainfall erosion model parameters"
@with_kw struct RainfallErosionAnswersParameters
    # Soil erodibility factor [-]
    usle_k::Vector{Float64}
    # Crop management factor [-]
    usle_c::Vector{Float64}
    # ANSWERS rainfall erosion factor [-]
    answers_rainfall_factor::Vector{Float64}
end

"Initialize ANSWERS rainfall erosion model parameters"
function RainfallErosionAnswersParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    usle_k =
        ncread(dataset, config, "soil_erosion__usle_k_factor", SoilLossModel; sel = indices)
    usle_c =
        ncread(dataset, config, "soil_erosion__usle_c_factor", SoilLossModel; sel = indices)
    answers_rainfall_factor = ncread(
        dataset,
        config,
        "soil_erosion__answers_rainfall_factor",
        SoilLossModel;
        sel = indices,
    )

    answers_parameters =
        RainfallErosionAnswersParameters(; usle_k, usle_c, answers_rainfall_factor)
    return answers_parameters
end

"ANSWERS rainfall erosion model"
@with_kw struct RainfallErosionAnswersModel <: AbstractRainfallErosionModel
    n_cells::Int
    boundary_conditions::RainfallErosionAnswersBC = RainfallErosionAnswersBC(; n_cells)
    parameters::RainfallErosionAnswersParameters
    variables::RainfallErosionModelVariables = RainfallErosionModelVariables(; n_cells)
end

"Initialize ANSWERS rainfall erosion model"
function RainfallErosionAnswersModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n_cells = length(indices)
    parameters = RainfallErosionAnswersParameters(dataset, config, indices)
    rainfall_erosion_model = RainfallErosionAnswersModel(; n_cells, parameters)
    return rainfall_erosion_model
end

"Update ANSWERS rainfall erosion model boundary conditions for a single timestep"
function update_bc_rainfall_erosion_model!(
    rainfall_erosion_model::RainfallErosionAnswersModel,
    atmospheric_forcing::AtmosphericForcing,
    ::HydrologicalForcing,
)
    (; precipitation) = rainfall_erosion_model.boundary_conditions
    @. precipitation = atmospheric_forcing.precipitation
end

"Update ANSWERS rainfall erosion model for a single timestep"
function update_rainfall_erosion_model!(
    rainfall_erosion_model::RainfallErosionAnswersModel,
    parameters::LandParameters,
    dt::Float64,
)
    (; precipitation) = rainfall_erosion_model.boundary_conditions
    (; usle_k, usle_c, answers_rainfall_factor) = rainfall_erosion_model.parameters
    (; soil_erosion_rate) = rainfall_erosion_model.variables

    n_cells = length(precipitation)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        soil_erosion_rate[cell_idx] = rainfall_erosion_answers(
            precipitation[cell_idx],
            usle_k[cell_idx],
            usle_c[cell_idx],
            answers_rainfall_factor[cell_idx],
            parameters.area[cell_idx],
        )
    end
end
