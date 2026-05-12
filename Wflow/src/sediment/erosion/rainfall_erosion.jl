abstract type AbstractRainfallErosionModel end

"Struct for storing rainfall erosion model variables"
@with_data_lookup struct RainfallErosionModelVariables
    n::Int
    # Total soil erosion rate [t dt-1] from rainfall (splash)
    "rainfall_soil_erosion__mass_flow_rate"
    soil_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing EUROSEM rainfall erosion model boundary conditions"
@kwdef struct RainfallErosionEurosemBC
    n::Int
    # precipitation [mm dt-1]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Interception [mm dt-1]
    interception::Vector{Float64} = fill(MISSING_VALUE, n)
    # Waterlevel on land [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing EUROSEM rainfall erosion model parameters"
@with_data_lookup struct RainfallErosionEurosemParameters
    # Soil detachability factor [g J-1]
    "soil_erosion__rainfall_soil_detachability_factor"
    soil_detachability::Vector{Float64}
    # Exponent EUROSEM [-]
    "soil_erosion__eurosem_exponent"
    eurosem_exponent::Vector{Float64}
    # Canopy height [m]
    "vegetation_canopy__height"
    canopyheight::Vector{Float64}
    # Canopy gap fraction [-]
    "vegetation_canopy__gap_fraction"
    canopygapfraction::Vector{Float64}
    # Fraction of the soil that is covered (eg paved, snow, etc) [-]
    "compacted_soil__area_fraction"
    soilcover_fraction::Vector{Float64}
end

"Initialize EUROSEM rainfall erosion model parameters"
function RainfallErosionEurosemParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
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
    canopyheight =
        ncread(dataset, config, "vegetation_canopy__height", SoilLossModel; sel = indices)
    canopygapfraction = ncread(
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

    eurosem_parameters = RainfallErosionEurosemParameters(
        data_lookup;
        soil_detachability,
        eurosem_exponent,
        canopyheight,
        canopygapfraction,
        soilcover_fraction,
    )
    return eurosem_parameters
end

"EUROSEM rainfall erosion model"
@kwdef struct RainfallErosionEurosemModel <: AbstractRainfallErosionModel
    n::Int
    boundary_conditions::RainfallErosionEurosemBC = RainfallErosionEurosemBC(; n)
    parameters::RainfallErosionEurosemParameters
    variables::RainfallErosionModelVariables = RainfallErosionModelVariables(; n)
end

"Initialize EUROSEM rainfall erosion model"
function RainfallErosionEurosemModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
)
    n = length(indices)
    parameters = RainfallErosionEurosemParameters(dataset, config, indices; data_lookup)
    variables = RainfallErosionModelVariables(data_lookup; n)
    rainfall_erosion_model = RainfallErosionEurosemModel(; n, parameters, variables)
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
        canopyheight,
        canopygapfraction,
        soilcover_fraction,
    ) = rainfall_erosion_model.parameters
    (; soil_erosion_rate) = rainfall_erosion_model.variables

    n = length(precipitation)
    threaded_foreach(1:n; basesize = 1000) do i
        soil_erosion_rate[i] = rainfall_erosion_eurosem(
            precipitation[i],
            interception[i],
            waterlevel[i],
            soil_detachability[i],
            eurosem_exponent[i],
            canopyheight[i],
            canopygapfraction[i],
            soilcover_fraction[i],
            parameters.area[i],
            dt,
        )
    end
end

"Struct for storing ANSWERS rainfall erosion model boundary conditions"
@kwdef struct RainfallErosionAnswersBC
    n::Int
    # precipitation [mm dt-1]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing ANSWERS rainfall erosion model parameters"
@with_data_lookup struct RainfallErosionAnswersParameters
    # Soil erodibility factor [-]
    usle_k::Vector{Float64}
    # Crop management factor [-]
    usle_c::Vector{Float64}
    # ANSWERS rainfall erosion factor [-]
    "soil_erosion__answers_rainfall_factor"
    answers_rainfall_factor::Vector{Float64}
end

"Initialize ANSWERS rainfall erosion model parameters"
function RainfallErosionAnswersParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
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

    answers_parameters = RainfallErosionAnswersParameters(
        data_lookup;
        usle_k,
        usle_c,
        answers_rainfall_factor,
    )
    return answers_parameters
end

"ANSWERS rainfall erosion model"
@kwdef struct RainfallErosionAnswersModel <: AbstractRainfallErosionModel
    n::Int
    boundary_conditions::RainfallErosionAnswersBC = RainfallErosionAnswersBC(; n)
    parameters::RainfallErosionAnswersParameters
    variables::RainfallErosionModelVariables = RainfallErosionModelVariables(; n)
end

"Initialize ANSWERS rainfall erosion model"
function RainfallErosionAnswersModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
)
    n = length(indices)
    parameters = RainfallErosionAnswersParameters(dataset, config, indices; data_lookup)
    variables = RainfallErosionModelVariables(data_lookup; n)
    rainfall_erosion_model = RainfallErosionAnswersModel(; n, parameters, variables)
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

    n = length(precipitation)
    threaded_foreach(1:n; basesize = 1000) do i
        soil_erosion_rate[i] = rainfall_erosion_answers(
            precipitation[i],
            usle_k[i],
            usle_c[i],
            answers_rainfall_factor[i],
            parameters.area[i],
            dt,
        )
    end
end
