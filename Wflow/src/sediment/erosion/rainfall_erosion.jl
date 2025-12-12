abstract type AbstractRainfallErosionModel end

"Struct for storing rainfall erosion model variables"
@with_kw struct RainfallErosionModelVariables
    n::Int
    # Total soil erosion rate [t dt-1] from rainfall (splash)
    soil_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing EUROSEM rainfall erosion model boundary conditions"
@with_kw struct RainfallErosionEurosemBC
    n::Int
    # precipitation [mm dt-1]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Interception [mm dt-1]
    interception::Vector{Float64} = fill(MISSING_VALUE, n)
    # Waterlevel on land [m]
    waterlevel::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing EUROSEM rainfall erosion model parameters"
@with_kw struct RainfallErosionEurosemParameters
    # Soil detachability factor [g J-1]
    soil_detachability::Vector{Float64}
    # Exponent EUROSEM [-]
    eurosem_exponent::Vector{Float64}
    # Canopy height [m]
    canopyheight::Vector{Float64}
    # Canopy gap fraction [-]
    canopygapfraction::Vector{Float64}
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
        "soil_erosion__rainfall_soil_detachability_factor";
        sel = indices,
        defaults = 0.6,
        type = Float64,
    )
    eurosem_exponent = ncread(
        dataset,
        config,
        "soil_erosion__eurosem_exponent";
        sel = indices,
        defaults = 2.0,
        type = Float64,
    )
    canopyheight = ncread(
        dataset,
        config,
        "vegetation_canopy__height";
        sel = indices,
        defaults = 0.5,
        type = Float64,
    )
    canopygapfraction = ncread(
        dataset,
        config,
        "vegetation_canopy__gap_fraction";
        sel = indices,
        defaults = 0.1,
        type = Float64,
    )
    soilcover_fraction = ncread(
        dataset,
        config,
        "compacted_soil__area_fraction";
        sel = indices,
        defaults = 0.01,
        type = Float64,
    )

    eurosem_parameters = RainfallErosionEurosemParameters(;
        soil_detachability,
        eurosem_exponent,
        canopyheight,
        canopygapfraction,
        soilcover_fraction,
    )
    return eurosem_parameters
end

"EUROSEM rainfall erosion model"
@with_kw struct RainfallErosionEurosemModel <: AbstractRainfallErosionModel
    n::Int
    boundary_conditions::RainfallErosionEurosemBC = RainfallErosionEurosemBC(; n)
    parameters::RainfallErosionEurosemParameters
    variables::RainfallErosionModelVariables = RainfallErosionModelVariables(; n)
end

"Initialize EUROSEM rainfall erosion model"
function RainfallErosionEurosemModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = RainfallErosionEurosemParameters(dataset, config, indices)
    model = RainfallErosionEurosemModel(; n, parameters)
    return model
end

"Update EUROSEM rainfall erosion model boundary conditions for a single timestep"
function update_boundary_conditions!(
    model::RainfallErosionEurosemModel,
    atmospheric_forcing::AtmosphericForcing,
    hydrological_forcing::HydrologicalForcing,
)
    (; precipitation, interception, waterlevel) = model.boundary_conditions
    @. precipitation = atmospheric_forcing.precipitation
    @. waterlevel = hydrological_forcing.waterlevel_land
    @. interception = hydrological_forcing.interception
end

"Update EUROSEM rainfall erosion model for a single timestep"
function update!(
    model::RainfallErosionEurosemModel,
    parameters::LandParameters,
    dt::Float64,
)
    (; precipitation, interception, waterlevel) = model.boundary_conditions
    (;
        soil_detachability,
        eurosem_exponent,
        canopyheight,
        canopygapfraction,
        soilcover_fraction,
    ) = model.parameters
    (; soil_erosion_rate) = model.variables

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
@with_kw struct RainfallErosionAnswersBC
    n::Int
    # precipitation [mm dt-1]
    precipitation::Vector{Float64} = fill(MISSING_VALUE, n)
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
    usle_k = ncread(
        dataset,
        config,
        "soil_erosion__usle_k_factor";
        sel = indices,
        defaults = 0.1,
        type = Float64,
    )
    usle_c = ncread(
        dataset,
        config,
        "soil_erosion__usle_c_factor";
        sel = indices,
        defaults = 0.01,
        type = Float64,
    )
    answers_rainfall_factor = ncread(
        dataset,
        config,
        "soil_erosion__answers_rainfall_factor";
        sel = indices,
        defaults = 0.108,
        type = Float64,
    )

    answers_parameters =
        RainfallErosionAnswersParameters(; usle_k, usle_c, answers_rainfall_factor)
    return answers_parameters
end

"ANSWERS rainfall erosion model"
@with_kw struct RainfallErosionAnswersModel <: AbstractRainfallErosionModel
    n::Int
    boundary_conditions::RainfallErosionAnswersBC = RainfallErosionAnswersBC(; n)
    parameters::RainfallErosionAnswersParameters
    variables::RainfallErosionModelVariables = RainfallErosionModelVariables(; n)
end

"Initialize ANSWERS rainfall erosion model"
function RainfallErosionAnswersModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    parameters = RainfallErosionAnswersParameters(dataset, config, indices)
    model = RainfallErosionAnswersModel(; n, parameters)
    return model
end

"Update ANSWERS rainfall erosion model boundary conditions for a single timestep"
function update_boundary_conditions!(
    model::RainfallErosionAnswersModel,
    atmospheric_forcing::AtmosphericForcing,
    hydrological_forcing::HydrologicalForcing,
)
    (; precipitation) = model.boundary_conditions
    @. precipitation = atmospheric_forcing.precipitation
end

"Update ANSWERS rainfall erosion model for a single timestep"
function update!(
    model::RainfallErosionAnswersModel,
    parameters::LandParameters,
    dt::Float64,
)
    (; precipitation) = model.boundary_conditions
    (; usle_k, usle_c, answers_rainfall_factor) = model.parameters
    (; soil_erosion_rate) = model.variables

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
