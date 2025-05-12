abstract type AbstractRainfallErosionModel end

"Struct for storing rainfall erosion model variables"
@with_kw struct RainfallErosionModelVariables
    # Total soil erosion rate [t dt-1] from rainfall (splash)
    amount::Vector{Float}
end

"Initialize rainfall erosion model variables"
function RainfallErosionModelVariables(
    n::Int;
    amount::Vector{Float} = fill(MISSING_VALUE, n),
)
    return RainfallErosionModelVariables(; amount = amount)
end

"Struct for storing EUROSEM rainfall erosion model boundary conditions"
@with_kw struct RainfallErosionEurosemBC
    # precipitation [mm dt-1]
    precipitation::Vector{Float}
    # Interception [mm dt-1]
    interception::Vector{Float}
    # Waterlevel on land [m]
    waterlevel::Vector{Float}
end

"Initialize EUROSEM rainfall erosion model boundary conditions"
function RainfallErosionEurosemBC(
    n::Int;
    precipitation::Vector{Float} = fill(MISSING_VALUE, n),
    interception::Vector{Float} = fill(MISSING_VALUE, n),
    waterlevel::Vector{Float} = fill(MISSING_VALUE, n),
)
    return RainfallErosionEurosemBC(;
        precipitation = precipitation,
        interception = interception,
        waterlevel = waterlevel,
    )
end

"Struct for storing EUROSEM rainfall erosion model parameters"
@with_kw struct RainfallErosionEurosemParameters
    # Soil detachability factor [g J-1]
    soil_detachability::Vector{Float}
    # Exponent EUROSEM [-]
    eurosem_exponent::Vector{Float}
    # Canopy height [m]
    canopyheight::Vector{Float}
    # Canopy gap fraction [-]
    canopygapfraction::Vector{Float}
    # Fraction of the soil that is covered (eg paved, snow, etc) [-]
    soilcover_fraction::Vector{Float}
end

"Initialize EUROSEM rainfall erosion model parameters"
function RainfallErosionEurosemParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    lens = lens_input_parameter(config, "soil_erosion__rainfall_soil_detachability_factor")
    soil_detachability =
        ncread(dataset, config, lens; sel = indices, defaults = 0.6, type = Float)
    lens = lens_input_parameter(config, "soil_erosion__eurosem_exponent")
    eurosem_exponent =
        ncread(dataset, config, lens; sel = indices, defaults = 2.0, type = Float)
    lens = lens_input_parameter(config, "vegetation_canopy__height")
    canopyheight =
        ncread(dataset, config, lens; sel = indices, defaults = 0.5, type = Float)
    lens = lens_input_parameter(config, "vegetation_canopy__gap_fraction")
    canopygapfraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)
    lens = lens_input_parameter(config, "soil~compacted__area_fraction")
    soilcover_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.01, type = Float)

    eurosem_parameters = RainfallErosionEurosemParameters(;
        soil_detachability = soil_detachability,
        eurosem_exponent = eurosem_exponent,
        canopyheight = canopyheight,
        canopygapfraction = canopygapfraction,
        soilcover_fraction = soilcover_fraction,
    )
    return eurosem_parameters
end

"EUROSEM rainfall erosion model"
@with_kw struct RainfallErosionEurosemModel <: AbstractRainfallErosionModel
    boundary_conditions::RainfallErosionEurosemBC
    parameters::RainfallErosionEurosemParameters
    variables::RainfallErosionModelVariables
end

"Initialize EUROSEM rainfall erosion model"
function RainfallErosionEurosemModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    vars = RainfallErosionModelVariables(n)
    params = RainfallErosionEurosemParameters(dataset, config, indices)
    bc = RainfallErosionEurosemBC(n)
    model = RainfallErosionEurosemModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
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
function update!(model::RainfallErosionEurosemModel, parameters::LandParameters, dt::Float)
    (; precipitation, interception, waterlevel) = model.boundary_conditions
    (;
        soil_detachability,
        eurosem_exponent,
        canopyheight,
        canopygapfraction,
        soilcover_fraction,
    ) = model.parameters
    (; amount) = model.variables

    n = length(precipitation)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = rainfall_erosion_eurosem(
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
    # precipitation [mm dt-1]
    precipitation::Vector{Float}
end

"Initialize ANSWERS rainfall erosion model boundary conditions"
function RainfallErosionAnswersBC(
    n::Int;
    precipitation::Vector{Float} = fill(MISSING_VALUE, n),
)
    return RainfallErosionAnswersBC(; precipitation = precipitation)
end

"Struct for storing ANSWERS rainfall erosion model parameters"
@with_kw struct RainfallErosionAnswersParameters
    # Soil erodibility factor [-]
    usle_k::Vector{Float}
    # Crop management factor [-]
    usle_c::Vector{Float}
    # ANSWERS rainfall erosion factor [-]
    answers_rainfall_factor::Vector{Float}
end

"Initialize ANSWERS rainfall erosion model parameters"
function RainfallErosionAnswersParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    lens = lens_input_parameter(config, "soil_erosion__usle_k_factor")
    usle_k = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)
    lens = lens_input_parameter(config, "soil_erosion__usle_c_factor")
    usle_c = ncread(dataset, config, lens; sel = indices, defaults = 0.01, type = Float)
    lens = lens_input_parameter(config, "soil_erosion__answers_rainfall_factor")
    answers_rainfall_factor =
        ncread(dataset, config, lens; sel = indices, defaults = 0.108, type = Float)

    answers_parameters = RainfallErosionAnswersParameters(;
        usle_k = usle_k,
        usle_c = usle_c,
        answers_rainfall_factor = answers_rainfall_factor,
    )
    return answers_parameters
end

"ANSWERS rainfall erosion model"
@with_kw struct RainfallErosionAnswersModel <: AbstractRainfallErosionModel
    boundary_conditions::RainfallErosionAnswersBC
    parameters::RainfallErosionAnswersParameters
    variables::RainfallErosionModelVariables
end

"Initialize ANSWERS rainfall erosion model"
function RainfallErosionAnswersModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    bc = RainfallErosionAnswersBC(n)
    vars = RainfallErosionModelVariables(n)
    params = RainfallErosionAnswersParameters(dataset, config, indices)
    model = RainfallErosionAnswersModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
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
function update!(model::RainfallErosionAnswersModel, parameters::LandParameters, dt::Float)
    (; precipitation) = model.boundary_conditions
    (; usle_k, usle_c, answers_rainfall_factor) = model.parameters
    (; amount) = model.variables

    n = length(precipitation)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = rainfall_erosion_answers(
            precipitation[i],
            usle_k[i],
            usle_c[i],
            answers_rainfall_factor[i],
            parameters.area[i],
            dt,
        )
    end
end