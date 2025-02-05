abstract type AbstractRainfallErosionModel{T} end

"Struct for storing rainfall erosion model variables"
@get_units @grid_loc @with_kw struct RainfallErosionModelVariables{T}
    # Total soil erosion from rainfall (splash)
    amount::Vector{T} | "t dt-1"
end

"Initialize rainfall erosion model variables"
function RainfallErosionModelVariables(n; amount::Vector{T} = fill(mv, n)) where {T}
    return RainfallErosionModelVariables{T}(; amount = amount)
end

"Struct for storing EUROSEM rainfall erosion model boundary conditions"
@get_units @grid_loc @with_kw struct RainfallErosionEurosemBC{T}
    # precipitation
    precipitation::Vector{T} | "mm dt-1"
    # Interception
    interception::Vector{T} | "mm dt-1"
    # Waterlevel on land
    waterlevel::Vector{T} | "m"
end

"Initialize EUROSEM rainfall erosion model boundary conditions"
function RainfallErosionEurosemBC(
    n;
    precipitation::Vector{T} = fill(mv, n),
    interception::Vector{T} = fill(mv, n),
    waterlevel::Vector{T} = fill(mv, n),
) where {T}
    return RainfallErosionEurosemBC{T}(;
        precipitation = precipitation,
        interception = interception,
        waterlevel = waterlevel,
    )
end

"Struct for storing EUROSEM rainfall erosion model parameters"
@get_units @grid_loc @with_kw struct RainfallErosionEurosemParameters{T}
    # Soil detachability factor
    soil_detachability::Vector{T} | "g J-1"
    # Exponent EUROSEM
    eurosem_exponent::Vector{T} | "-"
    # Canopy height
    canopyheight::Vector{T} | "m"
    # Canopy gap fraction
    canopygapfraction::Vector{T} | "-"
    # Fraction of the soil that is covered (eg paved, snow, etc)
    soilcover_fraction::Vector{T} | "-"
end

"Initialize EUROSEM rainfall erosion model parameters"
function RainfallErosionEurosemParameters(dataset, config, indices)
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
@with_kw struct RainfallErosionEurosemModel{T} <: AbstractRainfallErosionModel{T}
    boundary_conditions::RainfallErosionEurosemBC{T}
    parameters::RainfallErosionEurosemParameters{T}
    variables::RainfallErosionModelVariables{T}
end

"Initialize EUROSEM rainfall erosion model"
function RainfallErosionEurosemModel(dataset, config, indices)
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
function update!(model::RainfallErosionEurosemModel, geometry::LandGeometry, dt)
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
            geometry.area[i],
            dt,
        )
    end
end

"Struct for storing ANSWERS rainfall erosion model boundary conditions"
@get_units @grid_loc @with_kw struct RainfallErosionAnswersBC{T}
    # precipitation
    precipitation::Vector{T} | "mm dt-1"
end

"Initialize ANSWERS rainfall erosion model boundary conditions"
function RainfallErosionAnswersBC(n; precipitation::Vector{T} = fill(mv, n)) where {T}
    return RainfallErosionAnswersBC{T}(; precipitation = precipitation)
end

"Struct for storing ANSWERS rainfall erosion model parameters"
@get_units @grid_loc @with_kw struct RainfallErosionAnswersParameters{T}
    # Soil erodibility factor
    usle_k::Vector{T} | "-"
    # Crop management factor
    usle_c::Vector{T} | "-"
end

"Initialize ANSWERS rainfall erosion model parameters"
function RainfallErosionAnswersParameters(dataset, config, indices)
    lens = lens_input_parameter(config, "soil_erosion__usle_k_factor")
    usle_k = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)
    lens = lens_input_parameter(config, "soil_erosion__usle_c_factor")
    usle_c = ncread(dataset, config, lens; sel = indices, defaults = 0.01, type = Float)
    answers_parameters =
        RainfallErosionAnswersParameters(; usle_k = usle_k, usle_c = usle_c)
    return answers_parameters
end

"ANSWERS rainfall erosion model"
@with_kw struct RainfallErosionAnswersModel{T} <: AbstractRainfallErosionModel{T}
    boundary_conditions::RainfallErosionAnswersBC{T}
    parameters::RainfallErosionAnswersParameters{T}
    variables::RainfallErosionModelVariables{T}
end

"Initialize ANSWERS rainfall erosion model"
function RainfallErosionAnswersModel(dataset, config, indices)
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
function update!(model::RainfallErosionAnswersModel, geometry::LandGeometry, dt)
    (; precipitation) = model.boundary_conditions
    (; usle_k, usle_c) = model.parameters
    (; amount) = model.variables

    n = length(precipitation)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i] = rainfall_erosion_answers(
            precipitation[i],
            usle_k[i],
            usle_c[i],
            geometry.area[i],
            dt,
        )
    end
end