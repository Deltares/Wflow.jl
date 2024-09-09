abstract type AbstractRainfallErosionModel end

struct NoRainfallErosionModel <: RainfallErosionModel end

## General rainfall erosion functions and structs
@get_units @with_kw struct RainfallErosionModelVars{T}
    # Total soil erosion from rainfall (splash)
    rainfall_erosion::Vector{T} | "t dt-1"
end

function rainfall_erosion_model_vars(n)
    vars = RainfallErosionModelVars(; rainfall_erosion = fill(mv, n))
    return vars
end

@get_units @with_kw struct RainfallErosionBC{T}
    # Precipitation
    precip::Vector{T} | "mm dt-1"
end

function rainfall_erosion_bc(n)
    bc = RainfallErosionBC(; precip = fill(mv, n))
    return bc
end

## EUROSEM specific structs and functions for rainfall erosion
@get_units @with_kw struct RainfallErosionEurosemBC{T}
    # Precipitation
    precip::Vector{T} | "mm dt-1"
    # Interception
    interception::Vector{T} | "mm dt-1"
    # Water level
    waterlevel::Vector{T} | "m"
end

function rainfall_erosion_eurosem_bc(n)
    bc = RainfallErosionEurosemBC(;
        precip = fill(mv, n),
        interception = fill(mv, n),
        waterlevel = fill(mv, n),
    )
    return bc
end

@get_units @with_kw struct RainfallErosionEurosemParameters{T}
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

@get_units @with_kw struct RainfallErosionEurosemModel{T} <: AbstractRainfallErosionModel
    boundary_conditions::RainfallErosionEurosemBC{T} | "-"
    parameters::RainfallErosionEurosemParameters{T} | "-"
    variables::RainfallErosionModelVars{T} | "-"
end

function initialize_eurosem_params(nc, config, inds)
    n = length(inds)
    soil_detachability = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.soil_detachability";
        sel = inds,
        defaults = 0.6,
        type = Float,
    )
    eurosem_exponent = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.eurosem_exponent";
        sel = inds,
        defaults = 2.0,
        type = Float,
    )
    canopyheight = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.canopyheight";
        sel = inds,
        defaults = 0.5,
        type = Float,
    )
    canopygapfraction = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.canopygapfraction";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    soilcover_fraction =
        ncread(nc, config, "vertical.pathfrac"; sel = inds, defaults = 0.01, type = Float)
    eurosem_parameters = RainfallErosionEurosemParameters(;
        soil_detachability = soil_detachability,
        eurosem_exponent = eurosem_exponent,
        canopyheight = canopyheight,
        canopygapfraction = canopygapfraction,
        soilcover_fraction = soilcover_fraction,
    )
    return eurosem_parameters
end

function initialize_eurosem_rainfall_erosion_model(nc, config, inds)
    n = length(inds)
    vars = rainfall_erosion_model_vars(n)
    params = initialize_eurosem_params(nc, config, inds)
    bc = rainfall_erosion_eurosem_bc(n)
    model = RainfallErosionEurosemModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::RainfallErosionEurosemModel, area, dt)
    (; precip, interception, waterlevel) = model.boundary_conditions
    (;
        soil_detachability,
        eurosem_exponent,
        canopyheight,
        canopygapfraction,
        soilcover_fraction,
    ) = model.parameters
    (; rainfall_erosion) = model.variables

    n = length(precip)
    threaded_foreach(1:n; basesize = 1000) do i
        rainfall_erosion[i] = rainfall_erosion_eurosem(
            precip[i],
            interception[i],
            waterlevel[i],
            soil_detachability[i],
            eurosem_exponent[i],
            canopyheight[i],
            canopygapfraction[i],
            soilcover_fraction[i],
            area[i],
            dt,
        )
    end
end

# ANSWERS specific structs and functions for rainfall erosion
@get_units @with_kw struct RainfallErosionAnswersParameters{T}
    # Soil erodibility factor
    usle_k::Vector{T} | "-"
    # Crop management factor
    usle_c::Vector{T} | "-"
end

@get_units @with_kw struct RainfallErosionAnswersModel{T} <: AbstractRainfallErosionModel
    boundary_conditions::RainfallErosionBC{T} | "-"
    parameters::RainfallErosionAnswersParameters{T} | "-"
    variables::RainfallErosionModelVars{T} | "-"
end

function initialize_answers_params_rainfall(nc, config, inds)
    usle_k = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.usle_k";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    usle_c = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.usle_c";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )
    answers_parameters =
        RainfallErosionAnswersParameters(; usle_k = usle_k, usle_c = usle_c)
    return answers_parameters
end

function initialize_answers_rainfall_erosion_model(nc, config, inds)
    n = length(inds)
    vars = rainfall_erosion_model_vars(n)
    params = initialize_answers_params_rainfall(nc, config, inds)
    bc = rainfall_erosion_bc(n)
    model = RainfallErosionAnswersModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update!(model::RainfallErosionAnswersModel, area, dt)
    (; precip) = model.boundary_conditions
    (; usle_k, usle_c) = model.parameters
    (; rainfall_erosion) = model.variables

    n = length(precip)
    threaded_foreach(1:n; basesize = 1000) do i
        rainfall_erosion[i] =
            rainfall_erosion_answers(precip[i], usle_k[i], usle_c[i], area[i], dt)
    end
end

function update!(model::NoRainfallErosionModel)
    return nothing
end

get_rainfall_erosion(model::NoRainfallErosionModel) = 0.0
get_rainfall_erosion(model::RainfallErosionEurosemModel) = model.variables.rainfall_erosion
get_rainfall_erosion(model::RainfallErosionAnswersModel) = model.variables.rainfall_erosion