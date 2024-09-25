abstract type AbstractRainfallErosionModel{T} end

struct NoRainfallErosionModel{T} <: AbstractRainfallErosionModel{T} end

## General rainfall erosion functions and structs
@get_units @with_kw struct RainfallErosionModelVariables{T}
    # Total soil erosion from rainfall (splash)
    amount::Vector{T} | "t dt-1"
end

function RainfallErosionModelVariables(n; amount::Vector{T} = fill(mv, n)) where {T}
    return RainfallErosionModelVariables{T}(; amount = amount)
end

## EUROSEM specific structs and functions for rainfall erosion
@get_units @with_kw struct RainfallErosionEurosemBC{T}
    # precipitation
    precipitation::Vector{T} | "mm dt-1"
    # Interception
    interception::Vector{T} | "mm dt-1"
    # Waterlevel on land
    waterlevel::Vector{T} | "m"
end

function RainfallErosionEurosemBC(
    n;
    precipitation::Vector{T} = fill(mv, n),
    interception::Vector{T} = fill(0.0, n),
    waterlevel::Vector{T} = fill(mv, n),
) where {T}
    return RainfallErosionEurosemBC{T}(;
        precipitation = precipitation,
        interception = interception,
        waterlevel = waterlevel,
    )
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

function RainfallErosionEurosemParameters(nc, config, inds)
    soil_detachability = ncread(
        nc,
        config,
        "vertical.rainfall_erosion.parameters.soil_detachability";
        sel = inds,
        defaults = 0.6,
        type = Float,
    )
    eurosem_exponent = ncread(
        nc,
        config,
        "vertical.rainfall_erosion.parameters.eurosem_exponent";
        sel = inds,
        defaults = 2.0,
        type = Float,
    )
    canopyheight = ncread(
        nc,
        config,
        "vertical.rainfall_erosion.parameters.canopyheight";
        sel = inds,
        defaults = 0.5,
        type = Float,
    )
    canopygapfraction = ncread(
        nc,
        config,
        "vertical.rainfall_erosion.parameters.canopygapfraction";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    soilcover_fraction = ncread(
        nc,
        config,
        "vertical.rainfall_erosion.parameters.pathfrac";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )
    eurosem_parameters = RainfallErosionEurosemParameters(;
        soil_detachability = soil_detachability,
        eurosem_exponent = eurosem_exponent,
        canopyheight = canopyheight,
        canopygapfraction = canopygapfraction,
        soilcover_fraction = soilcover_fraction,
    )
    return eurosem_parameters
end

@with_kw struct RainfallErosionEurosemModel{T} <: AbstractRainfallErosionModel{T}
    boundary_conditions::RainfallErosionEurosemBC{T}
    parameters::RainfallErosionEurosemParameters{T}
    variables::RainfallErosionModelVariables{T}
end

function RainfallErosionEurosemModel(nc, config, inds)
    n = length(inds)
    vars = RainfallErosionModelVariables(n)
    params = RainfallErosionEurosemParameters(nc, config, inds)
    bc = RainfallErosionEurosemBC(n)
    model = RainfallErosionEurosemModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update_boundary_conditions!(
    model::RainfallErosionEurosemModel,
    hydrometeo_forcing::HydrometeoForcing,
)
    (; precipitation, waterlevel) = model.boundary_conditions
    @. precipitation = hydrometeo_forcing.precipitation
    @. waterlevel = model.boundary_conditions.waterlevel_land
end

function update!(model::RainfallErosionEurosemModel, geometry::LandGeometry, ts)
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
            ts,
        )
    end
end

### ANSWERS specific structs and functions for rainfall erosion
@get_units @with_kw struct RainfallErosionAnswersBC{T}
    # precipitation
    precipitation::Vector{T} | "mm dt-1"
end

function RainfallErosionAnswersBC(n; precipitation::Vector{T} = fill(mv, n)) where {T}
    return RainfallErosionAnswersBC{T}(; precipitation = precipitation)
end

@get_units @with_kw struct RainfallErosionAnswersParameters{T}
    # Soil erodibility factor
    usle_k::Vector{T} | "-"
    # Crop management factor
    usle_c::Vector{T} | "-"
end

function RainfallErosionAnswersParameters(nc, config, inds)
    usle_k = ncread(
        nc,
        config,
        "vertical.rainfall_erosion.parameters.usle_k";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )
    usle_c = ncread(
        nc,
        config,
        "vertical.rainfall_erosion.parameters.usle_c";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )
    answers_parameters =
        RainfallErosionAnswersParameters(; usle_k = usle_k, usle_c = usle_c)
    return answers_parameters
end

@with_kw struct RainfallErosionAnswersModel{T} <: AbstractRainfallErosionModel{T}
    boundary_conditions::RainfallErosionAnswersBC{T}
    parameters::RainfallErosionAnswersParameters{T}
    variables::RainfallErosionModelVariables{T}
end

function RainfallErosionAnswersModel(nc, config, inds)
    n = length(inds)
    bc = RainfallErosionAnswersBC(n)
    vars = RainfallErosionModelVariables(n)
    params = RainfallErosionAnswersParameters(nc, config, inds)
    model = RainfallErosionAnswersModel(;
        boundary_conditions = bc,
        parameters = params,
        variables = vars,
    )
    return model
end

function update_boundary_conditions!(
    model::RainfallErosionAnswersModel,
    hydrometeo_forcing::HydrometeoForcing,
)
    (; precipitation) = model.boundary_conditions
    @. precipitation = hydrometeo_forcing.precipitation
end

function update!(model::RainfallErosionAnswersModel, geometry::LandGeometry, ts)
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
            ts,
        )
    end
end

function update_boundary_conditions!(
    model::NoRainfallErosionModel,
    hydrometeo_forcing::HydrometeoForcing,
)
    return nothing
end

function update!(model::NoRainfallErosionModel)
    return nothing
end

get_rainfall_erosion(model::NoRainfallErosionModel) = 0.0
get_rainfall_erosion(model::AbstractRainfallErosionModel) = model.variables.amount