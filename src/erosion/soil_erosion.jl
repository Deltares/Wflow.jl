abstract type AbstractSoilErosionModel end

## Total soil erosion and differentiation structs and functions
@get_units @with_kw struct SoilErosionModelVars{T}
    # Total soil erosion
    soil_erosion::Vector{T} | "t dt-1"
    # Total clay erosion
    clay_erosion::Vector{T} | "t dt-1"
    # Total silt erosion
    silt_erosion::Vector{T} | "t dt-1"
    # Total sand erosion
    sand_erosion::Vector{T} | "t dt-1"
    # Total small aggregates erosion
    sagg_erosion::Vector{T} | "t dt-1"
    # Total large aggregates erosion
    lagg_erosion::Vector{T} | "t dt-1"
end

function soil_erosion_model_vars(n)
    vars = SoilErosionModelVars(;
        soil_erosion = fill(mv, n),
        clay_erosion = fill(mv, n),
        silt_erosion = fill(mv, n),
        sand_erosion = fill(mv, n),
        sagg_erosion = fill(mv, n),
        lagg_erosion = fill(mv, n),
    )
    return vars
end

@get_units @with_kw struct SoilErosionBC{T}
    # Rainfall erosion
    rainfall_erosion::Vector{T} | "t dt-1"
    # Overland flow erosion
    overland_flow_erosion::Vector{T} | "m dt-1"
end

function soil_erosion_bc(n)
    bc =
        SoilErosionBC(; rainfall_erosion = fill(mv, n), overland_flow_erosion = fill(mv, n))
    return bc
end

# Parameters for particle differentiation
@get_units @with_kw struct SoilErosionParameters{T}
    # Soil content clay
    clay_fraction::Vector{T} | "-"
    # Soil content silt
    silt_fraction::Vector{T} | "-"
    # Soil content sand
    sand_fraction::Vector{T} | "-"
    # Soil content small aggregates
    sagg_fraction::Vector{T} | "-"
    # Soil content large aggregates
    lagg_fraction::Vector{T} | "-"
end

@get_units @with_kw struct SoilErosionModel{T} <: AbstractSoilErosionModel
    boundary_conditions::SoilErosionBC{T} | "-"
    parameters::SoilErosionParameters{T} | "-"
    variables::SoilErosionModelVars{T} | "-"
end

function initialize_soil_erosion_params(nc, config, inds)
    clay_fraction = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.clay_fraction";
        sel = inds,
        defaults = 0.4,
        type = Float,
    )
    silt_fraction = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.silt_fraction";
        sel = inds,
        defaults = 0.3,
        type = Float,
    )
    sand_fraction = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.sand_fraction";
        sel = inds,
        defaults = 0.3,
        type = Float,
    )
    sagg_fraction = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.sagg_fraction";
        sel = inds,
        defaults = 0.0,
        type = Float,
    )
    lagg_fraction = ncread(
        nc,
        config,
        "vertical.soil_erosion.parameters.lagg_fraction";
        sel = inds,
        defaults = 0.0,
        type = Float,
    )
    # Check that soil fractions sum to 1
    soil_fractions =
        clay_fraction + silt_fraction + sand_fraction + sagg_fraction + lagg_fraction
    if any(abs.(soil_fractions .- 1.0) .> 1e-6)
        error("Particle fractions in the soil must sum to 1")
    end
    soil_parameters = SoilErosionParameters(;
        clay_fraction = clay_fraction,
        silt_fraction = silt_fraction,
        sand_fraction = sand_fraction,
        sagg_fraction = sagg_fraction,
        lagg_fraction = lagg_fraction,
    )

    return soil_parameters
end

function initialize_soil_erosion_model(nc, config, inds)
    n = length(inds)
    vars = soil_erosion_model_vars(n)
    params = initialize_soil_erosion_params(nc, config, inds)
    bc = soil_erosion_bc(n)
    model =
        SoilErosionModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

function update!(model::SoilErosionModel)
    (; rainfall_erosion, overland_flow_erosion) = model.boundary_conditions
    (; clay_fraction, silt_fraction, sand_fraction, sagg_fraction, lagg_fraction) =
        model.parameters
    (; soil_erosion, clay_erosion, silt_erosion, sand_erosion, sagg_erosion, lagg_erosion) =
        model.variables

    n = length(rainfall_erosion)
    threaded_foreach(1:n; basesize = 1000) do i
        soil_erosion[i],
        clay_erosion[i],
        silt_erosion[i],
        sand_erosion[i],
        sagg_erosion[i],
        lagg_erosion[i] = total_soil_erosion(
            rainfall_erosion[i],
            overland_flow_erosion[i],
            clay_fraction[i],
            silt_fraction[i],
            sand_fraction[i],
            sagg_fraction[i],
            lagg_fraction[i],
        )
    end
end