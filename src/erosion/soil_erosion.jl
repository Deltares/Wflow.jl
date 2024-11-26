abstract type AbstractSoilErosionModel{T} end

## Total soil erosion and differentiation structs and functions
@get_units @grid_loc @with_kw struct SoilErosionModelVariables{T}
    # Total soil erosion
    amount::Vector{T} | "t dt-1"
    # Total clay erosion
    clay::Vector{T} | "t dt-1"
    # Total silt erosion
    silt::Vector{T} | "t dt-1"
    # Total sand erosion
    sand::Vector{T} | "t dt-1"
    # Total small aggregates erosion
    sagg::Vector{T} | "t dt-1"
    # Total large aggregates erosion
    lagg::Vector{T} | "t dt-1"
end

function SoilErosionModelVariables(
    n;
    amount::Vector{T} = fill(mv, n),
    clay::Vector{T} = fill(mv, n),
    silt::Vector{T} = fill(mv, n),
    sand::Vector{T} = fill(mv, n),
    sagg::Vector{T} = fill(mv, n),
    lagg::Vector{T} = fill(mv, n),
) where {T}
    return SoilErosionModelVariables{T}(;
        amount = amount,
        clay = clay,
        silt = silt,
        sand = sand,
        sagg = sagg,
        lagg = lagg,
    )
end

@get_units @grid_loc @with_kw struct SoilErosionBC{T}
    # Rainfall erosion
    rainfall_erosion::Vector{T} | "t dt-1"
    # Overland flow erosion
    overland_flow_erosion::Vector{T} | "m dt-1"
end

function SoilErosionBC(
    n;
    rainfall_erosion::Vector{T} = fill(mv, n),
    overland_flow_erosion::Vector{T} = fill(mv, n),
) where {T}
    return SoilErosionBC{T}(;
        rainfall_erosion = rainfall_erosion,
        overland_flow_erosion = overland_flow_erosion,
    )
end

# Parameters for particle differentiation
@get_units @grid_loc @with_kw struct SoilErosionParameters{T}
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

function SoilErosionParameters(nc, config, inds)
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
    if any(abs.(soil_fractions .- 1.0) .> 1e-3)
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

@with_kw struct SoilErosionModel{T} <: AbstractSoilErosionModel{T}
    boundary_conditions::SoilErosionBC{T}
    parameters::SoilErosionParameters{T}
    variables::SoilErosionModelVariables{T}
end

function SoilErosionModel(nc, config, inds)
    n = length(inds)
    vars = SoilErosionModelVariables(n)
    params = SoilErosionParameters(nc, config, inds)
    bc = SoilErosionBC(n)
    model =
        SoilErosionModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

function update_boundary_conditions!(
    model::SoilErosionModel,
    rainfall_erosion::AbstractRainfallErosionModel,
    overland_flow_erosion::AbstractOverlandFlowErosionModel,
)
    re = get_rainfall_erosion(rainfall_erosion)
    ole = get_overland_flow_erosion(overland_flow_erosion)
    (; rainfall_erosion, overland_flow_erosion) = model.boundary_conditions
    @. rainfall_erosion = re
    @. overland_flow_erosion = ole
end

function update!(model::SoilErosionModel)
    (; rainfall_erosion, overland_flow_erosion) = model.boundary_conditions
    (; clay_fraction, silt_fraction, sand_fraction, sagg_fraction, lagg_fraction) =
        model.parameters
    (; amount, clay, silt, sand, sagg, lagg) = model.variables

    n = length(rainfall_erosion)
    threaded_foreach(1:n; basesize = 1000) do i
        amount[i], clay[i], silt[i], sand[i], sagg[i], lagg[i] = total_soil_erosion(
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