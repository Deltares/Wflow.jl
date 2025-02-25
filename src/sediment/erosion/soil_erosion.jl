abstract type AbstractSoilErosionModel{T} end

"Struct for storing total soil erosion with differentiation model variables"
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

"Initialize soil erosion model variables"
function SoilErosionModelVariables(
    n;
    amount::Vector{T} = fill(MISSING_VALUE, n),
    clay::Vector{T} = fill(MISSING_VALUE, n),
    silt::Vector{T} = fill(MISSING_VALUE, n),
    sand::Vector{T} = fill(MISSING_VALUE, n),
    sagg::Vector{T} = fill(MISSING_VALUE, n),
    lagg::Vector{T} = fill(MISSING_VALUE, n),
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

"Struct for storing soil erosion model boundary conditions"
@get_units @grid_loc @with_kw struct SoilErosionBC{T}
    # Rainfall erosion
    rainfall_erosion::Vector{T} | "t dt-1"
    # Overland flow erosion
    overland_flow_erosion::Vector{T} | "m dt-1"
end

"Initialize soil erosion model boundary conditions"
function SoilErosionBC(
    n;
    rainfall_erosion::Vector{T} = fill(MISSING_VALUE, n),
    overland_flow_erosion::Vector{T} = fill(MISSING_VALUE, n),
) where {T}
    return SoilErosionBC{T}(;
        rainfall_erosion = rainfall_erosion,
        overland_flow_erosion = overland_flow_erosion,
    )
end

"Struct for storing soil erosion model parameters"
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

"Initialize soil erosion model parameters"
function SoilErosionParameters(dataset, config, indices)
    lens = lens_input_parameter(config, "soil_clay__mass_fraction")
    clay_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.4, type = Float64)
    lens = lens_input_parameter(config, "soil_silt__mass_fraction")
    silt_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.3, type = Float64)
    lens = lens_input_parameter(config, "soil_sand__mass_fraction")
    sand_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.3, type = Float64)
    lens = lens_input_parameter(config, "soil_aggregates~small__mass_fraction")
    sagg_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float64)
    lens = lens_input_parameter(config, "soil_aggregates~large__mass_fraction")
    lagg_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float64)
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

"Total soil erosion with differentiation model"
@with_kw struct SoilErosionModel{T} <: AbstractSoilErosionModel{T}
    boundary_conditions::SoilErosionBC{T}
    parameters::SoilErosionParameters{T}
    variables::SoilErosionModelVariables{T}
end

"Initialize soil erosion model"
function SoilErosionModel(dataset, config, indices)
    n = length(indices)
    vars = SoilErosionModelVariables(n)
    params = SoilErosionParameters(dataset, config, indices)
    bc = SoilErosionBC(n)
    model =
        SoilErosionModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

"Update boundary conditions for soil erosion model"
function update_boundary_conditions!(
    model::SoilErosionModel,
    rainfall_erosion::AbstractRainfallErosionModel,
    overland_flow_erosion::AbstractOverlandFlowErosionModel,
)
    re = rainfall_erosion.variables.amount
    ole = overland_flow_erosion.variables.amount
    (; rainfall_erosion, overland_flow_erosion) = model.boundary_conditions
    @. rainfall_erosion = re
    @. overland_flow_erosion = ole
end

"Update soil erosion model for a single timestep"
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