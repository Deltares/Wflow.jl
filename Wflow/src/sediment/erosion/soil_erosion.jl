abstract type AbstractSoilErosionModel end

"Struct for storing total soil erosion with differentiation model variables"
@with_kw struct SoilErosionModelVariables
    n::Int
    # Total soil erosion rate [t dt-1]
    amount::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total clay erosion rate [t dt-1]
    clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total silt erosion rate [t dt-1]
    silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total sand erosion rate [t dt-1]
    sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total small aggregates erosion rate [t dt-1]
    sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total large aggregates erosion rate [t dt-1]
    lagg::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing soil erosion model boundary conditions"
@with_kw struct SoilErosionBC
    n::Int
    # Rainfall erosion rate [t dt-1]
    rainfall_erosion::Vector{Float64} = fill(MISSING_VALUE, n)
    # Overland flow erosion rate [t dt-1]
    overland_flow_erosion::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing soil erosion model parameters"
@with_kw struct SoilErosionParameters
    # Soil content clay [-]
    clay_fraction::Vector{Float64}
    # Soil content silt [-]
    silt_fraction::Vector{Float64}
    # Soil content sand [-]
    sand_fraction::Vector{Float64}
    # Soil content small aggregates [-]
    sagg_fraction::Vector{Float64}
    # Soil content large aggregates [-]
    lagg_fraction::Vector{Float64}
end

"Initialize soil erosion model parameters"
function SoilErosionParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    lens = lens_input_parameter(config, "soil_clay__mass_fraction")
    clay_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.4, type = Float64)
    lens = lens_input_parameter(config, "soil_silt__mass_fraction")
    silt_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.3, type = Float64)
    lens = lens_input_parameter(config, "soil_sand__mass_fraction")
    sand_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.3, type = Float64)
    lens = lens_input_parameter(config, "soil_small_aggregates__mass_fraction")
    sagg_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float64)
    lens = lens_input_parameter(config, "soil_large_aggregates__mass_fraction")
    lagg_fraction =
        ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float64)
    # Check that soil fractions sum to 1
    soil_fractions =
        clay_fraction + silt_fraction + sand_fraction + sagg_fraction + lagg_fraction
    if any(abs.(soil_fractions .- 1.0) .> 1e-3)
        error("Particle fractions in the soil must sum to 1")
    end
    soil_parameters = SoilErosionParameters(;
        clay_fraction,
        silt_fraction,
        sand_fraction,
        sagg_fraction,
        lagg_fraction,
    )

    return soil_parameters
end

"Total soil erosion with differentiation model"
@with_kw struct SoilErosionModel <: AbstractSoilErosionModel
    boundary_conditions::SoilErosionBC
    parameters::SoilErosionParameters
    variables::SoilErosionModelVariables
end

"Initialize soil erosion model"
function SoilErosionModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    vars = SoilErosionModelVariables(; n)
    params = SoilErosionParameters(dataset, config, indices)
    bc = SoilErosionBC(; n)
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
