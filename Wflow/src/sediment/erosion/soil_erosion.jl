abstract type AbstractSoilErosionModel end

"Struct for storing total soil erosion with differentiation model variables"
@with_kw struct SoilErosionModelVariables
    n_cells::Int
    # Total soil erosion rate [kg s⁻¹]
    soil_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Total clay erosion rate [kg s⁻¹]
    clay_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Total silt erosion rate [kg s⁻¹]
    silt_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Total sand erosion rate [kg s⁻¹]
    sand_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Total small aggregates erosion rate [kg s⁻¹]
    small_aggregates_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Total large aggregates erosion rate [kg s⁻¹]
    large_aggregates_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct for storing soil erosion model boundary conditions"
@with_kw struct SoilErosionBC
    n_cells::Int
    # Rainfall erosion rate [kg s⁻¹]
    rainfall_erosion::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Overland flow erosion rate [kg s⁻¹]
    overland_flow_erosion::Vector{Float64} = fill(MISSING_VALUE, n_cells)
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
    small_aggregates_fraction::Vector{Float64}
    # Soil content large aggregates [-]
    large_aggregates_fraction::Vector{Float64}
end

"Initialize soil erosion model parameters"
function SoilErosionParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    clay_fraction =
        ncread(dataset, config, "soil_clay__mass_fraction", SoilLossModel; sel = indices)
    silt_fraction =
        ncread(dataset, config, "soil_silt__mass_fraction", SoilLossModel; sel = indices)
    sand_fraction =
        ncread(dataset, config, "soil_sand__mass_fraction", SoilLossModel; sel = indices)
    small_aggregates_fraction = ncread(
        dataset,
        config,
        "soil_small_aggregates__mass_fraction",
        SoilLossModel;
        sel = indices,
    )
    large_aggregates_fraction = ncread(
        dataset,
        config,
        "soil_large_aggregates__mass_fraction",
        SoilLossModel;
        sel = indices,
    )
    # Check that soil fractions sum to 1
    soil_fractions =
        clay_fraction +
        silt_fraction +
        sand_fraction +
        small_aggregates_fraction +
        large_aggregates_fraction
    if !all(
        hydraulic_conductivity_scale_parameter ->
            isapprox(hydraulic_conductivity_scale_parameter, 1.0; rtol = 1e-3),
        soil_fractions,
    )
        error("Particle fractions in the soil must sum to 1.")
    end
    soil_parameters = SoilErosionParameters(;
        clay_fraction,
        silt_fraction,
        sand_fraction,
        small_aggregates_fraction,
        large_aggregates_fraction,
    )

    return soil_parameters
end

"Total soil erosion with differentiation model"
@with_kw struct SoilErosionModel <: AbstractSoilErosionModel
    n_cells::Int
    boundary_conditions::SoilErosionBC = SoilErosionBC(; n_cells)
    parameters::SoilErosionParameters
    variables::SoilErosionModelVariables = SoilErosionModelVariables(; n_cells)
end

"Initialize soil erosion model"
function SoilErosionModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n_cells = length(indices)
    parameters = SoilErosionParameters(dataset, config, indices)
    soil_erosion_model = SoilErosionModel(; n_cells, parameters)
    return soil_erosion_model
end

"Update boundary conditions for soil erosion model"
function update_bc_soil_erosion_model!(
    soil_erosion_model::SoilErosionModel,
    rainfall_erosion::AbstractRainfallErosionModel,
    overland_flow_erosion::OverlandFlowErosionAnswersModel,
)
    re = rainfall_erosion.variables.soil_erosion_rate
    ole = overland_flow_erosion.variables.soil_erosion_rate
    (; rainfall_erosion, overland_flow_erosion) = soil_erosion_model.boundary_conditions
    @. rainfall_erosion = re
    @. overland_flow_erosion = ole
end

"Update soil erosion model for a single timestep"
function update_soil_erosion_model!(soil_erosion_model::SoilErosionModel)
    (; rainfall_erosion, overland_flow_erosion) = soil_erosion_model.boundary_conditions
    (;
        clay_fraction,
        silt_fraction,
        sand_fraction,
        small_aggregates_fraction,
        large_aggregates_fraction,
    ) = soil_erosion_model.parameters
    (;
        soil_erosion_rate,
        clay_erosion_rate,
        silt_erosion_rate,
        sand_erosion_rate,
        small_aggregates_erosion_rate,
        large_aggregates_erosion_rate,
    ) = soil_erosion_model.variables

    n_cells = length(rainfall_erosion)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        soil_erosion_rate[cell_idx],
        clay_erosion_rate[cell_idx],
        silt_erosion_rate[cell_idx],
        sand_erosion_rate[cell_idx],
        small_aggregates_erosion_rate[cell_idx],
        large_aggregates_erosion_rate[cell_idx] = total_soil_erosion(
            rainfall_erosion[cell_idx],
            overland_flow_erosion[cell_idx],
            clay_fraction[cell_idx],
            silt_fraction[cell_idx],
            sand_fraction[cell_idx],
            small_aggregates_fraction[cell_idx],
            large_aggregates_fraction[cell_idx],
        )
    end
end
