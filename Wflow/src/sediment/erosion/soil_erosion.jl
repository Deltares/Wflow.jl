abstract type AbstractSoilErosionModel end

"Struct for storing total soil erosion with differentiation model variables"
@with_kw struct SoilErosionModelVariables
    n_land_cells::Int
    # Total soil erosion rate [t dt-1]
    soil_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Total clay erosion rate [t dt-1]
    clay_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Total silt erosion rate [t dt-1]
    silt_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Total sand erosion rate [t dt-1]
    sand_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Total small aggregates erosion rate [t dt-1]
    sagg_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Total large aggregates erosion rate [t dt-1]
    lagg_erosion_rate::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
end

"Struct for storing soil erosion model boundary conditions"
@with_kw struct SoilErosionBC
    n_land_cells::Int
    # Rainfall erosion rate [t dt-1]
    rainfall_erosion::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Overland flow erosion rate [t dt-1]
    overland_flow_erosion::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
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
    land_indices_2d::Vector{CartesianIndex{2}},
)
    clay_fraction = ncread(
        dataset,
        config,
        "soil_clay__mass_fraction",
        SoilLossModel;
        sel = land_indices_2d,
    )
    silt_fraction = ncread(
        dataset,
        config,
        "soil_silt__mass_fraction",
        SoilLossModel;
        sel = land_indices_2d,
    )
    sand_fraction = ncread(
        dataset,
        config,
        "soil_sand__mass_fraction",
        SoilLossModel;
        sel = land_indices_2d,
    )
    sagg_fraction = ncread(
        dataset,
        config,
        "soil_small_aggregates__mass_fraction",
        SoilLossModel;
        sel = land_indices_2d,
    )
    lagg_fraction = ncread(
        dataset,
        config,
        "soil_large_aggregates__mass_fraction",
        SoilLossModel;
        sel = land_indices_2d,
    )
    # Check that soil fractions sum to 1
    soil_fractions =
        clay_fraction + silt_fraction + sand_fraction + sagg_fraction + lagg_fraction
    if !all(f -> isapprox(f, 1.0; rtol = 1e-3), soil_fractions)
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
    n_land_cells::Int
    boundary_conditions::SoilErosionBC = SoilErosionBC(; n_land_cells)
    parameters::SoilErosionParameters
    variables::SoilErosionModelVariables = SoilErosionModelVariables(; n_land_cells)
end

"Initialize soil erosion model"
function SoilErosionModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
)
    n_land_cells = length(land_indices_2d)
    parameters = SoilErosionParameters(dataset, config, land_indices_2d)
    soil_erosion_model = SoilErosionModel(; n_land_cells, parameters)
    return soil_erosion_model
end

"Update boundary conditions for soil erosion model"
function update_bc_soil_erosion_model!(
    soil_erosion_model::SoilErosionModel,
    rainfall_erosion_model::AbstractRainfallErosionModel,
    overland_flow_erosion_model::AbstractOverlandFlowErosionModel,
)
    re = rainfall_erosion_model.variables.soil_erosion_rate
    ole = overland_flow_erosion_model.variables.soil_erosion_rate
    (; rainfall_erosion, overland_flow_erosion) = soil_erosion_model.boundary_conditions
    @. rainfall_erosion = re
    @. overland_flow_erosion = ole
end

"Update soil erosion model for a single timestep"
function update_soil_erosion_model!(soil_erosion_model::SoilErosionModel)
    (; n_land_cells) = soil_erosion_model
    (; rainfall_erosion, overland_flow_erosion) = soil_erosion_model.boundary_conditions
    (; clay_fraction, silt_fraction, sand_fraction, sagg_fraction, lagg_fraction) =
        soil_erosion_model.parameters
    (;
        soil_erosion_rate,
        clay_erosion_rate,
        silt_erosion_rate,
        sand_erosion_rate,
        sagg_erosion_rate,
        lagg_erosion_rate,
    ) = soil_erosion_model.variables

    threaded_foreach(1:n_land_cells; basesize = 1000) do land_cell_idx
        soil_erosion_rate[land_cell_idx],
        clay_erosion_rate[land_cell_idx],
        silt_erosion_rate[land_cell_idx],
        sand_erosion_rate[land_cell_idx],
        sagg_erosion_rate[land_cell_idx],
        lagg_erosion_rate[land_cell_idx] = total_soil_erosion(
            rainfall_erosion[land_cell_idx],
            overland_flow_erosion[land_cell_idx],
            clay_fraction[land_cell_idx],
            silt_fraction[land_cell_idx],
            sand_fraction[land_cell_idx],
            sagg_fraction[land_cell_idx],
            lagg_fraction[land_cell_idx],
        )
    end
end
