abstract type AbstractSedimentToRiverModel end

"Struct to store total sediment reaching the river model variables"
@with_kw struct SedimentToRiverVariables
    n_river_cells::Int
    # Total sediment rate to the river [t dt-1]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
end

"Struct to store total sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverBC
    n_river_cells::Int
    # Deposition material rate [t dt-1]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
end

"Struct to store total sediment reaching the river model"
@with_kw struct SedimentToRiverModel <: AbstractSedimentToRiverModel
    n_river_cells::Int
    boundary_conditions::SedimentToRiverBC = SedimentToRiverBC(; n_river_cells)
    variables::SedimentToRiverVariables = SedimentToRiverVariables(; n_river_cells)
end

"Initialize total sediment reaching the river model"
function SedimentToRiverModel(river_indices_2d::Vector{CartesianIndex{2}})
    n_river_cells = length(river_indices_2d)
    sediment_to_river_model = SedimentToRiverModel(; n_river_cells)
    return sediment_to_river_model
end

"Update total sediment reaching the river model boundary conditions"
function update_bc_sediment_to_river_model!(
    sediment_to_river_model::SedimentToRiverModel,
    sediment_transport_model::SedimentLandTransportModel,
)
    (; deposition) = sediment_to_river_model.boundary_conditions
    @. deposition = sediment_transport_model.variables.deposition
end

"Update total sediment reaching the river model for a single timestep"
function update_sediment_to_river_model!(
    sediment_to_river_model::SedimentToRiverModel,
    rivers::Vector{Bool},
)
    (; deposition) = sediment_to_river_model.boundary_conditions
    (; sediment_rate, n_river_cells) = sediment_to_river_model.variables

    map!(
        river_cell_idx -> rivers[river_cell_idx] ? deposition[river_cell_idx] : 0.0,
        sediment_rate,
        1:n_river_cells,
    )
end

"Struct to store differentiated sediment reaching the river model variables"
@with_kw struct SedimentToRiverDifferentiationVariables
    n_river_cells::Int
    # Total sediment rate [t dt-1]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Clay rate [t dt-1]
    clay_rate::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Silt rate [t dt-1]
    silt_rate::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Sand rate [t dt-1]
    sand_rate::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Small aggregates rate [t dt-1]
    sagg_rate::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Large aggregates rate [t dt-1]
    lagg_rate::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
end

"Struct to store differentiated sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverDifferentiationBC
    n_river_cells::Int
    # Clay deposition rate [t dt-1]
    deposition_clay::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Silt deposition rate [t dt-1]
    deposition_silt::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Sand deposition rate [t dt-1]
    deposition_sand::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Small aggregates deposition rate [t dt-1]
    deposition_sagg::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
    # Large aggregates deposition rate [t dt-1]
    deposition_lagg::Vector{Float64} = fill(MISSING_VALUE, n_river_cells)
end

"Struct to store differentiated sediment reaching the river model"
@with_kw struct SedimentToRiverDifferentiationModel <: AbstractSedimentToRiverModel
    n_river_cells::Int
    boundary_conditions::SedimentToRiverDifferentiationBC =
        SedimentToRiverDifferentiationBC(; n_river_cells)
    variables::SedimentToRiverDifferentiationVariables =
        SedimentToRiverDifferentiationVariables(; n_river_cells)
end

"Initialize differentiated sediment reaching the river model"
function SedimentToRiverDifferentiationModel(river_indices_2d::Vector{CartesianIndex{2}})
    n_river_cells = length(river_indices_2d)
    sediment_to_river_model = SedimentToRiverDifferentiationModel(; n_river_cells)
    return sediment_to_river_model
end

"Update differentiated sediment reaching the river model boundary conditions"
function update_bc_sediment_to_river_model!(
    sediment_to_river_model::SedimentToRiverDifferentiationModel,
    sediment_transport_model::SedimentLandTransportDifferentiationModel,
)
    (;
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_sagg,
        deposition_lagg,
    ) = sediment_to_river_model.boundary_conditions
    @. deposition_clay = sediment_transport_model.variables.deposition_clay
    @. deposition_silt = sediment_transport_model.variables.deposition_silt
    @. deposition_sand = sediment_transport_model.variables.deposition_sand
    @. deposition_sagg = sediment_transport_model.variables.deposition_sagg
    @. deposition_lagg = sediment_transport_model.variables.deposition_lagg
end

"Update differentiated sediment reaching the river model for a single timestep"
function update_sediment_to_river_model!(
    sediment_to_river_model::SedimentToRiverDifferentiationModel,
    rivers::Vector{Bool},
)
    (; n_river_cells) = sediment_to_river_model
    (;
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_sagg,
        deposition_lagg,
    ) = sediment_to_river_model.boundary_conditions
    (; sediment_rate, clay_rate, silt_rate, sand_rate, sagg_rate, lagg_rate) =
        sediment_to_river_model.variables

    for river_cell_idx in 1:n_river_cells
        if rivers[river_cell_idx]
            clay_rate[river_cell_idx] = deposition_clay[river_cell_idx]
            silt_rate[river_cell_idx] = deposition_silt[river_cell_idx]
            sand_rate[river_cell_idx] = deposition_sand[river_cell_idx]
            sagg_rate[river_cell_idx] = deposition_sagg[river_cell_idx]
            lagg_rate[river_cell_idx] = deposition_lagg[river_cell_idx]
        else
            clay_rate[river_cell_idx] = 0.0
            silt_rate[river_cell_idx] = 0.0
            sand_rate[river_cell_idx] = 0.0
            sagg_rate[river_cell_idx] = 0.0
            lagg_rate[river_cell_idx] = 0.0
        end
    end
    @. sediment_rate = clay_rate + silt_rate + sand_rate + sagg_rate + lagg_rate
end
