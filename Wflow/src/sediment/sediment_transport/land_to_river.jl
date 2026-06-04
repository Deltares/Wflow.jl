abstract type AbstractSedimentToRiverModel end

"Struct to store total sediment reaching the river model variables"
@with_kw struct SedimentToRiverVariables
    n::Int
    # Total sediment rate to the river [kg s⁻¹]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store total sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverBC
    n::Int
    # Deposition material rate [kg s⁻¹]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store total sediment reaching the river model"
@with_kw struct SedimentToRiverModel <: AbstractSedimentToRiverModel
    n::Int
    boundary_conditions::SedimentToRiverBC = SedimentToRiverBC(; n)
    variables::SedimentToRiverVariables = SedimentToRiverVariables(; n)
end

"Initialize total sediment reaching the river model"
function SedimentToRiverModel(indices::Vector{CartesianIndex{2}})
    n = length(indices)
    sediment_to_river_model = SedimentToRiverModel(; n)
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
    dt::Float64,
)
    (; deposition) = sediment_to_river_model.boundary_conditions
    (; sediment_rate) = sediment_to_river_model.variables

    for (cell_idx, river) in enumerate(rivers)
        sediment_rate[cell_idx] = river ? deposition[cell_idx] : 0.0
    end
end

"Struct to store differentiated sediment reaching the river model variables"
@with_kw struct SedimentToRiverDifferentiationVariables
    n::Int
    # Total sediment rate [kg s⁻¹]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Clay rate [kg s⁻¹]
    clay_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt rate [kg s⁻¹]
    silt_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand rate [kg s⁻¹]
    sand_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates rate [kg s⁻¹]
    small_aggregates_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates rate [kg s⁻¹]
    large_aggregates_rate::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store differentiated sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverDifferentiationBC
    n::Int
    # Clay deposition rate [kg s⁻¹]
    deposition_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt deposition rate [kg s⁻¹]
    deposition_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand deposition rate [kg s⁻¹]
    deposition_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates deposition rate [kg s⁻¹]
    deposition_small_aggregates::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates deposition rate [kg s⁻¹]
    deposition_large_aggregates::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store differentiated sediment reaching the river model"
@with_kw struct SedimentToRiverDifferentiationModel <: AbstractSedimentToRiverModel
    n::Int
    boundary_conditions::SedimentToRiverDifferentiationBC =
        SedimentToRiverDifferentiationBC(; n)
    variables::SedimentToRiverDifferentiationVariables =
        SedimentToRiverDifferentiationVariables(; n)
end

"Initialize differentiated sediment reaching the river model"
function SedimentToRiverDifferentiationModel(indices::Vector{CartesianIndex{2}})
    n = length(indices)
    sediment_to_river_model = SedimentToRiverDifferentiationModel(; n)
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
        deposition_small_aggregates,
        deposition_large_aggregates,
    ) = sediment_to_river_model.boundary_conditions
    @. deposition_clay = sediment_transport_model.variables.deposition_clay
    @. deposition_silt = sediment_transport_model.variables.deposition_silt
    @. deposition_sand = sediment_transport_model.variables.deposition_sand
    @. deposition_small_aggregates =
        sediment_transport_model.variables.deposition_small_aggregates
    @. deposition_large_aggregates =
        sediment_transport_model.variables.deposition_large_aggregates
end

"Update differentiated sediment reaching the river model for a single timestep"
function update_sediment_to_river_model!(
    sediment_to_river_model::SedimentToRiverDifferentiationModel,
    rivers::Vector{Bool},
    dt::Float64,
)
    (;
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_small_aggregates,
        deposition_large_aggregates,
    ) = sediment_to_river_model.boundary_conditions
    (;
        sediment_rate,
        clay_rate,
        silt_rate,
        sand_rate,
        small_aggregates_rate,
        large_aggregates_rate,
    ) = sediment_to_river_model.variables

    for (cell_idx, river) in enumerate(rivers)
        if river
            clay_rate[cell_idx] = deposition_clay[cell_idx]
            silt_rate[cell_idx] = deposition_silt[cell_idx]
            sand_rate[cell_idx] = deposition_sand[cell_idx]
            small_aggregates_rate[cell_idx] = deposition_small_aggregates[cell_idx]
            large_aggregates_rate[cell_idx] = deposition_large_aggregates[cell_idx]
        else
            clay_rate[cell_idx] = 0.0
            silt_rate[cell_idx] = 0.0
            sand_rate[cell_idx] = 0.0
            small_aggregates_rate[cell_idx] = 0.0
            large_aggregates_rate[cell_idx] = 0.0
        end
    end
    @. sediment_rate =
        clay_rate + silt_rate + sand_rate + small_aggregates_rate + large_aggregates_rate
end
