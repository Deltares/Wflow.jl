abstract type AbstractSedimentToRiverModel end

"Struct to store total sediment reaching the river model variables"
@with_kw struct SedimentToRiverVariables
    n_river::Int
    # Total sediment rate to the river [kg s⁻¹]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n_river)
end

"Struct to store total sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverBC
    n_river::Int
    # Deposition material rate [kg s⁻¹]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n_river)
end

"Struct to store total sediment reaching the river model"
@with_kw struct SedimentToRiverModel <: AbstractSedimentToRiverModel
    n_river::Int
    boundary_conditions::SedimentToRiverBC = SedimentToRiverBC(; n_river)
    variables::SedimentToRiverVariables = SedimentToRiverVariables(; n_river)
end

"Initialize total sediment reaching the river model"
function SedimentToRiverModel(indices::Vector{CartesianIndex{2}})
    n_river = length(indices)
    sediment_to_river_model = SedimentToRiverModel(; n_river)
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

    for (idx, river) in enumerate(rivers)
        sediment_rate[idx] = river ? deposition[idx] : 0.0
    end
end

"Struct to store differentiated sediment reaching the river model variables"
@with_kw struct SedimentToRiverDifferentiationVariables
    n_river::Int
    # Total sediment rate [kg s⁻¹]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n_river)
    # Clay rate [kg s⁻¹]
    clay_rate::Vector{Float64} = fill(MISSING_VALUE, n_river)
    # Silt rate [kg s⁻¹]
    silt_rate::Vector{Float64} = fill(MISSING_VALUE, n_river)
    # Sand rate [kg s⁻¹]
    sand_rate::Vector{Float64} = fill(MISSING_VALUE, n_river)
    # Small aggregates rate [kg s⁻¹]
    small_aggregates_rate::Vector{Float64} = fill(MISSING_VALUE, n_river)
    # Large aggregates rate [kg s⁻¹]
    large_aggregates_rate::Vector{Float64} = fill(MISSING_VALUE, n_river)
end

"Struct to store differentiated sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverDifferentiationBC
    n_river::Int
    # Clay deposition rate [kg s⁻¹]
    deposition_clay::Vector{Float64} = fill(MISSING_VALUE, n_river)
    # Silt deposition rate [kg s⁻¹]
    deposition_silt::Vector{Float64} = fill(MISSING_VALUE, n_river)
    # Sand deposition rate [kg s⁻¹]
    deposition_sand::Vector{Float64} = fill(MISSING_VALUE, n_river)
    # Small aggregates deposition rate [kg s⁻¹]
    deposition_small_aggregates::Vector{Float64} = fill(MISSING_VALUE, n_river)
    # Large aggregates deposition rate [kg s⁻¹]
    deposition_large_aggregates::Vector{Float64} = fill(MISSING_VALUE, n_river)
end

"Struct to store differentiated sediment reaching the river model"
@with_kw struct SedimentToRiverDifferentiationModel <: AbstractSedimentToRiverModel
    n_river::Int
    boundary_conditions::SedimentToRiverDifferentiationBC =
        SedimentToRiverDifferentiationBC(; n_river)
    variables::SedimentToRiverDifferentiationVariables =
        SedimentToRiverDifferentiationVariables(; n_river)
end

"Initialize differentiated sediment reaching the river model"
function SedimentToRiverDifferentiationModel(indices::Vector{CartesianIndex{2}})
    n_river = length(indices)
    sediment_to_river_model = SedimentToRiverDifferentiationModel(; n_river)
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

    for (idx, river) in enumerate(rivers)
        if river
            clay_rate[idx] = deposition_clay[idx]
            silt_rate[idx] = deposition_silt[idx]
            sand_rate[idx] = deposition_sand[idx]
            small_aggregates_rate[idx] = deposition_small_aggregates[idx]
            large_aggregates_rate[idx] = deposition_large_aggregates[idx]
        else
            clay_rate[idx] = 0.0
            silt_rate[idx] = 0.0
            sand_rate[idx] = 0.0
            small_aggregates_rate[idx] = 0.0
            large_aggregates_rate[idx] = 0.0
        end
    end
    @. sediment_rate =
        clay_rate + silt_rate + sand_rate + small_aggregates_rate + large_aggregates_rate
end
