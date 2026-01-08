abstract type AbstractSedimentToRiverModel end

"Struct to store total sediment reaching the river model variables"
@with_kw struct SedimentToRiverVariables
    n::Int
    # Total sediment rate to the river [t dt-1]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store total sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverBC
    n::Int
    # Deposition material rate [t dt-1]
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
    model = SedimentToRiverModel(; n)
    return model
end

"Update total sediment reaching the river model boundary conditions"
function update_boundary_conditions!(
    model::SedimentToRiverModel,
    transport_model::SedimentLandTransportModel,
)
    (; deposition) = model.boundary_conditions
    @. deposition = transport_model.variables.deposition
end

"Update total sediment reaching the river model for a single timestep"
function update!(model::SedimentToRiverModel, rivers::Vector{Bool})
    (; deposition) = model.boundary_conditions
    (; sediment_rate) = model.variables

    zeros = fill(0.0, length(sediment_rate))
    sediment_rate .= ifelse.(rivers, deposition, zeros)
end

"Struct to store differentiated sediment reaching the river model variables"
@with_kw struct SedimentToRiverDifferentiationVariables
    n::Int
    # Total sediment rate [t dt-1]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Clay rate [t dt-1]
    clay_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt rate [t dt-1]
    silt_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand rate [t dt-1]
    sand_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates rate [t dt-1]
    sagg_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates rate [t dt-1]
    lagg_rate::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store differentiated sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverDifferentiationBC
    n::Int
    # Clay deposition rate [t dt-1]
    deposition_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt deposition rate [t dt-1]
    deposition_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand deposition rate [t dt-1]
    deposition_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates deposition rate [t dt-1]
    deposition_sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates deposition rate [t dt-1]
    deposition_lagg::Vector{Float64} = fill(MISSING_VALUE, n)
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
    model = SedimentToRiverDifferentiationModel(; n)
    return model
end

"Update differentiated sediment reaching the river model boundary conditions"
function update_boundary_conditions!(
    model::SedimentToRiverDifferentiationModel,
    transport_model::SedimentLandTransportDifferentiationModel,
)
    (;
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_sagg,
        deposition_lagg,
    ) = model.boundary_conditions
    @. deposition_clay = transport_model.variables.deposition_clay
    @. deposition_silt = transport_model.variables.deposition_silt
    @. deposition_sand = transport_model.variables.deposition_sand
    @. deposition_sagg = transport_model.variables.deposition_sagg
    @. deposition_lagg = transport_model.variables.deposition_lagg
end

"Update differentiated sediment reaching the river model for a single timestep"
function update!(model::SedimentToRiverDifferentiationModel, rivers::Vector{Bool})
    (;
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_sagg,
        deposition_lagg,
    ) = model.boundary_conditions
    (; sediment_rate, clay_rate, silt_rate, sand_rate, sagg_rate, lagg_rate) =
        model.variables

    for (i, river) in enumerate(rivers)
        if river
            clay_rate[i] = deposition_clay[i]
            silt_rate[i] = deposition_silt[i]
            sand_rate[i] = deposition_sand[i]
            sagg_rate[i] = deposition_sagg[i]
            lagg_rate[i] = deposition_lagg[i]
        else
            clay_rate[i] = 0.0
            silt_rate[i] = 0.0
            sand_rate[i] = 0.0
            sagg_rate[i] = 0.0
            lagg_rate[i] = 0.0
        end
    end
    @. sediment_rate = clay_rate + silt_rate + sand_rate + sagg_rate + lagg_rate
end
