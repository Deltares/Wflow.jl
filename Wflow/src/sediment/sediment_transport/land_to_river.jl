abstract type AbstractSedimentToRiverModel end

"Struct to store total sediment reaching the river model variables"
@with_kw struct SedimentToRiverVariables
    # Total sediment rate to the river [t dt-1]
    sediment_rate::Vector{Float64}
end

"Initialize total sediment reaching the river model variables"
function SedimentToRiverVariables(
    n::Int;
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n),
)
    return SedimentToRiverVariables(; sediment_rate)
end

"Struct to store total sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverBC
    # Deposition material rate [t dt-1]
    deposition::Vector{Float64}
end

"Initialize total sediment reaching the river model boundary conditions"
function SedimentToRiverBC(n::Int; deposition::Vector{Float64} = fill(MISSING_VALUE, n))
    return SedimentToRiverBC(; deposition)
end

"Struct to store total sediment reaching the river model"
@with_kw struct SedimentToRiverModel <: AbstractSedimentToRiverModel
    boundary_conditions::SedimentToRiverBC
    variables::SedimentToRiverVariables
end

"Initialize total sediment reaching the river model"
function SedimentToRiverModel(indices::Vector{CartesianIndex{2}})
    n = length(indices)
    vars = SedimentToRiverVariables(n)
    bc = SedimentToRiverBC(n)
    model = SedimentToRiverModel(; boundary_conditions = bc, variables = vars)
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
    clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt rate [t dt-1]
    silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand rate [t dt-1]
    sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates rate [t dt-1]
    sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates rate [t dt-1]
    lagg::Vector{Float64} = fill(MISSING_VALUE, n)
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
    boundary_conditions::SedimentToRiverDifferentiationBC
    variables::SedimentToRiverDifferentiationVariables
end

"Initialize differentiated sediment reaching the river model"
function SedimentToRiverDifferentiationModel(indices::Vector{CartesianIndex{2}})
    n = length(indices)
    vars = SedimentToRiverDifferentiationVariables(; n)
    bc = SedimentToRiverDifferentiationBC(; n)
    model =
        SedimentToRiverDifferentiationModel(; boundary_conditions = bc, variables = vars)
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
    (; sediment_rate, clay, silt, sand, sagg, lagg) = model.variables

    for (i, river) in enumerate(rivers)
        if river
            clay[i] = deposition_clay[i]
            silt[i] = deposition_silt[i]
            sand[i] = deposition_sand[i]
            sagg[i] = deposition_sagg[i]
            lagg[i] = deposition_lagg[i]
        else
            clay[i] = 0.0
            silt[i] = 0.0
            sand[i] = 0.0
            sagg[i] = 0.0
            lagg[i] = 0.0
        end
    end
    @. sediment_rate = clay + silt + sand + sagg + lagg
end
