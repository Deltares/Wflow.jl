abstract type AbstractSedimentToRiverModel end

"Struct to store total sediment reaching the river model variables"
@with_kw struct SedimentToRiverVariables
    # Total sediment rate to the river [t dt-1]
    amount::Vector{Float}
end

"Initialize total sediment reaching the river model variables"
function SedimentToRiverVariables(n::Int; amount::Vector{Float} = fill(MISSING_VALUE, n))
    return SedimentToRiverVariables(; amount = amount)
end

"Struct to store total sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverBC
    # Deposition material rate [t dt-1]
    deposition::Vector{Float}
end

"Initialize total sediment reaching the river model boundary conditions"
function SedimentToRiverBC(n::Int; deposition::Vector{Float} = fill(MISSING_VALUE, n))
    return SedimentToRiverBC(; deposition = deposition)
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
    (; amount) = model.variables

    zeros = fill(Float(0.0), length(amount))
    amount .= ifelse.(rivers, deposition, zeros)
end

"Struct to store differentiated sediment reaching the river model variables"
@with_kw struct SedimentToRiverDifferentiationVariables
    # Total sediment rate [t dt-1]
    amount::Vector{Float}
    # Clay rate [t dt-1]
    clay::Vector{Float}
    # Silt rate [t dt-1]
    silt::Vector{Float}
    # Sand rate [t dt-1]
    sand::Vector{Float}
    # Small aggregates rate [t dt-1]
    sagg::Vector{Float}
    # Large aggregates rate [t dt-1]
    lagg::Vector{Float}
end

"Initialize differentiated sediment reaching the river model variables"
function SedimentToRiverDifferentiationVariables(
    n::Int;
    amount::Vector{Float} = fill(MISSING_VALUE, n),
    clay::Vector{Float} = fill(MISSING_VALUE, n),
    silt::Vector{Float} = fill(MISSING_VALUE, n),
    sand::Vector{Float} = fill(MISSING_VALUE, n),
    sagg::Vector{Float} = fill(MISSING_VALUE, n),
    lagg::Vector{Float} = fill(MISSING_VALUE, n),
)
    return SedimentToRiverDifferentiationVariables(;
        amount = amount,
        clay = clay,
        silt = silt,
        sand = sand,
        sagg = sagg,
        lagg = lagg,
    )
end

"Struct to store differentiated sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverDifferentiationBC
    # Clay deposition rate [t dt-1]
    deposition_clay::Vector{Float}
    # Silt deposition rate [t dt-1]
    deposition_silt::Vector{Float}
    # Sand deposition rate [t dt-1]
    deposition_sand::Vector{Float}
    # Small aggregates deposition rate [t dt-1]
    deposition_sagg::Vector{Float}
    # Large aggregates deposition rate [t dt-1]
    deposition_lagg::Vector{Float}
end

"Initialize differentiated sediment reaching the river model boundary conditions"
function SedimentToRiverDifferentiationBC(
    n::Int;
    deposition_clay::Vector{Float} = fill(MISSING_VALUE, n),
    deposition_silt::Vector{Float} = fill(MISSING_VALUE, n),
    deposition_sand::Vector{Float} = fill(MISSING_VALUE, n),
    deposition_sagg::Vector{Float} = fill(MISSING_VALUE, n),
    deposition_lagg::Vector{Float} = fill(MISSING_VALUE, n),
)
    return SedimentToRiverDifferentiationBC(;
        deposition_clay = deposition_clay,
        deposition_silt = deposition_silt,
        deposition_sand = deposition_sand,
        deposition_sagg = deposition_sagg,
        deposition_lagg = deposition_lagg,
    )
end

"Struct to store differentiated sediment reaching the river model"
@with_kw struct SedimentToRiverDifferentiationModel <: AbstractSedimentToRiverModel
    boundary_conditions::SedimentToRiverDifferentiationBC
    variables::SedimentToRiverDifferentiationVariables
end

"Initialize differentiated sediment reaching the river model"
function SedimentToRiverDifferentiationModel(indices::Vector{CartesianIndex{2}})
    n = length(indices)
    vars = SedimentToRiverDifferentiationVariables(n)
    bc = SedimentToRiverDifferentiationBC(n)
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
    (; amount, clay, silt, sand, sagg, lagg) = model.variables

    zeros = fill(0.0, length(amount))
    clay .= ifelse.(rivers .> 0, deposition_clay, zeros)
    silt .= ifelse.(rivers .> 0, deposition_silt, zeros)
    sand .= ifelse.(rivers .> 0, deposition_sand, zeros)
    sagg .= ifelse.(rivers .> 0, deposition_sagg, zeros)
    lagg .= ifelse.(rivers .> 0, deposition_lagg, zeros)

    amount .= clay .+ silt .+ sand .+ sagg .+ lagg
end
