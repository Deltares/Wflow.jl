abstract type AbstractSedimentToRiverModel{T} end

"Struct to store total sediment reaching the river model variables"
@with_kw struct SedimentToRiverVariables{T}
    # Total sediment rate to the river [t dt-1]
    amount::Vector{T}
end

"Initialize total sediment reaching the river model variables"
function SedimentToRiverVariables(n; amount::Vector{T} = fill(MISSING_VALUE, n)) where {T}
    return SedimentToRiverVariables{T}(; amount = amount)
end

"Struct to store total sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverBC{T}
    # Deposition material rate [t dt-1]
    deposition::Vector{T}
end

"Initialize total sediment reaching the river model boundary conditions"
function SedimentToRiverBC(n; deposition::Vector{T} = fill(MISSING_VALUE, n)) where {T}
    return SedimentToRiverBC{T}(; deposition = deposition)
end

"Struct to store total sediment reaching the river model"
@with_kw struct SedimentToRiverModel{T} <: AbstractSedimentToRiverModel{T}
    boundary_conditions::SedimentToRiverBC{T}
    variables::SedimentToRiverVariables{T}
end

"Initialize total sediment reaching the river model"
function SedimentToRiverModel(indices)
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
function update!(model::SedimentToRiverModel, rivers)
    (; deposition) = model.boundary_conditions
    (; amount) = model.variables

    zeros = fill(0.0, length(amount))
    amount .= ifelse.(rivers, deposition, zeros)
end

"Struct to store differentiated sediment reaching the river model variables"
@with_kw struct SedimentToRiverDifferentiationVariables{T}
    # Total sediment rate [t dt-1]
    amount::Vector{T}
    # Clay rate [t dt-1]
    clay::Vector{T}
    # Silt rate [t dt-1]
    silt::Vector{T}
    # Sand rate [t dt-1]
    sand::Vector{T}
    # Small aggregates rate [t dt-1]
    sagg::Vector{T}
    # Large aggregates rate [t dt-1]
    lagg::Vector{T}
end

"Initialize differentiated sediment reaching the river model variables"
function SedimentToRiverDifferentiationVariables(
    n;
    amount::Vector{T} = fill(MISSING_VALUE, n),
    clay::Vector{T} = fill(MISSING_VALUE, n),
    silt::Vector{T} = fill(MISSING_VALUE, n),
    sand::Vector{T} = fill(MISSING_VALUE, n),
    sagg::Vector{T} = fill(MISSING_VALUE, n),
    lagg::Vector{T} = fill(MISSING_VALUE, n),
) where {T}
    return SedimentToRiverDifferentiationVariables{T}(;
        amount = amount,
        clay = clay,
        silt = silt,
        sand = sand,
        sagg = sagg,
        lagg = lagg,
    )
end

"Struct to store differentiated sediment reaching the river model boundary conditions"
@with_kw struct SedimentToRiverDifferentiationBC{T}
    # Clay deposition rate [t dt-1]
    deposition_clay::Vector{T}
    # Silt deposition rate [t dt-1]
    deposition_silt::Vector{T}
    # Sand deposition rate [t dt-1]
    deposition_sand::Vector{T}
    # Small aggregates deposition rate [t dt-1]
    deposition_sagg::Vector{T}
    # Large aggregates deposition rate [t dt-1]
    deposition_lagg::Vector{T}
end

"Initialize differentiated sediment reaching the river model boundary conditions"
function SedimentToRiverDifferentiationBC(
    n;
    deposition_clay::Vector{T} = fill(MISSING_VALUE, n),
    deposition_silt::Vector{T} = fill(MISSING_VALUE, n),
    deposition_sand::Vector{T} = fill(MISSING_VALUE, n),
    deposition_sagg::Vector{T} = fill(MISSING_VALUE, n),
    deposition_lagg::Vector{T} = fill(MISSING_VALUE, n),
) where {T}
    return SedimentToRiverDifferentiationBC{T}(;
        deposition_clay = deposition_clay,
        deposition_silt = deposition_silt,
        deposition_sand = deposition_sand,
        deposition_sagg = deposition_sagg,
        deposition_lagg = deposition_lagg,
    )
end

"Struct to store differentiated sediment reaching the river model"
@with_kw struct SedimentToRiverDifferentiationModel{T} <: AbstractSedimentToRiverModel{T}
    boundary_conditions::SedimentToRiverDifferentiationBC{T}
    variables::SedimentToRiverDifferentiationVariables{T}
end

"Initialize differentiated sediment reaching the river model"
function SedimentToRiverDifferentiationModel(indices)
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
function update!(model::SedimentToRiverDifferentiationModel, rivers)
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
