abstract type AbstractSedimentToRiverModel{T} end

## Total sediment transport in overland flow structs and functions
@get_units @grid_loc @with_kw struct SedimentToRiverVariables{T}
    # Total sediment reaching the river
    amount::Vector{T} | "t dt-1"
end

function SedimentToRiverVariables(n; amount::Vector{T} = fill(mv, n)) where {T}
    return SedimentToRiverVariables{T}(; amount = amount)
end

@get_units @grid_loc @with_kw struct SedimentToRiverBC{T}
    # Deposited material
    deposition::Vector{T} | "t dt-1"
end

function SedimentToRiverBC(n; deposition::Vector{T} = fill(mv, n)) where {T}
    return SedimentToRiverBC{T}(; deposition = deposition)
end

@with_kw struct SedimentToRiverModel{T} <: AbstractSedimentToRiverModel{T}
    boundary_conditions::SedimentToRiverBC{T}
    variables::SedimentToRiverVariables{T}
end

function SedimentToRiverModel(indices)
    n = length(indices)
    vars = SedimentToRiverVariables(n)
    bc = SedimentToRiverBC(n)
    model = SedimentToRiverModel(; boundary_conditions = bc, variables = vars)
    return model
end

function update_boundary_conditions!(
    model::SedimentToRiverModel,
    transport_model::SedimentLandTransportModel,
)
    (; deposition) = model.boundary_conditions
    @. deposition = transport_model.variables.deposition
end

function update!(model::SedimentToRiverModel, rivers)
    (; deposition) = model.boundary_conditions
    (; amount) = model.variables

    zeros = fill(0.0, length(amount))
    amount .= ifelse.(rivers, deposition, zeros)
end

## Different particles reaching the river structs and functions
@get_units @grid_loc @with_kw struct SedimentToRiverDifferentiationVariables{T}
    # Total sediment flux
    amount::Vector{T} | "t dt-1"
    # Clay flux
    clay::Vector{T} | "t dt-1"
    # Silt
    silt::Vector{T} | "t dt-1"
    # Sand flux
    sand::Vector{T} | "t dt-1"
    # Small aggregates flux
    sagg::Vector{T} | "t dt-1"
    # Large aggregates flux
    lagg::Vector{T} | "t dt-1"
end

function SedimentToRiverDifferentiationVariables(
    n;
    amount::Vector{T} = fill(mv, n),
    clay::Vector{T} = fill(mv, n),
    silt::Vector{T} = fill(mv, n),
    sand::Vector{T} = fill(mv, n),
    sagg::Vector{T} = fill(mv, n),
    lagg::Vector{T} = fill(mv, n),
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

@get_units @grid_loc @with_kw struct SedimentToRiverDifferentiationBC{T}
    # Deposited clay
    deposition_clay::Vector{T} | "t dt-1"
    # Deposited silt
    deposition_silt::Vector{T} | "t dt-1"
    # Deposited sand
    deposition_sand::Vector{T} | "t dt-1"
    # Deposited small aggregates
    deposition_sagg::Vector{T} | "t dt-1"
    # Deposited large aggregates
    deposition_lagg::Vector{T} | "t dt-1"
end

function SedimentToRiverDifferentiationBC(
    n;
    deposition_clay::Vector{T} = fill(mv, n),
    deposition_silt::Vector{T} = fill(mv, n),
    deposition_sand::Vector{T} = fill(mv, n),
    deposition_sagg::Vector{T} = fill(mv, n),
    deposition_lagg::Vector{T} = fill(mv, n),
) where {T}
    return SedimentToRiverDifferentiationBC{T}(;
        deposition_clay = deposition_clay,
        deposition_silt = deposition_silt,
        deposition_sand = deposition_sand,
        deposition_sagg = deposition_sagg,
        deposition_lagg = deposition_lagg,
    )
end

@with_kw struct SedimentToRiverDifferentiationModel{T} <: AbstractSedimentToRiverModel{T}
    boundary_conditions::SedimentToRiverDifferentiationBC{T}
    variables::SedimentToRiverDifferentiationVariables{T}
end

function SedimentToRiverDifferentiationModel(indices)
    n = length(indices)
    vars = SedimentToRiverDifferentiationVariables(n)
    bc = SedimentToRiverDifferentiationBC(n)
    model =
        SedimentToRiverDifferentiationModel(; boundary_conditions = bc, variables = vars)
    return model
end

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
