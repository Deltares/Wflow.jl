abstract type AbstractSedimentToRiverModel end

## Total sediment transport in overland flow structs and functions
@get_units @with_kw struct SedimentToRiverVars{T}
    # Total sediment reaching the river
    amount::Vector{T} | "t dt-1"
end

function sediment_to_river_vars(n)
    vars = SedimentToRiverVars(; amount = fill(mv, n))
    return vars
end

@get_units @with_kw struct SedimentToRiverBC{T}
    # Deposited material
    deposition::Vector{T} | "t dt-1"
end

function sediment_to_river_bc(n)
    bc = SedimentToRiverBC(; deposition = fill(mv, n))
    return bc
end

@get_units @with_kw struct SedimentToRiverModel{T} <: AbstractSedimentToRiverModel
    boundary_conditions::SedimentToRiverBC{T} | "-"
    variables::SedimentToRiverVars{T} | "-"
end

function initialize_sediment_to_river_model(inds)
    n = length(inds)
    vars = sediment_to_river_vars(n)
    bc = sediment_to_river_bc(n)
    model = SedimentToRiverModel(; boundary_conditions = bc, variables = vars)
    return model
end

function update_bc(model::SedimentToRiverModel, transport_model::SedimentLandTransportModel)
    (; deposition) = model.boundary_conditions
    (; amount) = transport_model.variables
    @. deposition = amount
end

function update!(model::SedimentToRiverModel, rivers)
    (; deposition) = model.boundary_conditions
    (; amount) = model.variables

    zeros = fill(0.0, length(amount))
    amount .= ifelse.(rivers, deposition, zeros)
end

## Different particles reaching the river structs and functions
@get_units @with_kw struct SedimentToRiverDifferentiationVars{T}
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

function sediment_to_river_differentiation_vars(n)
    vars = SedimentToRiverDifferentiationVars(;
        amount = fill(mv, n),
        clay = fill(mv, n),
        silt = fill(mv, n),
        sand = fill(mv, n),
        sagg = fill(mv, n),
        lagg = fill(mv, n),
    )
    return vars
end

@get_units @with_kw struct SedimentToRiverDifferentiationBC{T}
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

function sediment_to_river_differentiation_bc(n)
    bc = SedimentToRiverDifferentiationBC(;
        deposition_clay = fill(mv, n),
        deposition_silt = fill(mv, n),
        deposition_sand = fill(mv, n),
        deposition_sagg = fill(mv, n),
        deposition_lagg = fill(mv, n),
    )
    return bc
end

@get_units @with_kw struct SedimentToRiverDifferentiationModel{T} <:
                           AbstractSedimentToRiverModel
    boundary_conditions::SedimentToRiverDifferentiationBC{T} | "-"
    variables::SedimentToRiverDifferentiationVars{T} | "-"
end

function initialize_sediment_to_river_differentiation_model(inds)
    n = length(inds)
    vars = sediment_to_river_differentiation_vars(n)
    bc = sediment_to_river_differentiation_bc(n)
    model =
        SedimentToRiverDifferentiationModel(; boundary_conditions = bc, variables = vars)
    return model
end

function update_bc(
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
    (; clay, silt, sand, sagg, lagg) = transport_model.variables
    @. deposition_clay = clay
    @. deposition_silt = silt
    @. deposition_sand = sand
    @. deposition_sagg = sagg
    @. deposition_lagg = lagg
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
    clay .= ifelse.(rivers, deposition_clay, zeros)
    silt .= ifelse.(rivers, deposition_silt, zeros)
    sand .= ifelse.(rivers, deposition_sand, zeros)
    sagg .= ifelse.(rivers, deposition_sagg, zeros)
    lagg .= ifelse.(rivers, deposition_lagg, zeros)

    amount .= clay .+ silt .+ sand .+ sagg .+ lagg
end
