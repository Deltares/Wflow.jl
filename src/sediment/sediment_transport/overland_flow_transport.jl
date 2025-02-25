abstract type AbstractSedimentLandTransportModel{T} end

"Struct to store total sediment flux in overland flow model variables"
@with_kw struct SedimentLandTransportVariables{T}
    # Total sediment rate [t dt-1]
    amount::Vector{T}
    # Total sediment deposition rate [t dt-1]
    deposition::Vector{T}
end

"Initialize total sediment flux in overland flow model variables"
function SedimentLandTransportVariables(
    n;
    amount::Vector{T} = fill(mv, n),
    deposition::Vector{T} = fill(mv, n),
) where {T}
    return SedimentLandTransportVariables{T}(; amount = amount, deposition = deposition)
end

"Struct to store total sediment flux in overland flow model boundary conditions"
@with_kw struct SedimentLandTransportBC{T}
    # Erosion rate material [t dt-1]
    erosion::Vector{T}
    # Transport capacity [t dt-1]
    transport_capacity::Vector{T}
end

"Initialize total sediment flux in overland flow model boundary conditions"
function SedimentLandTransportBC(
    n;
    erosion::Vector{T} = fill(mv, n),
    transport_capacity::Vector{T} = fill(mv, n),
) where {T}
    return SedimentLandTransportBC{T}(;
        erosion = erosion,
        transport_capacity = transport_capacity,
    )
end

"Struct to store total sediment flux in overland flow model"
@with_kw struct SedimentLandTransportModel{T} <: AbstractSedimentLandTransportModel{T}
    boundary_conditions::SedimentLandTransportBC{T}
    variables::SedimentLandTransportVariables{T}
end

"Initialize total sediment flux in overland flow model"
function SedimentLandTransportModel(indices)
    n = length(indices)
    vars = SedimentLandTransportVariables(n)
    bc = SedimentLandTransportBC(n)
    model = SedimentLandTransportModel(; boundary_conditions = bc, variables = vars)
    return model
end

"Update total sediment flux in overland flow model boundary conditions"
function update_boundary_conditions!(
    model::SedimentLandTransportModel,
    erosion_model::SoilErosionModel,
    transport_capacity_model::AbstractTransportCapacityModel,
)
    (; erosion, transport_capacity) = model.boundary_conditions
    (; amount) = erosion_model.variables
    @. erosion = amount

    (; amount) = transport_capacity_model.variables
    @. transport_capacity = amount
end

"Update total sediment flux in overland flow model for a single timestep"
function update!(model::SedimentLandTransportModel, network)
    (; erosion, transport_capacity) = model.boundary_conditions
    (; amount, deposition) = model.variables

    accucapacityflux!(amount, erosion, network, transport_capacity)
    deposition .= erosion
end

"Struct to store differentiated sediment flux in overland flow model variables"
@with_kw struct SedimentLandTransportDifferentiationVariables{T}
    # Total sediment rate [t dt-1]
    amount::Vector{T}
    # Deposition rate [t dt-1]
    deposition::Vector{T}
    # Clay rate [t dt-1]
    clay::Vector{T}
    # Deposition clay rate [t dt-1]
    deposition_clay::Vector{T}
    # Silt rate [t dt-1]
    silt::Vector{T}
    # Deposition silt rate [t dt-1]
    deposition_silt::Vector{T}
    # Sand rate [t dt-1]
    sand::Vector{T}
    # Deposition sand rate [t dt-1]
    deposition_sand::Vector{T}
    # Small aggregates rate [t dt-1]
    sagg::Vector{T}
    # Deposition rate small aggregates [t dt-1]
    deposition_sagg::Vector{T}
    # Large aggregates rate [t dt-1]
    lagg::Vector{T}
    # Deposition rate large aggregates [t dt-1]
    deposition_lagg::Vector{T}
end

"Initialize differentiated sediment flux in overland flow model variables"
function SedimentLandTransportDifferentiationVariables(
    n;
    amount::Vector{T} = fill(mv, n),
    deposition::Vector{T} = fill(mv, n),
    clay::Vector{T} = fill(mv, n),
    deposition_clay::Vector{T} = fill(mv, n),
    silt::Vector{T} = fill(mv, n),
    deposition_silt::Vector{T} = fill(mv, n),
    sand::Vector{T} = fill(mv, n),
    deposition_sand::Vector{T} = fill(mv, n),
    sagg::Vector{T} = fill(mv, n),
    deposition_sagg::Vector{T} = fill(mv, n),
    lagg::Vector{T} = fill(mv, n),
    deposition_lagg::Vector{T} = fill(mv, n),
) where {T}
    return SedimentLandTransportDifferentiationVariables{T}(;
        amount = amount,
        deposition = deposition,
        clay = clay,
        deposition_clay = deposition_clay,
        silt = silt,
        deposition_silt = deposition_silt,
        sand = sand,
        deposition_sand = deposition_sand,
        sagg = sagg,
        deposition_sagg = deposition_sagg,
        lagg = lagg,
        deposition_lagg = deposition_lagg,
    )
end

"Struct to store differentiated sediment flux in overland flow model boundary conditions"
@with_kw struct SedimentLandTransportDifferentiationBC{T}
    # Erosion rate clay [t dt-1]
    erosion_clay::Vector{T}
    # Erosion rate silt [t dt-1]
    erosion_silt::Vector{T}
    # Erosion rate sand [t dt-1]
    erosion_sand::Vector{T}
    # Erosion rate small aggregates [t dt-1]
    erosion_sagg::Vector{T}
    # Erosion large aggregates [t dt-1]
    erosion_lagg::Vector{T}
    # Transport capacity clay [t dt-1]
    transport_capacity_clay::Vector{T}
    # Transport capacity silt [t dt-1]
    transport_capacity_silt::Vector{T}
    # Transport capacity sand [t dt-1]
    transport_capacity_sand::Vector{T}
    # Transport capacity small aggregates [t dt-1]
    transport_capacity_sagg::Vector{T}
    # Transport capacity large aggregates [t dt-1]
    transport_capacity_lagg::Vector{T}
end

"Initialize differentiated sediment flux in overland flow model boundary conditions"
function SedimentLandTransportDifferentiationBC(
    n;
    erosion_clay::Vector{T} = fill(mv, n),
    erosion_silt::Vector{T} = fill(mv, n),
    erosion_sand::Vector{T} = fill(mv, n),
    erosion_sagg::Vector{T} = fill(mv, n),
    erosion_lagg::Vector{T} = fill(mv, n),
    transport_capacity_clay::Vector{T} = fill(mv, n),
    transport_capacity_silt::Vector{T} = fill(mv, n),
    transport_capacity_sand::Vector{T} = fill(mv, n),
    transport_capacity_sagg::Vector{T} = fill(mv, n),
    transport_capacity_lagg::Vector{T} = fill(mv, n),
) where {T}
    return SedimentLandTransportDifferentiationBC{T}(;
        erosion_clay = erosion_clay,
        erosion_silt = erosion_silt,
        erosion_sand = erosion_sand,
        erosion_sagg = erosion_sagg,
        erosion_lagg = erosion_lagg,
        transport_capacity_clay = transport_capacity_clay,
        transport_capacity_silt = transport_capacity_silt,
        transport_capacity_sand = transport_capacity_sand,
        transport_capacity_sagg = transport_capacity_sagg,
        transport_capacity_lagg = transport_capacity_lagg,
    )
end

"Struct to store differentiated sediment flux in overland flow model"
@with_kw struct SedimentLandTransportDifferentiationModel{T} <:
                AbstractSedimentLandTransportModel{T}
    boundary_conditions::SedimentLandTransportDifferentiationBC{T}
    variables::SedimentLandTransportDifferentiationVariables{T}
end

"Initialize differentiated sediment flux in overland flow model"
function SedimentLandTransportDifferentiationModel(indices)
    n = length(indices)
    vars = SedimentLandTransportDifferentiationVariables(n)
    bc = SedimentLandTransportDifferentiationBC(n)
    model = SedimentLandTransportDifferentiationModel(;
        boundary_conditions = bc,
        variables = vars,
    )
    return model
end

"Update differentiated sediment flux in overland flow model boundary conditions"
function update_boundary_conditions!(
    model::SedimentLandTransportDifferentiationModel,
    erosion_model::SoilErosionModel,
    transport_capacity_model::TransportCapacityYalinDifferentiationModel,
)
    (;
        erosion_clay,
        erosion_silt,
        erosion_sand,
        erosion_sagg,
        erosion_lagg,
        transport_capacity_clay,
        transport_capacity_silt,
        transport_capacity_sand,
        transport_capacity_sagg,
        transport_capacity_lagg,
    ) = model.boundary_conditions
    (; clay, silt, sand, sagg, lagg) = erosion_model.variables
    @. erosion_clay = clay
    @. erosion_silt = silt
    @. erosion_sand = sand
    @. erosion_sagg = sagg
    @. erosion_lagg = lagg

    (; clay, silt, sand, sagg, lagg) = transport_capacity_model.variables
    @. transport_capacity_clay = clay
    @. transport_capacity_silt = silt
    @. transport_capacity_sand = sand
    @. transport_capacity_sagg = sagg
    @. transport_capacity_lagg = lagg
end

"Update differentiated sediment flux in overland flow model for a single timestep"
function update!(model::SedimentLandTransportDifferentiationModel, network)
    (;
        erosion_clay,
        erosion_silt,
        erosion_sand,
        erosion_sagg,
        erosion_lagg,
        transport_capacity_clay,
        transport_capacity_silt,
        transport_capacity_sand,
        transport_capacity_sagg,
        transport_capacity_lagg,
    ) = model.boundary_conditions
    (;
        amount,
        deposition,
        clay,
        deposition_clay,
        silt,
        deposition_silt,
        sand,
        deposition_sand,
        sagg,
        deposition_sagg,
        lagg,
        deposition_lagg,
    ) = model.variables

    accucapacityflux!(clay, erosion_clay, network, transport_capacity_clay)
    deposition_clay .= erosion_clay
    accucapacityflux!(silt, erosion_silt, network, transport_capacity_silt)
    deposition_silt .= erosion_silt
    accucapacityflux!(sand, erosion_sand, network, transport_capacity_sand)
    deposition_sand .= erosion_sand
    accucapacityflux!(sagg, erosion_sagg, network, transport_capacity_sagg)
    deposition_sagg .= erosion_sagg
    accucapacityflux!(lagg, erosion_lagg, network, transport_capacity_lagg)
    deposition_lagg .= erosion_lagg
    amount .= clay .+ silt .+ sand .+ sagg .+ lagg
    deposition .=
        deposition_clay .+ deposition_silt .+ deposition_sand .+ deposition_sagg .+
        deposition_lagg
end