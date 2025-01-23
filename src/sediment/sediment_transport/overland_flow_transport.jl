abstract type AbstractSedimentLandTransportModel{T} end

"Struct to store total sediment flux in overland flow model variables"
@get_units @grid_loc @with_kw struct SedimentLandTransportVariables{T}
    # Total sediment flux
    amount::Vector{T} | "t dt-1"
    deposition::Vector{T} | "t dt-1"
end

"Initialize total sediment flux in overland flow model variables"
function SedimentLandTransportVariables(
    n;
    amount::Vector{T} = fill(MISSING_VALUE, n),
    deposition::Vector{T} = fill(MISSING_VALUE, n),
) where {T}
    return SedimentLandTransportVariables{T}(; amount = amount, deposition = deposition)
end

"Struct to store total sediment flux in overland flow model boundary conditions"
@get_units @grid_loc @with_kw struct SedimentLandTransportBC{T}
    # Eroded material
    erosion::Vector{T} | "t dt-1"
    # Transport capacity
    transport_capacity::Vector{T} | "t dt-1"
end

"Initialize total sediment flux in overland flow model boundary conditions"
function SedimentLandTransportBC(
    n;
    erosion::Vector{T} = fill(MISSING_VALUE, n),
    transport_capacity::Vector{T} = fill(MISSING_VALUE, n),
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
@get_units @grid_loc @with_kw struct SedimentLandTransportDifferentiationVariables{T}
    # Total sediment flux
    amount::Vector{T} | "t dt-1"
    # Deposition
    deposition::Vector{T} | "t dt-1"
    # Clay flux
    clay::Vector{T} | "t dt-1"
    # Deposition clay
    deposition_clay::Vector{T} | "t dt-1"
    # Silt
    silt::Vector{T} | "t dt-1"
    # Deposition silt
    deposition_silt::Vector{T} | "t dt-1"
    # Sand flux
    sand::Vector{T} | "t dt-1"
    # Deposition sand
    deposition_sand::Vector{T} | "t dt-1"
    # Small aggregates flux
    sagg::Vector{T} | "t dt-1"
    # Deposition small aggregates
    deposition_sagg::Vector{T} | "t dt-1"
    # Large aggregates flux
    lagg::Vector{T} | "t dt-1"
    # Deposition large aggregates
    deposition_lagg::Vector{T} | "t dt-1"
end

"Initialize differentiated sediment flux in overland flow model variables"
function SedimentLandTransportDifferentiationVariables(
    n;
    amount::Vector{T} = fill(MISSING_VALUE, n),
    deposition::Vector{T} = fill(MISSING_VALUE, n),
    clay::Vector{T} = fill(MISSING_VALUE, n),
    deposition_clay::Vector{T} = fill(MISSING_VALUE, n),
    silt::Vector{T} = fill(MISSING_VALUE, n),
    deposition_silt::Vector{T} = fill(MISSING_VALUE, n),
    sand::Vector{T} = fill(MISSING_VALUE, n),
    deposition_sand::Vector{T} = fill(MISSING_VALUE, n),
    sagg::Vector{T} = fill(MISSING_VALUE, n),
    deposition_sagg::Vector{T} = fill(MISSING_VALUE, n),
    lagg::Vector{T} = fill(MISSING_VALUE, n),
    deposition_lagg::Vector{T} = fill(MISSING_VALUE, n),
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
@get_units @grid_loc @with_kw struct SedimentLandTransportDifferentiationBC{T}
    # Eroded clay
    erosion_clay::Vector{T} | "t dt-1"
    # Eroded silt
    erosion_silt::Vector{T} | "t dt-1"
    # Eroded sand
    erosion_sand::Vector{T} | "t dt-1"
    # Eroded small aggregates
    erosion_sagg::Vector{T} | "t dt-1"
    # Eroded large aggregates
    erosion_lagg::Vector{T} | "t dt-1"
    # Transport capacity clay
    transport_capacity_clay::Vector{T} | "t dt-1"
    # Transport capacity silt
    transport_capacity_silt::Vector{T} | "t dt-1"
    # Transport capacity sand
    transport_capacity_sand::Vector{T} | "t dt-1"
    # Transport capacity small aggregates
    transport_capacity_sagg::Vector{T} | "t dt-1"
    # Transport capacity large aggregates
    transport_capacity_lagg::Vector{T} | "t dt-1"
end

"Initialize differentiated sediment flux in overland flow model boundary conditions"
function SedimentLandTransportDifferentiationBC(
    n;
    erosion_clay::Vector{T} = fill(MISSING_VALUE, n),
    erosion_silt::Vector{T} = fill(MISSING_VALUE, n),
    erosion_sand::Vector{T} = fill(MISSING_VALUE, n),
    erosion_sagg::Vector{T} = fill(MISSING_VALUE, n),
    erosion_lagg::Vector{T} = fill(MISSING_VALUE, n),
    transport_capacity_clay::Vector{T} = fill(MISSING_VALUE, n),
    transport_capacity_silt::Vector{T} = fill(MISSING_VALUE, n),
    transport_capacity_sand::Vector{T} = fill(MISSING_VALUE, n),
    transport_capacity_sagg::Vector{T} = fill(MISSING_VALUE, n),
    transport_capacity_lagg::Vector{T} = fill(MISSING_VALUE, n),
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