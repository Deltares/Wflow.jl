abstract type AbstractSedimentLandTransportModel end

"Struct to store total sediment flux in overland flow model variables"
@with_kw struct SedimentLandTransportVariables
    # Total sediment rate [t dt-1]
    amount::Vector{Float}
    # Total sediment deposition rate [t dt-1]
    deposition::Vector{Float}
end

"Initialize total sediment flux in overland flow model variables"
function SedimentLandTransportVariables(
    n::Int;
    amount::Vector{Float} = fill(MISSING_VALUE, n),
    deposition::Vector{Float} = fill(MISSING_VALUE, n),
)
    return SedimentLandTransportVariables(; amount = amount, deposition = deposition)
end

"Struct to store total sediment flux in overland flow model boundary conditions"
@with_kw struct SedimentLandTransportBC
    # Erosion rate material [t dt-1]
    erosion::Vector{Float}
    # Transport capacity [t dt-1]
    transport_capacity::Vector{Float}
end

"Initialize total sediment flux in overland flow model boundary conditions"
function SedimentLandTransportBC(
    n::Int;
    erosion::Vector{Float} = fill(MISSING_VALUE, n),
    transport_capacity::Vector{Float} = fill(MISSING_VALUE, n),
)
    return SedimentLandTransportBC(;
        erosion = erosion,
        transport_capacity = transport_capacity,
    )
end

"Struct to store total sediment flux in overland flow model"
@with_kw struct SedimentLandTransportModel <: AbstractSedimentLandTransportModel
    boundary_conditions::SedimentLandTransportBC
    variables::SedimentLandTransportVariables
end

"Initialize total sediment flux in overland flow model"
function SedimentLandTransportModel(indices::Vector{CartesianIndex{2}})
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
function update!(model::SedimentLandTransportModel, network::NetworkLand)
    (; erosion, transport_capacity) = model.boundary_conditions
    (; amount, deposition) = model.variables

    accucapacityflux!(amount, erosion, network, transport_capacity)
    deposition .= erosion
end

"Struct to store differentiated sediment flux in overland flow model variables"
@with_kw struct SedimentLandTransportDifferentiationVariables
    # Total sediment rate [t dt-1]
    amount::Vector{Float}
    # Deposition rate [t dt-1]
    deposition::Vector{Float}
    # Clay rate [t dt-1]
    clay::Vector{Float}
    # Deposition clay rate [t dt-1]
    deposition_clay::Vector{Float}
    # Silt rate [t dt-1]
    silt::Vector{Float}
    # Deposition silt rate [t dt-1]
    deposition_silt::Vector{Float}
    # Sand rate [t dt-1]
    sand::Vector{Float}
    # Deposition sand rate [t dt-1]
    deposition_sand::Vector{Float}
    # Small aggregates rate [t dt-1]
    sagg::Vector{Float}
    # Deposition rate small aggregates [t dt-1]
    deposition_sagg::Vector{Float}
    # Large aggregates rate [t dt-1]
    lagg::Vector{Float}
    # Deposition rate large aggregates [t dt-1]
    deposition_lagg::Vector{Float}
end

"Initialize differentiated sediment flux in overland flow model variables"
function SedimentLandTransportDifferentiationVariables(
    n::Int;
    amount::Vector{Float} = fill(MISSING_VALUE, n),
    deposition::Vector{Float} = fill(MISSING_VALUE, n),
    clay::Vector{Float} = fill(MISSING_VALUE, n),
    deposition_clay::Vector{Float} = fill(MISSING_VALUE, n),
    silt::Vector{Float} = fill(MISSING_VALUE, n),
    deposition_silt::Vector{Float} = fill(MISSING_VALUE, n),
    sand::Vector{Float} = fill(MISSING_VALUE, n),
    deposition_sand::Vector{Float} = fill(MISSING_VALUE, n),
    sagg::Vector{Float} = fill(MISSING_VALUE, n),
    deposition_sagg::Vector{Float} = fill(MISSING_VALUE, n),
    lagg::Vector{Float} = fill(MISSING_VALUE, n),
    deposition_lagg::Vector{Float} = fill(MISSING_VALUE, n),
)
    return SedimentLandTransportDifferentiationVariables(;
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
@with_kw struct SedimentLandTransportDifferentiationBC
    # Erosion rate clay [t dt-1]
    erosion_clay::Vector{Float}
    # Erosion rate silt [t dt-1]
    erosion_silt::Vector{Float}
    # Erosion rate sand [t dt-1]
    erosion_sand::Vector{Float}
    # Erosion rate small aggregates [t dt-1]
    erosion_sagg::Vector{Float}
    # Erosion large aggregates [t dt-1]
    erosion_lagg::Vector{Float}
    # Transport capacity clay [t dt-1]
    transport_capacity_clay::Vector{Float}
    # Transport capacity silt [t dt-1]
    transport_capacity_silt::Vector{Float}
    # Transport capacity sand [t dt-1]
    transport_capacity_sand::Vector{Float}
    # Transport capacity small aggregates [t dt-1]
    transport_capacity_sagg::Vector{Float}
    # Transport capacity large aggregates [t dt-1]
    transport_capacity_lagg::Vector{Float}
end

"Initialize differentiated sediment flux in overland flow model boundary conditions"
function SedimentLandTransportDifferentiationBC(
    n::Int;
    erosion_clay::Vector{Float} = fill(MISSING_VALUE, n),
    erosion_silt::Vector{Float} = fill(MISSING_VALUE, n),
    erosion_sand::Vector{Float} = fill(MISSING_VALUE, n),
    erosion_sagg::Vector{Float} = fill(MISSING_VALUE, n),
    erosion_lagg::Vector{Float} = fill(MISSING_VALUE, n),
    transport_capacity_clay::Vector{Float} = fill(MISSING_VALUE, n),
    transport_capacity_silt::Vector{Float} = fill(MISSING_VALUE, n),
    transport_capacity_sand::Vector{Float} = fill(MISSING_VALUE, n),
    transport_capacity_sagg::Vector{Float} = fill(MISSING_VALUE, n),
    transport_capacity_lagg::Vector{Float} = fill(MISSING_VALUE, n),
)
    return SedimentLandTransportDifferentiationBC(;
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
@with_kw struct SedimentLandTransportDifferentiationModel <:
                AbstractSedimentLandTransportModel
    boundary_conditions::SedimentLandTransportDifferentiationBC
    variables::SedimentLandTransportDifferentiationVariables
end

"Initialize differentiated sediment flux in overland flow model"
function SedimentLandTransportDifferentiationModel(indices::Vector{CartesianIndex{2}})
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
function update!(model::SedimentLandTransportDifferentiationModel, network::NetworkLand)
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