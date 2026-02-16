abstract type AbstractSedimentLandTransportModel end

"Struct to store total sediment flux in overland flow model variables"
@with_kw struct SedimentLandTransportVariables
    n::Int
    # Total sediment rate [t dt-1]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total sediment deposition rate [t dt-1]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store total sediment flux in overland flow model boundary conditions"
@with_kw struct SedimentLandTransportBC
    n::Int
    # Erosion rate material [t dt-1]
    erosion::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity [t dt-1]
    transport_capacity::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store total sediment flux in overland flow model"
@with_kw struct SedimentLandTransportModel <: AbstractSedimentLandTransportModel
    n::Int
    boundary_conditions::SedimentLandTransportBC = SedimentLandTransportBC(; n)
    variables::SedimentLandTransportVariables = SedimentLandTransportVariables(; n)
end

"Initialize total sediment flux in overland flow model"
function SedimentLandTransportModel(indices::Vector{CartesianIndex{2}})
    n = length(indices)
    model = SedimentLandTransportModel(; n)
    return model
end

"Update total sediment flux in overland flow model boundary conditions"
function update_boundary_conditions!(
    model::SedimentLandTransportModel,
    erosion_model::SoilErosionModel,
    transport_capacity_model::AbstractTransportCapacityModel,
)
    (; erosion, transport_capacity) = model.boundary_conditions
    (; soil_erosion_rate) = erosion_model.variables
    @. erosion = soil_erosion_rate

    (; sediment_transport_capacity) = transport_capacity_model.variables
    @. transport_capacity = sediment_transport_capacity
end

"Update total sediment flux in overland flow model for a single timestep"
function update!(model::SedimentLandTransportModel, network::NetworkLand)
    (; erosion, transport_capacity) = model.boundary_conditions
    (; sediment_rate, deposition) = model.variables

    accucapacityflux!(sediment_rate, erosion, network, transport_capacity)
    deposition .= erosion
end

"Struct to store differentiated sediment flux in overland flow model variables"
@with_kw struct SedimentLandTransportDifferentiationVariables
    n::Int
    # Total sediment rate [t dt-1]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition rate [t dt-1]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n)
    # Clay rate [t dt-1]
    clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition clay rate [t dt-1]
    deposition_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt rate [t dt-1]
    silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition silt rate [t dt-1]
    deposition_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand rate [t dt-1]
    sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition sand rate [t dt-1]
    deposition_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates rate [t dt-1]
    sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition rate small aggregates [t dt-1]
    deposition_sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates rate [t dt-1]
    lagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition rate large aggregates [t dt-1]
    deposition_lagg::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store differentiated sediment flux in overland flow model boundary conditions"
@with_kw struct SedimentLandTransportDifferentiationBC
    n::Int
    # Erosion rate clay [t dt-1]
    erosion_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Erosion rate silt [t dt-1]
    erosion_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Erosion rate sand [t dt-1]
    erosion_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Erosion rate small aggregates [t dt-1]
    erosion_sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Erosion large aggregates [t dt-1]
    erosion_lagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity clay [t dt-1]
    transport_capacity_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity silt [t dt-1]
    transport_capacity_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity sand [t dt-1]
    transport_capacity_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity small aggregates [t dt-1]
    transport_capacity_sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity large aggregates [t dt-1]
    transport_capacity_lagg::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store differentiated sediment flux in overland flow model"
@with_kw struct SedimentLandTransportDifferentiationModel <:
                AbstractSedimentLandTransportModel
    n::Int
    boundary_conditions::SedimentLandTransportDifferentiationBC =
        SedimentLandTransportDifferentiationBC(; n)
    variables::SedimentLandTransportDifferentiationVariables =
        SedimentLandTransportDifferentiationVariables(; n)
end

"Initialize differentiated sediment flux in overland flow model"
function SedimentLandTransportDifferentiationModel(indices::Vector{CartesianIndex{2}})
    n = length(indices)
    model = SedimentLandTransportDifferentiationModel(; n)
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
    (;
        clay_erosion_rate,
        silt_erosion_rate,
        sand_erosion_rate,
        sagg_erosion_rate,
        lagg_erosion_rate,
    ) = erosion_model.variables
    @. erosion_clay = clay_erosion_rate
    @. erosion_silt = silt_erosion_rate
    @. erosion_sand = sand_erosion_rate
    @. erosion_sagg = sagg_erosion_rate
    @. erosion_lagg = lagg_erosion_rate

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
        sediment_rate,
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
    @. sediment_rate = clay + silt + sand + sagg + lagg
    @. deposition =
        deposition_clay +
        deposition_silt +
        deposition_sand +
        deposition_sagg +
        deposition_lagg
end
