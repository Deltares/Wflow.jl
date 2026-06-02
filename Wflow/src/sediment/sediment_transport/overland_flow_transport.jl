abstract type AbstractSedimentLandTransportModel end

"Struct to store total sediment flux in overland flow model variables"
@with_kw struct SedimentLandTransportVariables
    n::Int
    # Total sediment rate [kg s⁻¹]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total sediment deposition rate [kg s⁻¹]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store total sediment flux in overland flow model boundary conditions"
@with_kw struct SedimentLandTransportBC
    n::Int
    # Erosion rate material [kg s⁻¹]
    erosion::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity [kg s⁻¹]
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
    sediment_transport_model = SedimentLandTransportModel(; n)
    return sediment_transport_model
end

"Update total sediment flux in overland flow model boundary conditions"
function update_bc_sediment_land_transport_model!(
    sediment_transport_model::SedimentLandTransportModel,
    erosion_model::SoilErosionModel,
    transport_capacity_model::AbstractTransportCapacityModel,
)
    (; erosion, transport_capacity) = sediment_transport_model.boundary_conditions
    (; soil_erosion_rate) = erosion_model.variables
    @. erosion = soil_erosion_rate

    (; sediment_transport_capacity) = transport_capacity_model.variables
    @. transport_capacity = sediment_transport_capacity
end

"Update total sediment flux in overland flow model for a single timestep"
function update_sediment_overland_model!(
    sediment_transport_model::SedimentLandTransportModel,
    network::NetworkLand,
    dt::Float64,
)
    (; erosion, transport_capacity) = sediment_transport_model.boundary_conditions
    (; sediment_rate, deposition) = sediment_transport_model.variables

    # All inputs and outputs are rates [kg s⁻¹]
    accucapacityflux_rate!(
        sediment_rate,
        deposition,
        erosion,
        network,
        transport_capacity,
    )
end

"Struct to store differentiated sediment flux in overland flow model variables"
@with_kw struct SedimentLandTransportDifferentiationVariables
    n::Int
    # Total sediment rate [kg s⁻¹]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition rate [kg s⁻¹]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n)
    # Clay rate [kg s⁻¹]
    clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition clay rate [kg s⁻¹]
    deposition_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt rate [kg s⁻¹]
    silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition silt rate [kg s⁻¹]
    deposition_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand rate [kg s⁻¹]
    sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition sand rate [kg s⁻¹]
    deposition_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates rate [kg s⁻¹]
    small_aggregates::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition rate small aggregates [kg s⁻¹]
    deposition_small_aggregates::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates rate [kg s⁻¹]
    large_aggregates::Vector{Float64} = fill(MISSING_VALUE, n)
    # Deposition rate large aggregates [kg s⁻¹]
    deposition_large_aggregates::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store differentiated sediment flux in overland flow model boundary conditions"
@with_kw struct SedimentLandTransportDifferentiationBC
    n::Int
    # Erosion rate clay [kg s⁻¹]
    erosion_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Erosion rate silt [kg s⁻¹]
    erosion_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Erosion rate sand [kg s⁻¹]
    erosion_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Erosion rate small aggregates [kg s⁻¹]
    erosion_small_aggregates::Vector{Float64} = fill(MISSING_VALUE, n)
    # Erosion large aggregates [kg s⁻¹]
    erosion_large_aggregates::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity clay [kg s⁻¹]
    transport_capacity_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity silt [kg s⁻¹]
    transport_capacity_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity sand [kg s⁻¹]
    transport_capacity_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity small aggregates [kg s⁻¹]
    transport_capacity_small_aggregates::Vector{Float64} = fill(MISSING_VALUE, n)
    # Transport capacity large aggregates [kg s⁻¹]
    transport_capacity_large_aggregates::Vector{Float64} = fill(MISSING_VALUE, n)
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
    sediment_transport_model = SedimentLandTransportDifferentiationModel(; n)
    return sediment_transport_model
end

"Update differentiated sediment flux in overland flow model boundary conditions"
function update_bc_sediment_land_transport_model!(
    sediment_transport_model::SedimentLandTransportDifferentiationModel,
    erosion_model::SoilErosionModel,
    transport_capacity_model::TransportCapacityYalinDifferentiationModel,
)
    (;
        erosion_clay,
        erosion_silt,
        erosion_sand,
        erosion_small_aggregates,
        erosion_large_aggregates,
        transport_capacity_clay,
        transport_capacity_silt,
        transport_capacity_sand,
        transport_capacity_small_aggregates,
        transport_capacity_large_aggregates,
    ) = sediment_transport_model.boundary_conditions
    (;
        clay_erosion_rate,
        silt_erosion_rate,
        sand_erosion_rate,
        small_aggregates_erosion_rate,
        large_aggregates_erosion_rate,
    ) = erosion_model.variables
    @. erosion_clay = clay_erosion_rate
    @. erosion_silt = silt_erosion_rate
    @. erosion_sand = sand_erosion_rate
    @. erosion_small_aggregates = small_aggregates_erosion_rate
    @. erosion_large_aggregates = large_aggregates_erosion_rate

    (; clay, silt, sand, small_aggregates, large_aggregates) = transport_capacity_model.variables
    @. transport_capacity_clay = clay
    @. transport_capacity_silt = silt
    @. transport_capacity_sand = sand
    @. transport_capacity_small_aggregates = small_aggregates
    @. transport_capacity_large_aggregates = large_aggregates
end

"Update differentiated sediment flux in overland flow model for a single timestep"
function update_sediment_overland_model!(
    sediment_transport_model::SedimentLandTransportDifferentiationModel,
    network::NetworkLand,
    dt::Float64,
)
    (;
        erosion_clay,
        erosion_silt,
        erosion_sand,
        erosion_small_aggregates,
        erosion_large_aggregates,
        transport_capacity_clay,
        transport_capacity_silt,
        transport_capacity_sand,
        transport_capacity_small_aggregates,
        transport_capacity_large_aggregates,
    ) = sediment_transport_model.boundary_conditions
    (;
        sediment_rate,
        deposition,
        clay,
        deposition_clay,
        silt,
        deposition_silt,
        sand,
        deposition_sand,
        small_aggregates,
        deposition_small_aggregates,
        large_aggregates,
        deposition_large_aggregates,
    ) = sediment_transport_model.variables

    # All inputs and outputs are rates [kg s⁻¹]
    do_accucapacityflux!(rate, dep, erosion, tc) =
        accucapacityflux_rate!(rate, dep, erosion, network, tc)

    do_accucapacityflux!(clay, deposition_clay, erosion_clay, transport_capacity_clay)
    do_accucapacityflux!(silt, deposition_silt, erosion_silt, transport_capacity_silt)
    do_accucapacityflux!(sand, deposition_sand, erosion_sand, transport_capacity_sand)
    do_accucapacityflux!(small_aggregates, deposition_small_aggregates, erosion_small_aggregates, transport_capacity_small_aggregates)
    do_accucapacityflux!(large_aggregates, deposition_large_aggregates, erosion_large_aggregates, transport_capacity_large_aggregates)

    @. sediment_rate = clay + silt + sand + small_aggregates + large_aggregates
    @. deposition =
        deposition_clay +
        deposition_silt +
        deposition_sand +
        deposition_small_aggregates +
        deposition_large_aggregates
end
