abstract type AbstractSedimentLandTransportModel end

"Struct to store total sediment flux in overland flow model variables"
@with_kw struct SedimentLandTransportVariables
    n_land_cells::Int
    # Total sediment rate [t dt-1]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Total sediment deposition rate [t dt-1]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
end

"Struct to store total sediment flux in overland flow model boundary conditions"
@with_kw struct SedimentLandTransportBC
    n_land_cells::Int
    # Erosion rate material [t dt-1]
    erosion::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Transport capacity [t dt-1]
    transport_capacity::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
end

"Struct to store total sediment flux in overland flow model"
@with_kw struct SedimentLandTransportModel <: AbstractSedimentLandTransportModel
    n_land_cells::Int
    boundary_conditions::SedimentLandTransportBC = SedimentLandTransportBC(; n_land_cells)
    variables::SedimentLandTransportVariables =
        SedimentLandTransportVariables(; n_land_cells)
end

"Initialize total sediment flux in overland flow model"
function SedimentLandTransportModel(land_indices_2d::Vector{CartesianIndex{2}})
    n_land_cells = length(land_indices_2d)
    sediment_transport_model = SedimentLandTransportModel(; n_land_cells)
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
)
    (; erosion, transport_capacity) = sediment_transport_model.boundary_conditions
    (; sediment_rate, deposition) = sediment_transport_model.variables

    accucapacityflux!(sediment_rate, erosion, network, transport_capacity)
    deposition .= erosion
end

"Struct to store differentiated sediment flux in overland flow model variables"
@with_kw struct SedimentLandTransportDifferentiationVariables
    n_land_cells::Int
    # Total sediment rate [t dt-1]
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Deposition rate [t dt-1]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Clay rate [t dt-1]
    clay::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Deposition clay rate [t dt-1]
    deposition_clay::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Silt rate [t dt-1]
    silt::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Deposition silt rate [t dt-1]
    deposition_silt::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Sand rate [t dt-1]
    sand::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Deposition sand rate [t dt-1]
    deposition_sand::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Small aggregates rate [t dt-1]
    sagg::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Deposition rate small aggregates [t dt-1]
    deposition_sagg::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Large aggregates rate [t dt-1]
    lagg::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Deposition rate large aggregates [t dt-1]
    deposition_lagg::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
end

"Struct to store differentiated sediment flux in overland flow model boundary conditions"
@with_kw struct SedimentLandTransportDifferentiationBC
    n_land_cells::Int
    # Erosion rate clay [t dt-1]
    erosion_clay::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Erosion rate silt [t dt-1]
    erosion_silt::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Erosion rate sand [t dt-1]
    erosion_sand::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Erosion rate small aggregates [t dt-1]
    erosion_sagg::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Erosion large aggregates [t dt-1]
    erosion_lagg::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Transport capacity clay [t dt-1]
    transport_capacity_clay::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Transport capacity silt [t dt-1]
    transport_capacity_silt::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Transport capacity sand [t dt-1]
    transport_capacity_sand::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Transport capacity small aggregates [t dt-1]
    transport_capacity_sagg::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Transport capacity large aggregates [t dt-1]
    transport_capacity_lagg::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
end

"Struct to store differentiated sediment flux in overland flow model"
@with_kw struct SedimentLandTransportDifferentiationModel <:
                AbstractSedimentLandTransportModel
    n_land_cells::Int
    boundary_conditions::SedimentLandTransportDifferentiationBC =
        SedimentLandTransportDifferentiationBC(; n_land_cells)
    variables::SedimentLandTransportDifferentiationVariables =
        SedimentLandTransportDifferentiationVariables(; n_land_cells)
end

"Initialize differentiated sediment flux in overland flow model"
function SedimentLandTransportDifferentiationModel(
    land_indices_2d::Vector{CartesianIndex{2}},
)
    n_land_cells = length(land_indices_2d)
    sediment_transport_model = SedimentLandTransportDifferentiationModel(; n_land_cells)
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
        erosion_sagg,
        erosion_lagg,
        transport_capacity_clay,
        transport_capacity_silt,
        transport_capacity_sand,
        transport_capacity_sagg,
        transport_capacity_lagg,
    ) = sediment_transport_model.boundary_conditions
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
function update_sediment_overland_model!(
    sediment_transport_model::SedimentLandTransportDifferentiationModel,
    network::NetworkLand,
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
        sagg,
        deposition_sagg,
        lagg,
        deposition_lagg,
    ) = sediment_transport_model.variables

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
