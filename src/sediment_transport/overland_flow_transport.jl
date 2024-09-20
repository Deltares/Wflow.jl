abstract type AbstractSedimentLandTransportModel end

## Total sediment transport in overland flow structs and functions
@get_units @with_kw struct SedimentLandTransportVars{T}
    # Total sediment flux
    amount::Vector{T} | "t dt-1"
    deposition::Vector{T} | "t dt-1"
end

function sediment_land_transport_vars(n)
    vars = SedimentLandTransportVars(; amount = fill(mv, n), deposition = fill(mv, n))
    return vars
end

@get_units @with_kw struct SedimentLandTransportBC{T}
    # Eroded material
    erosion::Vector{T} | "t dt-1"
    # Transport capacity
    transport_capacity::Vector{T} | "t dt-1"
end

function sediment_land_transport_bc(n)
    bc = SedimentLandTransportBC(; erosion = fill(mv, n), transport_capacity = fill(mv, n))
    return bc
end

@get_units @with_kw struct SedimentLandTransportModel{T} <:
                           AbstractSedimentLandTransportModel
    boundary_conditions::SedimentLandTransportBC{T} | "-"
    variables::SedimentLandTransportVars{T} | "-"
end

function initialize_sediment_land_transport_model(inds)
    n = length(inds)
    vars = sediment_land_transport_vars(n)
    bc = sediment_land_transport_bc(n)
    model = SedimentLandTransportModel(; boundary_conditions = bc, variables = vars)
    return model
end

function update_bc(
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

function update!(model::SedimentLandTransportModel, network)
    (; erosion, transport_capacity) = model.boundary_conditions
    (; amount, deposition) = model.variables

    accucapacityflux!(amount, erosion, network, transport_capacity)
    deposition .= erosion
end

## Total transport capacity with particle differentiation structs and functions
@get_units @with_kw struct SedimentLandTransportDifferentiationVars{T}
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

function sediment_land_transport_differentiation_vars(n)
    vars = SedimentLandTransportDifferentiationVars(;
        amount = fill(mv, n),
        deposition = fill(mv, n),
        clay = fill(mv, n),
        deposition_clay = fill(mv, n),
        silt = fill(mv, n),
        deposition_silt = fill(mv, n),
        sand = fill(mv, n),
        deposition_sand = fill(mv, n),
        sagg = fill(mv, n),
        deposition_sagg = fill(mv, n),
        lagg = fill(mv, n),
        deposition_lagg = fill(mv, n),
    )
    return vars
end

@get_units @with_kw struct SedimentLandTransportDifferentiationBC{T}
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

function sediment_land_transport_differentiation_bc(n)
    bc = SedimentLandTransportDifferentiationBC(;
        erosion_clay = fill(mv, n),
        erosion_silt = fill(mv, n),
        erosion_sand = fill(mv, n),
        erosion_sagg = fill(mv, n),
        erosion_lagg = fill(mv, n),
        transport_capacity_clay = fill(mv, n),
        transport_capacity_silt = fill(mv, n),
        transport_capacity_sand = fill(mv, n),
        transport_capacity_sagg = fill(mv, n),
        transport_capacity_lagg = fill(mv, n),
    )
    return bc
end

@get_units @with_kw struct SedimentLandTransportDifferentiationModel{T} <:
                           AbstractSedimentLandTransportModel
    boundary_conditions::SedimentLandTransportDifferentiationBC{T} | "-"
    variables::SedimentLandTransportDifferentiationVars{T} | "-"
end

function initialize_sediment_land_transport_differentiation_model(inds)
    n = length(inds)
    vars = sediment_land_transport_differentiation_vars(n)
    bc = sediment_land_transport_differentiation_bc(n)
    model = SedimentLandTransportDifferentiationModel(;
        boundary_conditions = bc,
        variables = vars,
    )
    return model
end

function update_bc(
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