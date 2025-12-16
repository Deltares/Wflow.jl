"Sediment transport in overland flow model"
@with_kw struct OverlandFlowSediment{
    TT <: AbstractTransportCapacityModel,
    SF <: AbstractSedimentLandTransportModel,
    TR <: AbstractSedimentToRiverModel,
} <: AbstractOverlandFlowModel
    hydrological_forcing::HydrologicalForcing
    transport_capacity::TT
    sediment_flux::SF
    to_river::TR
end

function get_transport_capacity(
    transport_methods::Dict{<:EnumX.Enum, Type{<:AbstractTransportCapacityModel}},
    transport_method::Union{LandTransportType.T, RiverTransportType.T},
    dataset::NCDataset,
    config::Config,
    indices,
)::AbstractTransportCapacityModel
    transport_capacity_constr = get(transport_methods, transport_method, nothing)
    @assert !isnothing(transport_capacity_constr)
    return transport_capacity_constr(dataset, config, indices)
end

const land_transport_method =
    Dict{LandTransportType.T, Type{<:AbstractTransportCapacityModel}}(
        LandTransportType.yalinpart => TransportCapacityYalinDifferentiationModel,
        LandTransportType.govers => TransportCapacityGoversModel,
        LandTransportType.yalin => TransportCapacityYalinModel,
    )

"Initialize the overland flow sediment transport model"
function OverlandFlowSediment(
    dataset::NCDataset,
    config::Config,
    domain::DomainLand,
    soilloss::SoilLoss,
)
    (; indices) = domain.network
    (; hydrological_forcing) = soilloss

    # Check what transport capacity equation will be used
    do_river = config.model.run_river_model__flag
    # Overland flow transport capacity method: ["yalinpart", "govers", "yalin"]
    (; land_transport) = config.model

    do_river && (land_transport = LandTransportType.yalinpart)
    transport_capacity = get_transport_capacity(
        land_transport_method,
        land_transport,
        dataset,
        config,
        indices,
    )

    if do_river || land_transport == LandTransportType.yalinpart
        sediment_flux = SedimentLandTransportDifferentiationModel(indices)
        to_river = SedimentToRiverDifferentiationModel(indices)
    else
        sediment_flux = SedimentLandTransportModel(indices)
        to_river = SedimentToRiverModel(indices)
    end

    overland_flow_sediment = OverlandFlowSediment{
        typeof(transport_capacity),
        typeof(sediment_flux),
        typeof(to_river),
    }(;
        hydrological_forcing,
        transport_capacity,
        sediment_flux,
        to_river,
    )
    return overland_flow_sediment
end

"Update the overland flow sediment transport model for a single timestep"
function update!(
    model::OverlandFlowSediment,
    erosion_model::SoilErosionModel,
    domain::DomainLand,
    dt::Float64,
)
    # Transport capacity
    update_boundary_conditions!(model.transport_capacity, model.hydrological_forcing, :land)
    update!(model.transport_capacity, domain.parameters, dt)

    # Update boundary conditions before transport
    update_boundary_conditions!(
        model.sediment_flux,
        erosion_model,
        model.transport_capacity,
    )
    # Compute transport
    update!(model.sediment_flux, domain.network)

    # Update boundary conditions before computing sediment reaching the river
    update_boundary_conditions!(model.to_river, model.sediment_flux)
    # Compute sediment reaching the river
    update!(model.to_river, domain.parameters.river_location)
end

### River ###
"Sediment transport in river model"
@with_kw struct RiverSediment{
    TTR <: AbstractTransportCapacityModel,
    ER <: AbstractRiverErosionModel,
    SFR <: AbstractSedimentRiverTransportModel,
    CR <: AbstractSedimentConcentrationsRiverModel,
} <: AbstractRiverFlowModel
    hydrological_forcing::HydrologicalForcing
    transport_capacity::TTR
    potential_erosion::ER
    sediment_flux::SFR
    concentrations::CR
end

const river_transport_method =
    Dict{RiverTransportType.T, Type{<:AbstractTransportCapacityModel}}(
        RiverTransportType.bagnold => TransportCapacityBagnoldModel,
        RiverTransportType.engelund => TransportCapacityEngelundModel,
        RiverTransportType.yang => TransportCapacityYangModel,
        RiverTransportType.kodatie => TransportCapacityKodatieModel,
        RiverTransportType.molinas => TransportCapacityMolinasModel,
    )

"Initialize the river sediment transport model"
function RiverSediment(dataset::NCDataset, config::Config, domain::DomainRiver)
    (; indices) = domain.network
    n = length(indices)
    hydrological_forcing = HydrologicalForcing(; n)

    # Check what transport capacity equation will be used
    # River flow transport capacity method: ["bagnold", "engelund", "yang", "kodatie", "molinas"]
    (; river_transport) = config.model
    transport_capacity = get_transport_capacity(
        river_transport_method,
        river_transport,
        dataset,
        config,
        indices,
    )

    # Potential river erosion
    potential_erosion = RiverErosionJulianTorresModel(dataset, config, indices)

    # Sediment flux in river / mass balance
    sediment_flux = SedimentRiverTransportModel(dataset, config, indices)

    # Concentrations
    concentrations = SedimentConcentrationsRiverModel(dataset, config, indices)

    river_sediment = RiverSediment(;
        hydrological_forcing,
        transport_capacity,
        potential_erosion,
        sediment_flux,
        concentrations,
    )
    return river_sediment
end

"Update the river sediment transport model for a single timestep"
function update!(
    model::RiverSediment,
    to_river_model::SedimentToRiverDifferentiationModel,
    domain::DomainRiver,
    dt::Float64,
)
    # Transport capacity
    update_boundary_conditions!(
        model.transport_capacity,
        model.hydrological_forcing,
        :river,
    )
    update!(model.transport_capacity, domain.parameters, dt)

    # Potential maximum river erosion
    update_boundary_conditions!(model.potential_erosion, model.hydrological_forcing)
    update!(model.potential_erosion, domain.parameters, dt)

    # River transport
    update_boundary_conditions!(
        model.sediment_flux,
        model.hydrological_forcing,
        model.transport_capacity,
        to_river_model,
        model.potential_erosion,
        domain.network.land_indices,
    )
    update!(model.sediment_flux, domain, dt)

    # Concentrations
    update_boundary_conditions!(
        model.concentrations,
        model.hydrological_forcing,
        model.sediment_flux,
    )
    update!(model.concentrations, domain.parameters, dt)
end
