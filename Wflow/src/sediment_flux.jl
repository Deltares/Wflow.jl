"Sediment transport in overland flow model"
@kwdef struct OverlandFlowSedimentModel{
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
    indices;
    data_lookup::DataLookup = DataLookup(),
)::AbstractTransportCapacityModel
    transport_capacity_constr = get(transport_methods, transport_method, nothing)
    @assert !isnothing(transport_capacity_constr)
    return transport_capacity_constr(dataset, config, indices; data_lookup)
end

const land_transport_method =
    Dict{LandTransportType.T, Type{<:AbstractTransportCapacityModel}}(
        LandTransportType.yalinpart => TransportCapacityYalinDifferentiationModel,
        LandTransportType.govers => TransportCapacityGoversModel,
        LandTransportType.yalin => TransportCapacityYalinModel,
    )

"Initialize the overland flow sediment transport model"
function OverlandFlowSedimentModel(
    dataset::NCDataset,
    config::Config,
    domain::DomainLand,
    soilloss::SoilLossModel;
    data_lookup::DataLookup = DataLookup(),
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
        indices;
        data_lookup,
    )

    if do_river || land_transport == LandTransportType.yalinpart
        sediment_flux = SedimentLandTransportDifferentiationModel(indices; data_lookup)
        to_river = SedimentToRiverDifferentiationModel(indices; data_lookup)
    else
        sediment_flux = SedimentLandTransportModel(indices; data_lookup)
        to_river = SedimentToRiverModel(indices; data_lookup)
    end

    overland_flow_sediment = OverlandFlowSedimentModel{
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
function update_overland_flow_model!(
    overland_flow_model::OverlandFlowSedimentModel,
    erosion_model::SoilErosionModel,
    domain::DomainLand,
    dt::Float64,
)
    # Transport capacity
    update_bc_transport_capacity_model!(
        overland_flow_model.transport_capacity,
        overland_flow_model.hydrological_forcing,
        :land,
    )
    update_transport_capacity_model!(
        overland_flow_model.transport_capacity,
        domain.parameters,
        dt,
    )

    # Update boundary conditions before transport
    update_bc_sediment_land_transport_model!(
        overland_flow_model.sediment_flux,
        erosion_model,
        overland_flow_model.transport_capacity,
    )
    # Compute transport
    update_sediment_overland_model!(overland_flow_model.sediment_flux, domain.network)

    # Update boundary conditions before computing sediment reaching the river
    update_bc_sediment_to_river_model!(
        overland_flow_model.to_river,
        overland_flow_model.sediment_flux,
    )
    # Compute sediment reaching the river
    update_sediment_to_river_model!(
        overland_flow_model.to_river,
        domain.parameters.river_location,
    )
end

### River ###
"Sediment transport in river model"
@kwdef struct RiverSedimentModel{
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
function RiverSedimentModel(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver;
    data_lookup::DataLookup = DataLookup(),
)
    (; indices) = domain.network
    n = length(indices)
    # Construct HydrologicalForcing without data_lookup to avoid overwriting land-domain
    # registrations (e.g. "land_surface_water__volume_flow_rate") with river-sized vectors.
    # Only register river-domain names explicitly.
    hydrological_forcing = HydrologicalForcing(; n)
    data_lookup["river_water__depth"] = hydrological_forcing.waterlevel_river
    data_lookup["river_water__volume_flow_rate"] = hydrological_forcing.q_river

    # Check what transport capacity equation will be used
    # River flow transport capacity method: ["bagnold", "engelund", "yang", "kodatie", "molinas"]
    (; river_transport) = config.model
    transport_capacity = get_transport_capacity(
        river_transport_method,
        river_transport,
        dataset,
        config,
        indices;
        data_lookup,
    )

    # Potential river erosion
    potential_erosion = RiverErosionJulianTorresModel(dataset, config, indices; data_lookup)

    # Sediment flux in river / mass balance
    sediment_flux = SedimentRiverTransportModel(dataset, config, indices; data_lookup)

    # Concentrations
    concentrations = SedimentConcentrationsRiverModel(dataset, config, indices; data_lookup)

    river_sediment = RiverSedimentModel(;
        hydrological_forcing,
        transport_capacity,
        potential_erosion,
        sediment_flux,
        concentrations,
    )
    return river_sediment
end

"Update the river sediment transport model for a single timestep"
function update_river_sediment_model!(
    river_flow_model::RiverSedimentModel,
    sediment_to_river_model::SedimentToRiverDifferentiationModel,
    domain::DomainRiver,
    dt::Float64,
)
    # Transport capacity
    update_bc_transport_capacity_model!(
        river_flow_model.transport_capacity,
        river_flow_model.hydrological_forcing,
        :river,
    )
    update_transport_capacity_model!(
        river_flow_model.transport_capacity,
        domain.parameters,
        dt,
    )

    # Potential maximum river erosion
    update_bc_river_erosion_model!(
        river_flow_model.potential_erosion,
        river_flow_model.hydrological_forcing,
    )
    update_river_erosion_model!(river_flow_model.potential_erosion, domain.parameters, dt)

    # River transport
    update_bc_river_sediment_transport_model!(
        river_flow_model.sediment_flux,
        river_flow_model.hydrological_forcing,
        river_flow_model.transport_capacity,
        sediment_to_river_model,
        river_flow_model.potential_erosion,
        domain.network.land_indices,
    )
    update_sediment_river_transport_model!(river_flow_model.sediment_flux, domain, dt)

    # Concentrations
    update_bc_river_sediment_concentration_model!(
        river_flow_model.concentrations,
        river_flow_model.hydrological_forcing,
        river_flow_model.sediment_flux,
    )
    update_river_sediment_concentration_model!(
        river_flow_model.concentrations,
        domain.parameters,
        dt,
    )
end
