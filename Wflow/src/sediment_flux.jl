"Sediment transport in overland flow model"
@with_kw struct OverlandFlowSediment{TT, SF, TR} <: AbstractOverlandFlowModel
    hydrological_forcing::HydrologicalForcing
    transport_capacity::TT
    sediment_flux::SF
    to_river::TR
end

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
    do_river = get(config.model, "run_river_model__flag", false)::Bool
    # Overland flow transport capacity method: ["yalinpart", "govers", "yalin"]
    landtransportmethod = get(config.model, "land_transport", "yalinpart")::String

    if do_river || landtransportmethod == "yalinpart"
        transport_capacity =
            TransportCapacityYalinDifferentiationModel(dataset, config, indices)
    elseif landtransportmethod == "govers"
        transport_capacity = TransportCapacityGoversModel(dataset, config, indices)
    elseif landtransportmethod == "yalin"
        transport_capacity = TransportCapacityYalinModel(dataset, config, indices)
    else
        error("Unknown land transport method: $landtransportmethod")
    end

    if do_river || landtransportmethod == "yalinpart"
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
    dt::Float,
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
@with_kw struct RiverSediment{TTR, ER, SFR, CR} <: AbstractRiverFlowModel
    hydrological_forcing::HydrologicalForcing
    transport_capacity::TTR
    potential_erosion::ER
    sediment_flux::SFR
    concentrations::CR
end

"Initialize the river sediment transport model"
function RiverSediment(dataset::NCDataset, config::Config, domain::DomainRiver)
    (; indices) = domain.network
    n = Int(length(indices))
    hydrological_forcing = HydrologicalForcing(n)

    # Check what transport capacity equation will be used
    # River flow transport capacity method: ["bagnold", "engelund", "yang", "kodatie", "molinas"]
    transport_method = get(config.model, "river_transport", "bagnold")::String
    if transport_method == "bagnold"
        transport_capacity = TransportCapacityBagnoldModel(dataset, config, indices)
    elseif transport_method == "engelund"
        transport_capacity = TransportCapacityEngelundModel(dataset, config, indices)
    elseif transport_method == "yang"
        transport_capacity = TransportCapacityYangModel(dataset, config, indices)
    elseif transport_method == "kodatie"
        transport_capacity = TransportCapacityKodatieModel(dataset, config, indices)
    elseif transport_method == "molinas"
        transport_capacity = TransportCapacityMolinasModel(dataset, config, indices)
    else
        error("Unknown river transport method: $transport_method")
    end

    # Potential river erosion
    potential_erosion = RiverErosionJulianTorresModel(dataset, config, indices)

    # Sediment flux in river / mass balance
    sediment_flux = SedimentRiverTransportModel(dataset, config, indices)

    # Concentrations
    concentrations = SedimentConcentrationsRiverModel(dataset, config, indices)

    river_sediment = RiverSediment{
        typeof(transport_capacity),
        typeof(potential_erosion),
        typeof(sediment_flux),
        typeof(concentrations),
    }(;
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
    dt::Float,
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