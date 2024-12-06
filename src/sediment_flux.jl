"Sediment transport in overland flow model"
@get_units @grid_loc @with_kw struct OverlandFlowSediment{TT, SF, TR}
    hydrological_forcing::HydrologicalForcing
    geometry::LandGeometry
    transport_capacity::TT
    sediment_flux::SF
    to_river::TR
    waterbodies::Vector{Bool} | "-"
    rivers::Vector{Bool} | "-"
end

"Initialize the overland flow sediment transport model"
function OverlandFlowSediment(dataset, config, indices, waterbodies, rivers)
    n = length(indices)
    hydrological_forcing = HydrologicalForcing(n)
    geometry = LandGeometry(dataset, config, indices)
    # Check what transport capacity equation will be used
    do_river = get(config.model, "run_river_model", false)::Bool
    # Overland flow transport capacity method: ["yalinpart", "govers", "yalin"]
    landtransportmethod = get(config.model, "land_transport", "yalinpart")::String

    if do_river || landtransportmethod == "yalinpart"
        transport_capacity_model =
            TransportCapacityYalinDifferentiationModel(dataset, config, indices)
    elseif landtransportmethod == "govers"
        transport_capacity_model = TransportCapacityGoversModel(dataset, config, indices)
    elseif landtransportmethod == "yalin"
        transport_capacity_model = TransportCapacityYalinModel(dataset, config, indices)
    else
        error("Unknown land transport method: $landtransportmethod")
    end

    if do_river || landtransportmethod == "yalinpart"
        sediment_flux_model = SedimentLandTransportDifferentiationModel(indices)
        to_river_model = SedimentToRiverDifferentiationModel(indices)
    else
        sediment_flux_model = SedimentLandTransportModel(indices)
        to_river_model = SedimentToRiverModel(indices)
    end

    overland_flow_sediment = OverlandFlowSediment{
        typeof(transport_capacity_model),
        typeof(sediment_flux_model),
        typeof(to_river_model),
    }(;
        hydrological_forcing = hydrological_forcing,
        geometry = geometry,
        transport_capacity = transport_capacity_model,
        sediment_flux = sediment_flux_model,
        to_river = to_river_model,
        waterbodies = waterbodies,
        rivers = rivers,
    )
    return overland_flow_sediment
end

"Update the overland flow sediment transport model for a single timestep"
function update!(model::OverlandFlowSediment, erosion_model::SoilErosionModel, network, dt)
    # Transport capacity
    update_boundary_conditions!(model.transport_capacity, model.hydrological_forcing, :land)
    update!(model.transport_capacity, model.geometry, model.waterbodies, model.rivers, dt)

    # Update boundary conditions before transport
    update_boundary_conditions!(
        model.sediment_flux,
        erosion_model,
        model.transport_capacity,
    )
    # Compute transport
    update!(model.sediment_flux, network)

    # Update boundary conditions before computing sediment reaching the river
    update_boundary_conditions!(model.to_river, model.sediment_flux)
    # Compute sediment reaching the river
    update!(model.to_river, model.rivers)
end

### River ###
"Sediment transport in river model"
@get_units @grid_loc @with_kw struct RiverSediment{TTR, ER, SFR, CR}
    hydrological_forcing::HydrologicalForcing
    geometry::RiverGeometry
    transport_capacity::TTR
    potential_erosion::ER
    sediment_flux::SFR
    concentrations::CR
    waterbodies::Vector{Bool} | "-"
end

"Initialize the river sediment transport model"
function RiverSediment(dataset, config, indices, waterbodies)
    n = length(indices)
    hydrological_forcing = HydrologicalForcing(n)
    geometry = RiverGeometry(dataset, config, indices)

    # Check what transport capacity equation will be used
    # River flow transport capacity method: ["bagnold", "engelund", "yang", "kodatie", "molinas"]
    transport_method = get(config.model, "river_transport", "bagnold")::String
    if transport_method == "bagnold"
        transport_capacity_model = TransportCapacityBagnoldModel(dataset, config, indices)
    elseif transport_method == "engelund"
        transport_capacity_model = TransportCapacityEngelundModel(dataset, config, indices)
    elseif transport_method == "yang"
        transport_capacity_model = TransportCapacityYangModel(dataset, config, indices)
    elseif transport_method == "kodatie"
        transport_capacity_model = TransportCapacityKodatieModel(dataset, config, indices)
    elseif transport_method == "molinas"
        transport_capacity_model = TransportCapacityMolinasModel(dataset, config, indices)
    else
        error("Unknown river transport method: $transport_method")
    end

    # Potential river erosion
    potential_erosion_model = RiverErosionJulianTorresModel(dataset, config, indices)

    # Sediment flux in river / mass balance
    sediment_flux_model = SedimentRiverTransportModel(dataset, config, indices)

    # Concentrations
    concentrations_model = SedimentConcentrationsRiverModel(dataset, config, indices)

    river_sediment = RiverSediment{
        typeof(transport_capacity_model),
        typeof(potential_erosion_model),
        typeof(sediment_flux_model),
        typeof(concentrations_model),
    }(;
        hydrological_forcing = hydrological_forcing,
        geometry = geometry,
        transport_capacity = transport_capacity_model,
        potential_erosion = potential_erosion_model,
        sediment_flux = sediment_flux_model,
        concentrations = concentrations_model,
        waterbodies = waterbodies,
    )
    return river_sediment
end

"Update the river sediment transport model for a single timestep"
function update!(
    model::RiverSediment,
    to_river_model::SedimentToRiverDifferentiationModel,
    network,
    indices_river,
    dt,
)
    # Transport capacity
    update_boundary_conditions!(
        model.transport_capacity,
        model.hydrological_forcing,
        :river,
    )
    update!(model.transport_capacity, model.geometry, dt)

    # Potential maximum river erosion
    update_boundary_conditions!(model.potential_erosion, model.hydrological_forcing)
    update!(model.potential_erosion, model.geometry, dt)

    # River transport
    update_boundary_conditions!(
        model.sediment_flux,
        model.hydrological_forcing,
        model.transport_capacity,
        to_river_model,
        model.potential_erosion,
        indices_river,
    )
    update!(model.sediment_flux, network, model.geometry, model.waterbodies, dt)

    # Concentrations
    update_boundary_conditions!(
        model.concentrations,
        model.hydrological_forcing,
        model.sediment_flux,
    )
    update!(model.concentrations, model.geometry, dt)
end