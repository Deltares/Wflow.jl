### Overland flow ###
@get_units @with_kw struct OverlandFlowSediment{TT, SF, TR, T}
    hydrometeo_forcing::HydrometeoForcing | "-"
    geometry::LandGeometry | "-"
    transport_capacity::TT | "-"
    sediment_flux::SF | "-"
    to_river::TR | "-"
    waterbodies::Vector{Bool} | "-"
    rivers::Vector{Bool} | "-"
end

function initialize_overland_flow_sediment(nc, config, inds, waterbodies, rivers)
    n = length(inds)
    hydrometeo_forcing = HydrometeoForcing(n)
    geometry = LandGeometry(nc, config, inds)
    # Check what transport capacity equation will be used
    do_river = get(config.model, "runrivermodel", false)::Bool
    # Overland flow transport capacity method: ["yalinpart", "govers", "yalin"]
    landtransportmethod = get(config.model, "landtransportmethod", "yalinpart")::String

    if do_river || landtransportmethod == "yalinpart"
        transport_capacity_model =
            TransportCapacityYalinDifferentiationModel(nc, config, inds)
    elseif landtransportmethod == "govers"
        transport_capacity_model = TransportCapacityGoversModel(nc, config, inds)
    elseif landtransportmethod == "yalin"
        transport_capacity_model = TransportCapacityYalinModel(nc, config, inds)
    else
        error("Unknown land transport method: $landtransportmethod")
    end

    if do_river || landtransportmethod == "yalinpart"
        sediment_flux_model = SedimentLandTransportDifferentiationModel(inds)
        to_river_model = SedimentToRiverDifferentiationModel(inds)
    else
        sediment_flux_model = SedimentLandTransportModel(inds)
        to_river_model = SedimentToRiverModel(inds)
    end

    overland_flow_sediment = OverlandFlowSediment{
        typeof(transport_capacity_model),
        typeof(sediment_flux_model),
        typeof(to_river_model),
        Float,
    }(;
        hydrometeo_forcing = hydrometeo_forcing,
        geometry = geometry,
        transport_capacity = transport_capacity_model,
        sediment_flux = sediment_flux_model,
        to_river = to_river_model,
        waterbodies = waterbodies,
        rivers = rivers,
    )
    return overland_flow_sediment
end

function update!(model::OverlandFlowSediment, erosion_model::SoilErosionModel, network, dt)
    # Convert dt to integer
    ts = tosecond(dt)

    # Transport capacity
    update_boundary_conditions!(model.transport_capacity, model.hydrometeo_forcing, :land)
    update!(model.transport_capacity, model.geometry, model.waterbodies, model.rivers, ts)

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
@get_units @with_kw struct RiverSediment{TTR, ER, SFR, CR, T}
    hydrometeo_forcing::HydrometeoForcing | "-"
    geometry::RiverGeometry | "-"
    transport_capacity::TTR | "-"
    potential_erosion::ER | "-"
    sediment_flux::SFR | "-"
    concentrations::CR | "-"
    waterbodies::Vector{Bool} | "-"
end

function initialize_river_flow_sediment(nc, config, inds, waterbodies)
    n = length(inds)
    hydrometeo_forcing = HydrometeoForcing(n)
    geometry = RiverGeometry(nc, config, inds)

    # Check what transport capacity equation will be used
    # River flow transport capacity method: ["bagnold", "engelund", "yang", "kodatie", "molinas"]
    transport_method = get(config.model, "rivtransportmethod", "bagnold")::String
    if transport_method == "bagnold"
        transport_capacity_model = TransportCapacityBagnoldModel(nc, config, inds)
    elseif transport_method == "engelund"
        transport_capacity_model = TransportCapacityEngelundModel(nc, config, inds)
    elseif transport_method == "yang"
        transport_capacity_model = TransportCapacityYangModel(nc, config, inds)
    elseif transport_method == "kodatie"
        transport_capacity_model = TransportCapacityKodatieModel(nc, config, inds)
    elseif transport_method == "molinas"
        transport_capacity_model = TransportCapacityMolinasModel(nc, config, inds)
    else
        error("Unknown river transport method: $transport_method")
    end

    # Potential river erosion
    potential_erosion_model = RiverErosionJulianTorresModel(nc, config, inds)

    # Sediment flux in river / mass balance
    sediment_flux_model = SedimentRiverTransportModel(nc, config, inds)

    # Concentrations
    concentrations_model = SedimentConcentrationsRiverModel(nc, config, inds)

    river_sediment = RiverSediment{
        typeof(transport_capacity_model),
        typeof(potential_erosion_model),
        typeof(sediment_flux_model),
        typeof(concentrations_model),
        Float,
    }(;
        hydrometeo_forcing = hydrometeo_forcing,
        geometry = geometry,
        transport_capacity = transport_capacity_model,
        potential_erosion = potential_erosion_model,
        sediment_flux = sediment_flux_model,
        concentrations = concentrations_model,
        waterbodies = waterbodies,
    )
    return river_sediment
end

function update!(
    model::RiverSediment,
    to_river_model::SedimentToRiverDifferentiationModel,
    network,
    inds_riv,
    dt,
)
    # Convert dt to integer
    ts = tosecond(dt)

    # Transport capacity
    update_boundary_conditions!(model.transport_capacity, model.hydrometeo_forcing, :river)
    update!(model.transport_capacity, model.geometry, ts)

    # Potential maximum river erosion
    update_boundary_conditions!(model.potential_erosion, model.hydrometeo_forcing)
    update!(model.potential_erosion, model.geometry, ts)

    # River transport
    update_boundary_conditions!(
        model.sediment_flux,
        model.hydrometeo_forcing,
        model.transport_capacity,
        to_river_model,
        model.potential_erosion,
        inds_riv,
    )
    update!(model.sediment_flux, network, model.geometry, model.waterbodies, ts)

    # Concentrations
    update_boundary_conditions!(
        model.concentrations,
        model.hydrometeo_forcing,
        model.sediment_flux,
    )
    update!(model.concentrations, model.geometry, ts)
end