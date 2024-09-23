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
    hydrometeo_forcing = initialize_hydrometeo_forcing(n)
    geometry = initialize_land_geometry(nc, config, inds)
    # Check what transport capacity equation will be used
    do_river = get(config.model, "runrivermodel", false)::Bool
    # Overland flow transport capacity method: ["yalinpart", "govers", "yalin"]
    landtransportmethod = get(config.model, "landtransportmethod", "yalinpart")::String

    if do_river || landtransportmethod == "yalinpart"
        transport_capacity_model =
            initialize_transport_capacity_yalin_diff_model(nc, config, inds)
    elseif landtransportmethod == "govers"
        transport_capacity_model =
            initialize_transport_capacity_govers_model(nc, config, inds)
    elseif landtransportmethod == "yalin"
        transport_capacity_model =
            initialize_transport_capacity_yalin_model(nc, config, inds)
    else
        error("Unknown land transport method: $landtransportmethod")
    end

    if do_river || landtransportmethod == "yalinpart"
        sediment_flux_model = initialize_sediment_land_transport_differentiation_model(inds)
        to_river_model = initialize_sediment_to_river_differentiation_model(inds)
    else
        sediment_flux_model = initialize_sediment_land_transport_model(inds)
        to_river_model = initialize_sediment_to_river_model(inds)
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
    # Update the boundary conditions of transport capacity
    (; q, waterlevel) = model.transport_capacity.boundary_conditions
    (; q_land, waterlevel_land) = model.hydrometeo_forcing
    @. q = q_land
    @. waterlevel = waterlevel_land
    # Transport capacity
    update!(model.transport_capacity, model.width, model.waterbodies, model.rivers, ts)

    # Update boundary conditions before transport
    update_bc(model.sediment_flux, erosion_model, model.transport_capacity)
    # Compute transport
    update!(model.sediment_flux, network)

    # Update boundary conditions before computing sediment reaching the river
    update_bc(model.to_river, model.sediment_flux)
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
    hydrometeo_forcing = initialize_hydrometeo_forcing(n)
    geometry = initialize_river_geometry(nc, config, inds)

    # Check what transport capacity equation will be used
    # River flow transport capacity method: ["bagnold", "engelund", "yang", "kodatie", "molinas"]
    transport_method = get(config.model, "rivtransportmethod", "bagnold")::String
    if transport_method == "bagnold"
        transport_capacity_model =
            initialize_transport_capacity_bagnold_model(nc, config, inds)
    elseif transport_method == "engelund"
        transport_capacity_model =
            initialize_transport_capacity_engelund_model(nc, config, inds)
    elseif transport_method == "yang"
        transport_capacity_model =
            initialize_transport_capacity_yang_model(nc, config, inds)
    elseif transport_method == "kodatie"
        transport_capacity_model =
            initialize_transport_capacity_kodatie_model(nc, config, inds)
    elseif transport_method == "molinas"
        transport_capacity_model =
            initialize_transport_capacity_molinas_model(nc, config, inds)
    else
        error("Unknown river transport method: $transport_method")
    end

    # Potential river erosion
    potential_erosion_model = initialize_river_erosion_julian_torres_model(nc, config, inds)

    # Sediment flux in river / mass balance
    sediment_flux_model = initialize_sediment_river_transport_model(nc, config, inds)

    # Concentrations
    concentrations_model = initialize_sediment_concentrations_river_model(nc, config, inds)

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
    # Update the boundary conditions of transport capacity
    (; q, waterlevel) = model.transport_capacity.boundary_conditions
    (; q_river, waterlevel_river) = model.hydrometeo_forcing
    @. q = q_river
    @. waterlevel = waterlevel_river

    # Transport capacity
    update!(model.transport_capacity, model.geometry, ts)

    # Potential maximum river erosion
    (; waterlevel) = model.potential_erosion.boundary_conditions
    @. waterlevel = waterlevel_river
    update!(model.potential_erosion, model.geometry, ts)

    # River transport
    (;
        waterlevel,
        q,
        transport_capacity,
        erosion_land_clay,
        erosion_land_silt,
        erosion_land_sand,
        erosion_land_sagg,
        erosion_land_lagg,
        potential_erosion_river,
    ) = model.boundary_conditions
    @. q = q_river
    @. waterlevel = waterlevel_river

    @. transport_capacity = model.transport_capacity.variables.amount

    (; clay, silt, sand, sagg, lagg) = to_river_model.variables
    @. erosion_land_clay = clay[inds_riv]
    @. erosion_land_silt = silt[inds_riv]
    @. erosion_land_sand = sand[inds_riv]
    @. erosion_land_sagg = sagg[inds_riv]
    @. erosion_land_lagg = lagg[inds_riv]

    @. potential_erosion_river = model.potential_erosion.variables.amount

    update!(model.sediment_flux, network, model.geometry, ts)

    # Concentrations
    (; q, waterlevel, clay, silt, sand, sagg, lagg, gravel) =
        model.concentrations.boundary_conditions
    @. q = q_river
    @. waterlevel = waterlevel_river
    @. clay = model.sediment_flux.variables.clay
    @. silt = model.sediment_flux.variables.silt
    @. sand = model.sediment_flux.variables.sand
    @. sagg = model.sediment_flux.variables.sagg
    @. lagg = model.sediment_flux.variables.lagg
    @. gravel = model.sediment_flux.variables.gravel

    update!(model.concentrations, model.geometry, ts)
end