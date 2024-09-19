@get_units @with_kw struct OverlandFlowSediment{TT, SF, TR, T}
    hydrometeo_forcing::HydrometeoForcing | "-"
    transport_capacity::TT | "-"
    sediment_flux::SF | "-"
    to_river::TR | "-"
    width::Vector{T} | "m"
    waterbodies::Vector{Bool} | "-"
    rivers::Vector{Bool} | "-"
end

function initialize_overland_flow_sediment(nc, config, inds, width, waterbodies, rivers)
    n = length(inds)
    hydrometeo_forcing = initialize_hydrometeo_forcing(n)
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
        transport_capacity = transport_capacity_model,
        sediment_flux = sediment_flux_model,
        to_river = to_river_model,
        width = width,
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