"""
    surface_routing!(model)

Run surface routing (land and river) for a single timestep. Kinematic wave for overland flow
and kinematic wave or local inertial model for river flow.
"""
function surface_routing!(model)
    (; land, routing, domain, config, clock) = model
    (; soil, runoff, allocation) = land
    (; overland_flow, river_flow, subsurface_flow) = routing
    (; reservoir) = river_flow.boundary_conditions

    dt = tosecond(clock.dt)
    # update lateral inflow for kinematic wave overland flow
    update_lateral_inflow!(
        overland_flow,
        (; soil, allocation, subsurface_flow),
        domain.land.parameters.area,
        config,
        dt,
    )
    # run kinematic wave overland flow
    update!(overland_flow, domain.land, dt)

    # update lateral inflow river flow
    update_lateral_inflow!(
        river_flow,
        (; allocation = river_flow.allocation, runoff, overland_flow, subsurface_flow),
        domain,
        dt,
    )
    # update reservoir boundary conditions external inflow and inflow (overland and
    # subsurface flow), inflow from river flow is added within the river routing scheme
    update_boundary_conditions!(
        reservoir,
        river_flow,
        (; overland_flow, subsurface_flow),
        domain.reservoir.network,
    )
    # update river flow
    update!(river_flow, domain, clock)
    return nothing
end

"""
    surface_routing!(
        model::Model{R}
    ) where {R <: Routing{<:LocalInertialOverlandFlow, <:LocalInertialRiverFlow}}

Run surface routing (land and river) for a model type that contains the routing components
`LocalInertialOverlandFlow` and `LocalInertialRiverFlow` for a single timestep.
"""
function surface_routing!(
    model::Model{R},
) where {R <: Routing{<:LocalInertialOverlandFlow, <:LocalInertialRiverFlow}}
    (; routing, land, domain, clock) = model
    (; soil, runoff) = land
    (; overland_flow, river_flow, subsurface_flow) = routing
    (; reservoir) = river_flow.boundary_conditions

    dt = tosecond(clock.dt)
    update_boundary_conditions!(
        overland_flow,
        (; soil, runoff, subsurface_flow),
        domain,
        dt,
    )
    # update reservoir boundary conditions external inflow and subsurface flow, inflow from
    # river and overland flow is added within the river and overland routing schemes
    update_boundary_conditions!(
        reservoir,
        river_flow,
        subsurface_flow,
        domain.reservoir.network,
    )
    # update overland and river flow
    update!(overland_flow, river_flow, domain, clock)

    return nothing
end