"""
    surface_routing!(model)

Run surface routing (land and river) for a single timestep. Kinematic wave for overland flow
and kinematic wave or local inertial model for river flow.
"""
function surface_routing!(model)
    (; land, routing, domain, config, clock) = model
    (; soil, runoff, allocation) = land
    (; overland_flow, river_flow, subsurface_flow) = routing

    dt = Float(tosecond(clock.dt))
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
    update_inflow_waterbody!(
        river_flow,
        (; overland_flow, subsurface_flow),
        domain.river.network.land_indices,
    )
    update!(river_flow, domain, julian_day(clock.time - clock.dt), dt)
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

    dt = Float(tosecond(clock.dt))
    update_boundary_conditions!(
        overland_flow,
        (; river_flow, subsurface_flow, soil, runoff),
        domain,
        dt,
    )
    update!(overland_flow, river_flow, domain, julian_day(clock.time - clock.dt), dt)
    return nothing
end