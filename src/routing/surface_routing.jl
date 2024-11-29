"""
    surface_routing!(model)

Run surface routing (land and river) for a single timestep. Kinematic wave for overland flow
and kinematic wave or local inertial model for river flow.
"""
function surface_routing!(model)
    (; vertical, lateral, network, config, clock) = model
    (; soil, runoff, allocation) = vertical
    (; land, river, subsurface) = lateral

    dt = tosecond(clock.dt)
    # update lateral inflow for kinematic wave overland flow
    update_lateral_inflow!(
        land,
        (; soil, allocation, subsurface),
        network.land.area,
        config,
        dt,
    )
    # run kinematic wave overland flow
    update!(land, network.land, dt)

    # update lateral inflow river flow
    update_lateral_inflow!(
        river,
        (; allocation = river.allocation, runoff, land, subsurface),
        network.river.cell_area,
        network.land.area,
        network.river.land_indices,
        dt,
    )
    update_inflow_waterbody!(river, (; land, subsurface), network.river.land_indices)
    update!(river, network, julian_day(clock.time - clock.dt), dt)
    return nothing
end

"""
    surface_routing!(
        model::Model{N,L,V,R,W,T}
    ) where {N,L<:NamedTuple{<:Any,<:Tuple{Any,LocalInertialOverlandFlow,LocalInertialRiverFlow}},V,R,W,T}

Run surface routing (land and river) for a model type that contains the lateral components
`LocalInertialOverlandFlow` and `LocalInertialRiverFlow` for a single timestep.
"""
function surface_routing!(
    model::Model{N, L, V, R, W, T},
) where {
    N,
    L <: NamedTuple{<:Any, <:Tuple{Any, LocalInertialOverlandFlow, LocalInertialRiverFlow}},
    V,
    R,
    W,
    T,
}
    (; lateral, vertical, network, clock) = model
    (; land, river, subsurface) = lateral
    (; soil, runoff) = vertical

    dt = tosecond(clock.dt)
    update_boundary_conditions!(land, (; river, subsurface, soil, runoff), network, dt)

    update!(land, river, network, julian_day(clock.time - clock.dt), dt)

    return nothing
end