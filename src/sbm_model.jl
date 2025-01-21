"update SBM model for a single timestep"
function update!(model::AbstractModel{<:SbmModel})
    (; routing, land, network, clock, config) = model
    dt = tosecond(clock.dt)
    do_water_demand = haskey(config.model, "water_demand")
    (; kv_profile) = land.soil.parameters

    update_until_recharge!(model)
    # exchange of recharge between SBM soil model and subsurface flow domain
    routing.subsurface_flow.boundary_conditions.recharge .=
        land.soil.variables.recharge ./ 1000.0
    if do_water_demand
        @. routing.subsurface_flow.boundary_conditions.recharge -=
            land.allocation.variables.act_groundwater_abst / 1000.0
    end
    routing.subsurface_flow.boundary_conditions.recharge .*=
        routing.subsurface_flow.parameters.flow_width
    routing.subsurface_flow.variables.zi .= land.soil.variables.zi ./ 1000.0
    # update lateral subsurface flow domain (kinematic wave)
    kh_layered_profile!(land.soil, routing.subsurface_flow, kv_profile, dt)
    update!(routing.subsurface_flow, network.land, clock.dt / basetimestep)
    update_after_subsurfaceflow!(model)
    update_total_water_storage!(model)
    return nothing
end

"""
    update_until_recharge!model::AbstractModel{<:SbmModel})

Update SBM model until recharge for a single timestep. This function is also accessible
through BMI, to couple the SBM model to an external groundwater model.
"""
function update_until_recharge!(model::AbstractModel{<:SbmModel})
    (; routing, land, network, clock, config) = model
    dt = tosecond(clock.dt)
    update!(land, routing, network, config, dt)
    return nothing
end

"""
    update_after_subsurfaceflow!(model::AbstractModel{<:SbmModel})

Update SBM model after subsurface flow for a single timestep. This function is also
accessible through BMI, to couple the SBM model to an external groundwater model.
"""
function update_after_subsurfaceflow!(model::AbstractModel{<:SbmModel})
    (; routing, land) = model
    (; soil, runoff, demand) = land
    (; subsurface_flow) = routing

    # update SBM soil model (runoff, ustorelayerdepth and satwaterdepth)
    update!(soil, (; runoff, demand, subsurface_flow))

    surface_routing!(model)

    return nothing
end

"""
Update of the total water storage at the end of each timestep per model cell.

This is done here at model level.
"""
function update_total_water_storage!(model::AbstractModel{<:SbmModel})
    (; routing, land, network) = model

    # Update the total water storage based on land states
    # TODO Maybe look at routing in the near future
    update_total_water_storage!(
        land,
        network.river.land_indices,
        network.land.area,
        routing.river_flow,
        routing.overland_flow,
    )
    return nothing
end

function set_states!(model::AbstractModel{<:Union{SbmModel, SbmGwfModel}})
    (; routing, land, network, config) = model
    land_v = routing.overland_flow.variables
    land_p = routing.overland_flow.parameters
    river_v = routing.river_flow.variables
    river_p = routing.river_flow.parameters

    reinit = get(config.model, "reinit", true)::Bool
    routing_options = ("kinematic-wave", "local-inertial")
    land_routing =
        get_options(config.model, "land_routing", routing_options, "kinematic-wave")::String
    do_lakes = get(config.model, "lakes", false)::Bool
    floodplain_1d = get(config.model, "floodplain_1d", false)::Bool

    # read and set states in model object if reinit=false
    if reinit == false
        nriv = length(network.river.indices)
        instate_path = input_path(config, config.state.path_input)
        @info "Set initial conditions from state file `$instate_path`."
        set_states!(instate_path, model; type = Float, dimname = :layer)
        # update zi for SBM soil model
        zi =
            max.(
                0.0,
                land.soil.parameters.soilthickness .-
                land.soil.variables.satwaterdepth ./
                (land.soil.parameters.theta_s .- land.soil.parameters.theta_r),
            )
        land.soil.variables.zi .= zi
        if land_routing == "kinematic-wave"
            # make sure land cells with zero flow width are set to zero q and h
            for i in eachindex(land_p.flow_width)
                if land_p.flow_width[i] <= 0.0
                    land_v.q[i] = 0.0
                    land_v.h[i] = 0.0
                end
            end
            land_v.volume .= land_v.h .* land_p.flow_width .* land_p.flow_length
        elseif land_routing == "local-inertial"
            for i in eachindex(routing.overland_flow.volume)
                if land_p.rivercells[i]
                    j = network.land.index_river[i]
                    if land_v.h[i] > 0.0
                        land_v.volume[i] =
                            land_v.h[i] * land_p.xl[i] * land_p.yl[i] +
                            land_p.bankfull_volume[j]
                    else
                        land_v.volume[i] =
                            river_v.h[j] * river_p.flow_width[j] * river_p.flow_length[j]
                    end
                else
                    routing.overland_flow.volume[i] =
                        routing.overland_flow.h[i] *
                        routing.overland_flow.xl[i] *
                        routing.overland_flow.yl[i]
                end
            end
        end
        # only set active cells for river (ignore boundary conditions/ghost points)
        river_v.volume[1:nriv] .=
            river_v.h[1:nriv] .* river_p.flow_width[1:nriv] .* river_p.flow_length[1:nriv]

        if floodplain_1d
            initialize_volume!(routing.river_flow, nriv)
        end

        if do_lakes
            # storage must be re-initialized after loading the state with the current
            # waterlevel otherwise the storage will be based on the initial water level
            lakes = routing.river_flow.boundary_conditions.lake
            lakes.storage .=
                initialize_storage(lakes.storfunc, lakes.area, lakes.waterlevel, lakes.sh)
        end
    else
        @info "Set initial conditions from default values."
    end
    return nothing
end
