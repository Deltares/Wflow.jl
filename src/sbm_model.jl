"""
    Model(config::Config, type::SbmModel)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function Model(config::Config, type::SbmModel)
    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)
    dataset = NCDataset(static_path)

    reader = prepare_reader(config)
    clock = Clock(config, reader)

    modelsettings = (;
        snow = get(config.model, "snow", false)::Bool,
        gravitational_snow_transport = get(
            config.model,
            "gravitational_snow_transport",
            false,
        )::Bool,
        glacier = get(config.model, "glacier", false)::Bool,
        lakes = get(config.model, "lakes", false)::Bool,
        reservoirs = get(config.model, "reservoirs", false)::Bool,
        pits = get(config.model, "pits", false)::Bool,
        subsurface_flow = get(config.model, "kinematic-wave_subsurface", true)::Bool,
        water_demand = haskey(config.model, "water_demand"),
        kh_profile_type = get(
            config.model,
            "saturated_hydraulic_conductivity_profile",
            "exponential",
        )::String,
        min_streamorder_river = get(config.model, "min_streamorder_river", 6),
        min_streamorder_land = get(config.model, "min_streamorder_land", 5),
    )

    @info "General model settings" modelsettings[keys(modelsettings)[1:8]]...

    routing_types = get_routing_types(config)
    network = Network(dataset, config, modelsettings, routing_types)
    parameters = SharedParameters(dataset, config, network)

    land_hydrology = LandHydrologySBM(dataset, config, parameters.land, network.land)
    routing =
        Routing(dataset, config, land_hydrology.soil, network, parameters, routing_types)

    (; maxlayers) = land_hydrology.soil.parameters
    model_components = (land = land_hydrology, routing)
    writer = prepare_writer(
        config,
        model_components,
        network,
        dataset;
        extra_dim = (name = "layer", value = Float64.(1:(maxlayers))),
    )
    close(dataset)

    model =
        Model(config, network, routing, land_hydrology, clock, reader, writer, SbmModel())

    set_states!(model)

    @info "Initialized model"
    return model
end

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
    update!(routing.subsurface_flow, network.land, clock.dt / BASETIMESTEP)
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
        set_states!(instate_path, model; type = Float64, dimname = :layer)
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
            land_v.storage .= land_v.h .* land_p.flow_width .* land_p.flow_length
        elseif land_routing == "local-inertial"
            for i in eachindex(routing.overland_flow.storage)
                if land_p.rivercells[i]
                    j = network.land.index_river[i]
                    if land_v.h[i] > 0.0
                        land_v.storage[i] =
                            land_v.h[i] * land_p.xl[i] * land_p.yl[i] +
                            land_p.bankfull_storage[j]
                    else
                        land_v.storage[i] =
                            river_v.h[j] * river_p.flow_width[j] * river_p.flow_length[j]
                    end
                else
                    routing.overland_flow.storage[i] =
                        routing.overland_flow.h[i] *
                        routing.overland_flow.xl[i] *
                        routing.overland_flow.yl[i]
                end
            end
        end
        # only set active cells for river (ignore boundary conditions/ghost points)
        river_v.storage[1:nriv] .=
            river_v.h[1:nriv] .* river_p.flow_width[1:nriv] .* river_p.flow_length[1:nriv]

        if floodplain_1d
            initialize_storage!(routing.river_flow, nriv)
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
