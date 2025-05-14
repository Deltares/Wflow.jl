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
        snow = get(config.model, "snow__flag", false)::Bool,
        gravitational_snow_transport = get(
            config.model,
            "snow_gravitional_transport__flag",
            false,
        )::Bool,
        glacier = get(config.model, "glacier__flag", false)::Bool,
        lakes = get(config.model, "lake__flag", false)::Bool,
        reservoirs = get(config.model, "reservoir__flag", false)::Bool,
        pits = get(config.model, "pit__flag", false)::Bool,
        water_demand = haskey(config.model, "water_demand"),
        drains = get(config.model, "drain__flag", false)::Bool,
        kh_profile_type = get(
            config.model,
            "saturated_hydraulic_conductivity_profile",
            "exponential",
        )::String,
        min_streamorder_river = get(config.model, "river_streamorder__min_count", 6),
        min_streamorder_land = get(config.model, "land_streamorder__min_count", 5),
    )

    @info "General model settings" modelsettings[keys(modelsettings)[1:8]]...

    routing_types = get_routing_types(config)
    domain = Domain(dataset, config, modelsettings, routing_types)

    land_hydrology = LandHydrologySBM(dataset, config, domain.land)
    routing = Routing(dataset, config, domain, land_hydrology.soil, routing_types, type)

    (; maxlayers) = land_hydrology.soil.parameters
    modelmap = (land = land_hydrology, routing)
    writer = prepare_writer(
        config,
        modelmap,
        domain,
        dataset;
        extra_dim = (name = "layer", value = Float.(1:(maxlayers))),
    )
    close(dataset)

    model =
        Model(config, domain, routing, land_hydrology, clock, reader, writer, SbmModel())

    set_states!(model)

    @info "Initialized model"
    return model
end

"update the `sbm` model type for a single timestep"
function update!(model::AbstractModel{<:SbmModel})
    (; routing, land, domain, clock, config) = model
    dt = Float(tosecond(clock.dt))
    do_water_demand = haskey(config.model, "water_demand")
    (; kv_profile) = land.soil.parameters

    update_until_recharge!(model)
    # exchange of recharge between SBM soil model and subsurface flow domain
    routing.subsurface_flow.boundary_conditions.recharge .=
        land.soil.variables.recharge ./ Float(1000.0)
    if do_water_demand
        @. routing.subsurface_flow.boundary_conditions.recharge -=
            land.allocation.variables.act_groundwater_abst / Float(1000.0)
    end
    routing.subsurface_flow.boundary_conditions.recharge .*=
        domain.land.parameters.flow_width
    routing.subsurface_flow.variables.zi .= land.soil.variables.zi ./ Float(1000.0)
    # update lateral subsurface flow domain (kinematic wave)
    kh_layered_profile!(land.soil, routing.subsurface_flow, kv_profile, dt)
    update!(routing.subsurface_flow, domain.land, Float(clock.dt / BASETIMESTEP))
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
    (; routing, land, domain, clock, config) = model
    dt = Float(tosecond(clock.dt))
    update!(land, routing, domain, config, dt)
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
"""
function update_total_water_storage!(model::AbstractModel{<:SbmModel})
    (; routing, land, domain) = model

    update_total_water_storage!(land, domain, routing)
    return nothing
end

function set_states!(model::AbstractModel{<:Union{SbmModel, SbmGwfModel}})
    (; routing, land, domain, config) = model
    land_v = routing.overland_flow.variables
    river_v = routing.river_flow.variables

    cold_start = get(config.model, "cold_start__flag", true)::Bool
    routing_options = ("kinematic-wave", "local-inertial")
    land_routing =
        get_options(config.model, "land_routing", routing_options, "kinematic-wave")::String
    do_lakes = get(config.model, "lake__flag", false)::Bool
    floodplain_1d = get(config.model, "floodplain_1d__flag", false)::Bool

    # read and set states in model object if cold_start=false
    if cold_start == false
        nriv = Int(length(domain.river.network.indices))
        instate_path = input_path(config, config.state.path_input)
        @info "Set initial conditions from state file `$instate_path`."
        set_states!(instate_path, model; type = Float, dimname = :layer)
        # update zi for SBM soil model
        zi =
            max.(
                Float(0.0),
                land.soil.parameters.soilthickness .-
                land.soil.variables.satwaterdepth ./
                (land.soil.parameters.theta_s .- land.soil.parameters.theta_r),
            )
        land.soil.variables.zi .= zi
        if land_routing == "kinematic-wave"
            (; surface_flow_width, flow_length) = domain.land.parameters
            # make sure land cells with zero flow width are set to zero q and h
            for i in eachindex(surface_flow_width)
                if surface_flow_width[i] <= Float(0.0)
                    land_v.q[i] = Float(0.0)
                    land_v.h[i] = Float(0.0)
                end
            end
            land_v.storage .= land_v.h .* surface_flow_width .* flow_length
        elseif land_routing == "local-inertial"
            (; river, x_length, y_length) = domain.land.parameters
            (; flow_width, flow_length) = domain.river.parameters
            (; bankfull_storage) = routing.overland_flow.parameters
            for i in eachindex(routing.overland_flow.storage)
                if river[i]
                    j = domain.land.network.river_indices[i]
                    if land_v.h[i] > Float(0.0)
                        land_v.storage[i] =
                            land_v.h[i] * x_length[i] * y_length[i] + bankfull_storage[j]
                    else
                        land_v.storage[i] = river_v.h[j] * flow_width[j] * flow_length[j]
                    end
                else
                    routing.overland_flow.storage[i] =
                        routing.overland_flow.h[i] * x_length[i] * y_length[i]
                end
            end
        end
        # only set active cells for river (ignore boundary conditions/ghost points)
        (; flow_width, flow_length) = domain.river.parameters
        river_v.storage[1:nriv] .=
            river_v.h[1:nriv] .* flow_width[1:nriv] .* flow_length[1:nriv]

        if floodplain_1d
            initialize_storage!(routing.river_flow, nriv)
        end

        if do_lakes
            # storage must be re-initialized after loading the state with the current
            # waterlevel otherwise the storage will be based on the initial water level
            lakes = routing.river_flow.boundary_conditions.lake
            lakes.variables.storage .= initialize_storage(
                lakes.parameters.storfunc,
                lakes.parameters.area,
                lakes.variables.waterlevel,
                lakes.parameters.sh,
            )
        end
    else
        @info "Set initial conditions from default values."
    end
    return nothing
end
