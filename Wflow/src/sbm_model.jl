"""
    Model(config::Config, type::SbmModel)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function Model(config::Config, type::SbmModel)
    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)
    dataset = NCDataset(static_path)

    reader = NCReader(config)
    clock = Clock(config, reader)

    @info "General model settings." (;
        snow = config.model.snow__flag,
        gravitational_snow_transport = config.model.snow_gravitational_transport__flag,
        glacier = config.model.glacier__flag,
        reservoirs = config.model.reservoir__flag,
        pits = config.model.pit__flag,
        water_demand = do_water_demand(config),
    )...

    domain = Domain(dataset, config, type)

    land_hydrology = LandHydrologySBM(dataset, config, domain.land)
    routing = Routing(dataset, config, domain, land_hydrology.soil, type)
    mass_balance = HydrologicalMassBalance(domain, config)

    (; maxlayers) = land_hydrology.soil.parameters
    modelmap = (land = land_hydrology, routing, mass_balance)
    writer = Writer(
        config,
        modelmap,
        domain,
        dataset;
        extra_dim = (name = "layer", value = Float64.(1:(maxlayers))),
    )
    close(dataset)

    model = Model(
        config,
        domain,
        routing,
        land_hydrology,
        mass_balance,
        clock,
        reader,
        writer,
        SbmModel(),
    )

    set_states!(model)

    @info "Initialized model"
    return model
end

"update the `sbm` model type for a single timestep"
function update!(model::AbstractModel{<:SbmModel})
    (; routing, land, domain, clock, config) = model
    dt = tosecond(clock.dt)
    (; kv_profile) = land.soil.parameters

    update_until_recharge!(model)
    # exchange of recharge [mm dt⁻¹] between SBM soil model and subsurface flow domain
    routing.subsurface_flow.boundary_conditions.recharge .= land.soil.variables.recharge
    if do_water_demand(config)
        @. routing.subsurface_flow.boundary_conditions.recharge -=
            land.allocation.variables.act_groundwater_abst
    end
    # unit conversions
    routing.subsurface_flow.boundary_conditions.recharge .*=
        domain.land.parameters.flow_width * 0.001 * (tosecond(BASETIMESTEP) / dt)
    routing.subsurface_flow.variables.zi .= land.soil.variables.zi ./ 1000.0
    # update lateral subsurface flow domain (kinematic wave)
    kh_layered_profile!(land.soil, routing.subsurface_flow, kv_profile, dt)
    update!(routing.subsurface_flow, domain.land, clock.dt / BASETIMESTEP)
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
    dt = tosecond(clock.dt)
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

    (; land_routing, cold_start__flag, reservoir__flag, floodplain_1d__flag) = config.model

    # read and set states in model object if cold_start=false
    if !cold_start__flag
        nriv = length(domain.river.network.indices)
        instate_path = input_path(config, config.state.path_input)
        @info "Set initial conditions from state file `$instate_path`."
        set_states!(instate_path, model; type = Float64, dimname = :layer)

        update_diagnostic_vars!(land.soil)

        if land_routing == RoutingType.kinematic_wave
            (; surface_flow_width, flow_length) = domain.land.parameters
            # make sure land cells with zero flow width are set to zero q and h
            for i in eachindex(surface_flow_width)
                if surface_flow_width[i] <= 0.0
                    land_v.q[i] = 0.0
                    land_v.h[i] = 0.0
                end
            end
            land_v.storage .= land_v.h .* surface_flow_width .* flow_length
        elseif land_routing == RoutingType.local_inertial
            (; river_location, x_length, y_length) = domain.land.parameters
            (; flow_width, flow_length) = domain.river.parameters
            (; bankfull_storage) = routing.river_flow.parameters
            for i in eachindex(land_v.storage)
                if river_location[i]
                    j = domain.land.network.river_indices[i]
                    if land_v.h[i] > 0.0
                        land_v.storage[i] =
                            land_v.h[i] * x_length[i] * y_length[i] + bankfull_storage[j]
                    else
                        land_v.storage[i] = river_v.h[j] * flow_width[j] * flow_length[j]
                    end
                else
                    land_v.storage[i] = land_v.h[i] * x_length[i] * y_length[i]
                end
            end
        end
        if config.model.type == ModelType.sbm
            (; zi, storage) = routing.subsurface_flow.variables
            (; theta_s, theta_r, soilthickness) = routing.subsurface_flow.parameters
            @. zi = 0.001 * land.soil.variables.zi # convert from unit [mm] to [m]
            @. storage =
                (theta_s - theta_r) * (soilthickness - zi) * domain.land.parameters.area
        elseif config.model.type == ModelType.sbm_gwf
            (; aquifer) = routing.subsurface_flow
            aquifer.variables.storage .=
                saturated_thickness(aquifer) .* aquifer.parameters.area .*
                storativity(aquifer)
        end
        # only set active cells for river (ignore boundary conditions/ghost points)
        (; flow_width, flow_length) = domain.river.parameters
        river_v.storage[1:nriv] .=
            river_v.h[1:nriv] .* flow_width[1:nriv] .* flow_length[1:nriv]

        if floodplain_1d__flag
            initialize_storage!(routing.river_flow, domain, nriv)
        end

        if reservoir__flag
            # storage must be re-initialized after loading the state with the current
            # waterlevel otherwise the storage will be based on the initial water level
            reservoirs = routing.river_flow.boundary_conditions.reservoir
            reservoirs.variables.storage .= initialize_storage(
                reservoirs.parameters.storfunc,
                reservoirs.parameters.area,
                reservoirs.variables.waterlevel,
                reservoirs.parameters.sh,
            )
        end
    else
        @info "Set initial conditions from default values."
    end
    return nothing
end
