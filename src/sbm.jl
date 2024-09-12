@get_units @grid_loc @with_kw struct LandHydrologySBM{IM, SM, GM, D, A, T}
    atmospheric_forcing::AtmosphericForcing | "-" | "none"
    vegetation_parameter_set::VegetationParameters | "-" | "none"
    interception::IM | "-" | "none"
    snow::SM | "-" | "none"
    glacier::GM | "-" | "none"
    bucket::SimpleBucketModel | "-" | "none"
    demand::D | "-" | "none"
    allocation::A | "-" | "none"
    # Model time step [s]
    dt::T | "s" | "none"
end

# used by Sediment model
function initialize_canopy(nc, config, inds)
    n = length(inds)
    # if leaf area index climatology provided use sl, swood and kext to calculate cmax, e_r and canopygapfraction
    if haskey(config.input.vertical, "leaf_area_index")
        # TODO confirm if leaf area index climatology is present in the netCDF
        sl = ncread(
            nc,
            config,
            "vertical.specific_leaf";
            optional = false,
            sel = inds,
            type = Float,
        )
        swood = ncread(
            nc,
            config,
            "vertical.storage_wood";
            optional = false,
            sel = inds,
            type = Float,
        )
        kext =
            ncread(nc, config, "vertical.kext"; optional = false, sel = inds, type = Float)
        cmax = fill(mv, n)
        e_r = fill(mv, n)
        canopygapfraction = fill(mv, n)
    else
        sl = fill(mv, n)
        swood = fill(mv, n)
        kext = fill(mv, n)
        # cmax, e_r, canopygapfraction only required when leaf area index climatology not provided
        cmax = ncread(nc, config, "vertical.cmax"; sel = inds, defaults = 1.0, type = Float)
        e_r =
            ncread(nc, config, "vertical.eoverr"; sel = inds, defaults = 0.1, type = Float)
        canopygapfraction = ncread(
            nc,
            config,
            "vertical.canopygapfraction";
            sel = inds,
            defaults = 0.1,
            type = Float,
        )
    end
    return cmax, e_r, canopygapfraction, sl, swood, kext
end

function initialize_land_hydrology_sbm(nc, config, riverfrac, inds)
    dt = Second(config.timestepsecs)
    n = length(inds)

    atmospheric_forcing = initialize_atmospheric_forcing(n)
    vegetation_parameter_set = initialize_vegetation_params(nc, config, inds)
    if dt >= Hour(23)
        interception_model =
            initialize_gash_interception_model(nc, config, inds, vegetation_parameter_set)
    else
        interception_model =
            initialize_rutter_interception_model(vegetation_parameter_set, n)
    end

    modelsnow = get(config.model, "snow", false)::Bool
    if modelsnow
        snow_model = initialize_snow_hbv_model(nc, config, inds, dt)
    else
        snow_model = NoSnowModel()
    end
    modelglacier = get(config.model, "glacier", false)::Bool
    if modelsnow && modelglacier
        glacier_bc = glacier_model_bc(snow_model.variables.snow_storage)
        glacier_model = initialize_glacier_hbv_model(nc, config, inds, dt, glacier_bc)
    else
        glacier_model = NoGlacierModel()
    end
    bucket_model = initialize_simple_bucket_model(
        nc,
        config,
        vegetation_parameter_set,
        riverfrac,
        inds,
        dt,
    )
    @. vegetation_parameter_set.rootingdepth = min(
        bucket_model.parameters.soilthickness * 0.99,
        vegetation_parameter_set.rootingdepth,
    )

    do_water_demand = haskey(config.model, "water_demand")
    allocation = do_water_demand ? initialize_allocation_land(nc, config, inds) : nothing
    demand = do_water_demand ? initialize_water_demand(nc, config, inds, dt) : nothing

    # TODO (part of refactor v1.0): simplify typeof arguments
    lsm = LandHydrologySBM{
        typeof(interception_model),
        typeof(snow_model),
        typeof(glacier_model),
        typeof(demand),
        typeof(allocation),
        Float,
    }(;
        atmospheric_forcing = atmospheric_forcing,
        vegetation_parameter_set = vegetation_parameter_set,
        interception = interception_model,
        snow = snow_model,
        glacier = glacier_model,
        bucket = bucket_model,
        demand = demand,
        allocation = allocation,
        dt = tosecond(dt),
    )
    return lsm
end

function update(lsm::LandHydrologySBM, lateral, network, config)
    do_water_demand = haskey(config.model, "water_demand")::Bool

    update!(lsm.interception, lsm.atmospheric_forcing)

    (; throughfall, stemflow) = lsm.interception.variables
    update_boundary_conditions!(lsm.snow, throughfall .+ stemflow)
    update!(lsm.snow, lsm.atmospheric_forcing)

    # lateral snow transport
    if get(config.model, "masswasting", false)::Bool
        lateral_snow_transport!(
            lsm.snow.variables.snow_storage,
            lsm.snow.variables.snow_water,
            network.land.slope,
            network.land,
        )
    end

    update!(lsm.glacier, lsm.atmospheric_forcing)

    (; potential_evaporation) = lsm.atmospheric_forcing
    update_boundary_conditions!(
        lsm.bucket,
        lsm.interception,
        lsm.snow,
        lsm.glacier,
        potential_evaporation,
    )

    if do_water_demand
        (; potential_transpiration) = lsm.bucket.boundary_conditions
        (; h3_high, h3_low) = lsm.bucket.parameters
        @. lsm.bucket.variables.h3 =
            feddes_h3(h3_high, h3_low, potential_transpiration, lsm.dt)
        update_water_demand(lsm)
        update_water_allocation(lsm, lateral, network)
    end

    update!(lsm.bucket, lsm.demand, lsm.allocation, lsm.atmospheric_forcing, config, lsm.dt)
    @. lsm.bucket.variables.actevap += lsm.interception.variables.interception_flux
    return lsm
end

"""
Update the total water storage per cell at the end of a timestep.

Takes the following parameters:
- sbm:
    The vertical concept (SBM struct)
- river_network:
    The indices of the river cells in relation to the active cells, i.e. model.network.index_river
- area:
    Area of the cells acquired from model.network.land.area
- river_routing:
    The river routing struct, i.e. model.lateral.river
- land_routing:
    The land routing struct, i.e. model.lateral.land
"""
function update_total_water_storage(
    model::LandHydrologySBM,
    river_network,
    area,
    river_routing,
    land_routing,
)
    (; interception, snow, glacier, bucket, demand) = model
    (; total_storage, ustoredepth, satwaterdepth) = bucket.variables
    (; riverfrac) = bucket.parameters

    # Set the total storage to zero
    fill!(total_storage, 0)

    # Burn the river routing values
    for (i, index_river) in enumerate(river_network)
        total_storage[index_river] = (
            (river_routing.h_av[i] * river_routing.width[i] * river_routing.dl[i]) /
            (area[index_river]) * 1000 # Convert to mm
        )
    end

    # Add storage from interception, snow and glacier models
    snow_storage = get_snow_storage(snow)
    snow_water = get_snow_water(snow)
    glacier_store = get_glacier_store(glacier)
    paddy_h = !isnothing(demand) && !isnothing(demand.paddy) ? demand.paddy.h : 0.0
    @. total_storage +=
        snow_storage +
        snow_water +
        glacier_store +
        interception.variables.canopy_storage +
        paddy_h

    # Chunk the data for parallel computing
    n = length(ustoredepth)
    threaded_foreach(1:n; basesize = 1000) do i
        sub_surface = ustoredepth[i] + satwaterdepth[i]
        lateral = (
            land_routing.h_av[i] * (1 - riverfrac[i]) * 1000 # convert to mm
        )

        # Add everything to the total water storage
        total_storage[i] += (sub_surface + lateral)
    end
end

"""
    update_water_demand(lsm::LandHydrologySBM)

Update water demand for `LandHydrologySBM` for a single timestep. Water demand is computed
for sectors `industry`, `domestic` and `livestock`, and `paddy` rice fields and `nonpaddy`
(other crop) fields.

Gross water demand for irrigation `irri_demand_gross` and non-irrigation
`nonirri_demand_gross`, and total gross water demand `total_gross_demand` are updated as
part of `LandHydrologySBM` water allocation (`allocation`).
"""
function update_water_demand(lsm::LandHydrologySBM)
    (; rootingdepth) = lsm.vegetation_parameter_set
    (; nonpaddy, paddy, domestic, industry, livestock) = lsm.demand
    (; hb, theta_s, theta_r, c, sumlayers, act_thickl, pathfrac, infiltcapsoil) =
        lsm.bucket.parameters
    (; h3, n_unsatlayers, zi, ustorelayerdepth, soilinfredu) = lsm.bucket.variables

    n = length(zi)
    for i in 1:n
        industry_dem = update_non_irrigation_demand(industry, i)
        domestic_dem = update_non_irrigation_demand(domestic, i)
        livestock_dem = update_non_irrigation_demand(livestock, i)

        irri_dem_gross = 0.0
        if !isnothing(nonpaddy) && nonpaddy.irrigation_areas[i]
            if nonpaddy.irrigation_trigger[i]
                usl = set_layerthickness(zi[i], sumlayers[i], act_thickl[i])
                for k in 1:n_unsatlayers[i]
                    # compute water demand only for root zone through root fraction per layer
                    rootfrac =
                        min(1.0, (max(0.0, rootingdepth[i] - sumlayers[i][k]) / usl[k]))
                    # vwc_f and vwc_h3 can be precalculated.
                    vwc_fc =
                        vwc_brooks_corey(-100.0, hb[i], theta_s[i], theta_r[i], c[i][k])
                    vwc_h3 = vwc_brooks_corey(h3[i], hb[i], theta_s[i], theta_r[i], c[i][k])
                    depletion =
                        (vwc_fc * usl[k]) - (ustorelayerdepth[i][k] + theta_r[i] * usl[k])
                    depletion *= rootfrac
                    raw = (vwc_fc - vwc_h3) * usl[k] # readily available water
                    raw *= rootfrac

                    # check if maximum irrigation rate has been applied at the previous time step.
                    max_irri_rate_applied =
                        nonpaddy.demand_gross[i] == nonpaddy.maximum_irrigation_rate[i]
                    if depletion >= raw # start irrigation
                        irri_dem_gross += depletion
                        # add depletion to irrigation gross demand when the maximum irrigation rate has been 
                        # applied at the previous time step (to get volumetric water content at field capacity)
                    elseif depletion > 0.0 && max_irri_rate_applied # continue irrigation
                        irri_dem_gross += depletion
                    end
                end
                # limit irrigation demand to infiltration capacity 
                infiltration_capacity =
                    soilinfredu[i] * (1.0 - pathfrac[i]) * infiltcapsoil[i]
                irri_dem_gross = min(irri_dem_gross, infiltration_capacity)
                irri_dem_gross /= nonpaddy.irrigation_efficiency[i]
                # limit irrigation demand to the maximum irrigation rate
                irri_dem_gross = min(irri_dem_gross, nonpaddy.maximum_irrigation_rate[i])
            else
                irri_dem_gross = 0.0
            end
            nonpaddy.demand_gross[i] = irri_dem_gross
        elseif !isnothing(paddy) && paddy.irrigation_areas[i]
            if paddy.irrigation_trigger[i]
                # check if maximum irrigation rate has been applied at the previous time step.
                max_irri_rate_applied =
                    paddy.demand_gross[i] == paddy.maximum_irrigation_rate[i]
                # start irrigation
                if paddy.h[i] < paddy.h_min[i]
                    irr_depth_paddy = paddy.h_opt[i] - paddy.h[i]
                elseif paddy.h[i] < paddy.h_opt[i] && max_irri_rate_applied # continue irrigation
                    irr_depth_paddy = paddy.h_opt[i] - paddy.h[i]
                else
                    irr_depth_paddy = 0.0
                end
                irri_dem_gross += irr_depth_paddy / paddy.irrigation_efficiency[i]
                # limit irrigation demand to the maximum irrigation rate
                irri_dem_gross = min(irri_dem_gross, paddy.maximum_irrigation_rate[i])
            end
            paddy.demand_gross[i] = irri_dem_gross
        end
        # update gross water demands 
        lsm.allocation.irri_demand_gross[i] = irri_dem_gross
        lsm.allocation.nonirri_demand_gross[i] = industry_dem + domestic_dem + livestock_dem
        lsm.allocation.total_gross_demand[i] =
            irri_dem_gross + industry_dem + domestic_dem + livestock_dem
    end
end
