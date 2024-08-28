@get_units @with_kw struct LandHydrologySBM{IM, SM, GM, T}
    atmospheric_forcing::AtmosphericForcing | "-"
    vegetation_parameters::VegetationParameters | "-"
    interception::IM | "-"
    snow::SM | "-"
    glacier::GM | "-"
    bucket::SimpleBucketModel | "-"
    # Model time step [s]
    dt::T | "s"
end

function initialize_land_hydrology_sbm(nc, config, riverfrac, inds)
    dt = Second(config.timestepsecs)
    n = length(inds)

    atmospheric_forcing = initialize_atmospheric_forcing(n)
    vegetation_parameters = initialize_vegetation_params(nc, config, inds)
    if dt >= Hour(23)
        interception_model =
            initialize_gash_interception_model(nc, config, inds, vegetation_parameters)
    else
        interception_model = initialize_rutter_interception_model(vegetation_parameters, n)
    end

    modelsnow = get(config.model, "snow", false)::Bool
    if modelsnow
        snow_model = initialize_snow_hbv_model(nc, config, inds, dt)
    else
        snow_model = NoSnowModel()
    end
    modelglacier = get(config.model, "glacier", false)::Bool
    if modelsnow && modelglacier
        glacier_bc = glacier_model_bc(snow_model.variables.snow)
        glacier_model = initialize_glacier_hbv_model(nc, config, inds, dt, glacier_bc)
    else
        glacier_model = NoGlacierModel()
    end
    bucket_model = initialize_simple_bucket_model(
        nc,
        config,
        vegetation_parameters,
        riverfrac,
        inds,
        dt,
    )

    # TODO (part of refactor v1.0): simplify typeof arguments
    lsm = LandHydrologySBM{
        typeof(interception_model),
        typeof(snow_model),
        typeof(glacier_model),
        Float,
    }(;
        atmospheric_forcing = atmospheric_forcing,
        vegetation_parameters = vegetation_parameters,
        interception = interception_model,
        snow = snow_model,
        glacier = glacier_model,
        bucket = bucket_model,
        dt = tosecond(dt),
    )
    return lsm
end

function update!(lsm::LandHydrologySBM, config)
    modelsnow = get(config.model, "snow", false)::Bool
    (; potential_evaporation) = lsm.atmospheric_forcing
    (; canopy_potevap, interception_flux, throughfall, stemflow) =
        lsm.interception.variables

    update!(lsm.interception, lsm.atmospheric_forcing)

    if modelsnow
        (; effective_precip) = lsm.snow.boundary_conditions
        @. effective_precip = throughfall + stemflow
    end
    update!(lsm.snow, lsm.atmospheric_forcing)

    update!(lsm.glacier, lsm.atmospheric_forcing)

    (; potential_transpiration, surface_water_flux, potential_soilevaporation) =
        lsm.bucket.boundary_conditions
    @. potential_transpiration = max(0.0, canopy_potevap - interception_flux)
    snow_runoff = get_runoff(lsm.snow)
    glacier_melt = get_glacier_melt(lsm.glacier)
    @. surface_water_flux = throughfall + stemflow + snow_runoff + glacier_melt
    (; riverfrac, waterfrac) = lsm.bucket.parameters
    (; canopygapfraction) = lsm.interception.parameters.vegetation_parameters
    glacier_frac = get_glacier_frac(lsm.glacier)
    @. potential_soilevaporation =
        max(canopygapfraction - riverfrac - waterfrac - glacier_frac, 0.0) *
        potential_evaporation

    update!(lsm.bucket, lsm.atmospheric_forcing, config)
    @. lsm.bucket.variables.actevap += interception_flux
end

"""
Update the total water storage per cell at the end of a timestep.

Takes the following parameters:
- sbm:
    The vertical concept (SBM struct)
- river_network:
    The indices of the river cells in relation to the active cells, i.e. model.network.index_river
- cell_xsize:
    Size in X direction of the cells acquired from model.network.land.xl
- cell_ysize:
    Size in Y direction of the cells acquired from model.network.land.yl
- river_routing:
    The river routing struct, i.e. model.lateral.river
- land_routing:
    The land routing struct, i.e. model.lateral.land
"""
function update_total_water_storage(
    model::LandHydrologySBM,
    river_network,
    cell_xsize,
    cell_ysize,
    river_routing,
    land_routing,
)
    (; interception, snow, glacier, bucket) = model
    (; total_storage, ustoredepth, satwaterdepth) = bucket.variables
    (; riverfrac) = bucket.parameters
    # Get length active river cells
    nriv = length(river_network)

    # Set the total storage to zero
    fill!(total_storage, 0)

    # Burn the river routing values
    total_storage[river_network] = (
        (
            river_routing.h_av[1:nriv] .* river_routing.width[1:nriv] .*
            river_routing.dl[1:nriv]
        ) ./ (cell_xsize[river_network] .* cell_ysize[river_network]) * 1000 # Convert to mm
    )

    # Add storage from interception, snow and glacier models
    snow_storage = get_snow_storage(snow)
    snow_water = get_snow_water(snow)
    glacier_store = get_glacier_store(glacier)
    @. total_storage =
        snow_storage + snow_water + glacier_store + interception.variables.canopy_storage

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
