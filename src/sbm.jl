@get_units @with_kw struct SBM{IM, SM, GM, T}
    atmospheric_forcing::AtmosphericForcing | "-"
    veg_param_set::VegetationParameters | "-"
    interception_model::IM | "-"
    snow_model::SM | "-"
    glacier_model::GM | "-"
    soil_model::SoilSbmModel | "-"
    # Model time step [s]
    dt::T | "s"
end

function initialize_sbm(nc, config, riverfrac, inds)
    dt = Second(config.timestepsecs)
    n = length(inds)

    atmospheric_forcing = initialize_atmospheric_forcing(n)
    veg_param_set = initialize_vegetation_params(nc, config, inds)
    if dt >= Hour(23)
        interception_model =
            initialize_gash_interception_model(nc, config, inds, veg_param_set)
    else
        interception_model = initialize_rutter_interception_model(veg_param_set, n)
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
    soil_model = initialize_soil_sbm_model(nc, config, riverfrac, inds, dt)

    # TODO (part of refactor v1.0): simplify typeof arguments
    sbm =
        SBM{typeof(interception_model), typeof(snow_model), typeof(glacier_model), Float}(;
            atmospheric_forcing = atmospheric_forcing,
            veg_param_set = veg_param_set,
            interception_model = interception_model,
            snow_model = snow_model,
            glacier_model = glacier_model,
            soil_model = soil_model,
            dt = tosecond(dt),
        )
    return sbm
end

function update_until_snow(sbm::SBM, config)
    modelsnow = get(config.model, "snow", false)::Bool
    (; potential_evaporation) = sbm.atmospheric_forcing
    (; canopy_potevap, interception, throughfall, stemflow) =
        sbm.interception_model.variables

    update(sbm.interception_model, sbm.atmospheric_forcing)

    if modelsnow
        (; effective_precip) = sbm.snow_model.boundary_conditions
        @. effective_precip = throughfall + stemflow
    end
    update(sbm.snow_model, sbm.atmospheric_forcing)

    update(sbm.glacier_model, sbm.atmospheric_forcing)

    (; potential_transpiration, surface_water_flux, potential_soilevaporation) =
        sbm.soil_model.boundary_conditions
    @. potential_transpiration = max(0.0, canopy_potevap - interception)
    snow_runoff = runoff(sbm.snow_model)
    glacier_melt = glaciermelt(sbm.glacier_model)
    @. surface_water_flux = throughfall + stemflow + snow_runoff + glacier_melt
    (; riverfrac, waterfrac) = sbm.soil_model.parameters
    (; canopygapfraction) = sbm.interception_model.parameters.veg_param_set
    glacier_frac = glacierfrac(sbm.glacier_model)
    @. potential_soilevaporation =
        max(canopygapfraction - riverfrac - waterfrac - glacier_frac, 0.0) *
        potential_evaporation

    update(sbm.soil_model, sbm.atmospheric_forcing, config)
    @. sbm.soil_model.variables.actevap += interception
end

function update_after_subsurfaceflow(sbm::SBM, zi, exfiltsatwater)
    threaded_foreach(1:(sbm.n); basesize = 1000) do i
        usl, n_usl = set_layerthickness(zi[i], sbm.sumlayers[i], sbm.act_thickl[i])
        # exfiltration from ustore
        usld = sbm.ustorelayerdepth[i]
        exfiltustore = 0.0
        for k in sbm.n_unsatlayers[i]:-1:1
            if k <= n_usl
                exfiltustore = max(0, usld[k] - usl[k] * (sbm.theta_s[i] - sbm.theta_r[i]))
            else
                exfiltustore = usld[k]
            end
            usld = setindex(usld, usld[k] - exfiltustore, k)
            if k > 1
                usld = setindex(usld, usld[k - 1] + exfiltustore, k - 1)
            end
        end

        ustoredepth = sum(@view usld[1:n_usl])

        runoff =
            exfiltustore +
            exfiltsatwater[i] +
            sbm.excesswater[i] +
            sbm.runoff_land[i] +
            sbm.infiltexcess[i]

        # volumetric water content per soil layer and root zone
        vwc = sbm.vwc[i]
        vwc_perc = sbm.vwc_perc[i]
        for k in 1:sbm.nlayers[i]
            if k <= n_usl
                vwc = setindex(
                    vwc,
                    (
                        usld[k] +
                        (sbm.act_thickl[i][k] - usl[k]) * (sbm.theta_s[i] - sbm.theta_r[i])
                    ) / sbm.act_thickl[i][k] + sbm.theta_r[i],
                    k,
                )
            else
                vwc = setindex(vwc, sbm.theta_s[i], k)
            end
            vwc_perc = setindex(vwc_perc, (vwc[k] / sbm.theta_s[i]) * 100.0, k)
        end

        rootstore_unsat = 0
        for k in 1:n_usl
            rootstore_unsat =
                rootstore_unsat +
                min(1.0, (max(0.0, sbm.rootingdepth[i] - sbm.sumlayers[i][k]) / usl[k])) *
                usld[k]
        end

        rootstore_sat =
            max(0.0, sbm.rootingdepth[i] - zi[i]) * (sbm.theta_s[i] - sbm.theta_r[i])
        rootstore = rootstore_sat + rootstore_unsat
        vwc_root = rootstore / sbm.rootingdepth[i] + sbm.theta_r[i]
        vwc_percroot = (vwc_root / sbm.theta_s[i]) * 100.0

        satwaterdepth = (sbm.soilthickness[i] - zi[i]) * (sbm.theta_s[i] - sbm.theta_r[i])

        # update the outputs and states
        sbm.n_unsatlayers[i] = n_usl
        sbm.ustorelayerdepth[i] = usld
        sbm.ustoredepth[i] = ustoredepth
        sbm.satwaterdepth[i] = satwaterdepth
        sbm.exfiltsatwater[i] = exfiltsatwater[i]
        sbm.exfiltustore[i] = exfiltustore
        sbm.runoff[i] = runoff
        sbm.net_runoff[i] = runoff - sbm.ae_openw_l[i]
        sbm.vwc[i] = vwc
        sbm.vwc_perc[i] = vwc_perc
        sbm.rootstore[i] = rootstore
        sbm.vwc_root[i] = vwc_root
        sbm.vwc_percroot[i] = vwc_percroot
        sbm.zi[i] = zi[i]
    end
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
    sbm::SBM,
    river_network,
    cell_xsize,
    cell_ysize,
    river_routing,
    land_routing,
)
    # Get length active river cells
    nriv = length(river_network)

    # Set the total storage to zero
    fill!(sbm.total_storage, 0)

    # Burn the river routing values
    sbm.total_storage[river_network] = (
        (
            river_routing.h_av[1:nriv] .* river_routing.width[1:nriv] .*
            river_routing.dl[1:nriv]
        ) ./ (cell_xsize[river_network] .* cell_ysize[river_network]) * 1000 # Convert to mm
    )

    # Chunk the data for parallel computing
    threaded_foreach(1:(sbm.n); basesize = 1000) do i

        # Cumulate per vertical type
        # Maybe re-categorize in the future
        surface = (
            sbm.glacierstore[i] * sbm.glacierfrac[i] +
            sbm.snow[i] +
            sbm.snowwater[i] +
            sbm.canopystorage[i]
        )
        sub_surface = sbm.ustoredepth[i] + sbm.satwaterdepth[i]
        lateral = (
            land_routing.h_av[i] * (1 - sbm.riverfrac[i]) * 1000 # convert to mm
        )

        # Add everything to the total water storage
        sbm.total_storage[i] += (surface + sub_surface + lateral)
    end
end
