@get_units @with_kw struct SBM{IM, SM, GM, T, N, M}
    atmospheric_forcing::AtmosphericForcing | "-"
    veg_param_set::VegetationParameters | "-"
    interception_model::IM | "-"
    snow_model::SM | "-"
    glacier_model::GM | "-"
    # Model time step [s]
    dt::T | "s"
    # Maximum number of soil layers
    maxlayers::Int | "-"
    # number of cells
    n::Int | "-"
    # Number of soil layers
    nlayers::Vector{Int} | "-"
    # Number of unsaturated soil layers
    n_unsatlayers::Vector{Int} | "-"
    # Number of soil layers with vertical hydraulic conductivity value `kv`
    nlayers_kv::Vector{Int} | "-"
    # Fraction of river [-]
    riverfrac::Vector{T} | "-"
    # Saturated water content (porosity) [-]
    theta_s::Vector{T} | "-"
    # Residual water content [-]
    theta_r::Vector{T} | "-"
    # Vertical hydraulic conductivity [mm Δt⁻¹] at soil surface
    kv_0::Vector{T}
    # Vertical hydraulic conductivity [mm Δt⁻¹] per soil layer
    kv::Vector{SVector{N, T}} | "-"
    # Muliplication factor [-] applied to kv_z (vertical flow)
    kvfrac::Vector{SVector{N, T}} | "-"
    # Air entry pressure [cm] of soil (Brooks-Corey)
    hb::Vector{T} | "cm"
    # Soil thickness [mm]
    soilthickness::Vector{T} | "mm"
    # Thickness of soil layers [mm]
    act_thickl::Vector{SVector{N, T}} | "mm"
    # Cumulative sum of soil layers [mm], starting at soil surface (0)
    sumlayers::Vector{SVector{M, T}} | "mm"
    # Infiltration capacity of the compacted areas [mm Δt⁻¹]
    infiltcappath::Vector{T}
    # Soil infiltration capacity [mm Δt⁻¹]
    infiltcapsoil::Vector{T}
    # Maximum leakage [mm Δt⁻¹] from saturated zone
    maxleakage::Vector{T}
    # Fraction of open water (excluding rivers) [-]
    waterfrac::Vector{T} | "-"
    # Fraction of compacted area  [-]
    pathfrac::Vector{T} | "-"
    # Controls how roots are linked to water table [-]
    rootdistpar::Vector{T} | "-"
    # Parameter [mm] controlling capillary rise
    cap_hmax::Vector{T} | "mm"
    # Coefficient [-] controlling capillary rise
    cap_n::Vector{T} | "-"
    # Brooks-Corey power coefﬁcient [-] for each soil layer
    c::Vector{SVector{N, T}} | "-"
    # A scaling parameter [mm⁻¹] (controls exponential decline of kv_0)
    f::Vector{T} | "mm-1"
    # Depth [mm] from soil surface for which exponential decline of kv_0 is valid
    z_exp::Vector{T} | "mm"
    # Depth [mm] from soil surface for which layered profile is valid
    z_layered::Vector{T} | "mm"
    # Amount of water in the unsaturated store, per layer [mm]
    ustorelayerdepth::Vector{SVector{N, T}} | "mm"
    # Saturated store [mm]
    satwaterdepth::Vector{T} | "mm"
    # Pseudo-water table depth [mm] (top of the saturated zone)
    zi::Vector{T} | "mm"
    # Soilwater capacity [mm]
    soilwatercapacity::Vector{T} | "mm"
    # Potential transpiration (after subtracting interception from potential_evaporation)
    pottrans::Vector{T}
    # Transpiration [mm Δt⁻¹]
    transpiration::Vector{T}
    # Actual evaporation from unsaturated store [mm Δt⁻¹]
    ae_ustore::Vector{T}
    # Soil evaporation from unsaturated and saturated store [mm Δt⁻¹]
    soilevap::Vector{T}
    # Soil evaporation from saturated store [mm Δt⁻¹]
    soilevapsat::Vector{T}
    # Actual capillary rise [mm Δt⁻¹]
    actcapflux::Vector{T}
    # Actual transpiration from saturated store [mm Δt⁻¹]
    actevapsat::Vector{T}
    # Total actual evapotranspiration [mm Δt⁻¹]
    actevap::Vector{T}
    # Runoff from river based on riverfrac [mm Δt⁻¹]
    runoff_river::Vector{T}
    # Runoff from land based on waterfrac [mm Δt⁻¹]
    runoff_land::Vector{T}
    # Actual evaporation from open water (land) [mm Δt⁻¹]
    ae_openw_l::Vector{T}
    # Actual evaporation from river [mm Δt⁻¹]
    ae_openw_r::Vector{T}
    # Net runoff from river [mm Δt⁻¹]
    net_runoff_river::Vector{T}
    # Water available for infiltration [mm Δt⁻¹]
    avail_forinfilt::Vector{T}
    # Actual infiltration into the unsaturated zone [mm Δt⁻¹]
    actinfilt::Vector{T}
    # Actual infiltration non-compacted fraction [mm Δt⁻¹]
    actinfiltsoil::Vector{T}
    # Actual infiltration compacted fraction [mm Δt⁻¹]
    actinfiltpath::Vector{T}
    # Actual infiltration (compacted and the non-compacted areas) [mm Δt⁻¹]
    infiltsoilpath::Vector{T}
    # Infiltration excess water [mm Δt⁻¹]
    infiltexcess::Vector{T}
    # Water that cannot infiltrate due to saturated soil (saturation excess) [mm Δt⁻¹]
    excesswater::Vector{T}
    # Water exfiltrating during saturation excess conditions [mm Δt⁻¹]
    exfiltsatwater::Vector{T}
    # Water exfiltrating from unsaturated store because of change in water table [mm Δt⁻¹]
    exfiltustore::Vector{T}
    # Excess water for non-compacted fraction [mm Δt⁻¹]
    excesswatersoil::Vector{T}
    # Excess water for compacted fraction [mm Δt⁻¹]
    excesswaterpath::Vector{T}
    # Total surface runoff from infiltration and saturation excess (excluding actual open water evaporation) [mm Δt⁻¹]
    runoff::Vector{T}
    # Net surface runoff (surface runoff - actual open water evaporation) [mm Δt⁻¹]
    net_runoff::Vector{T}
    # Volumetric water content [-] per soil layer (including theta_r and saturated zone)
    vwc::Vector{SVector{N, T}} | "-"
    # Volumetric water content [%] per soil layer (including theta_r and saturated zone)
    vwc_perc::Vector{SVector{N, T}} | "%"
    # Root water storage [mm] in unsaturated and saturated zone (excluding theta_r)
    rootstore::Vector{T} | "mm"
    # Volumetric water content [-] in root zone (including theta_r and saturated zone)
    vwc_root::Vector{T} | "-"
    # Volumetric water content [%] in root zone (including theta_r and saturated zone)
    vwc_percroot::Vector{T} | "%"
    # Amount of available water in the unsaturated zone [mm]
    ustoredepth::Vector{T} | "mm"
    # Downward flux from unsaturated to saturated zone [mm Δt⁻¹]
    transfer::Vector{T}
    # Net recharge to saturated store [mm Δt⁻¹]
    recharge::Vector{T}
    # Actual leakage from saturated store [mm Δt⁻¹]
    actleakage::Vector{T}
    # Soil temperature smooth factor [-]
    w_soil::Vector{T} | "-"
    # Controls soil infiltration reduction factor when soil is frozen [-]
    cf_soil::Vector{T} | "-"
    # Top soil temperature [ᵒC]
    tsoil::Vector{T} | "ᵒC"
    # Water level land [mm]
    waterlevel_land::Vector{T} | "mm"
    # Water level river [mm]
    waterlevel_river::Vector{T} | "mm"
    # Total water storage (excluding floodplain volume, lakes and reservoirs) [mm]
    total_storage::Vector{T} | "mm"
end

function initialize_sbm(nc, config, riverfrac, inds)
    dt = Second(config.timestepsecs)
    config_thicknesslayers = get(config.model, "thicknesslayers", Float[])
    ksat_profile = get(config.input.vertical, "ksat_profile", "exponential")::String
    if length(config_thicknesslayers) > 0
        thicknesslayers = SVector(Tuple(push!(Float.(config_thicknesslayers), mv)))
        sumlayers = pushfirst(cumsum(thicknesslayers), 0.0)
        maxlayers = length(thicknesslayers) # max number of soil layers
    else
        maxlayers = 1
    end

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
    w_soil =
        ncread(
            nc,
            config,
            "vertical.w_soil";
            sel = inds,
            defaults = 0.1125,
            type = Float,
        ) .* (dt / basetimestep)
    cf_soil =
        ncread(nc, config, "vertical.cf_soil"; sel = inds, defaults = 0.038, type = Float)
    # soil parameters
    theta_s =
        ncread(nc, config, "vertical.theta_s"; sel = inds, defaults = 0.6, type = Float)
    theta_r =
        ncread(nc, config, "vertical.theta_r"; sel = inds, defaults = 0.01, type = Float)
    kv_0 =
        ncread(nc, config, "vertical.kv_0"; sel = inds, defaults = 3000.0, type = Float) .*
        (dt / basetimestep)
    f = ncread(nc, config, "vertical.f"; sel = inds, defaults = 0.001, type = Float)
    hb = ncread(nc, config, "vertical.hb"; sel = inds, defaults = 10.0, type = Float)
    soilthickness = ncread(
        nc,
        config,
        "vertical.soilthickness";
        sel = inds,
        defaults = 2000.0,
        type = Float,
    )
    infiltcappath =
        ncread(
            nc,
            config,
            "vertical.infiltcappath";
            sel = inds,
            defaults = 10.0,
            type = Float,
        ) .* (dt / basetimestep)
    infiltcapsoil =
        ncread(
            nc,
            config,
            "vertical.infiltcapsoil";
            sel = inds,
            defaults = 100.0,
            type = Float,
        ) .* (dt / basetimestep)
    maxleakage =
        ncread(
            nc,
            config,
            "vertical.maxleakage";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)

    c = ncread(
        nc,
        config,
        "vertical.c";
        sel = inds,
        defaults = 10.0,
        type = Float,
        dimname = :layer,
    )
    if size(c, 1) != maxlayers
        parname = param(config.input.vertical, "c")
        size1 = size(c, 1)
        error("$parname needs a layer dimension of size $maxlayers, but is $size1")
    end
    kvfrac = ncread(
        nc,
        config,
        "vertical.kvfrac";
        sel = inds,
        defaults = 1.0,
        type = Float,
        dimname = :layer,
    )
    if size(kvfrac, 1) != maxlayers
        parname = param(config.input, "vertical.kvfrac")
        size1 = size(kvfrac, 1)
        error("$parname needs a layer dimension of size $maxlayers, but is $size1")
    end

    # fraction open water and compacted area (land cover)
    waterfrac =
        ncread(nc, config, "vertical.waterfrac"; sel = inds, defaults = 0.0, type = Float)
    pathfrac =
        ncread(nc, config, "vertical.pathfrac"; sel = inds, defaults = 0.01, type = Float)

    # vegetation parameters
    rootdistpar = ncread(
        nc,
        config,
        "vertical.rootdistpar";
        sel = inds,
        defaults = -500.0,
        type = Float,
    )
    cap_hmax =
        ncread(nc, config, "vertical.cap_hmax"; sel = inds, defaults = 2000.0, type = Float)
    cap_n = ncread(nc, config, "vertical.cap_n"; sel = inds, defaults = 2.0, type = Float)

    if haskey(config.input.vertical, "et_reftopot")
        @warn string(
            "The `et_reftopot` key in `[input.vertical]` is now called ",
            "`kc`. Please update your TOML file.",
        )
    end

    theta_e = theta_s .- theta_r
    soilwatercapacity = soilthickness .* theta_e
    satwaterdepth = 0.85 .* soilwatercapacity # cold state value for satwaterdepth
    zi = max.(0.0, soilthickness .- satwaterdepth ./ theta_e) # cold state value for zi

    # these are filled in the loop below
    # TODO see if we can replace this approach
    nlayers = zeros(Int, n)
    n_unsatlayers = zeros(Int, n)
    act_thickl = zeros(Float, maxlayers, n)
    s_layers = zeros(Float, maxlayers + 1, n)

    for i in 1:n
        if length(config_thicknesslayers) > 0
            act_thickl_, nlayers_ =
                set_layerthickness(soilthickness[i], sumlayers, thicknesslayers)
            s_layers_ = pushfirst(cumsum(act_thickl_), 0.0)

            nlayers[i] = nlayers_
            act_thickl[:, i] = act_thickl_
            s_layers[:, i] = s_layers_
            _, n_unsatlayers[i] = set_layerthickness(zi[i], sumlayers, thicknesslayers)
        else
            nlayers[i] = 1
            act_thickl[:, i] = SVector(soilthickness[i])
            s_layers[:, i] = pushfirst(cumsum(SVector(soilthickness[i])), 0.0)
        end
    end
    act_thickl = svectorscopy(act_thickl, Val{maxlayers}())

    # copied to array of sarray below
    vwc = fill(mv, maxlayers, n)
    vwc_perc = fill(mv, maxlayers, n)
    sumlayers = svectorscopy(s_layers, Val{maxlayers + 1}())

    # ksat profiles
    if ksat_profile == "exponential"
        z_exp = soilthickness
        z_layered = fill(mv, n)
        kv = fill(mv, (maxlayers, n))
        nlayers_kv = fill(0, n)
    elseif ksat_profile == "exponential_constant"
        z_exp =
            ncread(nc, config, "vertical.z_exp"; optional = false, sel = inds, type = Float)
        z_layered = fill(mv, n)
        kv = fill(mv, (maxlayers, n))
        nlayers_kv = fill(0, n)
    elseif ksat_profile == "layered" || ksat_profile == "layered_exponential"
        z_exp = fill(mv, n)
        kv =
            ncread(
                nc,
                config,
                "vertical.kv";
                sel = inds,
                defaults = 1000.0,
                type = Float,
                dimname = :layer,
            ) .* (dt / basetimestep)
        if size(kv, 1) != maxlayers
            parname = param(config.input.vertical, "kv")
            size1 = size(kv, 1)
            error("$parname needs a layer dimension of size $maxlayers, but is $size1")
        end
        if ksat_profile == "layered"
            z_layered = soilthickness
            nlayers_kv = nlayers
        else
            z_layered = ncread(
                nc,
                config,
                "vertical.z_layered";
                optional = false,
                sel = inds,
                type = Float,
            )
            nlayers_kv = fill(0, n)
            for i in eachindex(nlayers_kv)
                layers = @view sumlayers[i][2:nlayers[i]]
                _, k = findmin(abs.(z_layered[i] .- layers))
                nlayers_kv[i] = k
                z_layered[i] = layers[k]
            end
        end
    else
        error("""An unknown "ksat_profile" is specified in the TOML file ($ksat_profile).
              This should be "exponential", "exponential_constant", "layered" or
              "layered_exponential".
              """)
    end

    # TODO (part of refactor v1.0): simplify typeof arguments
    sbm = SBM{
        typeof(interception_model),
        typeof(snow_model),
        typeof(glacier_model),
        Float,
        maxlayers,
        maxlayers + 1,
    }(;
        atmospheric_forcing = atmospheric_forcing,
        veg_param_set = veg_param_set,
        interception_model = interception_model,
        snow_model = snow_model,
        glacier_model = glacier_model,
        dt = tosecond(dt),
        maxlayers = maxlayers,
        n = n,
        nlayers = nlayers,
        n_unsatlayers = n_unsatlayers,
        nlayers_kv = nlayers_kv,
        riverfrac = riverfrac,
        theta_s = theta_s,
        theta_r = theta_r,
        kv_0 = kv_0,
        kv = svectorscopy(kv, Val{maxlayers}()),
        kvfrac = svectorscopy(kvfrac, Val{maxlayers}()),
        hb = hb,
        soilthickness = soilthickness,
        act_thickl = act_thickl,
        sumlayers = sumlayers,
        infiltcappath = infiltcappath,
        infiltcapsoil = infiltcapsoil,
        maxleakage = maxleakage,
        waterfrac = max.(waterfrac .- riverfrac, Float(0.0)),
        pathfrac = pathfrac,
        rootdistpar = rootdistpar,
        cap_hmax = cap_hmax,
        cap_n = cap_n,
        c = svectorscopy(c, Val{maxlayers}()),
        f = f,
        z_exp = z_exp,
        z_layered = z_layered,
        ustorelayerdepth = zero(act_thickl),
        satwaterdepth = satwaterdepth,
        zi = zi,
        soilwatercapacity = soilwatercapacity,
        pottrans = fill(mv, n),
        transpiration = fill(mv, n),
        ae_ustore = fill(mv, n),
        soilevap = fill(mv, n),
        soilevapsat = fill(mv, n),
        actcapflux = fill(mv, n),
        actevapsat = fill(mv, n),
        actevap = fill(mv, n),
        runoff_river = fill(mv, n),
        runoff_land = fill(mv, n),
        ae_openw_l = fill(mv, n),
        ae_openw_r = fill(mv, n),
        avail_forinfilt = fill(mv, n),
        actinfilt = fill(mv, n),
        actinfiltsoil = fill(mv, n),
        actinfiltpath = fill(mv, n),
        infiltsoilpath = fill(mv, n),
        infiltexcess = fill(mv, n),
        excesswater = fill(mv, n),
        exfiltsatwater = fill(mv, n),
        exfiltustore = fill(mv, n),
        excesswatersoil = fill(mv, n),
        excesswaterpath = fill(mv, n),
        runoff = fill(mv, n),
        net_runoff = fill(mv, n),
        net_runoff_river = fill(mv, n),
        vwc = svectorscopy(vwc, Val{maxlayers}()),
        vwc_perc = svectorscopy(vwc_perc, Val{maxlayers}()),
        rootstore = fill(mv, n),
        vwc_root = fill(mv, n),
        vwc_percroot = fill(mv, n),
        ustoredepth = fill(mv, n),
        transfer = fill(mv, n),
        recharge = fill(mv, n),
        actleakage = fill(mv, n),
        w_soil = w_soil,
        cf_soil = cf_soil,
        tsoil = fill(Float(10.0), n),
        waterlevel_land = fill(mv, n),
        waterlevel_river = zeros(Float, n), #set to zero to account for cells outside river domain
        total_storage = zeros(Float, n), # Set the total water storage from initialized values
    )

    return sbm
end

function update_until_snow(sbm::SBM, config)
    (; canopy_potevap, interception, throughfall, stemflow) =
        sbm.interception_model.variables
    (; effective_precip) = sbm.snow_model.boundary_conditions

    update(sbm.interception_model, sbm.atmospheric_forcing)
    @. sbm.pottrans = max(0.0, canopy_potevap - interception)

    @. effective_precip = throughfall + stemflow
    update(sbm.snow_model, sbm.atmospheric_forcing)

    update(sbm.glacier_model, sbm.atmospheric_forcing)
end

function update_until_recharge(sbm::SBM, config)

    # start dummy variables (should be generated from model reader and from Config.jl TOML)
    soilinfreduction = get(config.model, "soilinfreduction", false)::Bool
    modelglacier = get(config.model, "glacier", false)::Bool
    modelsnow = get(config.model, "snow", false)::Bool
    transfermethod = get(config.model, "transfermethod", false)::Bool
    ust = get(config.model, "whole_ust_available", false)::Bool # should be removed from optional setting and code?
    ksat_profile = get(config.input.vertical, "ksat_profile", "exponential")::String

    threaded_foreach(1:(sbm.n); basesize = 250) do i
        if modelsnow
            rainfallplusmelt = sbm.rainfallplusmelt[i]
            if modelglacier
                # Run Glacier module and add the snowpack on-top of it.
                # Estimate the fraction of snow turned into ice (HBV-light).
                # Estimate glacier melt.

                snow, _, glacierstore, glaciermelt = glacier_hbv(
                    sbm.glacierfrac[i],
                    sbm.glacierstore[i],
                    sbm.snow[i],
                    sbm.temperature[i],
                    sbm.g_tt[i],
                    sbm.g_cfmax[i],
                    sbm.g_sifrac[i],
                    Second(sbm.dt),
                )
                # Convert to mm per grid cell and add to snowmelt
                glaciermelt = glaciermelt * sbm.glacierfrac[i]
                rainfallplusmelt = rainfallplusmelt + glaciermelt
            end
        else
            rainfallplusmelt = sbm.stemflow[i] + sbm.throughfall[i]
        end

        avail_forinfilt = rainfallplusmelt
        ustoredepth = sum(@view sbm.ustorelayerdepth[i][1:sbm.nlayers[i]])

        runoff_river = min(1.0, sbm.riverfrac[i]) * avail_forinfilt
        runoff_land = min(1.0, sbm.waterfrac[i]) * avail_forinfilt
        avail_forinfilt = max(avail_forinfilt - runoff_river - runoff_land, 0.0)

        rootingdepth = min(sbm.soilthickness[i] * 0.99, sbm.rootingdepth[i])

        ae_openw_r = min(
            sbm.waterlevel_river[i] * sbm.riverfrac[i],
            sbm.riverfrac[i] * sbm.potential_evaporation[i],
        )
        ae_openw_l = min(
            sbm.waterlevel_land[i] * sbm.waterfrac[i],
            sbm.waterfrac[i] * sbm.potential_evaporation[i],
        )

        # evap available for soil evaporation
        soilevap_fraction = max(
            sbm.canopygapfraction[i] - sbm.riverfrac[i] - sbm.waterfrac[i] -
            sbm.glacierfrac[i],
            0.0,
        )
        potsoilevap = soilevap_fraction * sbm.potential_evaporation[i]

        # Calculate the initial capacity of the unsaturated store
        ustorecapacity = sbm.soilwatercapacity[i] - sbm.satwaterdepth[i] - ustoredepth

        # Calculate the infiltration flux into the soil column
        infiltsoilpath, infiltsoil, infiltpath, soilinf, pathinf, infiltexcess =
            infiltration(
                avail_forinfilt,
                sbm.pathfrac[i],
                sbm.cf_soil[i],
                sbm.tsoil[i],
                sbm.infiltcapsoil[i],
                sbm.infiltcappath[i],
                ustorecapacity,
                modelsnow,
                soilinfreduction,
            )

        usl, n_usl = set_layerthickness(sbm.zi[i], sbm.sumlayers[i], sbm.act_thickl[i])
        z = cumsum(usl)
        usld = sbm.ustorelayerdepth[i]

        ast = 0.0
        soilevapunsat = 0.0
        if n_usl > 0
            # Using the surface infiltration rate, calculate the flow rate between the
            # different soil layers that contain unsaturated storage assuming gravity
            # based flow only, estimate the gravity based flux rate to the saturated zone
            # (ast) and the updated unsaturated storage for each soil layer.
            if transfermethod && sbm.maxlayers == 1
                ustorelayerdepth = sbm.ustorelayerdepth[i][1] + infiltsoilpath
                kv_z = hydraulic_conductivity_at_depth(sbm, sbm.zi[i], i, 1, ksat_profile)
                ustorelayerdepth, ast = unsatzone_flow_sbm(
                    ustorelayerdepth,
                    sbm.soilwatercapacity[i],
                    sbm.satwaterdepth[i],
                    kv_z,
                    usl[1],
                    sbm.theta_s[i],
                    sbm.theta_r[i],
                )
                usld = setindex(usld, ustorelayerdepth, 1)
            else
                for m in 1:n_usl
                    l_sat = usl[m] * (sbm.theta_s[i] - sbm.theta_r[i])
                    kv_z = hydraulic_conductivity_at_depth(sbm, z[m], i, m, ksat_profile)
                    ustorelayerdepth = if m == 1
                        sbm.ustorelayerdepth[i][m] + infiltsoilpath
                    else
                        sbm.ustorelayerdepth[i][m] + ast
                    end
                    ustorelayerdepth, ast =
                        unsatzone_flow_layer(ustorelayerdepth, kv_z, l_sat, sbm.c[i][m])
                    usld = setindex(usld, ustorelayerdepth, m)
                end
            end

            # then evapotranspiration from layers
            # Calculate saturation deficit
            saturationdeficit = sbm.soilwatercapacity[i] - sbm.satwaterdepth[i]

            # First calculate the evaporation of unsaturated storage into the
            # atmosphere from the upper layer.
            if sbm.maxlayers == 1
                soilevapunsat =
                    potsoilevap * min(1.0, saturationdeficit / sbm.soilwatercapacity[i])
            else
                # In case only the most upper soil layer contains unsaturated storage
                if n_usl == 1
                    # Check if groundwater level lies below the surface
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (sbm.zi[i] * (sbm.theta_s[i] - sbm.theta_r[i])))
                else
                    # In case first layer contains no saturated storage
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (usl[1] * ((sbm.theta_s[i] - sbm.theta_r[i]))))
                end
            end
            # Ensure that the unsaturated evaporation rate does not exceed the
            # available unsaturated moisture
            soilevapunsat = min(soilevapunsat, usld[1])
            # Update the additional atmospheric demand
            potsoilevap = potsoilevap - soilevapunsat
            usld = setindex(usld, usld[1] - soilevapunsat, 1)
        end
        transfer = ast

        if sbm.maxlayers == 1
            soilevapsat = 0.0
        else
            if n_usl == 0 || n_usl == 1
                soilevapsat =
                    potsoilevap *
                    min(1.0, (sbm.act_thickl[i][1] - sbm.zi[i]) / sbm.act_thickl[i][1])
                soilevapsat = min(
                    soilevapsat,
                    (sbm.act_thickl[i][1] - sbm.zi[i]) * (sbm.theta_s[i] - sbm.theta_r[i]),
                )
            else
                soilevapsat = 0.0
            end
        end
        soilevap = soilevapunsat + soilevapsat
        satwaterdepth = sbm.satwaterdepth[i] - soilevapsat

        # transpiration from saturated store
        wetroots = scurve(sbm.zi[i], rootingdepth, Float(1.0), sbm.rootdistpar[i])
        actevapsat = min(sbm.pottrans[i] * wetroots, satwaterdepth)
        satwaterdepth = satwaterdepth - actevapsat
        restpottrans = sbm.pottrans[i] - actevapsat

        # actual transpiration from ustore
        actevapustore = 0.0
        for k in 1:n_usl
            ustorelayerdepth, actevapustore, restpottrans = acttransp_unsat_sbm(
                rootingdepth,
                usld[k],
                sbm.sumlayers[i][k],
                restpottrans,
                actevapustore,
                sbm.c[i][k],
                usl[k],
                sbm.theta_s[i],
                sbm.theta_r[i],
                sbm.hb[i],
                ust,
            )
            usld = setindex(usld, ustorelayerdepth, k)
        end

        # check soil moisture balance per layer
        du = 0.0
        for k in n_usl:-1:1
            du = max(0.0, usld[k] - usl[k] * (sbm.theta_s[i] - sbm.theta_r[i]))
            usld = setindex(usld, usld[k] - du, k)
            if k > 1
                usld = setindex(usld, usld[k - 1] + du, k - 1)
            end
        end

        actinfilt = infiltsoilpath - du
        excesswater = avail_forinfilt - infiltsoilpath - infiltexcess + du

        # Separation between compacted and non compacted areas (correction with the satflow du)
        # This is required for D-Emission/Delwaq
        if infiltsoil + infiltpath > 0.0
            actinfiltsoil = infiltsoil - du * infiltsoil / (infiltpath + infiltsoil)
            actinfiltpath = infiltpath - du * infiltpath / (infiltpath + infiltsoil)
        else
            actinfiltsoil = 0.0
            actinfiltpath = 0.0
        end
        excesswatersoil = max(soilinf - actinfiltsoil, 0.0)
        excesswaterpath = max(pathinf - actinfiltpath, 0.0)

        actcapflux = 0.0
        if n_usl > 0
            ksat = hydraulic_conductivity_at_depth(sbm, sbm.zi[i], i, n_usl, ksat_profile)
            ustorecapacity =
                sbm.soilwatercapacity[i] - satwaterdepth - sum(@view usld[1:sbm.nlayers[i]])
            maxcapflux = max(0.0, min(ksat, actevapustore, ustorecapacity, satwaterdepth))

            if sbm.zi[i] > rootingdepth
                capflux =
                    maxcapflux * pow(
                        1.0 - min(sbm.zi[i], sbm.cap_hmax[i]) / (sbm.cap_hmax[i]),
                        sbm.cap_n[i],
                    )
            else
                capflux = 0.0
            end

            netcapflux = capflux
            for k in n_usl:-1:1
                toadd = min(
                    netcapflux,
                    max(usl[k] * (sbm.theta_s[i] - sbm.theta_r[i]) - usld[k], 0.0),
                )
                usld = setindex(usld, usld[k] + toadd, k)
                netcapflux = netcapflux - toadd
                actcapflux = actcapflux + toadd
            end
        end
        deepksat = hydraulic_conductivity_at_depth(
            sbm,
            sbm.soilthickness[i],
            i,
            sbm.nlayers[i],
            ksat_profile,
        )
        deeptransfer = min(satwaterdepth, deepksat)
        actleakage = max(0.0, min(sbm.maxleakage[i], deeptransfer))

        # recharge (mm) for saturated zone
        recharge = (transfer - actcapflux - actleakage - actevapsat - soilevapsat)
        transpiration = actevapsat + actevapustore
        actevap = soilevap + transpiration + ae_openw_r + ae_openw_l + sbm.interception[i]

        # update the outputs and states
        sbm.n_unsatlayers[i] = n_usl
        sbm.net_runoff_river[i] = runoff_river - ae_openw_r
        sbm.avail_forinfilt[i] = avail_forinfilt
        sbm.actinfilt[i] = actinfilt
        sbm.infiltexcess[i] = infiltexcess
        sbm.recharge[i] = recharge
        sbm.transpiration[i] = transpiration
        sbm.soilevap[i] = soilevap
        sbm.soilevapsat[i] = soilevapsat
        sbm.ae_openw_r[i] = ae_openw_r
        sbm.ae_openw_l[i] = ae_openw_l
        sbm.runoff_land[i] = runoff_land
        sbm.runoff_river[i] = runoff_river
        sbm.actevapsat[i] = actevapsat
        sbm.actevap[i] = actevap
        sbm.ae_ustore[i] = actevapustore
        sbm.ustorelayerdepth[i] = usld
        sbm.transfer[i] = transfer
        sbm.actcapflux[i] = actcapflux
        sbm.actleakage[i] = actleakage
        sbm.actinfiltsoil[i] = actinfiltsoil
        sbm.actinfiltpath[i] = actinfiltpath
        sbm.excesswater[i] = excesswater
        sbm.excesswatersoil[i] = excesswatersoil
        sbm.excesswaterpath[i] = excesswaterpath
        sbm.rainfallplusmelt[i] = rainfallplusmelt
        sbm.infiltsoilpath[i] = infiltsoilpath
        sbm.satwaterdepth[i] = satwaterdepth
        if modelsnow
            if modelglacier
                sbm.snow[i] = snow
                sbm.glacierstore[i] = glacierstore
            end
        end
    end
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
