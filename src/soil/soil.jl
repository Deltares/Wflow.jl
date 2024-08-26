
# TODO (part of refactor v1.0): check if current variables and parameters of SBM soil model
# are fine or if for example some parameters should be part of another struct (e.g.
# riverfrac and waterfrac parameters). For ksat_profile parameters it is not required to
# store all parameters (depends on profile).

@get_units @with_kw struct SoilSbmModelVars{T, N}
    # Amount of water in the unsaturated store, per layer [mm]
    ustorelayerdepth::Vector{SVector{N, T}} | "mm"
    # Saturated store [mm]
    satwaterdepth::Vector{T} | "mm"
    # Pseudo-water table depth [mm] (top of the saturated zone)
    zi::Vector{T} | "mm"
    # Number of unsaturated soil layers
    n_unsatlayers::Vector{Int} | "-"
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
    # Water level land [mm]
    waterlevel_land::Vector{T} | "mm"
    # Water level river [mm]
    waterlevel_river::Vector{T} | "mm"
    # Total water storage (excluding floodplain volume, lakes and reservoirs) [mm]
    total_storage::Vector{T} | "mm"
    # Top soil temperature [ᵒC]
    tsoil::Vector{T} | "ᵒC"
end

function soil_sbm_model_vars(n, parameters)
    (;
        soilthickness,
        maxlayers,
        act_thickl,
        sumlayers,
        soilwatercapacity,
        theta_s,
        theta_r,
    ) = parameters
    satwaterdepth = 0.85 .* soilwatercapacity # cold state value for satwaterdepth
    zi = @. max(0.0, soilthickness - satwaterdepth / (theta_s - theta_r))
    thicknesslayers = set_layerthickness.(zi, sumlayers, act_thickl)
    n_unsatlayers = number_of_active_layers.(thicknesslayers)

    vwc = fill(mv, maxlayers, n)
    vwc_perc = fill(mv, maxlayers, n)

    vars = SoilSbmModelVars(;
        ustorelayerdepth = (act_thickl),
        satwaterdepth = satwaterdepth,
        zi = zi,
        n_unsatlayers = n_unsatlayers,
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
        tsoil = fill(Float(10.0), n),
        waterlevel_land = fill(mv, n),
        waterlevel_river = zeros(Float, n), #set to zero to account for cells outside river domain
        total_storage = zeros(Float, n), # Set the total water storage from initialized values
    )
    return vars
end

@get_units @with_kw struct SoilBC{T}
    surface_water_flux::Vector{T}
    potential_transpiration::Vector{T}
    potential_soilevaporation::Vector{T}
end

function soil_model_bc(n)
    bc = SoilBC(;
        surface_water_flux = fill(mv, n),
        potential_transpiration = fill(mv, n),
        potential_soilevaporation = fill(mv, n),
    )
    return bc
end

@get_units @with_kw struct SoilSbmParameters{T, N, M}
    # Maximum number of soil layers
    maxlayers::Int | "-"
    # Number of soil layers
    nlayers::Vector{Int} | "-"
    # Number of soil layers with vertical hydraulic conductivity value `kv`
    nlayers_kv::Vector{Int} | "-"
    # Saturated water content (porosity) [-]
    theta_s::Vector{T} | "-"
    # Residual water content [-]
    theta_r::Vector{T} | "-"
    # Soilwater capacity [mm]
    soilwatercapacity::Vector{T} | "mm"
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
    # Soil temperature smooth factor [-]
    w_soil::Vector{T} | "-"
    # Controls soil infiltration reduction factor when soil is frozen [-]
    cf_soil::Vector{T} | "-"
    # Fraction of river [-]
    riverfrac::Vector{T} | "-"
    # Fraction of open water (excluding rivers) [-]
    waterfrac::Vector{T} | "-"
    # Fraction of compacted area  [-]
    pathfrac::Vector{T} | "-"
    # Controls how roots are linked to water table [-]
    rootdistpar::Vector{T} | "-"
    # Rooting depth [mm]
    rootingdepth::Vector{T} | "mm"
end

abstract type AbstractSoilModel{T} end

@get_units @with_kw struct SoilSbmModel{T} <: AbstractSoilModel{T}
    boundary_conditions::SoilBC{T} | "-"
    parameters::SoilSbmParameters{T} | "-"
    variables::SoilSbmModelVars{T} | "-"
end

function sbm_ksat_profiles(nc, config, inds, soilthickness, maxlayers)
    ksat_profile = get(config.input.vertical, "ksat_profile", "exponential")::String
    n = length(inds)
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
    return z_exp, z_layered, kv, nlayers_kv
end

function initialize_soil_sbm_params(nc, config, riverfrac, inds, dt)
    config_thicknesslayers = get(config.model, "thicknesslayers", Float[])
    if length(config_thicknesslayers) > 0
        thicknesslayers = SVector(Tuple(push!(Float.(config_thicknesslayers), mv)))
        cum_depth_layers = pushfirst(cumsum(thicknesslayers), 0.0)
        maxlayers = length(thicknesslayers) # max number of soil layers
    else
        thicknesslayers = SVector.(soilthickness)
        cum_depth_layers = pushfirst(cumsum(thicknesslayers), 0.0)
        maxlayers = 1
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
    rootingdepth = ncread(
        nc,
        config,
        "vertical.rootingdepth";
        sel = inds,
        defaults = 750.0,
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

    act_thickl = set_layerthickness.(soilthickness, (cum_depth_layers,), (thicknesslayers,))
    sumlayers = @. pushfirst(cumsum(act_thickl), 0.0)
    nlayers = number_of_active_layers.(act_thickl)

    z_exp, z_layered, kv, nlayers_kv =
        sbm_ksat_profiles(nc, config, inds, soilthickness, maxlayers)

    soilwatercapacity = @. soilthickness * (theta_s - theta_r)

    soil_sbm_params = SoilSbmParameters(;
        maxlayers = maxlayers,
        nlayers = nlayers,
        nlayers_kv = nlayers_kv,
        soilwatercapacity = soilwatercapacity,
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
        rootingdepth = rootingdepth,
        cap_hmax = cap_hmax,
        cap_n = cap_n,
        c = svectorscopy(c, Val{maxlayers}()),
        f = f,
        z_exp = z_exp,
        z_layered = z_layered,
        w_soil = w_soil,
        cf_soil = cf_soil,
    )
    return soil_sbm_params
end

function initialize_soil_sbm_model(nc, config, riverfrac, inds, dt)
    n = length(inds)
    params = initialize_soil_sbm_params(nc, config, riverfrac, inds, dt)
    vars = soil_sbm_model_vars(n, params)
    bc = soil_model_bc(n)
    model = SoilSbmModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

function update(model::SoilSbmModel, atmospheric_forcing::AtmosphericForcing, config)
    soilinfreduction = get(config.model, "soilinfreduction", false)::Bool
    modelsnow = get(config.model, "snow", false)::Bool
    transfermethod = get(config.model, "transfermethod", false)::Bool
    ust = get(config.model, "whole_ust_available", false)::Bool # should be removed from optional setting and code?
    ksat_profile = get(config.input.vertical, "ksat_profile", "exponential")::String

    (; potential_evaporation, temperature) = atmospheric_forcing
    (; surface_water_flux, potential_soilevaporation, potential_transpiration) =
        model.boundary_conditions
    v = model.variables
    p = model.parameters

    n = length(potential_evaporation)
    threaded_foreach(1:n; basesize = 250) do i
        ustoredepth = sum(@view v.ustorelayerdepth[i][1:p.nlayers[i]])

        runoff_river = min(1.0, p.riverfrac[i]) * surface_water_flux[i]
        runoff_land = min(1.0, p.waterfrac[i]) * surface_water_flux[i]
        avail_forinfilt = max(surface_water_flux[i] - runoff_river - runoff_land, 0.0)

        rootingdepth = min(p.soilthickness[i] * 0.99, p.rootingdepth[i])

        ae_openw_r = min(
            v.waterlevel_river[i] * p.riverfrac[i],
            p.riverfrac[i] * potential_evaporation[i],
        )
        ae_openw_l = min(
            v.waterlevel_land[i] * p.waterfrac[i],
            p.waterfrac[i] * potential_evaporation[i],
        )

        # Calculate the initial capacity of the unsaturated store
        ustorecapacity = p.soilwatercapacity[i] - v.satwaterdepth[i] - ustoredepth

        if modelsnow
            tsoil = v.tsoil[i] + p.w_soil[i] * (temperature[i] - v.tsoil[i])
        end

        # Calculate the infiltration flux into the soil column
        infiltsoilpath, infiltsoil, infiltpath, soilinf, pathinf, infiltexcess =
            infiltration(
                avail_forinfilt,
                p.pathfrac[i],
                p.cf_soil[i],
                tsoil,
                p.infiltcapsoil[i],
                p.infiltcappath[i],
                ustorecapacity,
                modelsnow,
                soilinfreduction,
            )

        usl = set_layerthickness(v.zi[i], p.sumlayers[i], p.act_thickl[i])
        n_usl = number_of_active_layers(usl)
        z = cumsum(usl)
        usld = v.ustorelayerdepth[i]

        ast = 0.0
        soilevapunsat = 0.0
        if n_usl > 0
            # Using the surface infiltration rate, calculate the flow rate between the
            # different soil layers that contain unsaturated storage assuming gravity
            # based flow only, estimate the gravity based flux rate to the saturated zone
            # (ast) and the updated unsaturated storage for each soil layer.
            if transfermethod && p.maxlayers == 1
                ustorelayerdepth = v.ustorelayerdepth[i][1] + infiltsoilpath
                kv_z = hydraulic_conductivity_at_depth(p, v.zi[i], i, 1, ksat_profile)
                ustorelayerdepth, ast = unsatzone_flow_sbm(
                    ustorelayerdepth,
                    p.soilwatercapacity[i],
                    v.satwaterdepth[i],
                    kv_z,
                    usl[1],
                    p.theta_s[i],
                    p.theta_r[i],
                )
                usld = setindex(usld, ustorelayerdepth, 1)
            else
                for m in 1:n_usl
                    l_sat = usl[m] * (p.theta_s[i] - p.theta_r[i])
                    kv_z = hydraulic_conductivity_at_depth(p, z[m], i, m, ksat_profile)
                    ustorelayerdepth = if m == 1
                        v.ustorelayerdepth[i][m] + infiltsoilpath
                    else
                        v.ustorelayerdepth[i][m] + ast
                    end
                    ustorelayerdepth, ast =
                        unsatzone_flow_layer(ustorelayerdepth, kv_z, l_sat, p.c[i][m])
                    usld = setindex(usld, ustorelayerdepth, m)
                end
            end

            # then evapotranspiration from layers
            # Calculate saturation deficit
            saturationdeficit = p.soilwatercapacity[i] - v.satwaterdepth[i]

            # First calculate the evaporation of unsaturated storage into the
            # atmosphere from the upper layer.
            if p.maxlayers == 1
                soilevapunsat =
                    potential_soilevaporation[i] *
                    min(1.0, saturationdeficit / p.soilwatercapacity[i])
            else
                # In case only the most upper soil layer contains unsaturated storage
                if n_usl == 1
                    # Check if groundwater level lies below the surface
                    soilevapunsat =
                        potential_soilevaporation[i] *
                        min(1.0, usld[1] / (v.zi[i] * (p.theta_s[i] - p.theta_r[i])))
                else
                    # In case first layer contains no saturated storage
                    soilevapunsat =
                        potential_soilevaporation[i] *
                        min(1.0, usld[1] / (usl[1] * ((p.theta_s[i] - p.theta_r[i]))))
                end
            end
            # Ensure that the unsaturated evaporation rate does not exceed the
            # available unsaturated moisture
            soilevapunsat = min(soilevapunsat, usld[1])
            # Update the additional atmospheric demand
            potsoilevap = potential_soilevaporation[i] - soilevapunsat
            usld = setindex(usld, usld[1] - soilevapunsat, 1)
        end
        transfer = ast

        if p.maxlayers == 1
            soilevapsat = 0.0
        else
            if n_usl == 0 || n_usl == 1
                soilevapsat =
                    potsoilevap *
                    min(1.0, (p.act_thickl[i][1] - v.zi[i]) / p.act_thickl[i][1])
                soilevapsat = min(
                    soilevapsat,
                    (p.act_thickl[i][1] - v.zi[i]) * (p.theta_s[i] - p.theta_r[i]),
                )
            else
                soilevapsat = 0.0
            end
        end
        soilevap = soilevapunsat + soilevapsat
        satwaterdepth = v.satwaterdepth[i] - soilevapsat

        # transpiration from saturated store
        wetroots = scurve(v.zi[i], rootingdepth, Float(1.0), p.rootdistpar[i])
        actevapsat = min(potential_transpiration[i] * wetroots, satwaterdepth)
        satwaterdepth = satwaterdepth - actevapsat
        restpottrans = potential_transpiration[i] - actevapsat

        # actual transpiration from ustore
        actevapustore = 0.0
        for k in 1:n_usl
            ustorelayerdepth, actevapustore, restpottrans = acttransp_unsat_sbm(
                rootingdepth,
                usld[k],
                p.sumlayers[i][k],
                restpottrans,
                actevapustore,
                p.c[i][k],
                usl[k],
                p.theta_s[i],
                p.theta_r[i],
                p.hb[i],
                ust,
            )
            usld = setindex(usld, ustorelayerdepth, k)
        end

        # check soil moisture balance per layer
        du = 0.0
        for k in n_usl:-1:1
            du = max(0.0, usld[k] - usl[k] * (p.theta_s[i] - p.theta_r[i]))
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
            ksat = hydraulic_conductivity_at_depth(p, v.zi[i], i, n_usl, ksat_profile)
            ustorecapacity =
                p.soilwatercapacity[i] - satwaterdepth - sum(@view usld[1:p.nlayers[i]])
            maxcapflux = max(0.0, min(ksat, actevapustore, ustorecapacity, satwaterdepth))

            if v.zi[i] > rootingdepth
                capflux =
                    maxcapflux *
                    pow(1.0 - min(v.zi[i], p.cap_hmax[i]) / (p.cap_hmax[i]), p.cap_n[i])
            else
                capflux = 0.0
            end

            netcapflux = capflux
            for k in n_usl:-1:1
                toadd = min(
                    netcapflux,
                    max(usl[k] * (p.theta_s[i] - p.theta_r[i]) - usld[k], 0.0),
                )
                usld = setindex(usld, usld[k] + toadd, k)
                netcapflux = netcapflux - toadd
                actcapflux = actcapflux + toadd
            end
        end
        deepksat = hydraulic_conductivity_at_depth(
            p,
            p.soilthickness[i],
            i,
            p.nlayers[i],
            ksat_profile,
        )
        deeptransfer = min(satwaterdepth, deepksat)
        actleakage = max(0.0, min(p.maxleakage[i], deeptransfer))

        # recharge (mm) for saturated zone
        recharge = (transfer - actcapflux - actleakage - actevapsat - soilevapsat)
        transpiration = actevapsat + actevapustore
        actevap = soilevap + transpiration + ae_openw_r + ae_openw_l

        # update the outputs and states
        v.n_unsatlayers[i] = n_usl
        v.net_runoff_river[i] = runoff_river - ae_openw_r
        v.avail_forinfilt[i] = avail_forinfilt
        v.actinfilt[i] = actinfilt
        v.infiltexcess[i] = infiltexcess
        v.recharge[i] = recharge
        v.transpiration[i] = transpiration
        v.soilevap[i] = soilevap
        v.soilevapsat[i] = soilevapsat
        v.ae_openw_r[i] = ae_openw_r
        v.ae_openw_l[i] = ae_openw_l
        v.runoff_land[i] = runoff_land
        v.runoff_river[i] = runoff_river
        v.actevapsat[i] = actevapsat
        v.actevap[i] = actevap
        v.ae_ustore[i] = actevapustore
        v.ustorelayerdepth[i] = usld
        v.transfer[i] = transfer
        v.actcapflux[i] = actcapflux
        v.actleakage[i] = actleakage
        v.actinfiltsoil[i] = actinfiltsoil
        v.actinfiltpath[i] = actinfiltpath
        v.excesswater[i] = excesswater
        v.excesswatersoil[i] = excesswatersoil
        v.excesswaterpath[i] = excesswaterpath
        v.infiltsoilpath[i] = infiltsoilpath
        v.satwaterdepth[i] = satwaterdepth
    end
end