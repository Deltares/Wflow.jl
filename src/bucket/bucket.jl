
# TODO (part of refactor v1.0): check if current variables and parameters of SBM soil model
# are fine or if for example some parameters should be part of another struct (e.g.
# riverfrac and waterfrac parameters). For ksat_profile parameters it is not required to
# store all parameters (depends on profile).

@get_units @with_kw struct SimpleBucketModelVars{T, N}
    # Calculated soil water pressure head h3 of the root water uptake reduction function (Feddes) [cm]
    h3::Vector{T} | "cm"
    # Amount of water in the unsaturated store, per layer [mm]
    ustorelayerdepth::Vector{SVector{N, T}} | "mm"
    # Saturated store [mm]
    satwaterdepth::Vector{T} | "mm"
    # Pseudo-water table depth [mm] (top of the saturated zone)
    zi::Vector{T} | "mm"
    # Number of unsaturated soil layers
    n_unsatlayers::Vector{Int} | "-"
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

function simple_bucket_model_vars(n, parameters)
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

    vars = SimpleBucketModelVars(;
        ustorelayerdepth = zero(act_thickl),
        satwaterdepth = satwaterdepth,
        zi = zi,
        n_unsatlayers = n_unsatlayers,
        transpiration = fill(mv, n),
        h3 = fill(mv, n),
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

@get_units @with_kw struct SimpleBucketModelBC{T}
    surface_water_flux::Vector{T}
    potential_transpiration::Vector{T}
    potential_soilevaporation::Vector{T}
end

function simple_bucket_model_bc(n)
    bc = SimpleBucketModelBC(;
        surface_water_flux = fill(mv, n),
        potential_transpiration = fill(mv, n),
        potential_soilevaporation = fill(mv, n),
    )
    return bc
end

@get_units @with_kw struct SimpleBucketModelParameters{T, N, M}
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
    # Fraction of the root length density in each soil layer [-]
    rootfraction::Vector{SVector{N, T}} | "-"
    # Soil water pressure head h1 of the root water uptake reduction function (Feddes) [cm]
    h1::Vector{T} | "cm"
    # Soil water pressure head h2 of the root water uptake reduction function (Feddes) [cm]
    h2::Vector{T} | "cm"
    # Soil water pressure head h3_high of the root water uptake reduction function (Feddes) [cm]
    h3_high::Vector{T} | "cm"
    # Soil water pressure head h3_low of the root water uptake reduction function (Feddes) [cm]
    h3_low::Vector{T} | "cm"
    # Soil water pressure head h4 of the root water uptake reduction function (Feddes) [cm]
    h4::Vector{T} | "cm"
    # Root water uptake reduction at soil water pressure head h1 (0.0 or 1.0) [-]
    alpha_h1::Vector{T} | "-"
    vegetation_parameters::VegetationParameters{T} | "-"
end

abstract type AbstractBucketModel{T} end

@get_units @with_kw struct SimpleBucketModel{T} <: AbstractBucketModel{T}
    boundary_conditions::SimpleBucketModelBC{T} | "-"
    parameters::SimpleBucketModelParameters{T} | "-"
    variables::SimpleBucketModelVars{T} | "-"
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

function initialize_simple_bucket_model_params(
    nc,
    config,
    vegetation_parameters,
    riverfrac,
    inds,
    dt,
)
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
            "vertical.bucket.parameters.w_soil";
            sel = inds,
            defaults = 0.1125,
            type = Float,
        ) .* (dt / basetimestep)
    cf_soil = ncread(
        nc,
        config,
        "vertical.bucket.parameters.cf_soil";
        sel = inds,
        defaults = 0.038,
        type = Float,
    )
    # soil parameters
    theta_s = ncread(
        nc,
        config,
        "vertical.bucket.parameters.theta_s";
        sel = inds,
        defaults = 0.6,
        type = Float,
    )
    theta_r = ncread(
        nc,
        config,
        "vertical.bucket.parameters.theta_r";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )
    kv_0 =
        ncread(
            nc,
            config,
            "vertical.bucket.parameters.kv_0";
            sel = inds,
            defaults = 3000.0,
            type = Float,
        ) .* (dt / basetimestep)
    f = ncread(
        nc,
        config,
        "vertical.bucket.parameters.f";
        sel = inds,
        defaults = 0.001,
        type = Float,
    )
    hb = ncread(
        nc,
        config,
        "vertical.bucket.parameters.hb";
        sel = inds,
        defaults = 10.0,
        type = Float,
    )
    h1 = ncread(
        nc,
        config,
        "vertical.bucket.parameters.h1";
        sel = inds,
        defaults = 0.0,
        type = Float,
    )
    h2 = ncread(
        nc,
        config,
        "vertical.bucket.parameters.h2";
        sel = inds,
        defaults = -100.0,
        type = Float,
    )
    h3_high = ncread(
        nc,
        config,
        "vertical.bucket.parameters.h3_high";
        sel = inds,
        defaults = -400.0,
        type = Float,
    )
    h3_low = ncread(
        nc,
        config,
        "vertical.bucket.parameters.h3_low";
        sel = inds,
        defaults = -1000.0,
        type = Float,
    )
    h4 = ncread(
        nc,
        config,
        "vertical.bucket.parameters.h4";
        sel = inds,
        defaults = -15849.0,
        type = Float,
    )
    alpha_h1 = ncread(
        nc,
        config,
        "vertical.bucket.parameters.alpha_h1";
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    soilthickness = ncread(
        nc,
        config,
        "vertical.bucket.parameters.soilthickness";
        sel = inds,
        defaults = 2000.0,
        type = Float,
    )
    infiltcappath =
        ncread(
            nc,
            config,
            "vertical.bucket.parameters.infiltcappath";
            sel = inds,
            defaults = 10.0,
            type = Float,
        ) .* (dt / basetimestep)
    infiltcapsoil =
        ncread(
            nc,
            config,
            "vertical.bucket.parameters.infiltcapsoil";
            sel = inds,
            defaults = 100.0,
            type = Float,
        ) .* (dt / basetimestep)
    maxleakage =
        ncread(
            nc,
            config,
            "vertical.bucket.parameters.maxleakage";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)

    c = ncread(
        nc,
        config,
        "vertical.parameters.bucket.c";
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
        "vertical.bucket.parameters.kvfrac";
        sel = inds,
        defaults = 1.0,
        type = Float,
        dimname = :layer,
    )
    if size(kvfrac, 1) != maxlayers
        parname = param(config.input, "vertical.bucket.parameters.kvfrac")
        size1 = size(kvfrac, 1)
        error("$parname needs a layer dimension of size $maxlayers, but is $size1")
    end

    # fraction open water and compacted area (land cover)
    waterfrac = ncread(
        nc,
        config,
        "vertical.bucket.parameters.waterfrac";
        sel = inds,
        defaults = 0.0,
        type = Float,
    )
    pathfrac = ncread(
        nc,
        config,
        "vertical.bucket.parameters.pathfrac";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )

    # vegetation parameters
    rootdistpar = ncread(
        nc,
        config,
        "vertical.bucket.parameters.rootdistpar";
        sel = inds,
        defaults = -500.0,
        type = Float,
    )
    cap_hmax = ncread(
        nc,
        config,
        "vertical.bucket.parameters.cap_hmax";
        sel = inds,
        defaults = 2000.0,
        type = Float,
    )
    cap_n = ncread(
        nc,
        config,
        "vertical.bucket.parameters.cap_n";
        sel = inds,
        defaults = 2.0,
        type = Float,
    )

    if haskey(config.input.vertical.bucket, "et_reftopot")
        @warn string(
            "The `et_reftopot` key in `[input.vertical.bucket.parameters]` is now called ",
            "`kc`. Please update your TOML file.",
        )
    end

    act_thickl = set_layerthickness.(soilthickness, (cum_depth_layers,), (thicknesslayers,))
    sumlayers = @. pushfirst(cumsum(act_thickl), 0.0)
    nlayers = number_of_active_layers.(act_thickl)

    if length(config_thicknesslayers) > 0
        # root fraction read from nc file, in case of multiple soil layers and TOML file
        # includes "vertical.rootfraction"
        if haskey(config.input.vertical.bucket.parameters, "rootfraction")
            rootfraction = ncread(
                nc,
                config,
                "vertical.bucket.parameters.rootfraction";
                sel = inds,
                optional = false,
                type = Float,
                dimname = :layer,
            )
        else
            n = length(inds)
            (; rootingdepth) = vegetation_parameters
            # default root fraction in case of multiple soil layers
            rootfraction = zeros(Float, maxlayers, n)
            for i in 1:n
                if rootingdepth[i] > 0.0
                    for k in 1:maxlayers
                        if (rootingdepth[i] - sumlayers[i][k]) >= act_thickl[i][k]
                            rootfraction[k, i] = act_thickl[i][k] / rootingdepth[i]
                        else
                            rootfraction[k, i] =
                                max(rootingdepth[i] - sumlayers[i][k], 0.0) /
                                rootingdepth[i]
                        end
                    end
                end
            end
        end
    else
        # for the case of 1 soil layer
        rootfraction = ones(Float, maxlayers, n)
    end

    z_exp, z_layered, kv, nlayers_kv =
        sbm_ksat_profiles(nc, config, inds, soilthickness, maxlayers)

    soilwatercapacity = @. soilthickness * (theta_s - theta_r)

    sbm_params = SimpleBucketModelParameters(;
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
        h1 = h1,
        h2 = h2,
        h3_high = h3_high,
        h3_low = h3_low,
        h4 = h4,
        alpha_h1 = alpha_h1,
        soilthickness = soilthickness,
        act_thickl = act_thickl,
        sumlayers = sumlayers,
        infiltcappath = infiltcappath,
        infiltcapsoil = infiltcapsoil,
        maxleakage = maxleakage,
        waterfrac = max.(waterfrac .- riverfrac, Float(0.0)),
        pathfrac = pathfrac,
        rootdistpar = rootdistpar,
        rootfraction = svectorscopy(rootfraction, Val{maxlayers}()),
        cap_hmax = cap_hmax,
        cap_n = cap_n,
        c = svectorscopy(c, Val{maxlayers}()),
        f = f,
        z_exp = z_exp,
        z_layered = z_layered,
        w_soil = w_soil,
        cf_soil = cf_soil,
        vegetation_parameters = vegetation_parameters,
    )
    return sbm_params
end

function initialize_simple_bucket_model(
    nc,
    config,
    vegetation_parameters,
    riverfrac,
    inds,
    dt,
)
    n = length(inds)
    params = initialize_simple_bucket_model_params(
        nc,
        config,
        vegetation_parameters,
        riverfrac,
        inds,
        dt,
    )
    vars = simple_bucket_model_vars(n, params)
    bc = simple_bucket_model_bc(n)
    model =
        SimpleBucketModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

#TODO (v1.0): 
# add water demand functionality from Master branch
# check if it's feasible to make the SimpleBucketModel (SBM) modular:
#   - soil evaporation
#   - transpiration (root water uptake Feddes)
#   - unsaturated flow (soil_model)
#   - capillary rise
#   - deep transfer 

function update!(
    model::SimpleBucketModel,
    demand,
    allocation,
    atmospheric_forcing::AtmosphericForcing,
    config,
    dt,
)
    soilinfreduction = get(config.model, "soilinfreduction", false)::Bool
    modelsnow = get(config.model, "snow", false)::Bool
    transfermethod = get(config.model, "transfermethod", false)::Bool
    ust = get(config.model, "whole_ust_available", false)::Bool # should be removed from optional setting and code?
    ksat_profile = get(config.input.vertical, "ksat_profile", "exponential")::String

    (; rootingdepth) = model.parameters.vegetation_parameters
    (; potential_evaporation, temperature) = atmospheric_forcing
    (; surface_water_flux, potential_soilevaporation, potential_transpiration) =
        model.boundary_conditions
    v = model.variables
    p = model.parameters

    n = length(potential_evaporation)
    threaded_foreach(1:n; basesize = 250) do i
        h3 = feddes_h3(p.h3_high[i], p.h3_low[i], potential_transpiration[i], dt)
        ustoredepth = sum(@view v.ustorelayerdepth[i][1:p.nlayers[i]])

        runoff_river = min(1.0, p.riverfrac[i]) * surface_water_flux[i]
        runoff_land = min(1.0, p.waterfrac[i]) * surface_water_flux[i]
        if !isnothing(demand) && (!isnothing(demand.paddy) || !isnothing(demand.nonpaddy))
            avail_forinfilt = avail_forinfilt + allocation.irri_alloc[i]
        end
        avail_forinfilt = max(surface_water_flux[i] - runoff_river - runoff_land, 0.0)

        ae_openw_r = min(
            v.waterlevel_river[i] * p.riverfrac[i],
            p.riverfrac[i] * potential_evaporation[i],
        )
        ae_openw_l = min(
            v.waterlevel_land[i] * p.waterfrac[i],
            p.waterfrac[i] * potential_evaporation[i],
        )

        potsoilevap = potential_soilevaporation[i]
        if !isnothing(demand) &&
           !isnothing(demand.paddy) &&
           demand.paddy.irrigation_areas[i]
            evap_paddy_water = min(paddy.h[i], potsoilevap)
            demand.paddy.h[i] -= evap_paddy_water
            potsoilevap -= evap_paddy_water
            avail_forinfilt += paddy.h[i] # allow infiltration of paddy water
        else
            evap_paddy_water = 0.0
        end

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
                    potsoilevap * min(1.0, saturationdeficit / p.soilwatercapacity[i])
            else
                # In case only the most upper soil layer contains unsaturated storage
                if n_usl == 1
                    # Check if groundwater level lies below the surface
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (v.zi[i] * (p.theta_s[i] - p.theta_r[i])))
                else
                    # In case first layer contains no saturated storage
                    soilevapunsat =
                        potsoilevap *
                        min(1.0, usld[1] / (usl[1] * ((p.theta_s[i] - p.theta_r[i]))))
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

        # actual transpiration, first from ustore, potential transpiration is partitioned over
        # depth based on the rootfraction
        actevapustore = 0.0
        rootfraction_unsat = 0.0
        for k in 1:n_usl
            vwc = max(usld[k] / usl[k], Float(0.0000001))
            head = head_brooks_corey(vwc, p.theta_s[i], p.theta_r[i], p.c[i][k], p.hb[i])
            alpha = rwu_reduction_feddes(head, p.h1[i], p.h2[i], h3, p.h4[i], p.alpha_h1[i])
            # availcap is fraction of soil layer containing roots
            # if `ust` is `true`, the whole unsaturated store is available for transpiration
            if ust
                availcap = usld[k] * 0.99
            else
                availcap =
                    min(1.0, max(0.0, (rootingdepth[i] - p.sumlayers[i][k]) / usl[k]))
            end
            maxextr = usld[k] * availcap
            # the rootfraction is valid for the root length in a soil layer, if zi decreases the root length
            # the rootfraction needs to be adapted
            if k == n_usl && v.zi[i] < rootingdepth[i]
                rootlength = min(p.act_thickl[i][k], rootingdepth[i] - p.sumlayers[i][k])
                rootfraction_act = p.rootfraction[i][k] * (usl[k] / rootlength)
            else
                rootfraction_act = p.rootfraction[i][k]
            end
            actevapustore_layer =
                min(alpha * rootfraction_act * potential_transpiration[i], maxextr)
            rootfraction_unsat = rootfraction_unsat + rootfraction_act
            ustorelayerdepth = usld[k] - actevapustore_layer
            actevapustore = actevapustore + actevapustore_layer
            usld = setindex(usld, ustorelayerdepth, k)
        end

        # transpiration from saturated store
        wetroots = scurve(v.zi[i], rootingdepth[i], Float(1.0), p.rootdistpar[i])
        alpha =
            rwu_reduction_feddes(Float(0.0), p.h1[i], p.h2[i], h3, p.h4[i], p.alpha_h1[i])
        # include remaining root fraction if rooting depth is below water table zi
        if v.zi[i] >= rootingdepth[i]
            f_roots = wetroots
            restevap = potential_transpiration[i] - actevapustore
        else
            f_roots = wetroots * (1.0 - rootfraction_unsat)
            restevap = potential_transpiration[i]
        end
        actevapsat = min(restevap * f_roots * alpha, satwaterdepth)
        satwaterdepth = satwaterdepth - actevapsat
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

            if v.zi[i] > rootingdepth[i]
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
        actevap = soilevap + transpiration + ae_openw_r + ae_openw_l + evap_paddy_water

        # update the outputs and states
        v.n_unsatlayers[i] = n_usl
        v.h3[i] = h3
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

function update!(model::SimpleBucketModel, demand, zi, exfiltsatwater)
    v = model.variables
    (;
        ustorelayerdepth,
        n_unsatlayers,
        excesswater,
        runoff_land,
        infiltexcess,
        vwc,
        vwc_perc,
    ) = v
    (; theta_s, theta_r, act_thickl, soilthickness, sumlayers, nlayers) = model.parameters
    (; rootingdepth) = model.parameters.vegetation_parameters

    n = length(zi)
    threaded_foreach(1:n; basesize = 1000) do i
        usl = set_layerthickness(zi[i], sumlayers[i], act_thickl[i])
        n_usl = number_of_active_layers(usl)
        # exfiltration from ustore
        usld = ustorelayerdepth[i]
        exfiltustore = 0.0
        for k in n_unsatlayers[i]:-1:1
            if k <= n_usl
                exfiltustore = max(0, usld[k] - usl[k] * (theta_s[i] - theta_r[i]))
            else
                exfiltustore = usld[k]
            end
            usld = setindex(usld, usld[k] - exfiltustore, k)
            if k > 1
                usld = setindex(usld, usld[k - 1] + exfiltustore, k - 1)
            end
        end

        ustoredepth = sum(@view usld[1:n_usl])

        if !isnothing(demand) &&
           !isnothing(demand.paddy) &&
           demand.paddy.irrigation_areas[i]
            paddy_h_add =
                exfiltustore +
                exfiltsatwater[i] +
                excesswater[i] +
                runoff_land[i] +
                infiltexcess[i]
            runoff = max(paddy_h_add - paddy.h_max[i], 0.0)
            paddy.h[i] = paddy_h_add - runoff
        else
            runoff =
                exfiltustore +
                exfiltsatwater[i] +
                excesswater[i] +
                runoff_land[i] +
                infiltexcess[i]
        end

        # volumetric water content per soil layer and root zone
        vwc_ = vwc[i]
        vwc_perc_ = vwc_perc[i]
        for k in 1:nlayers[i]
            if k <= n_usl
                vwc_ = setindex(
                    vwc_,
                    (usld[k] + (act_thickl[i][k] - usl[k]) * (theta_s[i] - theta_r[i])) / act_thickl[i][k] + theta_r[i],
                    k,
                )
            else
                vwc_ = setindex(vwc_, theta_s[i], k)
            end
            vwc_perc_ = setindex(vwc_perc_, (vwc_[k] / theta_s[i]) * 100.0, k)
        end

        rootstore_unsat = 0
        for k in 1:n_usl
            rootstore_unsat =
                rootstore_unsat +
                min(1.0, (max(0.0, rootingdepth[i] - sumlayers[i][k]) / usl[k])) * usld[k]
        end

        rootstore_sat = max(0.0, rootingdepth[i] - zi[i]) * (theta_s[i] - theta_r[i])
        rootstore = rootstore_sat + rootstore_unsat
        vwc_root = rootstore / rootingdepth[i] + theta_r[i]
        vwc_percroot = (vwc_root / theta_s[i]) * 100.0

        satwaterdepth = (soilthickness[i] - zi[i]) * (theta_s[i] - theta_r[i])

        # update the outputs and states
        v.n_unsatlayers[i] = n_usl
        v.ustorelayerdepth[i] = usld
        v.ustoredepth[i] = ustoredepth
        v.satwaterdepth[i] = satwaterdepth
        v.exfiltsatwater[i] = exfiltsatwater[i]
        v.exfiltustore[i] = exfiltustore
        v.runoff[i] = runoff
        v.net_runoff[i] = runoff - v.ae_openw_l[i]
        v.vwc[i] = vwc_
        v.vwc_perc[i] = vwc_perc_
        v.rootstore[i] = rootstore
        v.vwc_root[i] = vwc_root
        v.vwc_percroot[i] = vwc_percroot
        v.zi[i] = zi[i]
    end
end