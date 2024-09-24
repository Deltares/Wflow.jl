
# TODO (part of refactor v1.0): check if current variables and parameters of SBM soil model
# are fine or if for example some parameters should be part of another struct (e.g.
# riverfrac and waterfrac parameters). For ksat_profile parameters it is not required to
# store all parameters (depends on profile).

abstract type AbstractSoilModel{T} end

@get_units @grid_loc @with_kw struct SbmSoilVariables{T, N}
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
    # Total water storage (excluding floodplain volume, lakes and reservoirs) [mm]
    total_storage::Vector{T} | "mm"
    # Top soil temperature [ᵒC]
    tsoil::Vector{T} | "ᵒC"
    # Soil infiltration reduction factor (when soil is frozen) [-]
    soilinfredu::Vector{T} | "-"
end

function SbmSoilVariables(n, parameters)
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

    vars = SbmSoilVariables(;
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
        soilinfredu = fill(Float(1), n),
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
        total_storage = zeros(Float, n), # Set the total water storage from initialized values
    )
    return vars
end

@get_units @grid_loc @with_kw struct SbmSoilBC{T}
    surface_water_flux::Vector{T}
    potential_transpiration::Vector{T}
    potential_soilevaporation::Vector{T}
end

function SbmSoilBC(
    n;
    surface_water_flux::Vector{T} = fill(mv, n),
    potential_transpiration::Vector{T} = fill(mv, n),
    potential_soilevaporation::Vector{T} = fill(mv, n),
) where {T}
    return SbmSoilBC{T}(;
        surface_water_flux = surface_water_flux,
        potential_transpiration = potential_transpiration,
        potential_soilevaporation = potential_soilevaporation,
    )
end

@get_units @grid_loc @with_kw struct SbmSoilParameters{T, N, M}
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
    # Soil fraction [-]
    soil_fraction::Vector{T} | "-"
    # Vegetation parameter set
    vegetation_parameter_set::VegetationParameters{T} | "-" | "none"
end

@get_units @with_kw struct SbmSoilModel{T} <: AbstractSoilModel{T}
    boundary_conditions::SbmSoilBC{T} | "-"
    parameters::SbmSoilParameters{T} | "-"
    variables::SbmSoilVariables{T} | "-"
end

get_rootingdepth(model::SbmSoilModel) =
    model.parameters.vegetation_parameter_set.rootingdepth

function sbm_ksat_profiles(
    nc,
    config,
    inds,
    soilthickness,
    maxlayers,
    nlayers,
    sumlayers,
    dt,
)
    ksat_profile = get(config.input.vertical, "ksat_profile", "exponential")::String
    n = length(inds)
    if ksat_profile == "exponential"
        z_exp = soilthickness
        z_layered = fill(mv, n)
        kv = fill(mv, (maxlayers, n))
        nlayers_kv = fill(0, n)
    elseif ksat_profile == "exponential_constant"
        z_exp = ncread(
            nc,
            config,
            "vertical.soil.parameters.z_exp";
            optional = false,
            sel = inds,
            type = Float,
        )
        z_layered = fill(mv, n)
        kv = fill(mv, (maxlayers, n))
        nlayers_kv = fill(0, n)
    elseif ksat_profile == "layered" || ksat_profile == "layered_exponential"
        z_exp = fill(mv, n)
        kv =
            ncread(
                nc,
                config,
                "vertical.soil.parameters.kv";
                sel = inds,
                defaults = 1000.0,
                type = Float,
                dimname = :layer,
            ) .* (dt / basetimestep)
        if size(kv, 1) != maxlayers
            parname = param(config.input.vertical.soil.parameter, "kv")
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
                "vertical.soil.parameters.z_layered";
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

function SbmSoilParameters(nc, config, vegetation_parameter_set, inds, dt)
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
            "vertical.soil.parameters.w_soil";
            sel = inds,
            defaults = 0.1125,
            type = Float,
        ) .* (dt / basetimestep)
    cf_soil = ncread(
        nc,
        config,
        "vertical.soil.parameters.cf_soil";
        sel = inds,
        defaults = 0.038,
        type = Float,
    )
    # soil parameters
    theta_s = ncread(
        nc,
        config,
        "vertical.soil.parameters.theta_s";
        sel = inds,
        defaults = 0.6,
        type = Float,
    )
    theta_r = ncread(
        nc,
        config,
        "vertical.soil.parameters.theta_r";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )
    kv_0 =
        ncread(
            nc,
            config,
            "vertical.soil.parameters.kv_0";
            sel = inds,
            defaults = 3000.0,
            type = Float,
        ) .* (dt / basetimestep)
    f = ncread(
        nc,
        config,
        "vertical.soil.parameters.f";
        sel = inds,
        defaults = 0.001,
        type = Float,
    )
    hb = ncread(
        nc,
        config,
        "vertical.soil.parameters.hb";
        sel = inds,
        defaults = -10.0,
        type = Float,
    )
    h1 = ncread(
        nc,
        config,
        "vertical.soil.parameters.h1";
        sel = inds,
        defaults = 0.0,
        type = Float,
    )
    h2 = ncread(
        nc,
        config,
        "vertical.soil.parameters.h2";
        sel = inds,
        defaults = -100.0,
        type = Float,
    )
    h3_high = ncread(
        nc,
        config,
        "vertical.soil.parameters.h3_high";
        sel = inds,
        defaults = -400.0,
        type = Float,
    )
    h3_low = ncread(
        nc,
        config,
        "vertical.soil.parameters.h3_low";
        sel = inds,
        defaults = -1000.0,
        type = Float,
    )
    h4 = ncread(
        nc,
        config,
        "vertical.soil.parameters.h4";
        sel = inds,
        defaults = -15849.0,
        type = Float,
    )
    alpha_h1 = ncread(
        nc,
        config,
        "vertical.soil.parameters.alpha_h1";
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    soilthickness = ncread(
        nc,
        config,
        "vertical.soil.parameters.soilthickness";
        sel = inds,
        defaults = 2000.0,
        type = Float,
    )
    infiltcappath =
        ncread(
            nc,
            config,
            "vertical.soil.parameters.infiltcappath";
            sel = inds,
            defaults = 10.0,
            type = Float,
        ) .* (dt / basetimestep)
    infiltcapsoil =
        ncread(
            nc,
            config,
            "vertical.soil.parameters.infiltcapsoil";
            sel = inds,
            defaults = 100.0,
            type = Float,
        ) .* (dt / basetimestep)
    maxleakage =
        ncread(
            nc,
            config,
            "vertical.soil.parameters.maxleakage";
            sel = inds,
            defaults = 0.0,
            type = Float,
        ) .* (dt / basetimestep)

    c = ncread(
        nc,
        config,
        "vertical.soil.parameters.c";
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
        "vertical.soil.parameters.kvfrac";
        sel = inds,
        defaults = 1.0,
        type = Float,
        dimname = :layer,
    )
    if size(kvfrac, 1) != maxlayers
        parname = param(config.input, "vertical.soil.parameters.kvfrac")
        size1 = size(kvfrac, 1)
        error("$parname needs a layer dimension of size $maxlayers, but is $size1")
    end
    # fraction compacted area
    pathfrac = ncread(
        nc,
        config,
        "vertical.soil.parameters.pathfrac";
        sel = inds,
        defaults = 0.01,
        type = Float,
    )

    # vegetation parameters
    rootdistpar = ncread(
        nc,
        config,
        "vertical.soil.parameters.rootdistpar";
        sel = inds,
        defaults = -500.0,
        type = Float,
    )
    cap_hmax = ncread(
        nc,
        config,
        "vertical.soil.parameters.cap_hmax";
        sel = inds,
        defaults = 2000.0,
        type = Float,
    )
    cap_n = ncread(
        nc,
        config,
        "vertical.soil.parameters.cap_n";
        sel = inds,
        defaults = 2.0,
        type = Float,
    )

    if haskey(config.input.vertical.soil.parameters, "et_reftopot")
        @warn string(
            "The `et_reftopot` key in `[input.vertical.soil.parameters]` is now called ",
            "`kc`. Please update your TOML file.",
        )
    end

    act_thickl = set_layerthickness.(soilthickness, (cum_depth_layers,), (thicknesslayers,))
    sumlayers = @. pushfirst(cumsum(act_thickl), 0.0)
    nlayers = number_of_active_layers.(act_thickl)

    if length(config_thicknesslayers) > 0
        # root fraction read from nc file, in case of multiple soil layers and TOML file
        # includes "vertical.rootfraction"
        if haskey(config.input.vertical.soil.parameters, "rootfraction")
            rootfraction = ncread(
                nc,
                config,
                "vertical.soil.parameters.rootfraction";
                sel = inds,
                optional = false,
                type = Float,
                dimname = :layer,
            )
        else
            n = length(inds)
            (; rootingdepth) = vegetation_parameter_set
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

    z_exp, z_layered, kv, nlayers_kv = sbm_ksat_profiles(
        nc,
        config,
        inds,
        soilthickness,
        maxlayers,
        nlayers,
        sumlayers,
        dt,
    )

    soilwatercapacity = @. soilthickness * (theta_s - theta_r)

    n = length(inds)
    sbm_params = SbmSoilParameters(;
        maxlayers = maxlayers,
        nlayers = nlayers,
        nlayers_kv = nlayers_kv,
        soilwatercapacity = soilwatercapacity,
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
        soil_fraction = fill(mv, n),
        vegetation_parameter_set = vegetation_parameter_set,
    )
    return sbm_params
end

function SbmSoilModel(nc, config, vegetation_parameter_set, inds, dt)
    n = length(inds)
    params = SbmSoilParameters(nc, config, vegetation_parameter_set, inds, dt)
    vars = SbmSoilVariables(n, params)
    bc = SbmSoilBC(n)
    model = SbmSoilModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

function soil_fraction!(soil, runoff, glacier)
    (; canopygapfraction) = soil.parameters.vegetation_parameter_set
    (; soil_fraction) = soil.parameters
    (; waterfrac, riverfrac) = runoff.parameters
    glacier_fraction = get_glacier_fraction(glacier)

    @. soil_fraction =
        max(canopygapfraction - waterfrac - riverfrac - glacier_fraction, 0.0)
end

function update_boundary_conditions!(
    model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    external_models::NamedTuple,
)
    (; interception, runoff, allocation) = external_models
    (; potential_transpiration, surface_water_flux, potential_soilevaporation) =
        model.boundary_conditions

    potential_transpiration .= get_potential_transpiration(interception)

    irrigation = get_irrigation_allocated(allocation)
    @. surface_water_flux = max(
        runoff.boundary_conditions.surface_water_flux + irrigation -
        runoff.variables.runoff_river - runoff.variables.runoff_land,
        0.0,
    )

    @. potential_soilevaporation =
        model.parameters.soil_fraction * atmospheric_forcing.potential_evaporation
end

#TODO (v1.0):
# check if it's feasible to make the SbmSoilModel modular:
#   - soil evaporation
#   - transpiration (root water uptake Feddes)
#   - unsaturated flow (soil_model)
#   - capillary rise
#   - deep transfer

function update!(
    model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    external_models::NamedTuple,
    config,
    dt,
)
    soilinfreduction = get(config.model, "soilinfreduction", false)::Bool
    modelsnow = get(config.model, "snow", false)::Bool
    transfermethod = get(config.model, "transfermethod", false)::Bool
    ust = get(config.model, "whole_ust_available", false)::Bool # should be removed from optional setting and code?
    ksat_profile = get(config.input.vertical, "ksat_profile", "exponential")::String

    (; runoff, demand) = external_models
    (; potential_evaporation, temperature) = atmospheric_forcing
    (; ae_openw_r, ae_openw_l) = runoff.variables

    (; surface_water_flux, potential_soilevaporation, potential_transpiration) =
        model.boundary_conditions
    rootingdepth = get_rootingdepth(model)
    v = model.variables
    p = model.parameters

    evaporation!(demand.paddy, potential_soilevaporation)
    waterdepth_paddy = get_water_depth(demand.paddy)
    @. v.avail_forinfilt = surface_water_flux + waterdepth_paddy # allow infiltration of paddy water
    evap_paddy = get_evaporation(demand.paddy)
    @. potential_soilevaporation = potential_soilevaporation - evap_paddy

    n = length(potential_evaporation)
    threaded_foreach(1:n; basesize = 250) do i
        h3 = feddes_h3(p.h3_high[i], p.h3_low[i], potential_transpiration[i], dt)
        ustoredepth = sum(@view v.ustorelayerdepth[i][1:p.nlayers[i]])
        potsoilevap = potential_soilevaporation[i]

        # Calculate the initial capacity of the unsaturated store
        ustorecapacity = p.soilwatercapacity[i] - v.satwaterdepth[i] - ustoredepth

        if modelsnow
            tsoil = v.tsoil[i] + p.w_soil[i] * (temperature[i] - v.tsoil[i])
        else
            tsoil = v.tsoil[i]
        end

        # Calculate the infiltration flux into the soil column
        infiltsoilpath,
        infiltsoil,
        infiltpath,
        soilinf,
        pathinf,
        infiltexcess,
        soilinfredu = infiltration(
            v.avail_forinfilt[i],
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
        excesswater = v.avail_forinfilt[i] - infiltsoilpath - infiltexcess + du

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
        actevap =
            soilevap +
            transpiration +
            ae_openw_r[i] +
            ae_openw_l[i] +
            get_evaporation(demand.paddy, i)

        # update the outputs and states
        v.n_unsatlayers[i] = n_usl
        v.h3[i] = h3
        #v.avail_forinfilt[i] = avail_forinfilt
        v.actinfilt[i] = actinfilt
        v.infiltexcess[i] = infiltexcess
        v.recharge[i] = recharge
        v.transpiration[i] = transpiration
        v.soilevap[i] = soilevap
        v.soilevapsat[i] = soilevapsat
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
        v.soilinfredu[i] = soilinfredu
    end
end

function update!(
    model::SbmSoilModel,
    external_models::NamedTuple,
    boundary_conditions::NamedTuple,
)
    (; runoff, demand) = external_models
    (; ustorelayerdepth, n_unsatlayers, excesswater, infiltexcess, vwc, vwc_perc) =
        model.variables
    (; runoff_land, ae_openw_l) = runoff.variables
    (; theta_s, theta_r, act_thickl, soilthickness, sumlayers, nlayers) = model.parameters
    rootingdepth = get_rootingdepth(model)

    v = model.variables

    n = length(boundary_conditions.zi)
    threaded_foreach(1:n; basesize = 1000) do i
        usl = set_layerthickness(boundary_conditions.zi[i], sumlayers[i], act_thickl[i])
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

        runoff =
            exfiltustore +
            boundary_conditions.exfiltsatwater[i] +
            excesswater[i] +
            runoff_land[i] +
            infiltexcess[i]

        runoff = update_runoff(demand.paddy, runoff, i)

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

        rootstore_sat =
            max(0.0, rootingdepth[i] - boundary_conditions.zi[i]) *
            (theta_s[i] - theta_r[i])
        rootstore = rootstore_sat + rootstore_unsat
        vwc_root = rootstore / rootingdepth[i] + theta_r[i]
        vwc_percroot = (vwc_root / theta_s[i]) * 100.0

        satwaterdepth =
            (soilthickness[i] - boundary_conditions.zi[i]) * (theta_s[i] - theta_r[i])

        # update the outputs and states
        v.n_unsatlayers[i] = n_usl
        v.ustorelayerdepth[i] = usld
        v.ustoredepth[i] = ustoredepth
        v.satwaterdepth[i] = satwaterdepth
        v.exfiltsatwater[i] = boundary_conditions.exfiltsatwater[i]
        v.exfiltustore[i] = exfiltustore
        v.runoff[i] = runoff
        v.net_runoff[i] = runoff - ae_openw_l[i]
        v.vwc[i] = vwc_
        v.vwc_perc[i] = vwc_perc_
        v.rootstore[i] = rootstore
        v.vwc_root[i] = vwc_root
        v.vwc_percroot[i] = vwc_percroot
        v.zi[i] = boundary_conditions.zi[i]
    end
end