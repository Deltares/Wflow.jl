@get_units @exchange @grid_type @grid_location @with_kw struct SBM{
    T,
    N,
    M,
    P<:Union{Paddy,Nothing},
    NP<:Union{NonPaddy,Nothing},
    D<:Union{NonIrrigationDemand,Nothing},
    L<:Union{NonIrrigationDemand,Nothing},
    I<:Union{NonIrrigationDemand,Nothing},
    A<:Union{AllocationLand,Nothing},
}
    # Model time step [s]
    dt::T | "s" | 0 | "none" | "none"
    # Maximum number of soil layers
    maxlayers::Int | "-" | 0 | "none" | "none"
    # number of cells
    n::Int | "-" | 0 | "none" | "none"
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
    kv::Vector{SVector{N,T}} | "-"
    # Muliplication factor [-] applied to kv_z (vertical flow)
    kvfrac::Vector{SVector{N,T}} | "-"
    # Air entry pressure [cm] of soil (Brooks-Corey)
    hb::Vector{T} | "cm"
    # Soil thickness [mm]
    soilthickness::Vector{T} | "mm"
    # Thickness of soil layers [mm]
    act_thickl::Vector{SVector{N,T}} | "mm"
    # Cumulative sum of soil layers [mm], starting at soil surface (0)
    sumlayers::Vector{SVector{M,T}} | "mm"
    # Infiltration capacity of the compacted areas [mm Δt⁻¹]
    infiltcappath::Vector{T}
    # Soil infiltration capacity [mm Δt⁻¹]
    infiltcapsoil::Vector{T}
    # Soil infiltration reduction factor (when soil is frozen) [-]
    soilinfredu::Vector{T} | "-"
    # Maximum leakage [mm Δt⁻¹] from saturated zone
    maxleakage::Vector{T}
    # Fraction of open water (excluding rivers) [-]
    waterfrac::Vector{T} | "-"
    # Fraction of compacted area  [-]
    pathfrac::Vector{T} | "-"
    # Rooting depth [mm]
    rootingdepth::Vector{T} | "mm"
    # Fraction of the root length density in each soil layer [-]
    rootfraction::Vector{SVector{N,T}} | "-"
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
    # Calculated soil water pressure head h3 of the root water uptake reduction function (Feddes) [cm]
    h3::Vector{T} | "cm"
    # Root water uptake reduction at soil water pressure head h1 (0.0 or 1.0) [-]
    alpha_h1::Vector{T} | "-"
    # Controls how roots are linked to water table [-]
    rootdistpar::Vector{T} | "-"
    # Parameter [mm] controlling capillary rise
    cap_hmax::Vector{T} | "mm"
    # Coefficient [-] controlling capillary rise
    cap_n::Vector{T} | "-"
    # Crop coefficient Kc [-]
    kc::Vector{T} | "-"
    # Brooks-Corey power coefﬁcient [-] for each soil layer
    c::Vector{SVector{N,T}} | "-"
    # Stemflow [mm Δt⁻¹]
    stemflow::Vector{T}
    # Throughfall [mm Δt⁻¹]
    throughfall::Vector{T}
    # A scaling parameter [mm⁻¹] (controls exponential decline of kv_0)
    f::Vector{T} | "mm-1"
    # Depth [mm] from soil surface for which exponential decline of kv_0 is valid
    z_exp::Vector{T} | "mm"
    # Depth [mm] from soil surface for which layered profile is valid
    z_layered::Vector{T} | "mm"
    # Amount of water in the unsaturated store, per layer [mm]
    ustorelayerdepth::Vector{SVector{N,T}} | "mm"
    # Saturated store [mm]
    satwaterdepth::Vector{T} | "mm"
    # Pseudo-water table depth [mm] (top of the saturated zone)
    zi::Vector{T} | "mm"
    # Soilwater capacity [mm]
    soilwatercapacity::Vector{T} | "mm"
    # Canopy storage [mm]
    canopystorage::Vector{T} | "mm"
    # Maximum canopy storage [mm]
    cmax::Vector{T} | "mm"
    # Canopy gap fraction [-]
    canopygapfraction::Vector{T} | "-"
    # Gash interception model parameter, ratio of the average evaporation from the
    # wet canopy [mm Δt⁻¹] and the average precipitation intensity [mm Δt⁻¹] on a saturated canopy
    e_r::Vector{T} | "-"
    # Precipitation [mm Δt⁻¹]
    precipitation::Vector{T}
    # Temperature [ᵒC]
    temperature::Vector{T} | "°C"
    # Potential reference evapotranspiration [mm Δt⁻¹]
    potential_evaporation::Vector{T}
    # Potential transpiration (after subtracting interception from potential_evaporation)
    pottrans::Vector{T}
    # Transpiration [mm Δt⁻¹]
    transpiration::Vector{T}
    # Actual evaporation from unsaturated store [mm Δt⁻¹]
    ae_ustore::Vector{T}
    # Interception loss by evaporation [mm Δt⁻¹]
    interception::Vector{T}
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
    # Infiltration from surface water [mm Δt⁻¹]
    infilt_surfacewater::Vector{T}
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
    vwc::Vector{SVector{N,T}} | "-"
    # Volumetric water content [%] per soil layer (including theta_r and saturated zone)
    vwc_perc::Vector{SVector{N,T}} | "%"
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
    ### Snow parameters ###
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹]
    cfmax::Vector{T} | "mm ᵒC-1 dt-1"
    # Threshold temperature for snowfall [ᵒC]
    tt::Vector{T} | "ᵒC"
    # Threshold temperature interval length [ᵒC]
    tti::Vector{T} | "ᵒC"
    # Threshold temperature for snowmelt [ᵒC]
    ttm::Vector{T} | "ᵒC"
    # Water holding capacity as fraction of current snow pack [-]
    whc::Vector{T} | "-"
    # Soil temperature smooth factor [-]
    w_soil::Vector{T} | "-"
    # Controls soil infiltration reduction factor when soil is frozen [-]
    cf_soil::Vector{T} | "-"
    # Snow storage [mm]
    snow::Vector{T} | "mm"
    # Liquid water content in the snow pack [mm]
    snowwater::Vector{T} | "mm"
    # Snow melt + precipitation as rainfall [mm]
    rainfallplusmelt::Vector{T} | "mm"
    # Threshold temperature for snowfall above glacier [ᵒC]
    g_tt::Vector{T} | "ᵒC"
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹] for glacier
    g_cfmax::Vector{T} | "mm ᵒC-1 dt-1"
    # Fraction of the snowpack on top of the glacier converted into ice [Δt⁻¹]
    g_sifrac::Vector{T} | "dt-1"
    # Water within the glacier [mm]
    glacierstore::Vector{T} | "mm"
    # Fraction covered by a glacier [-]
    glacierfrac::Vector{T} | "-"
    # Top soil temperature [ᵒC]
    tsoil::Vector{T} | "ᵒC"
    ## Interception related to leaf_area_index climatology ###
    # Specific leaf storage [mm]
    sl::Vector{T} | "mm"
    # Storage woody part of vegetation [mm]
    swood::Vector{T} | "mm"
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    kext::Vector{T} | "-"
    # Leaf area index [m² m⁻²]
    leaf_area_index::Vector{T} | "m2 m-2"
    # Water level land [mm]
    waterlevel_land::Vector{T} | "mm"
    # Water level river [mm]
    waterlevel_river::Vector{T} | "mm"
    # Total water storage (excluding floodplain volume, lakes and reservoirs) [mm]
    total_storage::Vector{T} | "mm"
    # Water demand structs (of arrays)
    paddy::P | "-" | 0
    nonpaddy::NP | "-" | 0
    domestic::D | "-" | 0
    livestock::L | "-" | 0
    industry::I | "-" | 0
    allocation::A | "-" | 0


    function SBM{T,N,M,P,NP,D,L,I,A}(args...) where {T,N,M,P,NP,D,L,I,A}
        equal_size_vectors(args)
        return new(args...)
    end
end


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

    cfmax =
        ncread(
            nc,
            config,
            "vertical.cfmax";
            sel = inds,
            defaults = 3.75653,
            type = Float,
        ) .* (dt / basetimestep)
    tt = ncread(nc, config, "vertical.tt"; sel = inds, defaults = 0.0, type = Float)
    tti = ncread(nc, config, "vertical.tti"; sel = inds, defaults = 1.0, type = Float)
    ttm = ncread(nc, config, "vertical.ttm"; sel = inds, defaults = 0.0, type = Float)
    whc = ncread(nc, config, "vertical.whc"; sel = inds, defaults = 0.1, type = Float)
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
    # glacier parameters
    g_tt = ncread(
        nc,
        config,
        "vertical.g_tt";
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    g_cfmax =
        ncread(
            nc,
            config,
            "vertical.g_cfmax";
            sel = inds,
            defaults = 3.0,
            type = Float,
            fill = 0.0,
        ) .* (dt / basetimestep)
    g_sifrac =
        ncread(
            nc,
            config,
            "vertical.g_sifrac";
            sel = inds,
            defaults = 0.001,
            type = Float,
            fill = 0.0,
        ) .* (dt / basetimestep)
    glacierfrac = ncread(
        nc,
        config,
        "vertical.glacierfrac";
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    glacierstore = ncread(
        nc,
        config,
        "vertical.glacierstore";
        sel = inds,
        defaults = 5500.0,
        type = Float,
        fill = 0.0,
    )
    # soil parameters
    theta_s =
        ncread(nc, config, "vertical.theta_s"; sel = inds, defaults = 0.6, type = Float)
    theta_r =
        ncread(nc, config, "vertical.theta_r"; sel = inds, defaults = 0.01, type = Float)
    kv_0 =
        ncread(nc, config, "vertical.kv_0"; sel = inds, defaults = 3000.0, type = Float) .*
        (dt / basetimestep)
    f = ncread(nc, config, "vertical.f"; sel = inds, defaults = 0.001, type = Float)
    hb = ncread(nc, config, "vertical.hb"; sel = inds, defaults = -10.0, type = Float)
    h1 = ncread(nc, config, "vertical.h1"; sel = inds, defaults = 0.0, type = Float)
    h2 = ncread(nc, config, "vertical.h2"; sel = inds, defaults = -100.0, type = Float)
    h3_high =
        ncread(nc, config, "vertical.h3_high"; sel = inds, defaults = -400.0, type = Float)
    h3_low =
        ncread(nc, config, "vertical.h3_low"; sel = inds, defaults = -1000.0, type = Float)
    h4 = ncread(nc, config, "vertical.h4"; sel = inds, defaults = -15849.0, type = Float)
    alpha_h1 =
        ncread(nc, config, "vertical.alpha_h1"; sel = inds, defaults = 1.0, type = Float)
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
    rootingdepth = ncread(
        nc,
        config,
        "vertical.rootingdepth";
        sel = inds,
        defaults = 750.0,
        type = Float,
    )
    # correct rooting depth for soilthickness
    rootingdepth = @. min(0.99 * soilthickness, rootingdepth)
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
    kc = ncread(
        nc,
        config,
        "vertical.kc";
        alias = "vertical.et_reftopot",
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    if haskey(config.input.vertical, "et_reftopot")
        @warn string(
            "The `et_reftopot` key in `[input.vertical]` is now called ",
            "`kc`. Please update your TOML file.",
        )
    end

    cmax, e_r, canopygapfraction, sl, swood, kext = initialize_canopy(nc, config, inds)

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

    for i = 1:n
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

    if length(config_thicknesslayers) > 0
        # root fraction read from nc file, in case of multiple soil layers and TOML file
        # includes "vertical.rootfraction"
        if haskey(config.input.vertical, "rootfraction")
            rootfraction = ncread(
                nc,
                config,
                "vertical.rootfraction";
                sel = inds,
                optional = false,
                type = Float,
                dimname = :layer,
            )
        else
            # default root fraction in case of multiple soil layers
            rootfraction = zeros(Float, maxlayers, n)
            for i = 1:n
                if rootingdepth[i] > 0.0
                    for k = 1:maxlayers
                        if (rootingdepth[i] - s_layers[k, i]) >= act_thickl[k, i]
                            rootfraction[k, i] = act_thickl[k, i] / rootingdepth[i]
                        else
                            rootfraction[k, i] =
                                max(rootingdepth[i] - s_layers[k, i], 0.0) / rootingdepth[i]
                        end
                    end
                end
            end
        end
    else
        # for the case of 1 soil layer
        rootfraction = ones(Float, maxlayers, n)
    end

    # needed for derived parameters below
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

    # water demand and irrigation options
    do_water_demand = haskey(config.model, "water_demand")
    domestic = do_water_demand ? get(config.model.water_demand, "domestic", false) : false
    industry = do_water_demand ? get(config.model.water_demand, "industry", false) : false
    livestock = do_water_demand ? get(config.model.water_demand, "livestock", false) : false
    paddy = do_water_demand ? get(config.model.water_demand, "paddy", false) : false
    nonpaddy = do_water_demand ? get(config.model.water_demand, "nonpaddy", false) : false

    sbm = SBM(
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
        h1 = h1,
        h2 = h2,
        h3_high = h3_high,
        h3_low = h3_low,
        h4 = h4,
        h3 = fill(mv, n),
        alpha_h1 = alpha_h1,
        soilthickness = soilthickness,
        act_thickl = act_thickl,
        sumlayers = sumlayers,
        infiltcappath = infiltcappath,
        infiltcapsoil = infiltcapsoil,
        soilinfredu = fill(Float(1), n),
        maxleakage = maxleakage,
        waterfrac = max.(waterfrac .- riverfrac, Float(0.0)),
        pathfrac = pathfrac,
        rootingdepth = rootingdepth,
        rootfraction = svectorscopy(rootfraction, Val{maxlayers}()),
        rootdistpar = rootdistpar,
        cap_hmax = cap_hmax,
        cap_n = cap_n,
        kc = kc,
        c = svectorscopy(c, Val{maxlayers}()),
        stemflow = fill(mv, n),
        throughfall = fill(mv, n),
        f = f,
        z_exp = z_exp,
        z_layered = z_layered,
        ustorelayerdepth = zero(act_thickl),
        satwaterdepth = satwaterdepth,
        zi = zi,
        soilwatercapacity = soilwatercapacity,
        canopystorage = zeros(Float, n),
        cmax = cmax,
        canopygapfraction = canopygapfraction,
        e_r = e_r,
        precipitation = fill(mv, n),
        temperature = fill(mv, n),
        potential_evaporation = fill(mv, n),
        pottrans = fill(mv, n),
        transpiration = fill(mv, n),
        ae_ustore = fill(mv, n),
        interception = fill(mv, n),
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
        infilt_surfacewater = fill(mv, n),
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
        # snow parameters
        cfmax = cfmax,
        tt = tt,
        tti = tti,
        ttm = ttm,
        whc = whc,
        w_soil = w_soil,
        cf_soil = cf_soil,
        snow = zeros(Float, n),
        snowwater = zeros(Float, n),
        rainfallplusmelt = fill(mv, n),
        tsoil = fill(Float(10.0), n),
        # glacier parameters
        g_tt = g_tt,
        g_sifrac = g_sifrac,
        g_cfmax = g_cfmax,
        glacierstore = glacierstore,
        glacierfrac = glacierfrac,
        # Interception related to climatology (leaf_area_index)
        sl = sl,
        swood = swood,
        kext = kext,
        leaf_area_index = fill(mv, n),
        # water level land and river domain
        waterlevel_land = fill(mv, n),
        waterlevel_river = zeros(Float, n), #set to zero to account for cells outside river domain
        total_storage = zeros(Float, n), # Set the total water storage from initialized values
        # water demand
        paddy = paddy ? initialize_paddy(nc, config, inds, dt) : nothing,
        nonpaddy = nonpaddy ? initialize_nonpaddy(nc, config, inds, dt) : nothing,
        domestic = domestic ? initialize_domestic_demand(nc, config, inds, dt) : nothing,
        industry = industry ? initialize_industry_demand(nc, config, inds, dt) : nothing,
        livestock = livestock ? initialize_livestock_demand(nc, config, inds, dt) : nothing,
        allocation = do_water_demand ? initialize_allocation_land(nc, config, inds) :
                     nothing,
    )

    return sbm

end


function update_until_snow(sbm::SBM, config)

    do_lai = haskey(config.input.vertical, "leaf_area_index")
    modelglacier = get(config.model, "glacier", false)::Bool
    modelsnow = get(config.model, "snow", false)::Bool

    threaded_foreach(1:sbm.n, basesize = 1000) do i
        if do_lai
            cmax = sbm.sl[i] * sbm.leaf_area_index[i] + sbm.swood[i]
            canopygapfraction = exp(-sbm.kext[i] * sbm.leaf_area_index[i])
            canopyfraction = 1.0 - canopygapfraction
            ewet = canopyfraction * sbm.potential_evaporation[i] * sbm.kc[i]
            e_r =
                sbm.precipitation[i] > 0.0 ?
                min(0.25, ewet / max(0.0001, canopyfraction * sbm.precipitation[i])) : 0.0
        else
            cmax = sbm.cmax[i]
            canopygapfraction = sbm.canopygapfraction[i]
            e_r = sbm.e_r[i]
        end

        canopy_potevap =
            sbm.kc[i] * sbm.potential_evaporation[i] * (1.0 - canopygapfraction)
        if Second(sbm.dt) >= Hour(23)
            throughfall, interception, stemflow, canopystorage = rainfall_interception_gash(
                cmax,
                e_r,
                canopygapfraction,
                sbm.precipitation[i],
                sbm.canopystorage[i],
                canopy_potevap,
            )
            pottrans = max(0.0, canopy_potevap - interception) # now in mm
        else
            netinterception, throughfall, stemflow, leftover, interception, canopystorage =
                rainfall_interception_modrut(
                    sbm.precipitation[i],
                    canopy_potevap,
                    sbm.canopystorage[i],
                    canopygapfraction,
                    cmax,
                )
            pottrans = max(0.0, leftover)  # now in mm
        end

        if modelsnow
            tsoil = sbm.tsoil[i] + sbm.w_soil[i] * (sbm.temperature[i] - sbm.tsoil[i])
            snow, snowwater, snowmelt, rainfallplusmelt, snowfall = snowpack_hbv(
                sbm.snow[i],
                sbm.snowwater[i],
                throughfall + stemflow,
                sbm.temperature[i],
                sbm.tti[i],
                sbm.tt[i],
                sbm.ttm[i],
                sbm.cfmax[i],
                sbm.whc[i],
            )
        end

        h3 = feddes_h3(sbm.h3_high[i], sbm.h3_low[i], pottrans, Second(sbm.dt))

        # update the outputs and states
        sbm.e_r[i] = e_r
        sbm.cmax[i] = cmax
        sbm.canopygapfraction[i] = canopygapfraction
        sbm.canopystorage[i] = canopystorage
        sbm.interception[i] = interception
        sbm.stemflow[i] = stemflow
        sbm.throughfall[i] = throughfall
        sbm.pottrans[i] = pottrans
        sbm.h3[i] = h3
        if modelsnow
            sbm.snow[i] = snow
            sbm.snowwater[i] = snowwater
            sbm.tsoil[i] = tsoil
            sbm.rainfallplusmelt[i] = rainfallplusmelt
        end
    end
end

function update_until_recharge(sbm::SBM, config)

    # start dummy variables (should be generated from model reader and from Config.jl TOML)
    soilinfreduction = get(config.model, "soilinfreduction", false)::Bool
    modelglacier = get(config.model, "glacier", false)::Bool
    modelsnow = get(config.model, "snow", false)::Bool
    transfermethod = get(config.model, "transfermethod", false)::Bool
    ust = get(config.model, "whole_ust_available", false)::Bool # should be removed from optional setting and code?
    ksat_profile = get(config.input.vertical, "ksat_profile", "exponential")::String
    do_surface_water_infiltration =
        get(config.model, "surface_water_infiltration", false)::Bool

    threaded_foreach(1:sbm.n, basesize = 250) do i
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
        if !isnothing(sbm.paddy) || !isnothing(sbm.nonpaddy)
            avail_forinfilt = avail_forinfilt + sbm.allocation.irri_alloc[i]
        end
        avail_forinfilt = max(avail_forinfilt - runoff_river - runoff_land, 0.0)

        ae_openw_r = min(
            sbm.waterlevel_river[i] * sbm.riverfrac[i],
            sbm.riverfrac[i] * sbm.potential_evaporation[i],
        )
        ae_openw_l = min(
            sbm.waterlevel_land[i] * sbm.waterfrac[i],
            sbm.waterfrac[i] * sbm.potential_evaporation[i],
        )

        # Update land waterlevel
        waterlevel_land = sbm.waterlevel_land[i] - ae_openw_l

        if do_surface_water_infiltration
            # Add land waterlevel to infiltration
            avail_forinfilt += waterlevel_land
        end

        # evap available for soil evaporation
        soilevap_fraction = max(
            sbm.canopygapfraction[i] - sbm.riverfrac[i] - sbm.waterfrac[i] -
            sbm.glacierfrac[i],
            0.0,
        )
        potsoilevap = soilevap_fraction * sbm.potential_evaporation[i]

        if !isnothing(sbm.paddy) && sbm.paddy.irrigation_areas[i]
            evap_paddy_water = min(sbm.paddy.h[i], potsoilevap)
            sbm.paddy.h[i] -= evap_paddy_water
            potsoilevap -= evap_paddy_water
            avail_forinfilt += sbm.paddy.h[i] # allow infiltration of paddy water
        else
            evap_paddy_water = 0.0
        end

        # Calculate the initial capacity of the unsaturated store
        ustorecapacity = sbm.soilwatercapacity[i] - sbm.satwaterdepth[i] - ustoredepth

        # Calculate the infiltration flux into the soil column
        infiltsoilpath,
        infiltsoil,
        infiltpath,
        soilinf,
        pathinf,
        infiltexcess,
        soilinfredu = infiltration(
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
                for m = 1:n_usl
                    l_sat = usl[m] * (sbm.theta_s[i] - sbm.theta_r[i])
                    kv_z = hydraulic_conductivity_at_depth(sbm, z[m], i, m, ksat_profile)
                    ustorelayerdepth =
                        m == 1 ? sbm.ustorelayerdepth[i][m] + infiltsoilpath :
                        sbm.ustorelayerdepth[i][m] + ast
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

        # actual transpiration, first from ustore, potential transpiration is partitioned over
        # depth based on the rootfraction
        actevapustore = 0.0
        rootfraction_unsat = 0.0
        for k = 1:n_usl
            vwc = max(usld[k] / usl[k], Float(0.0000001))
            head = head_brooks_corey(
                vwc,
                sbm.theta_s[i],
                sbm.theta_r[i],
                sbm.c[i][k],
                sbm.hb[i],
            )
            alpha = rwu_reduction_feddes(
                head,
                sbm.h1[i],
                sbm.h2[i],
                sbm.h3[i],
                sbm.h4[i],
                sbm.alpha_h1[i],
            )
            # availcap is fraction of soil layer containing roots
            # if `ust` is `true`, the whole unsaturated store is available for transpiration
            if ust
                availcap = usld[k] * 0.99
            else
                availcap =
                    min(1.0, max(0.0, (sbm.rootingdepth[i] - sbm.sumlayers[i][k]) / usl[k]))
            end
            maxextr = usld[k] * availcap
            # the rootfraction is valid for the root length in a soil layer, if zi decreases the root length
            # the rootfraction needs to be adapted
            if k == n_usl && sbm.zi[i] < sbm.rootingdepth[i]
                rootlength =
                    min(sbm.act_thickl[i][k], sbm.rootingdepth[i] - sbm.sumlayers[i][k])
                rootfraction_act = sbm.rootfraction[i][k] * (usl[k] / rootlength)
            else
                rootfraction_act = sbm.rootfraction[i][k]
            end
            actevapustore_layer = min(alpha * rootfraction_act * sbm.pottrans[i], maxextr)
            rootfraction_unsat = rootfraction_unsat + rootfraction_act
            ustorelayerdepth = usld[k] - actevapustore_layer
            actevapustore = actevapustore + actevapustore_layer
            usld = setindex(usld, ustorelayerdepth, k)
        end

        # transpiration from saturated store
        wetroots = scurve(sbm.zi[i], sbm.rootingdepth[i], Float(1.0), sbm.rootdistpar[i])
        alpha = rwu_reduction_feddes(
            Float(0.0),
            sbm.h1[i],
            sbm.h2[i],
            sbm.h3[i],
            sbm.h4[i],
            sbm.alpha_h1[i],
        )
        # include remaining root fraction if rooting depth is below water table zi
        if sbm.zi[i] >= sbm.rootingdepth[i]
            f_roots = wetroots
            restevap = sbm.pottrans[i] - actevapustore
        else
            f_roots = wetroots * (1.0 - rootfraction_unsat)
            restevap = sbm.pottrans[i]
        end
        actevapsat = min(restevap * f_roots * alpha, satwaterdepth)
        satwaterdepth = satwaterdepth - actevapsat

        # check soil moisture balance per layer
        du = 0.0
        for k = n_usl:-1:1
            du = max(0.0, usld[k] - usl[k] * (sbm.theta_s[i] - sbm.theta_r[i]))
            usld = setindex(usld, usld[k] - du, k)
            if k > 1
                usld = setindex(usld, usld[k-1] + du, k - 1)
            end
        end

        actinfilt = infiltsoilpath - du

        # Scale infiltration from surface water based on the ratio between actinfil and
        # infiltsoilpath (and prevent division by zero)
        if do_surface_water_infiltration
            infilt_ratio = iszero(avail_forinfilt) ? 0.0 : actinfilt / avail_forinfilt
            infilt_surfacewater = max(0.0, waterlevel_land * infilt_ratio)
        else
            infilt_surfacewater = 0.0
        end

        excesswater = avail_forinfilt - actinfilt - infiltexcess

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

            if sbm.zi[i] > sbm.rootingdepth[i]
                capflux =
                    maxcapflux * pow(
                        1.0 - min(sbm.zi[i], sbm.cap_hmax[i]) / (sbm.cap_hmax[i]),
                        sbm.cap_n[i],
                    )
            else
                capflux = 0.0
            end

            netcapflux = capflux
            for k = n_usl:-1:1
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
        actevap =
            soilevap +
            transpiration +
            ae_openw_r +
            ae_openw_l +
            sbm.interception[i] +
            evap_paddy_water

        # update the outputs and states
        sbm.n_unsatlayers[i] = n_usl
        sbm.net_runoff_river[i] = runoff_river - ae_openw_r
        sbm.avail_forinfilt[i] = avail_forinfilt
        sbm.actinfilt[i] = actinfilt
        sbm.infiltexcess[i] = infiltexcess
        sbm.infilt_surfacewater[i] = infilt_surfacewater
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
        sbm.soilinfredu[i] = soilinfredu
        if modelsnow
            if modelglacier
                sbm.snow[i] = snow
                sbm.glacierstore[i] = glacierstore
            end
        end
    end
end

function update_after_subsurfaceflow(sbm::SBM, zi, exfiltsatwater, config)

    threaded_foreach(1:sbm.n, basesize = 1000) do i
        usl, n_usl = set_layerthickness(zi[i], sbm.sumlayers[i], sbm.act_thickl[i])
        # exfiltration from ustore
        usld = sbm.ustorelayerdepth[i]
        exfiltustore = 0.0
        for k = sbm.n_unsatlayers[i]:-1:1
            if k <= n_usl
                exfiltustore = max(0, usld[k] - usl[k] * (sbm.theta_s[i] - sbm.theta_r[i]))
            else
                exfiltustore = usld[k]
            end
            usld = setindex(usld, usld[k] - exfiltustore, k)
            if k > 1
                usld = setindex(usld, usld[k-1] + exfiltustore, k - 1)
            end
        end

        ustoredepth = sum(@view usld[1:n_usl])

        # Correct the excess water to only include the water that is not (yet) surface water,
        # otherwise this water would be accounted for twice. This is only relevant if
        # do_surface_water_infiltration is true
        do_surface_water_infiltration =
            get(config.model, "surface_water_infiltration", false)::Bool
        if do_surface_water_infiltration
            contribution_surfacewater = sbm.waterlevel_land[i] - sbm.ae_openw_l[i]
            correction_surfacewater =
                iszero(sbm.avail_forinfilt[i]) ? 1.0 :
                1.0 - (contribution_surfacewater / sbm.avail_forinfilt[i])
        else
            correction_surfacewater = 1.0
        end

        if !isnothing(sbm.paddy) && sbm.paddy.irrigation_areas[i]
            paddy_h_add =
                exfiltustore +
                exfiltsatwater[i] +
                sbm.excesswater[i] * correction_surfacewater +
                sbm.runoff_land[i] +
                sbm.infiltexcess[i] * correction_surfacewater
            runoff = max(paddy_h_add - sbm.paddy.h_max[i], 0.0)
            sbm.paddy.h[i] = paddy_h_add - runoff
        else
            runoff =
                exfiltustore +
                exfiltsatwater[i] +
                sbm.excesswater[i] * correction_surfacewater +
                sbm.runoff_land[i] +
                sbm.infiltexcess[i] * correction_surfacewater
        end


        # volumetric water content per soil layer and root zone
        vwc = sbm.vwc[i]
        vwc_perc = sbm.vwc_perc[i]
        for k = 1:sbm.nlayers[i]
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
        for k = 1:n_usl
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
        sbm.net_runoff[i] = runoff - sbm.ae_openw_l[i] - sbm.infilt_surfacewater[i]
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
- area:
    Area of the cells acquired from model.network.land.area
- river_routing:
    The river routing struct, i.e. model.lateral.river
- land_routing:
    The land routing struct, i.e. model.lateral.land
"""
function update_total_water_storage(
    sbm::SBM,
    river_network,
    area,
    river_routing,
    land_routing,
)
    # Set the total storage to zero
    fill!(sbm.total_storage, 0)

    # Burn the river routing values
    for (i, index_river) in enumerate(river_network)
        sbm.total_storage[index_river] = (
            (river_routing.h_av[i] * river_routing.width[i] * river_routing.dl[i]) /
            (area[index_river]) * 1000 # Convert to mm
        )
    end

    # Chunk the data for parallel computing
    threaded_foreach(1:sbm.n, basesize = 1000) do i

        # Paddy water depth
        paddy_h = isnothing(sbm.paddy) ? 0.0 : sbm.paddy.h[i]

        # Cumulate per vertical type
        # Maybe re-categorize in the future
        surface = (
            sbm.glacierstore[i] * sbm.glacierfrac[i] +
            sbm.snow[i] +
            sbm.snowwater[i] +
            sbm.canopystorage[i] +
            paddy_h
        )
        sub_surface = sbm.ustoredepth[i] + sbm.satwaterdepth[i]
        lateral = (
            land_routing.h_av[i] * (1 - sbm.riverfrac[i]) * 1000 # convert to mm
        )

        # Add everything to the total water storage
        sbm.total_storage[i] += (surface + sub_surface + lateral)
    end
end

"""
    update_water_demand(sbm::SBM)

Update water demand for vertical `SBM` concept for a single timestep. Water demand is
computed for sectors `industry`, `domestic` and `livestock`, and `paddy` rice fields and
`nonpaddy` (other crop) fields.

Gross water demand for irrigation `irri_demand_gross` and non-irrigation
`nonirri_demand_gross`, and total gross water demand `total_gross_demand` are updated as
part of `SBM` water allocation (`allocation`)
"""
function update_water_demand(sbm::SBM)
    for i = 1:sbm.n

        industry_dem = update_non_irrigation_demand(sbm.industry, i)
        domestic_dem = update_non_irrigation_demand(sbm.domestic, i)
        livestock_dem = update_non_irrigation_demand(sbm.livestock, i)

        irri_dem_gross = 0.0
        if !isnothing(sbm.nonpaddy) && sbm.nonpaddy.irrigation_areas[i]
            if sbm.nonpaddy.irrigation_trigger[i]
                usl, _ = set_layerthickness(sbm.zi[i], sbm.sumlayers[i], sbm.act_thickl[i])
                for k = 1:sbm.n_unsatlayers[i]
                    # compute water demand only for root zone through root fraction per layer
                    rootfrac = min(
                        1.0,
                        (max(0.0, sbm.rootingdepth[i] - sbm.sumlayers[i][k]) / usl[k]),
                    )
                    # vwc_f and vwc_h3 can be precalculated.
                    vwc_fc = vwc_brooks_corey(
                        -100.0,
                        sbm.hb[i],
                        sbm.theta_s[i],
                        sbm.theta_r[i],
                        sbm.c[i][k],
                    )
                    vwc_h3 = vwc_brooks_corey(
                        sbm.h3[i],
                        sbm.hb[i],
                        sbm.theta_s[i],
                        sbm.theta_r[i],
                        sbm.c[i][k],
                    )
                    depletion =
                        (vwc_fc * usl[k]) -
                        (sbm.ustorelayerdepth[i][k] + sbm.theta_r[i] * usl[k])
                    depletion *= rootfrac
                    raw = (vwc_fc - vwc_h3) * usl[k] # readily available water
                    raw *= rootfrac

                    # check if maximum irrigation rate has been applied at the previous time step.
                    max_irri_rate_applied =
                        sbm.nonpaddy.demand_gross[i] ==
                        sbm.nonpaddy.maximum_irrigation_rate[i]
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
                    sbm.soilinfredu[i] * (1.0 - sbm.pathfrac[i]) * sbm.infiltcapsoil[i]
                irri_dem_gross = min(irri_dem_gross, infiltration_capacity)
                irri_dem_gross /= sbm.nonpaddy.irrigation_efficiency[i]
                # limit irrigation demand to the maximum irrigation rate
                irri_dem_gross =
                    min(irri_dem_gross, sbm.nonpaddy.maximum_irrigation_rate[i])
            else
                irri_dem_gross = 0.0
            end
            sbm.nonpaddy.demand_gross[i] = irri_dem_gross
        elseif !isnothing(sbm.paddy) && sbm.paddy.irrigation_areas[i]
            if sbm.paddy.irrigation_trigger[i]
                # check if maximum irrigation rate has been applied at the previous time step.
                max_irri_rate_applied =
                    sbm.paddy.demand_gross[i] == sbm.paddy.maximum_irrigation_rate[i]
                # start irrigation
                if sbm.paddy.h[i] < sbm.paddy.h_min[i]
                    irr_depth_paddy = sbm.paddy.h_opt[i] - sbm.paddy.h[i]
                elseif sbm.paddy.h[i] < sbm.paddy.h_opt[i] && max_irri_rate_applied # continue irrigation
                    irr_depth_paddy = sbm.paddy.h_opt[i] - sbm.paddy.h[i]
                else
                    irr_depth_paddy = 0.0
                end
                irri_dem_gross += irr_depth_paddy / sbm.paddy.irrigation_efficiency[i]
                # limit irrigation demand to the maximum irrigation rate
                irri_dem_gross = min(irri_dem_gross, sbm.paddy.maximum_irrigation_rate[i])
            end
            sbm.paddy.demand_gross[i] = irri_dem_gross
        end
        # update gross water demands
        sbm.allocation.irri_demand_gross[i] = irri_dem_gross
        sbm.allocation.nonirri_demand_gross[i] = industry_dem + domestic_dem + livestock_dem
        sbm.allocation.total_gross_demand[i] =
            irri_dem_gross + industry_dem + domestic_dem + livestock_dem
    end
end
