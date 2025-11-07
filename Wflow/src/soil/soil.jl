abstract type AbstractSoilModel end

"Struct for storing SBM soil model variables"
@with_kw struct SbmSoilVariables{N}
    # Calculated soil water pressure head h3 of the root water uptake reduction function (Feddes) [cm]
    h3::Vector{Float64}
    # Unsaturated store capacity [mm]
    ustorecapacity::Vector{Float64}
    # Amount of water in the unsaturated store, per layer [mm]
    ustorelayerdepth::Vector{SVector{N, Float64}}
    # Thickness of unsaturated zone, per layer [mm]
    ustorelayerthickness::Vector{SVector{N, Float64}}
    # Saturated store [mm]
    satwaterdepth::Vector{Float64}
    # Pseudo-water table depth [mm] (top of the saturated zone)
    zi::Vector{Float64}
    # Number of unsaturated soil layers
    n_unsatlayers::Vector{Int}
    # Transpiration [mm Δt⁻¹]
    transpiration::Vector{Float64}
    # Actual evaporation from unsaturated store [mm Δt⁻¹]
    ae_ustore::Vector{Float64}
    # Soil evaporation from unsaturated and saturated store [mm Δt⁻¹]
    soilevap::Vector{Float64}
    # Soil evaporation from saturated store [mm Δt⁻¹]
    soilevapsat::Vector{Float64}
    # Actual capillary rise [mm Δt⁻¹]
    actcapflux::Vector{Float64}
    # Actual transpiration from saturated store [mm Δt⁻¹]
    actevapsat::Vector{Float64}
    # Total actual evapotranspiration [mm Δt⁻¹]
    actevap::Vector{Float64}
    # Actual infiltration into the unsaturated zone [mm Δt⁻¹]
    actinfilt::Vector{Float64}
    # Actual infiltration non-compacted fraction [mm Δt⁻¹]
    actinfiltsoil::Vector{Float64}
    # Actual infiltration compacted fraction [mm Δt⁻¹]
    actinfiltpath::Vector{Float64}
    # Actual infiltration (compacted and the non-compacted areas) [mm Δt⁻¹]
    infiltsoilpath::Vector{Float64}
    # Infiltration excess water [mm Δt⁻¹]
    infiltexcess::Vector{Float64}
    # Water that cannot infiltrate due to saturated soil (saturation excess) [mm Δt⁻¹]
    excesswater::Vector{Float64}
    # Water exfiltrating during saturation excess conditions [mm Δt⁻¹]
    exfiltsatwater::Vector{Float64}
    # Water exfiltrating from unsaturated store because of change in water table [mm Δt⁻¹]
    exfiltustore::Vector{Float64}
    # Excess water for non-compacted fraction [mm Δt⁻¹]
    excesswatersoil::Vector{Float64}
    # Excess water for compacted fraction [mm Δt⁻¹]
    excesswaterpath::Vector{Float64}
    # Total surface runoff from infiltration and saturation excess (excluding actual open water evaporation) [mm Δt⁻¹]
    runoff::Vector{Float64}
    # Net surface runoff (surface runoff - actual open water evaporation) [mm Δt⁻¹]
    net_runoff::Vector{Float64}
    # Volumetric water content [-] per soil layer (including theta_r and saturated zone)
    vwc::Vector{SVector{N, Float64}}
    # Volumetric water content [%] per soil layer (including theta_r and saturated zone)
    vwc_perc::Vector{SVector{N, Float64}}
    # Root water storage [mm] in unsaturated and saturated zone (excluding theta_r)
    rootstore::Vector{Float64}
    # Volumetric water content [-] in root zone (including theta_r and saturated zone)
    vwc_root::Vector{Float64}
    # Volumetric water content [%] in root zone (including theta_r and saturated zone)
    vwc_percroot::Vector{Float64}
    # Amount of available water in the unsaturated zone [mm]
    ustoredepth::Vector{Float64}
    # Downward flux from unsaturated to saturated zone [mm Δt⁻¹]
    transfer::Vector{Float64}
    # Net recharge to saturated store [mm Δt⁻¹]
    recharge::Vector{Float64}
    # Actual leakage from saturated store [mm Δt⁻¹]
    actleakage::Vector{Float64}
    # Total water storage (excluding floodplain volume and reservoirs) [mm]
    total_storage::Vector{Float64}
    # Total soil water storage [mm]
    total_soilwater_storage::Vector{Float64}
    # Top soil temperature [ᵒC]
    tsoil::Vector{Float64}
    # Soil infiltration reduction factor (when soil is frozen) [-]
    f_infiltration_reduction::Vector{Float64}
end

"Struct for storing SBM soil model parameters"
@with_kw struct SbmSoilParameters{N, M, Kv}
    # Maximum number of soil layers [-]
    maxlayers::Int
    # Number of soil layers [-]
    nlayers::Vector{Int}
    # Saturated water content (porosity) [-]
    theta_s::Vector{Float64}
    # Residual water content [-]
    theta_r::Vector{Float64}
    # Soilwater capacity [mm]
    soilwatercapacity::Vector{Float64}
    # Muliplication factor [-] applied to kv_z (vertical flow)
    kvfrac::Vector{SVector{N, Float64}}
    # Air entry pressure [cm] of soil (Brooks-Corey)
    hb::Vector{Float64}
    # Soil thickness [mm]
    soilthickness::Vector{Float64}
    # Thickness of soil layers [mm]
    act_thickl::Vector{SVector{N, Float64}}
    # Cumulative sum of soil layers [mm], starting at soil surface (0)
    sumlayers::Vector{SVector{M, Float64}}
    # Infiltration capacity of the compacted areas [mm Δt⁻¹]
    infiltcappath::Vector{Float64}
    # Soil infiltration capacity [mm Δt⁻¹]
    infiltcapsoil::Vector{Float64}
    # Maximum leakage [mm Δt⁻¹] from saturated zone
    maxleakage::Vector{Float64}
    # Parameter [mm] controlling capillary rise
    cap_hmax::Vector{Float64}
    # Coefficient [-] controlling capillary rise
    cap_n::Vector{Float64}
    # Brooks-Corey power coefﬁcient [-] for each soil layer
    c::Vector{SVector{N, Float64}}
    # Soil temperature smooth factor [-]
    w_soil::Vector{Float64}
    # Controls soil infiltration reduction factor when soil is frozen [-]
    cf_soil::Vector{Float64}
    # Fraction of compacted area  [-]
    pathfrac::Vector{Float64}
    # Controls how roots are linked to water table [-]
    rootdistpar::Vector{Float64}
    # Fraction of the root length density in each soil layer [-]
    rootfraction::Vector{SVector{N, Float64}}
    # Soil water pressure head h1 of the root water uptake reduction function (Feddes) [cm]
    h1::Vector{Float64}
    # Soil water pressure head h2 of the root water uptake reduction function (Feddes) [cm]
    h2::Vector{Float64}
    # Soil water pressure head h3_high of the root water uptake reduction function (Feddes) [cm]
    h3_high::Vector{Float64}
    # Soil water pressure head h3_low of the root water uptake reduction function (Feddes) [cm]
    h3_low::Vector{Float64}
    # Soil water pressure head h4 of the root water uptake reduction function (Feddes) [cm]
    h4::Vector{Float64}
    # Root water uptake reduction at soil water pressure head h1 (0.0 or 1.0) [-]
    alpha_h1::Vector{Float64}
    # Soil fraction [-]
    soil_fraction::Vector{Float64}
    # Vertical hydraulic conductivity profile type
    kv_profile::Kv
    # Vegetation parameter set
    vegetation_parameter_set::VegetationParameters
end

"Initialize SBM soil model variables"
function SbmSoilVariables(n::Int, parameters::SbmSoilParameters)
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
    ustoredepth = zeros(Float64, n)
    zi = @. max(0.0, soilthickness - satwaterdepth / (theta_s - theta_r))
    ustorelayerthickness = set_layerthickness.(zi, sumlayers, act_thickl)
    n_unsatlayers = number_of_active_layers.(ustorelayerthickness)

    vwc = fill(MISSING_VALUE, maxlayers, n)
    vwc_perc = fill(MISSING_VALUE, maxlayers, n)
    total_soilwater_storage = satwaterdepth .+ ustoredepth

    vars = SbmSoilVariables(;
        ustorelayerdepth = zero(act_thickl),
        ustorecapacity = soilwatercapacity .- satwaterdepth,
        ustorelayerthickness,
        ustoredepth = zeros(n),
        satwaterdepth,
        zi,
        n_unsatlayers,
        transpiration = fill(MISSING_VALUE, n),
        h3 = fill(MISSING_VALUE, n),
        ae_ustore = fill(MISSING_VALUE, n),
        soilevap = fill(MISSING_VALUE, n),
        soilevapsat = fill(MISSING_VALUE, n),
        actcapflux = fill(MISSING_VALUE, n),
        actevapsat = fill(MISSING_VALUE, n),
        actevap = fill(MISSING_VALUE, n),
        f_infiltration_reduction = fill(1.0, n),
        actinfilt = fill(MISSING_VALUE, n),
        actinfiltsoil = fill(MISSING_VALUE, n),
        actinfiltpath = fill(MISSING_VALUE, n),
        infiltsoilpath = fill(MISSING_VALUE, n),
        infiltexcess = fill(MISSING_VALUE, n),
        excesswater = fill(MISSING_VALUE, n),
        exfiltsatwater = fill(MISSING_VALUE, n),
        exfiltustore = fill(MISSING_VALUE, n),
        excesswatersoil = fill(MISSING_VALUE, n),
        excesswaterpath = fill(MISSING_VALUE, n),
        runoff = fill(MISSING_VALUE, n),
        net_runoff = fill(MISSING_VALUE, n),
        vwc = svectorscopy(vwc, Val{maxlayers}()),
        vwc_perc = svectorscopy(vwc_perc, Val{maxlayers}()),
        rootstore = fill(MISSING_VALUE, n),
        vwc_root = fill(MISSING_VALUE, n),
        vwc_percroot = fill(MISSING_VALUE, n),
        transfer = fill(MISSING_VALUE, n),
        recharge = fill(MISSING_VALUE, n),
        actleakage = fill(MISSING_VALUE, n),
        tsoil = fill(10.0, n),
        total_storage = zeros(Float64, n),
        total_soilwater_storage,
    )
    return vars
end

"Struct for storing SBM soil model boundary conditions"
@with_kw struct SbmSoilBC
    n::Int
    # Water flux at the soil surface [mm Δt⁻¹]
    water_flux_surface::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential transpiration rate [mm Δt⁻¹]
    potential_transpiration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential soil evaporation rate [mm Δt⁻¹]
    potential_soilevaporation::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Exponential depth profile of vertical hydraulic conductivity at the soil surface"
struct KvExponential
    # Vertical hydraulic conductivity [mm Δt⁻¹] at soil surface
    kv_0::Vector{Float64}
    # A scaling parameter [mm⁻¹] (controls exponential decline of kv_0)
    f::Vector{Float64}
end

"Exponential constant depth profile of vertical hydraulic conductivity"
struct KvExponentialConstant
    exponential::KvExponential
    # Depth [mm] from soil surface for which exponential decline of kv_0 is valid
    z_exp::Vector{Float64}
end

"Layered depth profile of vertical hydraulic conductivity"
struct KvLayered{N}
    # Vertical hydraulic conductivity [mm Δt⁻¹] per soil layer
    kv::Vector{SVector{N, Float64}}
end

"Layered exponential depth profile of vertical hydraulic conductivity"
struct KvLayeredExponential{N}
    # A scaling parameter [mm⁻¹] (controls exponential decline of kv_0)
    f::Vector{Float64}
    # Vertical hydraulic conductivity [mm Δt⁻¹] per soil layer
    kv::Vector{SVector{N, Float64}}
    # Number of soil layers [-] with vertical hydraulic conductivity value `kv`
    nlayers_kv::Vector{Int}
    # Depth [mm] from soil surface for which layered profile is valid
    z_layered::Vector{Float64}
end

"Initialize SBM soil model hydraulic conductivity depth profile"
function sbm_kv_profiles(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    kv_0::Vector{Float64},
    f::Vector{Float64},
    maxlayers::Int,
    nlayers::Vector{Int},
    sumlayers::Vector,
    dt::Second,
)
    kv_profile_type = config.model.saturated_hydraulic_conductivity_profile
    n = length(indices)
    if kv_profile_type == VerticalConductivityProfile.exponential
        kv_profile = KvExponential(kv_0, f)
    elseif kv_profile_type == VerticalConductivityProfile.exponential_constant
        z_exp = ncread(
            dataset,
            config,
            "soil_exponential_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth";
            optional = false,
            sel = indices,
            type = Float64,
        )
        exp_profile = KvExponential(kv_0, f)
        kv_profile = KvExponentialConstant(exp_profile, z_exp)
    elseif kv_profile_type == VerticalConductivityProfile.layered ||
           kv_profile_type == VerticalConductivityProfile.layered_exponential
        kv = ncread(
            dataset,
            config,
            "soil_layer_water__vertical_saturated_hydraulic_conductivity";
            optional = false,
            sel = indices,
            type = Float64,
            dimname = :layer,
        )
        if size(kv, 1) != maxlayers
            parname = param(
                config.input.static,
                "soil_layer_water__vertical_saturated_hydraulic_conductivity",
            )
            size1 = size(kv, 1)
            error("$parname needs a layer dimension of size $maxlayers, but is $size1")
        end
        if kv_profile_type == VerticalConductivityProfile.layered
            kv_profile = KvLayered(svectorscopy(kv, Val{maxlayers}()))
        else
            z_layered = ncread(
                dataset,
                config,
                "soil_layered_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth";
                optional = false,
                sel = indices,
                type = Float64,
            )
            nlayers_kv = fill(0, n)
            for i in eachindex(nlayers_kv)
                layers = @view sumlayers[i][2:nlayers[i]]
                _, k = findmin(abs.(z_layered[i] .- layers))
                nlayers_kv[i] = k
                z_layered[i] = layers[k]
            end
            kv_profile = KvLayeredExponential(
                f,
                svectorscopy(kv, Val{maxlayers}()),
                nlayers_kv,
                z_layered,
            )
        end
    end
    return kv_profile
end

"Initialize SBM soil model parameters"
function SbmSoilParameters(
    dataset::NCDataset,
    config::Config,
    vegetation_parameter_set::VegetationParameters,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    config_soil_layer_thickness = config.model.soil_layer__thickness

    soil_layer_thickness =
        SVector(Tuple(push!(Float64.(config_soil_layer_thickness), MISSING_VALUE)))
    cum_depth_layers = pushfirst(cumsum(soil_layer_thickness), 0.0)
    maxlayers = length(soil_layer_thickness) # max number of soil layers

    @info "Using `$(maxlayers - 1)` soil layers with the following thickness: `$config_soil_layer_thickness`"

    w_soil = ncread(
        dataset,
        config,
        "soil_surface_temperature__weight_coefficient",
        LandHydrologySBM;
        sel = indices,
        defaults = 0.1125,
        type = Float64,
    )
    cf_soil = ncread(
        dataset,
        config,
        "soil_surface_water__infiltration_reduction_parameter",
        LandHydrologySBM;
        sel = indices,
        defaults = 0.038,
        type = Float64,
    )

    # soil parameters
    theta_s = ncread(
        dataset,
        config,
        "soil_water__saturated_volume_fraction",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Float64,
    )
    theta_r = ncread(
        dataset,
        config,
        "soil_water__residual_volume_fraction",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Float64,
    )
    kv_0 = ncread(
        dataset,
        config,
        "soil_surface_water__vertical_saturated_hydraulic_conductivity",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Float64,
    )
    f = ncread(
        dataset,
        config,
        "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Float64,
    )
    hb = ncread(
        dataset,
        config,
        "soil_water__air_entry_pressure_head",
        LandHydrologySBM;
        sel = indices,
        defaults = -10.0,
        type = Float64,
    )
    h1 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h1",
        LandHydrologySBM;
        sel = indices,
        defaults = 0.0,
        type = Float64,
    )
    h2 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h2",
        LandHydrologySBM;
        sel = indices,
        defaults = -100.0,
        type = Float64,
    )
    h3_high = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h3_high",
        LandHydrologySBM;
        sel = indices,
        defaults = -400.0,
        type = Float64,
    )
    h3_low = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h3_low",
        LandHydrologySBM;
        sel = indices,
        defaults = -1000.0,
        type = Float64,
    )
    h4 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h4",
        LandHydrologySBM;
        sel = indices,
        defaults = -16000.0,
        type = Float64,
    )
    alpha_h1 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h1_reduction_coefficient",
        LandHydrologySBM;
        sel = indices,
        defaults = 1.0,
        type = Float64,
    )
    soilthickness = ncread(
        dataset,
        config,
        "soil__thickness",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Float64,
    )
    infiltcappath = ncread(
        dataset,
        config,
        "compacted_soil_surface_water__infiltration_capacity",
        LandHydrologySBM;
        sel = indices,
        defaults = 10.0,
        type = Float64,
    )
    maxleakage = ncread(
        dataset,
        config,
        "soil_water_saturated_zone_bottom__max_leakage_volume_flux",
        LandHydrologySBM;
        sel = indices,
        defaults = 0.0,
        type = Float64,
    )
    c = ncread(
        dataset,
        config,
        "soil_layer_water__brooks_corey_exponent",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Float64,
        dimname = :layer,
    )
    if size(c, 1) != maxlayers
        parname = param(config.input.static, "soil_layer_water__brooks_corey_exponent")
        size1 = size(c, 1)
        error("$parname needs a layer dimension of size $maxlayers, but is $size1")
    end
    kvfrac = ncread(
        dataset,
        config,
        "soil_layer_water__vertical_saturated_hydraulic_conductivity_factor",
        LandHydrologySBM;
        sel = indices,
        defaults = 1.0,
        type = Float64,
        dimname = :layer,
    )
    if size(kvfrac, 1) != maxlayers
        parname = param(
            config.input.static,
            "soil_layer_water__vertical_saturated_hydraulic_conductivity_factor",
        )
        size1 = size(kvfrac, 1)
        error("$parname needs a layer dimension of size $maxlayers, but is $size1")
    end

    # soil infiltration capacity based on kv_0 and kvfrac upper soil layer
    infiltcapsoil = kv_0 .* @view kvfrac[1, :]
    # fraction compacted area
    pathfrac = ncread(
        dataset,
        config,
        "compacted_soil__area_fraction",
        LandHydrologySBM;
        optional = false,
        sel = indices,
        type = Float64,
    )

    # vegetation parameters
    rootdistpar = ncread(
        dataset,
        config,
        "soil_wet_root__sigmoid_function_shape_parameter",
        LandHydrologySBM;
        sel = indices,
        defaults = -500.0,
        type = Float64,
    )
    cap_hmax = ncread(
        dataset,
        config,
        "soil_water_saturated_zone_top__capillary_rise_max_water_table_depth",
        LandHydrologySBM;
        sel = indices,
        defaults = 2000.0,
        type = Float64,
    )
    cap_n = ncread(
        dataset,
        config,
        "soil_water_saturated_zone_top__capillary_rise_averianov_exponent",
        LandHydrologySBM;
        sel = indices,
        defaults = 2.0,
        type = Float64,
    )

    act_thickl =
        set_layerthickness.(soilthickness, (cum_depth_layers,), (soil_layer_thickness,))
    sumlayers = @. pushfirst(cumsum(act_thickl), 0.0)
    nlayers = number_of_active_layers.(act_thickl)

    # root fraction read from dataset file, in case of multiple soil layers and TOML file
    # includes "soil_root__length_density_fraction"
    par_name = "soil_root__length_density_fraction"
    do_root_fraction =
        do_cyclic(config) ? haskey(config.input.cyclic, par_name) :
        haskey(config.input.static, par_name)
    if do_root_fraction
        rootfraction = ncread(
            dataset,
            config,
            par_name,
            LandHydrologySBM;
            optional = false,
            sel = indices,
            type = Float64,
            dimname = :layer,
        )
    else
        n = length(indices)
        (; rootingdepth) = vegetation_parameter_set
        # default root fraction in case of multiple soil layers
        rootfraction = zeros(Float64, maxlayers, n)
        for i in 1:n
            if rootingdepth[i] > 0.0
                for k in 1:maxlayers
                    if (rootingdepth[i] - sumlayers[i][k]) >= act_thickl[i][k]
                        rootfraction[k, i] = act_thickl[i][k] / rootingdepth[i]
                    else
                        rootfraction[k, i] =
                            max(rootingdepth[i] - sumlayers[i][k], 0.0) / rootingdepth[i]
                    end
                end
            end
        end
    end

    kv_profile = sbm_kv_profiles(
        dataset,
        config,
        indices,
        kv_0,
        f,
        maxlayers,
        nlayers,
        sumlayers,
        dt,
    )

    soilwatercapacity = @. soilthickness * (theta_s - theta_r)

    n = length(indices)
    sbm_params = SbmSoilParameters(;
        maxlayers,
        nlayers,
        soilwatercapacity,
        theta_s,
        theta_r,
        kvfrac = svectorscopy(kvfrac, Val{maxlayers}()),
        hb,
        h1,
        h2,
        h3_high,
        h3_low,
        h4,
        alpha_h1,
        soilthickness,
        act_thickl,
        sumlayers,
        infiltcappath,
        infiltcapsoil,
        maxleakage,
        pathfrac,
        rootdistpar,
        rootfraction = svectorscopy(rootfraction, Val{maxlayers}()),
        cap_hmax,
        cap_n,
        c = svectorscopy(c, Val{maxlayers}()),
        w_soil,
        cf_soil,
        soil_fraction = fill(MISSING_VALUE, n),
        kv_profile,
        vegetation_parameter_set,
    )
    return sbm_params
end

"SBM soil model"
@with_kw struct SbmSoilModel{N, M, Kv} <: AbstractSoilModel
    boundary_conditions::SbmSoilBC
    parameters::SbmSoilParameters{N, M, Kv}
    variables::SbmSoilVariables{N}
end

"Initialize SBM soil model"
function SbmSoilModel(
    dataset::NCDataset,
    config::Config,
    vegetation_parameter_set::VegetationParameters,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    n = length(indices)
    params = SbmSoilParameters(dataset, config, vegetation_parameter_set, indices, dt)
    vars = SbmSoilVariables(n, params)
    bc = SbmSoilBC(; n)
    model = SbmSoilModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

"Return soil fraction"
function soil_fraction!(
    soil::AbstractSoilModel,
    glacier::AbstractGlacierModel,
    parameters::LandParameters,
)
    (; canopygapfraction) = soil.parameters.vegetation_parameter_set
    (; soil_fraction) = soil.parameters
    (; water_fraction, river_fraction) = parameters
    glacier_fraction = get_glacier_fraction(glacier)

    @. soil_fraction =
        max(canopygapfraction - water_fraction - river_fraction - glacier_fraction, 0.0)
    return nothing
end

"Update boundary conditions of the SBM soil model for a single timestep"
function update_boundary_conditions!(
    model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    external_models::NamedTuple,
)
    (; interception, runoff, demand, allocation) = external_models
    (; potential_transpiration, water_flux_surface, potential_soilevaporation) =
        model.boundary_conditions

    potential_transpiration .= get_potential_transpiration(interception)

    @. potential_soilevaporation =
        model.parameters.soil_fraction * atmospheric_forcing.potential_evaporation
    evaporation!(demand.paddy, potential_soilevaporation)
    potential_soilevaporation .= potential_soilevaporation .- get_evaporation(demand.paddy)

    water_flux_surface .=
        max.(
            runoff.boundary_conditions.water_flux_surface .+
            get_irrigation_allocated(allocation) .- runoff.variables.runoff_river .-
            runoff.variables.runoff_land .+ get_water_depth(demand.paddy),
            0.0,
        )
    return nothing
end

"Update soil temperature of the SBM soil model for a single timestep"
function soil_temperature!(
    model::SbmSoilModel,
    snow::AbstractSnowModel,
    temperature::Vector{Float64},
)
    v = model.variables
    p = model.parameters
    @. v.tsoil = soil_temperature(v.tsoil, p.w_soil, temperature)
    return nothing
end

soil_temperature!(model::SbmSoilModel, snow::NoSnowModel, temperature::Vector{Float64}) =
    nothing

"Update total available water in the unsaturated zone of the SBM soil model for a single timestep"
function ustoredepth!(model::SbmSoilModel)
    v = model.variables
    p = model.parameters
    for i in eachindex(v.ustorelayerdepth)
        v.ustoredepth[i] = sum(@view v.ustorelayerdepth[i][1:p.nlayers[i]])
    end
    return nothing
end

"Update the infiltration reduction factor of the SBM soil model for a single timestep"
function infiltration_reduction_factor!(
    model::SbmSoilModel;
    modelsnow = false,
    soil_infiltration_reduction = false,
)
    v = model.variables
    p = model.parameters

    n = length(v.tsoil)
    threaded_foreach(1:n; basesize = 1000) do i
        v.f_infiltration_reduction[i] = infiltration_reduction_factor(
            v.tsoil[i],
            p.cf_soil[i];
            modelsnow,
            soil_infiltration_reduction,
        )
    end
    return nothing
end

"""
    infiltration!(model::SbmSoilModel)

Update the infiltration rate `infiltsoilpath` and infiltration excess water rate
`infiltexcess` of the SBM soil model for a single timestep.
"""
function infiltration!(model::SbmSoilModel)
    v = model.variables
    p = model.parameters
    (; water_flux_surface) = model.boundary_conditions

    n = length(v.infiltsoilpath)
    threaded_foreach(1:n; basesize = 1000) do i
        v.infiltsoilpath[i], v.infiltexcess[i] = infiltration(
            water_flux_surface[i],
            p.pathfrac[i],
            p.infiltcapsoil[i],
            p.infiltcappath[i],
            v.ustorecapacity[i],
            v.f_infiltration_reduction[i],
        )
    end
    return nothing
end

"""
    unsaturated_zone_flow!(model::SbmSoilModel)

Update unsaturated storage `ustorelayerdepth` and the `transfer` of water from the unsaturated
to the saturated store of the SBM soil model for a single timestep, based on the Brooks-Corey
approach.
"""
function unsaturated_zone_flow!(model::SbmSoilModel)
    v = model.variables
    p = model.parameters

    n = length(v.transfer)
    threaded_foreach(1:n; basesize = 250) do i
        if v.n_unsatlayers[i] > 0
            # Brooks-Corey approach
            z = cumsum(v.ustorelayerthickness[i])
            flow_rate = 0.0
            for m in 1:v.n_unsatlayers[i]
                l_sat = v.ustorelayerthickness[i][m] * (p.theta_s[i] - p.theta_r[i])
                kv_z = hydraulic_conductivity_at_depth(p.kv_profile, p.kvfrac, z[m], i, m)
                ustorelayerdepth = if m == 1
                    v.ustorelayerdepth[i][m] + v.infiltsoilpath[i]
                else
                    v.ustorelayerdepth[i][m] + flow_rate
                end
                ustorelayerdepth, flow_rate =
                    unsatzone_flow_layer(ustorelayerdepth, kv_z, l_sat, p.c[i][m])
                v.ustorelayerdepth[i] = setindex(v.ustorelayerdepth[i], ustorelayerdepth, m)
            end
            v.transfer[i] = flow_rate
        else
            v.transfer[i] = 0.0
        end
    end
    return nothing
end

"""
    soil_evaporation!(model::SbmSoilModel)

Update soil evaporation from the saturated store `soilevapsat` and the total soil
evaporation from the unsaturated and saturated store `soilevap` of the SBM soil model for a
single timestep. Also unsaturated storage `ustorelayerdepth` and the saturated store
`satwaterdepth` are updated.
"""
function soil_evaporation!(model::SbmSoilModel)
    (; potential_soilevaporation) = model.boundary_conditions
    v = model.variables
    p = model.parameters

    n = length(potential_soilevaporation)
    threaded_foreach(1:n; basesize = 1000) do i
        potsoilevap = potential_soilevaporation[i]
        # First calculate the evaporation of unsaturated storage into the
        # atmosphere from the upper layer.
        soilevapunsat = soil_evaporation_unsatured_store(
            potsoilevap,
            v.ustorelayerdepth[i][1],
            v.ustorelayerthickness[i][1],
            v.n_unsatlayers[i],
            v.zi[i],
            p.theta_s[i] - p.theta_r[i],
        )
        # Ensure that the unsaturated evaporation rate does not exceed the
        # available unsaturated moisture
        soilevapunsat = min(soilevapunsat, v.ustorelayerdepth[i][1])
        # Update the additional atmospheric demand
        potsoilevap -= soilevapunsat
        v.ustorelayerdepth[i] =
            setindex(v.ustorelayerdepth[i], v.ustorelayerdepth[i][1] - soilevapunsat, 1)

        soilevapsat = soil_evaporation_satured_store(
            potsoilevap,
            v.n_unsatlayers[i],
            p.act_thickl[i][1],
            v.zi[i],
            p.theta_s[i] - p.theta_r[i],
        )

        v.soilevapsat[i] = soilevapsat
        v.soilevap[i] = soilevapunsat + soilevapsat
        v.satwaterdepth[i] = v.satwaterdepth[i] - soilevapsat
    end
    return nothing
end

"""
    transpiration!(model::SbmSoilModel, dt)

Update total `transpiration`, transpiration from the unsaturated store `ae_ustore` and
saturated store `actevapsat` of the SBM soil model for a single timestep. Also unsaturated
storage `ustorelayerdepth` and the saturated store `satwaterdepth` are updated.
"""
function transpiration!(model::SbmSoilModel, dt::Float64)
    (; potential_transpiration) = model.boundary_conditions
    v = model.variables
    p = model.parameters

    rootingdepth = get_rootingdepth(model)
    n = length(rootingdepth)

    threaded_foreach(1:n; basesize = 250) do i
        v.h3[i] = feddes_h3(p.h3_high[i], p.h3_low[i], potential_transpiration[i], dt)

        # compute sum of root fraction in unsaturated soil layers and adapt root fraction
        # lowest unsaturated soil layer if water table depth intersects the unsaturated root
        # zone
        sum_rootfraction_unsat = 0.0
        rootfraction_unsat_lowest = 0.0
        for k in 1:v.n_unsatlayers[i]
            # the root fraction is valid for the root length in a soil layer, if zi decreases
            # the root length the root fraction needs to be adapted
            if k == v.n_unsatlayers[i] && v.zi[i] < rootingdepth[i]
                rootlength = min(p.act_thickl[i][k], rootingdepth[i] - p.sumlayers[i][k])
                rootfraction_unsat =
                    p.rootfraction[i][k] * (v.ustorelayerthickness[i][k] / rootlength)
                sum_rootfraction_unsat += rootfraction_unsat
            else
                rootfraction_unsat = p.rootfraction[i][k]
                sum_rootfraction_unsat += rootfraction_unsat
            end
            # rootfraction lowest unsaturated layer
            rootfraction_unsat_lowest = rootfraction_unsat
        end

        actevapustore = 0.0
        for k in 1:v.n_unsatlayers[i]
            # scale rootfraction soil layer unsaturated zone based on sum of rootfraction in
            # unsaturated zone
            if k < v.n_unsatlayers[i]
                rootfraction_unsat = p.rootfraction[i][k]
            else
                rootfraction_unsat = rootfraction_unsat_lowest
            end
            rootfraction_unsat_scaled =
                rootingdepth[i] > 0.0 ?
                max((1.0 / sum_rootfraction_unsat), 1.0) * rootfraction_unsat : 0.0

            vwc = max(
                v.ustorelayerdepth[i][k] / v.ustorelayerthickness[i][k],
                Float64(0.0000001),
            )
            head = head_brooks_corey(vwc, p.theta_s[i], p.theta_r[i], p.c[i][k], p.hb[i])
            alpha = rwu_reduction_feddes(
                head,
                p.h1[i],
                p.h2[i],
                v.h3[i],
                p.h4[i],
                p.alpha_h1[i],
            )

            availcap = min(
                1.0,
                max(
                    0.0,
                    (rootingdepth[i] - p.sumlayers[i][k]) / v.ustorelayerthickness[i][k],
                ),
            )
            maxextr = v.ustorelayerdepth[i][k] * availcap
            actevapustore_layer =
                min(alpha * rootfraction_unsat_scaled * potential_transpiration[i], maxextr)
            ustorelayerdepth = v.ustorelayerdepth[i][k] - actevapustore_layer
            actevapustore += actevapustore_layer
            v.ustorelayerdepth[i] = setindex(v.ustorelayerdepth[i], ustorelayerdepth, k)
        end

        # transpiration from saturated store
        wetroots = scurve(v.zi[i], rootingdepth[i], Float64(1.0), p.rootdistpar[i])
        alpha = rwu_reduction_feddes(
            Float64(0.0),
            p.h1[i],
            p.h2[i],
            v.h3[i],
            p.h4[i],
            p.alpha_h1[i],
        )
        restpottrans = potential_transpiration[i] - actevapustore
        actevapsat = min(restpottrans * wetroots * alpha, v.satwaterdepth[i])

        v.ae_ustore[i] = actevapustore
        v.actevapsat[i] = actevapsat
        v.satwaterdepth[i] = v.satwaterdepth[i] - actevapsat
        v.transpiration[i] = actevapustore + actevapsat
    end
    return nothing
end

"""
    actual_infiltration!(model::SbmSoilModel)

Update the actual infiltration rate `actinfilt` of the SBM soil model for a single timestep.

A soil water balance check is performed. Unsaturated storage that exceeds the maximum
storage per unsaturated soil layer is transferred to the layer above (or surface), from the
bottom to the top unsaturated soil layer. The resulting excess water `ustoredepth_excess` is
subtracted from the infiltration rate `infiltsoilpath`.
"""
function actual_infiltration!(model::SbmSoilModel)
    v = model.variables
    p = model.parameters

    n = length(v.actinfilt)
    threaded_foreach(1:n; basesize = 1000) do i
        # check soil moisture balance per layer
        ustoredepth_excess = 0.0
        for k in v.n_unsatlayers[i]:-1:1
            ustoredepth_excess = max(
                0.0,
                v.ustorelayerdepth[i][k] -
                v.ustorelayerthickness[i][k] * (p.theta_s[i] - p.theta_r[i]),
            )
            v.ustorelayerdepth[i] = setindex(
                v.ustorelayerdepth[i],
                v.ustorelayerdepth[i][k] - ustoredepth_excess,
                k,
            )
            if k > 1
                v.ustorelayerdepth[i] = setindex(
                    v.ustorelayerdepth[i],
                    v.ustorelayerdepth[i][k - 1] + ustoredepth_excess,
                    k - 1,
                )
            end
        end

        v.actinfilt[i] = v.infiltsoilpath[i] - ustoredepth_excess
    end
    return nothing
end

"""
    actual_infiltration_soil_path!(model::SbmSoilModel)

Update the actual infiltration rate for soil `actinfiltsoil` and paved area `actinfiltpath`
of the SBM soil model for a single timestep.
"""
function actual_infiltration_soil_path!(model::SbmSoilModel)
    v = model.variables
    p = model.parameters
    (; water_flux_surface) = model.boundary_conditions

    n = length(water_flux_surface)
    threaded_foreach(1:n; basesize = 1000) do i
        v.actinfiltsoil[i], v.actinfiltpath[i] = actual_infiltration_soil_path(
            water_flux_surface[i],
            v.actinfilt[i],
            p.pathfrac[i],
            p.infiltcapsoil[i],
            p.infiltcappath[i],
            v.f_infiltration_reduction[i],
        )
    end
    return nothing
end

"""
    capillary_flux!(model::SbmSoilModel)

Update the capillary flux `actcapflux` of the SBM soil model for a single timestep.
"""
function capillary_flux!(model::SbmSoilModel)
    v = model.variables
    p = model.parameters
    rootingdepth = get_rootingdepth(model)

    n = length(rootingdepth)
    threaded_foreach(1:n; basesize = 1000) do i
        if v.n_unsatlayers[i] > 0
            ksat = hydraulic_conductivity_at_depth(
                p.kv_profile,
                p.kvfrac,
                v.zi[i],
                i,
                v.n_unsatlayers[i],
            )
            maxcapflux =
                max(0.0, min(ksat, v.ae_ustore[i], v.ustorecapacity[i], v.satwaterdepth[i]))

            if v.zi[i] > rootingdepth[i]
                capflux =
                    maxcapflux *
                    pow(1.0 - min(v.zi[i], p.cap_hmax[i]) / (p.cap_hmax[i]), p.cap_n[i])
            else
                capflux = 0.0
            end

            netcapflux = capflux
            actcapflux = 0.0
            for k in v.n_unsatlayers[i]:-1:1
                toadd = min(
                    netcapflux,
                    max(
                        v.ustorelayerthickness[i][k] * (p.theta_s[i] - p.theta_r[i]) -
                        v.ustorelayerdepth[i][k],
                        0.0,
                    ),
                )
                v.ustorelayerdepth[i] =
                    setindex(v.ustorelayerdepth[i], v.ustorelayerdepth[i][k] + toadd, k)
                netcapflux -= toadd
                actcapflux += toadd
            end
            v.actcapflux[i] = actcapflux
        else
            v.actcapflux[i] = 0.0
        end
    end
    return nothing
end

"""
    leakage!(model::SbmSoilModel)

Update the actual leakage rate `actleakage` of the SBM soil model for a single timestep.
"""
function leakage!(model::SbmSoilModel)
    v = model.variables
    p = model.parameters

    n = length(v.actleakage)
    threaded_foreach(1:n; basesize = 1000) do i
        deepksat = hydraulic_conductivity_at_depth(
            p.kv_profile,
            p.kvfrac,
            p.soilthickness[i],
            i,
            p.nlayers[i],
        )
        deeptransfer = min(v.satwaterdepth[i], deepksat)
        v.actleakage[i] = max(0.0, min(p.maxleakage[i], deeptransfer))
    end
    return nothing
end

"""
    update!(
        model::SbmSoilModel,
        atmospheric_forcing::AtmosphericForcing,
        external_models::NamedTuple,
        config::Config,
        dt::Float64,
    )

Update the SBM soil model (infiltration, unsaturated zone flow, soil evaporation and
transpiration, capillary flux and leakage) for a single timestep.
"""
function update!(
    model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    external_models::NamedTuple,
    config::Config,
    dt::Float64,
)
    (; snow, runoff, demand) = external_models
    (; temperature) = atmospheric_forcing
    (; water_flux_surface) = model.boundary_conditions
    v = model.variables
    p = model.parameters

    # mainly required for external state changes (e.g. through BMI)
    update_diagnostic_vars!(model)
    # infiltration
    soil_temperature!(model, snow, temperature)
    infiltration_reduction_factor!(
        model;
        modelsnow = config.model.snow__flag,
        soil_infiltration_reduction = config.model.soil_infiltration_reduction__flag,
    )
    infiltration!(model)
    # unsaturated zone flow
    unsaturated_zone_flow!(model)
    # soil evaporation and transpiration
    soil_evaporation!(model)
    transpiration!(model, dt)
    # actual infiltration and excess water
    actual_infiltration!(model)
    @. v.excesswater = water_flux_surface - v.actinfilt - v.infiltexcess
    actual_infiltration_soil_path!(model)
    @. v.excesswatersoil =
        max(water_flux_surface * (1.0 - p.pathfrac) - v.actinfiltsoil, 0.0)
    @. v.excesswaterpath = max(water_flux_surface * p.pathfrac - v.actinfiltpath, 0.0)
    # recompute the unsaturated store and ustorecapacity (for capillary flux)
    ustoredepth!(model)
    @. v.ustorecapacity = p.soilwatercapacity - v.satwaterdepth - v.ustoredepth
    # capillary flux and leakage
    capillary_flux!(model)
    leakage!(model)
    # recharge rate to the saturated store
    @. v.recharge =
        (v.transfer - v.actcapflux - v.actleakage - v.actevapsat - v.soilevapsat)
    # total actual evapotranspiration
    v.actevap .=
        v.soilevap .+ v.transpiration .+ runoff.variables.ae_openw_r .+
        runoff.variables.ae_openw_l .+ get_evaporation(demand.paddy)
    return nothing
end

"""
    update!(model::SbmSoilModel, external_models::NamedTuple)

Update the SBM soil model for a single timestep based on the update of a subsurface flow
model, resulting in a change in water table depth and an exfiltration rate `exfiltwater`.

The water table depth `zi`, unsaturated storage `ustorelayerdepth`, water exfiltrating from
the unsaturated store `exfiltustore`, land `runoff` and `net_runoff`, the saturated store
`satwaterdepth` and the water exfiltrating during saturation excess conditions
`exfiltsatwater` are updated. Addionally, volumetric water content per soil layer and for
the root zone are updated.
"""
function update!(model::SbmSoilModel, external_models::NamedTuple)
    (; runoff, demand, subsurface_flow) = external_models
    (; runoff_land, ae_openw_l) = runoff.variables
    p = model.parameters
    v = model.variables

    zi = get_water_depth(subsurface_flow)
    exfiltsatwater = get_exfiltwater(subsurface_flow)
    rootingdepth = get_rootingdepth(model)

    n = length(model.variables.zi)
    threaded_foreach(1:n; basesize = 1000) do i
        ustorelayerthickness = set_layerthickness(zi[i], p.sumlayers[i], p.act_thickl[i])
        n_unsatlayers = number_of_active_layers(ustorelayerthickness)
        # exfiltration from ustore
        ustorelayerdepth = v.ustorelayerdepth[i]
        exfiltustore = 0.0
        for k in v.n_unsatlayers[i]:-1:1
            if k <= n_unsatlayers
                exfiltustore = max(
                    0,
                    ustorelayerdepth[k] -
                    ustorelayerthickness[k] * (p.theta_s[i] - p.theta_r[i]),
                )
            else
                exfiltustore = ustorelayerdepth[k]
            end
            ustorelayerdepth =
                setindex(ustorelayerdepth, ustorelayerdepth[k] - exfiltustore, k)
            if k > 1
                ustorelayerdepth = setindex(
                    ustorelayerdepth,
                    ustorelayerdepth[k - 1] + exfiltustore,
                    k - 1,
                )
            end
        end

        ustoredepth = sum(@view ustorelayerdepth[1:n_unsatlayers])
        sbm_runoff = max(
            0.0,
            exfiltustore +
            exfiltsatwater[i] +
            v.excesswater[i] +
            runoff_land[i] +
            v.infiltexcess[i],
        )

        # volumetric water content per soil layer and root zone
        vwc = v.vwc[i]
        vwc_perc = v.vwc_perc[i]
        for k in 1:p.nlayers[i]
            if k <= n_unsatlayers
                vwc = setindex(
                    vwc,
                    (
                        ustorelayerdepth[k] +
                        (p.act_thickl[i][k] - ustorelayerthickness[k]) *
                        (p.theta_s[i] - p.theta_r[i])
                    ) / p.act_thickl[i][k] + p.theta_r[i],
                    k,
                )
            else
                vwc = setindex(vwc, p.theta_s[i], k)
            end
            vwc_perc = setindex(vwc_perc, (vwc[k] / p.theta_s[i]) * 100.0, k)
        end

        rootstore_unsat = 0
        for k in 1:n_unsatlayers
            rootstore_unsat =
                rootstore_unsat +
                min(
                    1.0,
                    (
                        max(0.0, rootingdepth[i] - p.sumlayers[i][k]) /
                        ustorelayerthickness[k]
                    ),
                ) * ustorelayerdepth[k]
        end

        rootstore_sat = max(0.0, rootingdepth[i] - zi[i]) * (p.theta_s[i] - p.theta_r[i])
        rootstore = rootstore_sat + rootstore_unsat
        vwc_root = rootstore / rootingdepth[i] + p.theta_r[i]
        vwc_percroot = (vwc_root / p.theta_s[i]) * 100.0

        satwaterdepth = (p.soilthickness[i] - zi[i]) * (p.theta_s[i] - p.theta_r[i])
        ustorecapacity = p.soilwatercapacity[i] - satwaterdepth - ustoredepth

        # update the outputs and states
        v.n_unsatlayers[i] = n_unsatlayers
        v.ustorelayerdepth[i] = ustorelayerdepth
        v.ustorelayerthickness[i] = ustorelayerthickness
        v.ustorecapacity[i] = ustorecapacity
        v.n_unsatlayers[i] = n_unsatlayers
        v.ustoredepth[i] = ustoredepth
        v.satwaterdepth[i] = satwaterdepth
        v.exfiltsatwater[i] = exfiltsatwater[i]
        v.exfiltustore[i] = exfiltustore
        v.runoff[i] = sbm_runoff
        v.vwc[i] = vwc
        v.vwc_perc[i] = vwc_perc
        v.rootstore[i] = rootstore
        v.vwc_root[i] = vwc_root
        v.vwc_percroot[i] = vwc_percroot
        v.zi[i] = zi[i]
        v.total_soilwater_storage[i] = satwaterdepth + ustoredepth
    end
    # update runoff and net_runoff (the runoff rate depends on the presence of paddy fields
    # and the h_max parameter of a paddy field)
    update_runoff!(demand.paddy, v.runoff)
    @. v.net_runoff = v.runoff - ae_openw_l
    return nothing
end

"""
    update_diagnostic_vars!(model::SbmSoilModel)

Update diagnostic variables of `SbmSoilModel` that are critical for subsequent computations
and depend on state variables `satwaterdepth` and `ustorelayerdepth`.
"""
function update_diagnostic_vars!(model::SbmSoilModel)
    (;
        zi,
        satwaterdepth,
        ustorelayerthickness,
        ustorecapacity,
        ustoredepth,
        total_soilwater_storage,
        n_unsatlayers,
    ) = model.variables
    (; soilthickness, theta_s, theta_r, soilwatercapacity, sumlayers, act_thickl) =
        model.parameters

    ustoredepth!(model)
    @. zi .= max(0.0, soilthickness - satwaterdepth / (theta_s .- theta_r))
    @. ustorecapacity = soilwatercapacity - satwaterdepth - ustoredepth
    @. ustorelayerthickness = set_layerthickness(zi, sumlayers, act_thickl)
    @. n_unsatlayers = number_of_active_layers(ustorelayerthickness)
    @. total_soilwater_storage = satwaterdepth + ustoredepth
end

# wrapper method
get_rootingdepth(model::SbmSoilModel) =
    model.parameters.vegetation_parameter_set.rootingdepth
