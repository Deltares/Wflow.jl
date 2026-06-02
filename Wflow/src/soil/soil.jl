abstract type AbstractSoilModel end

"Struct for storing SBM soil model variables"
@with_data_lookup struct SbmSoilVariables{N}
    n::Int
    # Calculated soil water pressure head h3 of the root water uptake reduction function (Feddes) [m]
    h3::Vector{Float64} = fill(MISSING_VALUE, n)
    # Unsaturated store capacity [m]
    ustorecapacity::Vector{Float64}
    # Amount of water in the unsaturated store, per layer [m]
    "soil_layer_water_unsaturated_zone__depth"
    ustorelayerdepth::Vector{SVector{N, Float64}}
    # Thickness of unsaturated zone, per layer [m]
    ustorelayerthickness::Vector{SVector{N, Float64}}
    # Saturated store [m]
    "soil_water_saturated_zone__depth"
    satwaterdepth::Vector{Float64}
    # Drainable water store [m]
    drainable_waterdepth::Vector{Float64}
    # Pseudo-water table depth [m] (top of the saturated zone)
    "soil_water_saturated_zone_top__depth"
    zi::Vector{Float64}
    # Number of unsaturated soil layers
    n_unsatlayers::Vector{Int}
    # Transpiration [m s⁻¹]
    "soil_water__transpiration_volume_flux"
    transpiration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual evaporation from unsaturated store [m s⁻¹]
    ae_ustore::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil evaporation from unsaturated and saturated store [m s⁻¹]
    soilevap::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil evaporation from saturated store [m s⁻¹]
    soilevapsat::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual capillary rise [m s⁻¹]
    "soil_water_saturated_zone_top__capillary_volume_flux"
    actcapflux::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual transpiration from saturated store [m s⁻¹]
    actevapsat::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total actual evapotranspiration [m s⁻¹]\
    "land_surface__evapotranspiration_volume_flux"
    actevap::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration into the unsaturated zone [m s⁻¹]
    "soil_water__infiltration_volume_flux"
    actinfilt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration non-compacted fraction [m s⁻¹]
    actinfiltsoil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration compacted fraction [m s⁻¹]
    actinfiltpath::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration (compacted and the non-compacted areas) [m s⁻¹]
    infiltsoilpath::Vector{Float64} = fill(MISSING_VALUE, n)
    # Infiltration excess water [m s⁻¹]
    infiltexcess::Vector{Float64} = fill(MISSING_VALUE, n)
    # Water that cannot infiltrate due to saturated soil (saturation excess) [m s⁻¹]
    excesswater::Vector{Float64} = fill(MISSING_VALUE, n)
    # Water exfiltrating during saturation excess conditions [m s⁻¹]
    "soil_surface_water_saturated_zone__exfiltration_volume_flux"
    exfiltsatwater::Vector{Float64} = fill(MISSING_VALUE, n)
    # Excess water for non-compacted fraction [m s⁻¹]
    "compacted_soil_surface_water__excess_volume_flux"
    excesswatersoil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Excess water for compacted fraction [m s⁻¹]
    "non_compacted_soil_surface_water__excess_volume_flux"
    excesswaterpath::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total surface runoff from infiltration and saturation excess (excluding actual open water evaporation) [m s⁻¹]
    "soil_surface_water__runoff_volume_flux"
    runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Net surface runoff (surface runoff - actual open water evaporation) [m s⁻¹]
    "soil_surface_water__net_runoff_volume_flux"
    net_runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Volumetric water content [-] per soil layer (including theta_r and saturated zone)
    "soil_layer_water__volume_fraction"
    vwc::Vector{SVector{N, Float64}}
    # Volumetric water content [%] per soil layer (including theta_r and saturated zone)
    "soil_layer_water__volume_percentage"
    vwc_perc::Vector{SVector{N, Float64}}
    # Root water storage [m] in unsaturated and saturated zone (excluding theta_r)
    "soil_water_root_zone__depth"
    rootstore::Vector{Float64} = fill(MISSING_VALUE, n)
    # Volumetric water content [-] in root zone (including theta_r and saturated zone)
    "soil_water_root_zone__volume_fraction"
    vwc_root::Vector{Float64} = fill(MISSING_VALUE, n)
    # Volumetric water content [%] in root zone (including theta_r and saturated zone)
    "soil_water_root_zone__volume_percentage"
    vwc_percroot::Vector{Float64} = fill(MISSING_VALUE, n)
    # Amount of available water in the unsaturated zone [m]
    "soil_water_unsaturated_zone__depth"
    ustoredepth::Vector{Float64} = zeros(n)
    # Downward flux from unsaturated to saturated zone [m s⁻¹]
    "soil_water_saturated_zone_top__recharge_volume_flux"
    transfer::Vector{Float64} = fill(MISSING_VALUE, n)
    # Net recharge to saturated store [m s⁻¹]
    "soil_water_saturated_zone_top__net_recharge_volume_flux"
    recharge::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual leakage from saturated store [m s⁻¹]
    "soil_water_saturated_zone_bottom__leakage_volume_flux"
    actleakage::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total water storage (excluding floodplain volume and reservoirs) [m]
    "land_water_storage__total_depth"
    total_storage::Vector{Float64} = zeros(Float64, n)
    # Total soil water storage [m]
    total_soilwater_storage::Vector{Float64}
    # Top soil temperature [K]
    "soil_surface__temperature"
    tsoil::Vector{Float64} = to_SI.(fill(10.0, n), Ref(ABSOLUTE_DEGREES))
    # Soil infiltration reduction factor (when soil is frozen) [-]
    f_infiltration_reduction::Vector{Float64} = ones(n)
end

"Struct for storing SBM soil model parameters"
@with_data_lookup struct SbmSoilParameters{N, M, Kv}
    # Maximum number of soil layers [-]
    maxlayers::Int
    # Number of soil layers [-]
    nlayers::Vector{Int}
    # Saturated water content (porosity) [-]
    "soil_water__saturated_volume_fraction"
    theta_s::Vector{Float64}
    # Residual water content [-]
    "soil_water__residual_volume_fraction"
    theta_r::Vector{Float64}
    # Field capacity water content [-]
    theta_fc::Vector{Float64}
    # Soilwater capacity [m]
    soilwatercapacity::Vector{Float64}
    # Multiplication factor [-] applied to kv_z (vertical flow)
    "soil_layer_water__vertical_saturated_hydraulic_conductivity_factor"
    kvfrac::Vector{SVector{N, Float64}}
    # Air entry pressure [m] of soil (Brooks-Corey)
    "soil_water__air_entry_pressure_head"
    hb::Vector{Float64}
    # Soil thickness [m]
    "soil__thickness"
    soilthickness::Vector{Float64}
    # Thickness of soil layers [m]
    act_thickl::Vector{SVector{N, Float64}}
    # Cumulative sum of soil layers [m], starting at soil surface (0)
    sumlayers::Vector{SVector{M, Float64}}
    # Infiltration capacity of the compacted areas [m s⁻¹]
    "compacted_soil_surface_water__infiltration_capacity"
    infiltcappath::Vector{Float64}
    # Soil infiltration capacity [m s⁻¹]
    infiltcapsoil::Vector{Float64}
    # Maximum leakage [m s⁻¹] from saturated zone
    "soil_water_saturated_zone_bottom__max_leakage_volume_flux"
    maxleakage::Vector{Float64}
    # Parameter [m] controlling capillary rise
    "soil_water_saturated_zone_top__capillary_rise_max_water_table_depth"
    cap_hmax::Vector{Float64}
    # Coefficient [-] controlling capillary rise
    "soil_water_saturated_zone_top__capillary_rise_averianov_exponent"
    cap_n::Vector{Float64}
    # Brooks-Corey power coefficient [-] for each soil layer
    "soil_layer_water__brooks_corey_exponent"
    c::Vector{SVector{N, Float64}}
    # Soil temperature smooth factor [-]
    "soil_surface_temperature__weight_coefficient"
    w_soil::Vector{Float64}
    # Controls soil infiltration reduction factor when soil is frozen [-]
    "soil_surface_water__infiltration_reduction_parameter"
    cf_soil::Vector{Float64}
    # Fraction of compacted area  [-]
    "compacted_soil__area_fraction"
    pathfrac::Vector{Float64}
    # Controls how roots are linked to water table [m⁻¹]
    "soil_wet_root__sigmoid_function_shape_parameter"
    rootdistpar::Vector{Float64}
    # Fraction of the root length density in each soil layer [-]
    "soil_root__length_density_fraction"
    rootfraction::Vector{SVector{N, Float64}}
    # Soil water pressure head h1 of the root water uptake reduction function (Feddes) [m]
    "vegetation_root__feddes_critical_pressure_head_h1"
    h1::Vector{Float64}
    # Soil water pressure head h2 of the root water uptake reduction function (Feddes) [m]
    "vegetation_root__feddes_critical_pressure_head_h2"
    h2::Vector{Float64}
    # Soil water pressure head h3_high of the root water uptake reduction function (Feddes) [m]
    "vegetation_root__feddes_critical_pressure_head_h3_high"
    h3_high::Vector{Float64}
    # Soil water pressure head h3_low of the root water uptake reduction function (Feddes) [m]
    "vegetation_root__feddes_critical_pressure_head_h3_low"
    h3_low::Vector{Float64}
    # Soil water pressure head h4 of the root water uptake reduction function (Feddes) [m]
    "vegetation_root__feddes_critical_pressure_head_h4"
    h4::Vector{Float64}
    # Root water uptake reduction at soil water pressure head h1 (0.0 or 1.0) [-]
    "vegetation_root__feddes_critical_pressure_head_h1_reduction_coefficient"
    alpha_h1::Vector{Float64}
    # Soil fraction [-]
    soil_fraction::Vector{Float64}
    # Vertical hydraulic conductivity profile type
    kv_profile::Kv
    # Vegetation parameter set
    vegetation_parameter_set::VegetationParameters
end

"Initialize SBM soil model variables"
function SbmSoilVariables(
    n::Int,
    parameters::SbmSoilParameters;
    data_lookup::DataLookup = DataLookup(),
)
    (;
        soilthickness,
        maxlayers,
        act_thickl,
        sumlayers,
        soilwatercapacity,
        theta_s,
        theta_r,
        theta_fc,
    ) = parameters
    satwaterdepth = 0.85 .* soilwatercapacity # cold state value for satwaterdepth
    ustoredepth = zeros(n)
    zi = @. max(0.0, soilthickness - satwaterdepth / (theta_s - theta_r))
    drainable_waterdepth =
        @. (soilthickness - zi) * lower_bound_drainable_porosity(theta_s, theta_fc)
    ustorelayerthickness = set_layerthickness.(zi, sumlayers, act_thickl)
    n_unsatlayers = number_of_active_layers.(ustorelayerthickness)
    vwc = fill(MISSING_VALUE, maxlayers, n)
    vwc_perc = fill(MISSING_VALUE, maxlayers, n)
    total_soilwater_storage = satwaterdepth .+ ustoredepth

    vars = SbmSoilVariables(
        data_lookup;
        n,
        ustorelayerdepth = zero(act_thickl),
        ustorecapacity = soilwatercapacity .- satwaterdepth,
        ustorelayerthickness,
        satwaterdepth,
        drainable_waterdepth,
        zi,
        n_unsatlayers,
        vwc = svectorscopy(vwc, Val{maxlayers}()),
        vwc_perc = svectorscopy(vwc_perc, Val{maxlayers}()),
        total_soilwater_storage,
    )
    return vars
end

"Struct for storing SBM soil model boundary conditions"
@kwdef struct SbmSoilBC
    n::Int
    # Water flux at the soil surface [m s⁻¹]
    water_flux_surface::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential transpiration rate [m s⁻¹]
    potential_transpiration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential soil evaporation rate [m s⁻¹]
    potential_soilevaporation::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Exponential depth profile of vertical hydraulic conductivity at the soil surface"
struct KvExponential
    # Vertical hydraulic conductivity [m s⁻¹] at soil surface
    kv_0::Vector{Float64}
    # A scaling parameter [m⁻¹] (controls exponential decline of kv_0)
    f::Vector{Float64}
end

"Exponential constant depth profile of vertical hydraulic conductivity"
struct KvExponentialConstant
    exponential::KvExponential
    # Depth [m] from soil surface for which exponential decline of kv_0 is valid
    z_exp::Vector{Float64}
end

"Layered depth profile of vertical hydraulic conductivity"
struct KvLayered{N}
    # Vertical hydraulic conductivity [m s⁻¹] per soil layer
    kv::Vector{SVector{N, Float64}}
end

"Layered exponential depth profile of vertical hydraulic conductivity"
struct KvLayeredExponential{N}
    # A scaling parameter [m⁻¹] (controls exponential decline of kv_0)
    f::Vector{Float64}
    # Vertical hydraulic conductivity [m s⁻¹] per soil layer
    kv::Vector{SVector{N, Float64}}
    # Number of soil layers [-] with vertical hydraulic conductivity value `kv`
    nlayers_kv::Vector{Int}
    # Depth [m] from soil surface for which layered profile is valid
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
            "soil_exponential_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth",
            LandHydrologySBM;
            sel = indices,
        )
        exp_profile = KvExponential(kv_0, f)
        kv_profile = KvExponentialConstant(exp_profile, z_exp)
    elseif kv_profile_type == VerticalConductivityProfile.layered ||
           kv_profile_type == VerticalConductivityProfile.layered_exponential
        kv = ncread(
            dataset,
            config,
            "soil_layer_water__vertical_saturated_hydraulic_conductivity",
            LandHydrologySBM;
            sel = indices,
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
                "soil_layered_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth",
                LandHydrologySBM;
                sel = indices,
            )
            nlayers_kv = fill(0, n)
            for i in eachindex(nlayers_kv)
                layers = @view sumlayers[i][2:nlayers[i]]
                k = argmin(abs.(z_layered[i] .- layers))
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
    dt::Second;
    data_lookup::DataLookup = DataLookup(),
)
    config_soil_layer_thickness =
        to_SI.(Float64.(config.model.soil_layer__thickness), Ref(MM))

    soil_layer_thickness = SVector(Tuple(push!(config_soil_layer_thickness, MISSING_VALUE)))
    cum_depth_layers = pushfirst(cumsum(soil_layer_thickness), 0.0)
    maxlayers = length(soil_layer_thickness) # max number of soil layers

    @info "Using `$(maxlayers - 1)` soil layers with the following thickness: `$config_soil_layer_thickness`"

    w_soil = ncread(
        dataset,
        config,
        "soil_surface_temperature__weight_coefficient",
        LandHydrologySBM;
        sel = indices,
    )
    cf_soil = ncread(
        dataset,
        config,
        "soil_surface_water__infiltration_reduction_parameter",
        LandHydrologySBM;
        sel = indices,
    )

    # soil parameters
    theta_s = ncread(
        dataset,
        config,
        "soil_water__saturated_volume_fraction",
        LandHydrologySBM;
        sel = indices,
    )
    theta_r = ncread(
        dataset,
        config,
        "soil_water__residual_volume_fraction",
        LandHydrologySBM;
        sel = indices,
    )
    kv_0 = ncread(
        dataset,
        config,
        "soil_surface_water__vertical_saturated_hydraulic_conductivity",
        LandHydrologySBM;
        sel = indices,
    )
    f = ncread(
        dataset,
        config,
        "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
        LandHydrologySBM;
        sel = indices,
    )
    hb = ncread(
        dataset,
        config,
        "soil_water__air_entry_pressure_head",
        LandHydrologySBM;
        sel = indices,
    )
    h1 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h1",
        LandHydrologySBM;
        sel = indices,
    )
    h2 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h2",
        LandHydrologySBM;
        sel = indices,
    )
    h3_high = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h3_high",
        LandHydrologySBM;
        sel = indices,
    )
    h3_low = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h3_low",
        LandHydrologySBM;
        sel = indices,
    )
    h4 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h4",
        LandHydrologySBM;
        sel = indices,
    )
    alpha_h1 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h1_reduction_coefficient",
        LandHydrologySBM;
        sel = indices,
    )
    soilthickness =
        ncread(dataset, config, "soil__thickness", LandHydrologySBM; sel = indices)
    infiltcappath = ncread(
        dataset,
        config,
        "compacted_soil_surface_water__infiltration_capacity",
        LandHydrologySBM;
        sel = indices,
    )
    maxleakage = ncread(
        dataset,
        config,
        "soil_water_saturated_zone_bottom__max_leakage_volume_flux",
        LandHydrologySBM;
        sel = indices,
    )
    c = ncread(
        dataset,
        config,
        "soil_layer_water__brooks_corey_exponent",
        LandHydrologySBM;
        sel = indices,
    )
    if size(c, 1) != maxlayers
        parname = param(config.input.static, "soil_layer_water__brooks_corey_exponent")
        size1 = size(c, 1)
        error("$parname needs a layer dimension of size $maxlayers, but is $size1")
    end
    c = svectorscopy(c, Val{maxlayers}())
    kvfrac = ncread(
        dataset,
        config,
        "soil_layer_water__vertical_saturated_hydraulic_conductivity_factor",
        LandHydrologySBM;
        sel = indices,
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
        sel = indices,
    )

    # vegetation parameters
    rootdistpar = ncread(
        dataset,
        config,
        "soil_wet_root__sigmoid_function_shape_parameter",
        LandHydrologySBM;
        sel = indices,
    )
    cap_hmax = ncread(
        dataset,
        config,
        "soil_water_saturated_zone_top__capillary_rise_max_water_table_depth",
        LandHydrologySBM;
        sel = indices,
    )
    cap_n = ncread(
        dataset,
        config,
        "soil_water_saturated_zone_top__capillary_rise_averianov_exponent",
        LandHydrologySBM;
        sel = indices,
    )

    act_thickl =
        set_layerthickness.(soilthickness, (cum_depth_layers,), (soil_layer_thickness,))
    sumlayers = @. pushfirst(cumsum(act_thickl), 0.0)
    nlayers = number_of_active_layers.(act_thickl)

    if haskey(config.input.static, "soil_water__field_capacity_volume_fraction")
        theta_fc = ncread(
            dataset,
            config,
            "soil_water__field_capacity_volume_fraction",
            LandHydrologySBM;
            sel = indices,
        )
    else
        theta_fc = field_capacity.(act_thickl, nlayers, theta_s, theta_r, c, hb)
    end

    # optional root fraction
    rootfraction_name = "soil_root__length_density_fraction"
    if haskey(config.input.static, rootfraction_name)
        rootfraction =
            ncread(dataset, config, rootfraction_name, LandHydrologySBM; sel = indices)
    else
        n = length(indices)
        (; rootingdepth) = vegetation_parameter_set
        # default root fraction
        rootfraction = zeros(maxlayers, n)
        for i in 1:n
            if rootingdepth[i] > 0.0
                for k in 1:maxlayers
                    if (rootingdepth[i] - sumlayers[i][k]) >= act_thickl[i][k]
                        rootfraction[k, i] = act_thickl[i][k] / rootingdepth[i]
                    else
                        rootfraction[k, i] =
                            max((rootingdepth[i] - sumlayers[i][k]) / rootingdepth[i], 0.0)
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
    sbm_params = SbmSoilParameters(
        data_lookup;
        maxlayers,
        nlayers,
        soilwatercapacity,
        theta_s,
        theta_r,
        theta_fc,
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
        c,
        w_soil,
        cf_soil,
        soil_fraction = fill(MISSING_VALUE, n),
        kv_profile,
        vegetation_parameter_set,
    )
    return sbm_params
end

"SBM soil model"
@kwdef struct SbmSoilModel{N, M, Kv} <: AbstractSoilModel
    n::Int
    boundary_conditions::SbmSoilBC = SbmSoilBC(; n)
    parameters::SbmSoilParameters{N, M, Kv}
    variables::SbmSoilVariables{N}
end

"Initialize SBM soil model"
function SbmSoilModel(
    dataset::NCDataset,
    config::Config,
    vegetation_parameter_set::VegetationParameters,
    indices::Vector{CartesianIndex{2}},
    dt::Second;
    data_lookup::DataLookup = DataLookup(),
)
    n = length(indices)
    parameters = SbmSoilParameters(
        dataset,
        config,
        vegetation_parameter_set,
        indices,
        dt;
        data_lookup,
    )
    variables = SbmSoilVariables(n, parameters; data_lookup)
    soil_model = SbmSoilModel(; n, parameters, variables)
    return soil_model
end

"Return soil fraction"
function soil_fraction!(
    soil_model::AbstractSoilModel,
    glacier_model::AbstractGlacierModel,
    parameters::LandParameters,
)
    (; canopygapfraction) = soil_model.parameters.vegetation_parameter_set
    (; soil_fraction) = soil_model.parameters
    (; water_fraction, river_fraction) = parameters
    glacier_fraction = get_glacier_fraction(glacier_model)
    @. soil_fraction =
        max(canopygapfraction - water_fraction - river_fraction - glacier_fraction, 0.0)
    return nothing
end

"Update boundary conditions of the SBM soil model for a single timestep"
function update_bc_soil_model!(
    soil_model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    external_models::NamedTuple,
    dt::Float64,
)
    (; interception, runoff, demand, allocation) = external_models
    (; potential_transpiration, water_flux_surface, potential_soilevaporation) =
        soil_model.boundary_conditions

    potential_transpiration .= get_potential_transpiration(interception)

    @. potential_soilevaporation =
        soil_model.parameters.soil_fraction * atmospheric_forcing.potential_evaporation

    evaporation!(demand.paddy, potential_soilevaporation, dt)
    potential_soilevaporation .-= get_evaporation(demand.paddy)
    water_flux_surface .=
        max.(
            runoff.boundary_conditions.water_flux_surface .+
            get_irrigation_allocated(allocation) .- runoff.variables.runoff_river .-
            runoff.variables.runoff_land .+ get_water_depth(demand.paddy) / dt,
            0.0,
        )
    return nothing
end

"Update soil temperature of the SBM soil model for a single timestep"
function soil_temperature!(
    soil_model::SbmSoilModel,
    ::AbstractSnowModel,
    temperature::Vector{Float64},
)
    v = soil_model.variables
    p = soil_model.parameters
    @. v.tsoil = soil_temperature(v.tsoil, p.w_soil, temperature)
    return nothing
end

soil_temperature!(::SbmSoilModel, ::NoSnowModel, ::Vector{Float64}) = nothing

"Update total available water in the unsaturated zone of the SBM soil model for a single timestep"
function ustoredepth!(soil_model::SbmSoilModel)
    v = soil_model.variables
    p = soil_model.parameters
    for i in eachindex(v.ustorelayerdepth)
        v.ustoredepth[i] = sum(@view v.ustorelayerdepth[i][1:p.nlayers[i]])
    end
    return nothing
end

"Update the infiltration reduction factor of the SBM soil model for a single timestep"
function infiltration_reduction_factor!(
    soil_model::SbmSoilModel;
    modelsnow = false,
    soil_infiltration_reduction = false,
)
    v = soil_model.variables
    p = soil_model.parameters

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
    infiltration!(soil_model::SbmSoilModel, dt::Float64)

Update the infiltration rate `infiltsoilpath` and infiltration excess water rate
`infiltexcess` of the SBM soil model for a single timestep.
"""
function infiltration!(soil_model::SbmSoilModel, dt::Float64)
    v = soil_model.variables
    p = soil_model.parameters
    (; water_flux_surface) = soil_model.boundary_conditions

    n = length(v.infiltsoilpath)
    threaded_foreach(1:n; basesize = 1000) do i
        v.infiltsoilpath[i], v.infiltexcess[i] = infiltration(
            water_flux_surface[i],
            p.pathfrac[i],
            p.infiltcapsoil[i],
            p.infiltcappath[i],
            v.ustorecapacity[i],
            v.f_infiltration_reduction[i],
            dt,
        )
    end
    return nothing
end

"""
    unsaturated_zone_flow!(soil_model::SbmSoilModel, dt::Float64)

Update unsaturated storage `ustorelayerdepth` and the `transfer` of water from the unsaturated
to the saturated store of the SBM soil model for a single timestep, based on the Brooks-Corey
approach.
"""
function unsaturated_zone_flow!(soil_model::SbmSoilModel, dt::Float64)
    v = soil_model.variables
    p = soil_model.parameters

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
                    v.ustorelayerdepth[i][m] + v.infiltsoilpath[i] * dt
                else
                    v.ustorelayerdepth[i][m] + flow_rate * dt
                end
                ustorelayerdepth, flow_rate =
                    unsatzone_flow_layer(ustorelayerdepth, kv_z, l_sat, p.c[i][m], dt)
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
    soil_evaporation!(soil_model::SbmSoilModel)

Update soil evaporation from the saturated store `soilevapsat` and the total soil
evaporation from the unsaturated and saturated store `soilevap` of the SBM soil model for a
single timestep. Also unsaturated storage `ustorelayerdepth` and the saturated store
`satwaterdepth` are updated.
"""
function soil_evaporation!(soil_model::SbmSoilModel, dt::Float64)
    (; potential_soilevaporation) = soil_model.boundary_conditions
    v = soil_model.variables
    p = soil_model.parameters

    n = length(potential_soilevaporation)
    threaded_foreach(1:n; basesize = 1000) do i
        potsoilevap = potential_soilevaporation[i]
        # First calculate the evaporation of unsaturated storage into the
        # atmosphere from the upper layer.
        soilevapunsat = soil_evaporation_unsaturated_store(
            potsoilevap,
            v.ustorelayerdepth[i][1],
            v.ustorelayerthickness[i][1],
            v.n_unsatlayers[i],
            v.zi[i],
            p.theta_s[i] - p.theta_r[i],
        )
        # Ensure that the unsaturated evaporation rate does not exceed the
        # available unsaturated moisture
        soilevapunsat = min(soilevapunsat, v.ustorelayerdepth[i][1] / dt)
        # Update the additional atmospheric demand
        potsoilevap -= soilevapunsat
        v.ustorelayerdepth[i] = setindex(
            v.ustorelayerdepth[i],
            v.ustorelayerdepth[i][1] - soilevapunsat * dt,
            1,
        )
        theta_drainable = lower_bound_drainable_porosity(p.theta_s[i], p.theta_fc[i])
        soilevapsat = soil_evaporation_saturated_store(
            potsoilevap,
            v.n_unsatlayers[i],
            p.act_thickl[i][1],
            v.zi[i],
            theta_drainable,
            dt,
        )
        v.soilevapsat[i] = soilevapsat
        v.soilevap[i] = soilevapunsat + soilevapsat
        v.drainable_waterdepth[i] -= soilevapsat * dt
    end
    return nothing
end

"""
    transpiration!(soil_model::SbmSoilModel, dt::Float64)

Update total `transpiration`, transpiration from the unsaturated store `ae_ustore` and
saturated store `actevapsat` of the SBM soil model for a single timestep. Also unsaturated
storage `ustorelayerdepth` and the saturated store `satwaterdepth` are updated.
"""
function transpiration!(soil_model::SbmSoilModel, dt::Float64)
    (; potential_transpiration) = soil_model.boundary_conditions
    v = soil_model.variables
    p = soil_model.parameters

    rootingdepth = get_rootingdepth(soil_model)
    n = length(rootingdepth)

    threaded_foreach(1:n; basesize = 250) do i
        v.h3[i] = feddes_h3(p.h3_high[i], p.h3_low[i], potential_transpiration[i])

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
            else
                rootfraction_unsat = p.rootfraction[i][k]
            end
            sum_rootfraction_unsat += rootfraction_unsat

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
            vwc = max(v.ustorelayerdepth[i][k] / v.ustorelayerthickness[i][k], 1e-7)
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
            maxextr = v.ustorelayerdepth[i][k] * availcap / dt
            actevapustore_layer =
                min(alpha * rootfraction_unsat_scaled * potential_transpiration[i], maxextr)
            ustorelayerdepth = v.ustorelayerdepth[i][k] - actevapustore_layer * dt
            actevapustore += actevapustore_layer
            v.ustorelayerdepth[i] = setindex(v.ustorelayerdepth[i], ustorelayerdepth, k)
        end

        # transpiration from saturated store
        wetroots = scurve(v.zi[i], rootingdepth[i], 1.0, p.rootdistpar[i])
        alpha = rwu_reduction_feddes(
            Float64(0.0),
            p.h1[i],
            p.h2[i],
            v.h3[i],
            p.h4[i],
            p.alpha_h1[i],
        )
        restpottrans = potential_transpiration[i] - actevapustore
        actevapsat = min(restpottrans * wetroots * alpha, v.drainable_waterdepth[i] / dt)

        v.ae_ustore[i] = actevapustore
        v.actevapsat[i] = actevapsat
        v.drainable_waterdepth[i] -= actevapsat * dt
        v.transpiration[i] = actevapustore + actevapsat
    end
    return nothing
end

"""
    actual_infiltration!(soil_model::SbmSoilModel)

Update the actual infiltration rate `actinfilt` of the SBM soil model for a single timestep.

A soil water balance check is performed. Unsaturated storage that exceeds the maximum
storage per unsaturated soil layer is transferred to the layer above (or surface), from the
bottom to the top unsaturated soil layer. The resulting excess water `ustoredepth_excess` is
subtracted from the infiltration rate `infiltsoilpath`.
"""
function actual_infiltration!(soil_model::SbmSoilModel, dt::Float64)
    v = soil_model.variables
    p = soil_model.parameters

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
        v.actinfilt[i] = v.infiltsoilpath[i] - ustoredepth_excess / dt
    end
    return nothing
end

"""
    actual_infiltration_soil_path!(soil_model::SbmSoilModel)

Update the actual infiltration rate for soil `actinfiltsoil` and paved area `actinfiltpath`
of the SBM soil model for a single timestep.
"""
function actual_infiltration_soil_path!(soil_model::SbmSoilModel)
    v = soil_model.variables
    p = soil_model.parameters
    (; water_flux_surface) = soil_model.boundary_conditions

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
    capillary_flux!(soil_model::SbmSoilModel)

Update the capillary flux `actcapflux` of the SBM soil model for a single timestep.
"""
function capillary_flux!(soil_model::SbmSoilModel, dt::Float64)
    v = soil_model.variables
    p = soil_model.parameters
    rootingdepth = get_rootingdepth(soil_model)

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
            maxcapflux = max(
                0.0,
                min(
                    ksat,
                    v.ae_ustore[i],
                    v.ustorecapacity[i] / dt,
                    v.drainable_waterdepth[i] / dt,
                ),
            )

            capflux = if v.zi[i] > rootingdepth[i]
                maxcapflux *
                pow(1.0 - min(v.zi[i], p.cap_hmax[i]) / (p.cap_hmax[i]), p.cap_n[i])
            else
                0.0
            end
            netcapflux = capflux
            actcapflux = 0.0
            for k in v.n_unsatlayers[i]:-1:1
                toadd = min(
                    netcapflux,
                    max(
                        (
                            v.ustorelayerthickness[i][k] * (p.theta_s[i] - p.theta_r[i]) - v.ustorelayerdepth[i][k]
                        ) / dt,
                        0.0,
                    ),
                )
                v.ustorelayerdepth[i] = setindex(
                    v.ustorelayerdepth[i],
                    v.ustorelayerdepth[i][k] + toadd * dt,
                    k,
                )
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
    leakage!(soil_model::SbmSoilModel, dt::Float64)

Update the actual leakage rate `actleakage` of the SBM soil model for a single timestep.
"""
function leakage!(soil_model::SbmSoilModel, dt::Float64)
    v = soil_model.variables
    p = soil_model.parameters

    n = length(v.actleakage)
    threaded_foreach(1:n; basesize = 1000) do i
        deepksat = hydraulic_conductivity_at_depth(
            p.kv_profile,
            p.kvfrac,
            p.soilthickness[i],
            i,
            p.nlayers[i],
        )

        deeptransfer = min(v.drainable_waterdepth[i] / dt, deepksat)
        v.actleakage[i] = max(0.0, min(p.maxleakage[i], deeptransfer))
    end
    return nothing
end

"""
    update_soil_water_flow!(
        soil_model::SbmSoilModel,
        atmospheric_forcing::AtmosphericForcing,
        external_models::NamedTuple,
        config::Config,
        dt::Float64,
    )

Update the SBM soil model (infiltration, unsaturated zone flow, soil evaporation and
transpiration, capillary flux and leakage) for a single timestep.
"""
function update_soil_water_flow!(
    soil_model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    external_models::NamedTuple,
    config::Config,
    dt::Float64,
)
    (; snow, runoff, demand) = external_models
    (; temperature) = atmospheric_forcing
    (; water_flux_surface) = soil_model.boundary_conditions
    v = soil_model.variables
    p = soil_model.parameters

    # mainly required for external state changes (e.g. through BMI)
    update_diagnostic_vars!(soil_model)
    # infiltration
    soil_temperature!(soil_model, snow, temperature)
    infiltration_reduction_factor!(
        soil_model;
        modelsnow = config.model.snow__flag,
        soil_infiltration_reduction = config.model.soil_infiltration_reduction__flag,
    )
    infiltration!(soil_model, dt)
    # unsaturated zone flow
    unsaturated_zone_flow!(soil_model, dt)
    # soil evaporation and transpiration
    soil_evaporation!(soil_model, dt)
    transpiration!(soil_model, dt)
    # actual infiltration and excess water
    actual_infiltration!(soil_model, dt)
    @. v.excesswater = (water_flux_surface - v.actinfilt) - v.infiltexcess
    actual_infiltration_soil_path!(soil_model)
    @. v.excesswatersoil =
        max(water_flux_surface * (1.0 - p.pathfrac) - v.actinfiltsoil, 0.0)
    @. v.excesswaterpath = max(water_flux_surface * p.pathfrac - v.actinfiltpath, 0.0)
    # recompute the unsaturated store and ustorecapacity (for capillary flux)
    ustoredepth!(soil_model)
    @. v.ustorecapacity = p.soilwatercapacity - v.satwaterdepth - v.ustoredepth
    # capillary flux and leakage
    capillary_flux!(soil_model, dt)
    leakage!(soil_model, dt)
    # recharge rate to the saturated store
    @. v.recharge =
        (v.transfer - v.actcapflux - v.actleakage - v.actevapsat - v.soilevapsat)
    # total actual evapotranspiration
    v.actevap .=
        v.soilevap .+ v.transpiration .+ runoff.variables.ae_openw_r .+
        runoff.variables.ae_openw_l .+ get_evaporation(demand.paddy)
    return nothing
end

function update_ustorelayerdepth!(soil, zi_prev, zi, i)
    v = soil.variables
    p = soil.parameters

    n_unsatlayers_prev = v.n_unsatlayers[i]
    ustorelayerthickness_prev = v.ustorelayerthickness[i]
    ustorelayerthickness = set_layerthickness(zi, p.sumlayers[i], p.act_thickl[i])
    n_unsatlayers = number_of_active_layers(ustorelayerthickness)
    ustorelayerdepth = v.ustorelayerdepth[i]
    if zi < zi_prev
        for k in n_unsatlayers:n_unsatlayers_prev
            k == 0 && continue
            if isnan(ustorelayerthickness[k])
                ustorelayerdepth = setindex(ustorelayerdepth, 0.0, k)
            else
                f = ustorelayerthickness[k] / ustorelayerthickness_prev[k]
                ustorelayerdepth = setindex(ustorelayerdepth, f * ustorelayerdepth[k], k)
            end
        end
    else
        for k in n_unsatlayers_prev:n_unsatlayers
            k == 0 && continue
            thickness_prev =
                isnan(ustorelayerthickness_prev[k]) ? 0.0 : ustorelayerthickness_prev[k]
            delta_thickness = ustorelayerthickness[k] - thickness_prev
            ustorelayerdepth = setindex(
                ustorelayerdepth,
                ustorelayerdepth[k] + delta_thickness * (p.theta_fc[i] - p.theta_r[i]),
                k,
            )
        end
    end
    v.n_unsatlayers[i] = n_unsatlayers
    v.ustorelayerdepth[i] = ustorelayerdepth
    v.ustorelayerthickness[i] = ustorelayerthickness
    v.zi[i] = zi
end

"""
    update_ustorelayerdepth!(soil_model::SbmSoilModel, subsurface_flow)

Update the `SbmSoilModel` variables unsaturated store depth of soil layers
`ustorelayerdepth`, number of unsaturated zone soil layers `n_unsatlayers`, thickness of
unsaturated zone soil layers `ustorelayerthickness` and water table depth `zi`, based on the
water table change computed by a subsurface flow model.
"""
function update_ustorelayerdepth!(soil_model::SbmSoilModel, subsurface_flow)
    p = soil_model.parameters
    v = soil_model.variables

    zi = get_water_depth(subsurface_flow)

    n = length(v.zi)
    threaded_foreach(1:n; basesize = 1000) do i
        zi_prev = v.zi[i]
        update_ustorelayerdepth!(soil_model, zi_prev, zi[i], i)
    end
end

"""
    update_soil_water_storage!(soil_model::SbmSoilModel, external_models::NamedTuple)

Update the SBM soil model for a single timestep based on the update of a subsurface flow
model, resulting in a change in water table depth and an exfiltration rate `exfiltwater_average`.

The available water in unsaturated zone `ustoredepth`, unsaturated store capacity
`ustorecapacity`, `total_soilwater_storage`, land `runoff` and `net_runoff`, the saturated
store `satwaterdepth` and the water exfiltrating during saturation excess conditions
`exfiltsatwater` are updated. Additionally, volumetric water content per soil layer and for
the root zone are updated.
"""
function update_soil_water_storage!(
    soil_model::SbmSoilModel,
    external_models::NamedTuple,
    dt::Float64,
)
    (; runoff, demand, subsurface_flow) = external_models
    (; runoff_land, ae_openw_l) = runoff.variables
    p = soil_model.parameters
    v = soil_model.variables
    exfiltsatwater = subsurface_flow.variables.exfiltwater_average

    rootingdepth = get_rootingdepth(soil_model)

    n = length(v.zi)
    threaded_foreach(1:n; basesize = 1000) do i
        ustorelayerdepth = v.ustorelayerdepth[i]
        ustorelayerthickness = v.ustorelayerthickness[i]
        ustoredepth = sum(@view ustorelayerdepth[1:(v.n_unsatlayers[i])])
        sbm_runoff = max(
            0.0,
            exfiltsatwater[i] + v.excesswater[i] + runoff_land[i] + v.infiltexcess[i],
        )

        # volumetric water content per soil layer and root zone
        vwc = v.vwc[i]
        vwc_perc = v.vwc_perc[i]
        for k in 1:p.nlayers[i]
            if k <= v.n_unsatlayers[i]
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
            vwc_perc = setindex(vwc_perc, from_SI(vwc[k] / p.theta_s[i], PERCENTAGE), k)
        end

        rootstore_unsat = 0
        for k in 1:(v.n_unsatlayers[i])
            rootstore_unsat +=
                min(
                    1.0,
                    (
                        max(0.0, rootingdepth[i] - p.sumlayers[i][k]) /
                        ustorelayerthickness[k]
                    ),
                ) * ustorelayerdepth[k]
        end

        rootstore_sat = max(0.0, rootingdepth[i] - v.zi[i]) * (p.theta_s[i] - p.theta_r[i])
        rootstore = rootstore_sat + rootstore_unsat
        vwc_root = rootstore / rootingdepth[i] + p.theta_r[i]
        vwc_percroot = from_SI(vwc_root / p.theta_s[i], PERCENTAGE)
        satwaterdepth = (p.soilthickness[i] - v.zi[i]) * (p.theta_s[i] - p.theta_r[i])
        drainable_waterdepth =
            (p.soilthickness[i] - v.zi[i]) *
            lower_bound_drainable_porosity(p.theta_s[i], p.theta_fc[i])
        ustorecapacity = p.soilwatercapacity[i] - satwaterdepth - ustoredepth

        # update the outputs and states
        v.ustorecapacity[i] = ustorecapacity
        v.ustoredepth[i] = ustoredepth
        v.satwaterdepth[i] = satwaterdepth
        v.drainable_waterdepth[i] = drainable_waterdepth
        v.exfiltsatwater[i] = exfiltsatwater[i]
        v.runoff[i] = sbm_runoff
        v.vwc[i] = vwc
        v.vwc_perc[i] = vwc_perc
        v.rootstore[i] = rootstore
        v.vwc_root[i] = vwc_root
        v.vwc_percroot[i] = vwc_percroot
        v.total_soilwater_storage[i] = satwaterdepth + ustoredepth
    end
    # update runoff and net_runoff (the runoff rate depends on the presence of paddy fields
    # and the h_max parameter of a paddy field)
    update_runoff!(demand.paddy, v.runoff, dt)
    @. v.net_runoff = v.runoff - ae_openw_l
    return nothing
end

"""
    update_diagnostic_vars!(soil_model::SbmSoilModel)

Update diagnostic variables of `SbmSoilModel` that are critical for subsequent computations
and depend on state variables `satwaterdepth` and `ustorelayerdepth`.
"""
function update_diagnostic_vars!(soil_model::SbmSoilModel)
    (;
        zi,
        satwaterdepth,
        drainable_waterdepth,
        ustorelayerthickness,
        ustorecapacity,
        ustoredepth,
        total_soilwater_storage,
        n_unsatlayers,
    ) = soil_model.variables
    (;
        soilthickness,
        theta_s,
        theta_r,
        theta_fc,
        soilwatercapacity,
        sumlayers,
        act_thickl,
    ) = soil_model.parameters

    ustoredepth!(soil_model)
    @. zi = max(0.0, soilthickness - satwaterdepth / (theta_s - theta_r))
    @. drainable_waterdepth =
        (soilthickness - zi) * lower_bound_drainable_porosity(theta_s, theta_fc)
    @. ustorecapacity = soilwatercapacity - satwaterdepth - ustoredepth
    @. ustorelayerthickness = set_layerthickness(zi, sumlayers, act_thickl)
    @. n_unsatlayers = number_of_active_layers(ustorelayerthickness)
    @. total_soilwater_storage = satwaterdepth + ustoredepth
end

# wrapper method
get_rootingdepth(soil_model::SbmSoilModel) =
    soil_model.parameters.vegetation_parameter_set.rootingdepth
