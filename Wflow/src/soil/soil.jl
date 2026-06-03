abstract type AbstractSoilModel end

"Struct for storing SBM soil model variables"
@with_kw struct SbmSoilVariables{N}
    n::Int
    # Calculated soil water pressure head h3 of the root water uptake reduction function (Feddes) [m]
    h3::Vector{Float64} = fill(MISSING_VALUE, n)
    # Unsaturated store capacity [m]
    unsaturated_store_capacity::Vector{Float64}
    # Amount of water in the unsaturated store, per layer [m]
    unsaturated_layer_depth::Vector{SVector{N, Float64}}
    # Thickness of unsaturated zone, per layer [m]
    unsaturated_layer_thickness::Vector{SVector{N, Float64}}
    # Saturated store [m]
    saturated_water_depth::Vector{Float64}
    # Drainable water store [m]
    drainable_water_depth::Vector{Float64}
    # Pseudo-water table depth [m] (top of the saturated zone)
    water_table_depth::Vector{Float64}
    # Number of unsaturated soil layers
    n_unsatlayers::Vector{Int}
    # Transpiration [m s⁻¹]
    transpiration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual evaporation from unsaturated store [m s⁻¹]
    actual_evaporation_unsaturated_store::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil evaporation from unsaturated and saturated store [m s⁻¹]
    soil_evaporation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil evaporation from saturated store [m s⁻¹]
    soil_evaporation_saturated_zone::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual capillary rise [m s⁻¹]
    actual_capillary_flux::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual transpiration from saturated store [m s⁻¹]
    actual_evaporation_saturated_zone::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total actual evapotranspiration [m s⁻¹]
    actual_evapotranspiration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration into the unsaturated zone [m s⁻¹]
    actual_infiltration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration non-compacted fraction [m s⁻¹]
    actual_infiltration_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration compacted fraction [m s⁻¹]
    actual_infiltration_compacted_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration (compacted and the non-compacted areas) [m s⁻¹]
    infiltration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Infiltration excess water [m s⁻¹]
    infiltration_excess::Vector{Float64} = fill(MISSING_VALUE, n)
    # Water that cannot infiltrate due to saturated soil (saturation excess) [m s⁻¹]
    saturation_excess_water::Vector{Float64} = fill(MISSING_VALUE, n)
    # Water exfiltrating during saturation excess conditions [m s⁻¹]
    exfiltration_saturated_water::Vector{Float64} = fill(MISSING_VALUE, n)
    # Excess water for non-compacted fraction [m s⁻¹]
    excess_water_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Excess water for compacted fraction [m s⁻¹]
    excess_water_compacted_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total surface runoff from infiltration and saturation excess (excluding actual open water evaporation) [m s⁻¹]
    runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Net surface runoff (surface runoff - actual open water evaporation) [m s⁻¹]
    net_runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Volumetric water content [-] per soil layer (including theta_r and saturated zone)
    volumetric_water_content::Vector{SVector{N, Float64}}
    # Volumetric water content [%] per soil layer (including theta_r and saturated zone)
    relative_volumetric_water_content::Vector{SVector{N, Float64}}
    # Root water storage [m] in unsaturated and saturated zone (excluding theta_r)
    root_zone_storage::Vector{Float64} = fill(MISSING_VALUE, n)
    # Volumetric water content [-] in root zone (including theta_r and saturated zone)
    volumetric_water_content_root_zone::Vector{Float64} = fill(MISSING_VALUE, n)
    # Volumetric water content [%] in root zone (including theta_r and saturated zone)
    relative_volumetric_water_content_root_zone::Vector{Float64} = fill(MISSING_VALUE, n)
    # Amount of available water in the unsaturated zone [m]
    unsaturated_store_depth::Vector{Float64} = zeros(n)
    # Downward flux from unsaturated to saturated zone [m s⁻¹]
    transfer::Vector{Float64} = fill(MISSING_VALUE, n)
    # Net recharge to saturated store [m s⁻¹]
    recharge::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual leakage from saturated store [m s⁻¹]
    actual_leakage::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total water storage (excluding floodplain volume and reservoirs) [m]
    total_storage::Vector{Float64} = zeros(Float64, n)
    # Total soil water storage [m]
    total_soil_water_storage::Vector{Float64}
    # Top soil temperature [K]
    soil_surface_temperature::Vector{Float64} = to_SI.(fill(10.0, n), Ref(ABSOLUTE_DEGREES))
    # Soil infiltration reduction factor (when soil is frozen) [-]
    f_infiltration_reduction::Vector{Float64} = ones(n)
end

"Struct for storing SBM soil model parameters"
@with_kw struct SbmSoilParameters{N, M, Kv}
    # Maximum number of soil layers [-]
    maximum_number_of_layers::Int
    # Number of soil layers [-]
    number_of_layers::Vector{Int}
    # Saturated water content (porosity) [-]
    theta_s::Vector{Float64}
    # Residual water content [-]
    theta_r::Vector{Float64}
    # Field capacity water content [-]
    theta_fc::Vector{Float64}
    # Soilwater capacity [m]
    soil_water_capacity::Vector{Float64}
    # Multiplication factor [-] applied to kv_z (vertical flow)
    vertical_hydraulic_conductivity_factor::Vector{SVector{N, Float64}}
    # Air entry pressure [m] of soil (Brooks-Corey)
    air_entry_pressure::Vector{Float64}
    # Soil thickness [m]
    soil_thickness::Vector{Float64}
    # Thickness of soil layers [m]
    actual_layer_thickness::Vector{SVector{N, Float64}}
    # Cumulative sum of soil layers [m], starting at soil surface (0)
    cumulative_layer_depth::Vector{SVector{M, Float64}}
    # Infiltration capacity of the compacted areas [m s⁻¹]
    infiltration_capacity_compacted_soil::Vector{Float64}
    # Soil infiltration capacity [m s⁻¹]
    infiltration_capacity_soil::Vector{Float64}
    # Maximum leakage [m s⁻¹] from saturated zone
    maximum_leakage::Vector{Float64}
    # Parameter [m] controlling capillary rise
    cap_hmax::Vector{Float64}
    # Coefficient [-] controlling capillary rise
    cap_n::Vector{Float64}
    # Brooks-Corey power coefficient [-] for each soil layer
    brooks_corey_exponent::Vector{SVector{N, Float64}}
    # Soil temperature smooth factor [-]
    w_soil::Vector{Float64}
    # Controls soil infiltration reduction factor when soil is frozen [-]
    cf_soil::Vector{Float64}
    # Fraction of compacted area  [-]
    compacted_soil_area_fraction::Vector{Float64}
    # Controls how roots are linked to water table [m⁻¹]
    wet_root_distribution_parameter::Vector{Float64}
    # Fraction of the root length density in each soil layer [-]
    rootfraction::Vector{SVector{N, Float64}}
    # Soil water pressure head h1 of the root water uptake reduction function (Feddes) [m]
    h1::Vector{Float64}
    # Soil water pressure head h2 of the root water uptake reduction function (Feddes) [m]
    h2::Vector{Float64}
    # Soil water pressure head h3_high of the root water uptake reduction function (Feddes) [m]
    h3_high::Vector{Float64}
    # Soil water pressure head h3_low of the root water uptake reduction function (Feddes) [m]
    h3_low::Vector{Float64}
    # Soil water pressure head h4 of the root water uptake reduction function (Feddes) [m]
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
        soil_thickness,
        maximum_number_of_layers,
        actual_layer_thickness,
        cumulative_layer_depth,
        soil_water_capacity,
        theta_s,
        theta_r,
        theta_fc,
    ) = parameters
    saturated_water_depth = 0.85 .* soil_water_capacity # cold state value for saturated_water_depth
    unsaturated_store_depth = zeros(n)
    water_table_depth = @. max(0.0, soil_thickness - saturated_water_depth / (theta_s - theta_r))
    drainable_water_depth =
        @. (soil_thickness - water_table_depth) * lower_bound_drainable_porosity(theta_s, theta_fc)
    unsaturated_layer_thickness = set_layerthickness.(water_table_depth, cumulative_layer_depth, actual_layer_thickness)
    n_unsatlayers = number_of_active_layers.(unsaturated_layer_thickness)
    volumetric_water_content = fill(MISSING_VALUE, maximum_number_of_layers, n)
    relative_volumetric_water_content = fill(MISSING_VALUE, maximum_number_of_layers, n)
    total_soil_water_storage = saturated_water_depth .+ unsaturated_store_depth

    vars = SbmSoilVariables(;
        n,
        unsaturated_layer_depth = zero(actual_layer_thickness),
        unsaturated_store_capacity = soil_water_capacity .- saturated_water_depth,
        unsaturated_layer_thickness,
        saturated_water_depth,
        drainable_water_depth,
        water_table_depth,
        n_unsatlayers,
        volumetric_water_content = svectorscopy(volumetric_water_content, Val{maximum_number_of_layers}()),
        relative_volumetric_water_content = svectorscopy(relative_volumetric_water_content, Val{maximum_number_of_layers}()),
        total_soil_water_storage,
    )
    return vars
end

"Struct for storing SBM soil model boundary conditions"
@with_kw struct SbmSoilBC
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
    hydraulic_conductivity_scale_parameter::Vector{Float64}
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
    hydraulic_conductivity_scale_parameter::Vector{Float64}
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
    hydraulic_conductivity_scale_parameter::Vector{Float64},
    maximum_number_of_layers::Int,
    number_of_layers::Vector{Int},
    cumulative_layer_depth::Vector,
    dt::Second,
)
    kv_profile_type = config.model.saturated_hydraulic_conductivity_profile
    n = length(indices)
    if kv_profile_type == VerticalConductivityProfile.exponential
        kv_profile = KvExponential(kv_0, hydraulic_conductivity_scale_parameter)
    elseif kv_profile_type == VerticalConductivityProfile.exponential_constant
        z_exp = ncread(
            dataset,
            config,
            "soil_exponential_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth",
            LandHydrologySBM;
            sel = indices,
        )
        exp_profile = KvExponential(kv_0, hydraulic_conductivity_scale_parameter)
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
        if size(kv, 1) != maximum_number_of_layers
            parname = param(
                config.input.static,
                "soil_layer_water__vertical_saturated_hydraulic_conductivity",
            )
            size1 = size(kv, 1)
            error("$parname needs a layer dimension of size $maximum_number_of_layers, but is $size1")
        end
        if kv_profile_type == VerticalConductivityProfile.layered
            kv_profile = KvLayered(svectorscopy(kv, Val{maximum_number_of_layers}()))
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
                layers = @view cumulative_layer_depth[i][2:number_of_layers[i]]
                k = argmin(abs.(z_layered[i] .- layers))
                nlayers_kv[i] = k
                z_layered[i] = layers[k]
            end
            kv_profile = KvLayeredExponential(
                hydraulic_conductivity_scale_parameter,
                svectorscopy(kv, Val{maximum_number_of_layers}()),
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
    config_soil_layer_thickness =
        to_SI.(Float64.(config.model.soil_layer__thickness), Ref(MM))

    soil_layer_thickness = SVector(Tuple(push!(config_soil_layer_thickness, MISSING_VALUE)))
    cum_depth_layers = pushfirst(cumsum(soil_layer_thickness), 0.0)
    maximum_number_of_layers = length(soil_layer_thickness) # max number of soil layers

    @info "Using `$(maximum_number_of_layers - 1)` soil layers with the following thickness: `$config_soil_layer_thickness`"

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
    hydraulic_conductivity_scale_parameter = ncread(
        dataset,
        config,
        "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
        LandHydrologySBM;
        sel = indices,
    )
    air_entry_pressure = ncread(
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
    soil_thickness =
        ncread(dataset, config, "soil__thickness", LandHydrologySBM; sel = indices)
    infiltration_capacity_compacted_soil = ncread(
        dataset,
        config,
        "compacted_soil_surface_water__infiltration_capacity",
        LandHydrologySBM;
        sel = indices,
    )
    maximum_leakage = ncread(
        dataset,
        config,
        "soil_water_saturated_zone_bottom__max_leakage_volume_flux",
        LandHydrologySBM;
        sel = indices,
    )
    brooks_corey_exponent = ncread(
        dataset,
        config,
        "soil_layer_water__brooks_corey_exponent",
        LandHydrologySBM;
        sel = indices,
    )
    if size(brooks_corey_exponent, 1) != maximum_number_of_layers
        parname = param(config.input.static, "soil_layer_water__brooks_corey_exponent")
        size1 = size(brooks_corey_exponent, 1)
        error("$parname needs a layer dimension of size $maximum_number_of_layers, but is $size1")
    end
    brooks_corey_exponent = svectorscopy(brooks_corey_exponent, Val{maximum_number_of_layers}())
    vertical_hydraulic_conductivity_factor = ncread(
        dataset,
        config,
        "soil_layer_water__vertical_saturated_hydraulic_conductivity_factor",
        LandHydrologySBM;
        sel = indices,
    )
    if size(vertical_hydraulic_conductivity_factor, 1) != maximum_number_of_layers
        parname = param(
            config.input.static,
            "soil_layer_water__vertical_saturated_hydraulic_conductivity_factor",
        )
        size1 = size(vertical_hydraulic_conductivity_factor, 1)
        error("$parname needs a layer dimension of size $maximum_number_of_layers, but is $size1")
    end

    # soil infiltration capacity based on kv_0 and vertical_hydraulic_conductivity_factor upper soil layer
    infiltration_capacity_soil = kv_0 .* @view vertical_hydraulic_conductivity_factor[1, :]
    # fraction compacted area
    compacted_soil_area_fraction = ncread(
        dataset,
        config,
        "compacted_soil__area_fraction",
        LandHydrologySBM;
        sel = indices,
    )

    # vegetation parameters
    wet_root_distribution_parameter = ncread(
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

    actual_layer_thickness =
        set_layerthickness.(soil_thickness, (cum_depth_layers,), (soil_layer_thickness,))
    cumulative_layer_depth = @. pushfirst(cumsum(actual_layer_thickness), 0.0)
    number_of_layers = number_of_active_layers.(actual_layer_thickness)

    if haskey(config.input.static, "soil_water__field_capacity_volume_fraction")
        theta_fc = ncread(
            dataset,
            config,
            "soil_water__field_capacity_volume_fraction",
            LandHydrologySBM;
            sel = indices,
        )
    else
        theta_fc = field_capacity.(actual_layer_thickness, number_of_layers, theta_s, theta_r, brooks_corey_exponent, air_entry_pressure)
    end

    # optional root fraction
    rootfraction_name = "soil_root__length_density_fraction"
    if haskey(config.input.static, rootfraction_name)
        rootfraction =
            ncread(dataset, config, rootfraction_name, LandHydrologySBM; sel = indices)
    else
        n = length(indices)
        (; rooting_depth) = vegetation_parameter_set
        # default root fraction
        rootfraction = zeros(maximum_number_of_layers, n)
        for i in 1:n
            if rooting_depth[i] > 0.0
                for k in 1:maximum_number_of_layers
                    if (rooting_depth[i] - cumulative_layer_depth[i][k]) >= actual_layer_thickness[i][k]
                        rootfraction[k, i] = actual_layer_thickness[i][k] / rooting_depth[i]
                    else
                        rootfraction[k, i] =
                            max((rooting_depth[i] - cumulative_layer_depth[i][k]) / rooting_depth[i], 0.0)
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
        hydraulic_conductivity_scale_parameter,
        maximum_number_of_layers,
        number_of_layers,
        cumulative_layer_depth,
        dt,
    )

    soil_water_capacity = @. soil_thickness * (theta_s - theta_r)

    n = length(indices)
    sbm_params = SbmSoilParameters(;
        maximum_number_of_layers,
        number_of_layers,
        soil_water_capacity,
        theta_s,
        theta_r,
        theta_fc,
        vertical_hydraulic_conductivity_factor = svectorscopy(vertical_hydraulic_conductivity_factor, Val{maximum_number_of_layers}()),
        air_entry_pressure,
        h1,
        h2,
        h3_high,
        h3_low,
        h4,
        alpha_h1,
        soil_thickness,
        actual_layer_thickness,
        cumulative_layer_depth,
        infiltration_capacity_compacted_soil,
        infiltration_capacity_soil,
        maximum_leakage,
        compacted_soil_area_fraction,
        wet_root_distribution_parameter,
        rootfraction = svectorscopy(rootfraction, Val{maximum_number_of_layers}()),
        cap_hmax,
        cap_n,
        brooks_corey_exponent,
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
    dt::Second,
)
    n = length(indices)
    parameters = SbmSoilParameters(dataset, config, vegetation_parameter_set, indices, dt)
    variables = SbmSoilVariables(n, parameters)
    soil_model = SbmSoilModel(; n, parameters, variables)
    return soil_model
end

"Return soil fraction"
function soil_fraction!(
    soil_model::AbstractSoilModel,
    glacier_model::AbstractGlacierModel,
    parameters::LandParameters,
)
    (; canopy_gap_fraction) = soil_model.parameters.vegetation_parameter_set
    (; soil_fraction) = soil_model.parameters
    (; water_fraction, river_fraction) = parameters
    glacier_fraction = get_glacier_fraction(glacier_model)
    @. soil_fraction =
        max(canopy_gap_fraction - water_fraction - river_fraction - glacier_fraction, 0.0)
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
    @. v.soil_surface_temperature = soil_temperature(v.soil_surface_temperature, p.w_soil, temperature)
    return nothing
end

soil_temperature!(::SbmSoilModel, ::NoSnowModel, ::Vector{Float64}) = nothing

"Update total available water in the unsaturated zone of the SBM soil model for a single timestep"
function unsaturated_store_depth!(soil_model::SbmSoilModel)
    v = soil_model.variables
    p = soil_model.parameters
    for i in eachindex(v.unsaturated_layer_depth)
        v.unsaturated_store_depth[i] = sum(@view v.unsaturated_layer_depth[i][1:p.number_of_layers[i]])
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

    n = length(v.soil_surface_temperature)
    threaded_foreach(1:n; basesize = 1000) do i
        v.f_infiltration_reduction[i] = infiltration_reduction_factor(
            v.soil_surface_temperature[i],
            p.cf_soil[i];
            modelsnow,
            soil_infiltration_reduction,
        )
    end
    return nothing
end

"""
    infiltration!(soil_model::SbmSoilModel, dt::Float64)

Update the infiltration rate `infiltration` and infiltration excess water rate
`infiltration_excess` of the SBM soil model for a single timestep.
"""
function infiltration!(soil_model::SbmSoilModel, dt::Float64)
    v = soil_model.variables
    p = soil_model.parameters
    (; water_flux_surface) = soil_model.boundary_conditions

    n = length(v.infiltration)
    threaded_foreach(1:n; basesize = 1000) do i
        v.infiltration[i], v.infiltration_excess[i] = infiltration(
            water_flux_surface[i],
            p.compacted_soil_area_fraction[i],
            p.infiltration_capacity_soil[i],
            p.infiltration_capacity_compacted_soil[i],
            v.unsaturated_store_capacity[i],
            v.f_infiltration_reduction[i],
            dt,
        )
    end
    return nothing
end

"""
    unsaturated_zone_flow!(soil_model::SbmSoilModel, dt::Float64)

Update unsaturated storage `unsaturated_layer_depth` and the `transfer` of water from the unsaturated
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
            z = cumsum(v.unsaturated_layer_thickness[i])
            flow_rate = 0.0
            for m in 1:v.n_unsatlayers[i]
                l_sat = v.unsaturated_layer_thickness[i][m] * (p.theta_s[i] - p.theta_r[i])
                kv_z = hydraulic_conductivity_at_depth(p.kv_profile, p.vertical_hydraulic_conductivity_factor, z[m], i, m)
                unsaturated_layer_depth = if m == 1
                    v.unsaturated_layer_depth[i][m] + v.infiltration[i] * dt
                else
                    v.unsaturated_layer_depth[i][m] + flow_rate * dt
                end
                unsaturated_layer_depth, flow_rate =
                    unsatzone_flow_layer(unsaturated_layer_depth, kv_z, l_sat, p.brooks_corey_exponent[i][m], dt)
                v.unsaturated_layer_depth[i] = setindex(v.unsaturated_layer_depth[i], unsaturated_layer_depth, m)
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

Update soil evaporation from the saturated store `soil_evaporation_saturated_zone` and the total soil
evaporation from the unsaturated and saturated store `soil_evaporation` of the SBM soil model for a
single timestep. Also unsaturated storage `unsaturated_layer_depth` and the saturated store
`saturated_water_depth` are updated.
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
            v.unsaturated_layer_depth[i][1],
            v.unsaturated_layer_thickness[i][1],
            v.n_unsatlayers[i],
            v.water_table_depth[i],
            p.theta_s[i] - p.theta_r[i],
        )
        # Ensure that the unsaturated evaporation rate does not exceed the
        # available unsaturated moisture
        soilevapunsat = min(soilevapunsat, v.unsaturated_layer_depth[i][1] / dt)
        # Update the additional atmospheric demand
        potsoilevap -= soilevapunsat
        v.unsaturated_layer_depth[i] = setindex(
            v.unsaturated_layer_depth[i],
            v.unsaturated_layer_depth[i][1] - soilevapunsat * dt,
            1,
        )
        theta_drainable = lower_bound_drainable_porosity(p.theta_s[i], p.theta_fc[i])
        soil_evaporation_saturated_zone = soil_evaporation_saturated_store(
            potsoilevap,
            v.n_unsatlayers[i],
            p.actual_layer_thickness[i][1],
            v.water_table_depth[i],
            theta_drainable,
            dt,
        )
        v.soil_evaporation_saturated_zone[i] = soil_evaporation_saturated_zone
        v.soil_evaporation[i] = soilevapunsat + soil_evaporation_saturated_zone
        v.drainable_water_depth[i] -= soil_evaporation_saturated_zone * dt
    end
    return nothing
end

"""
    transpiration!(soil_model::SbmSoilModel, dt::Float64)

Update total `transpiration`, transpiration from the unsaturated store `actual_evaporation_unsaturated_store` and
saturated store `actual_evaporation_saturated_zone` of the SBM soil model for a single timestep. Also unsaturated
storage `unsaturated_layer_depth` and the saturated store `saturated_water_depth` are updated.
"""
function transpiration!(soil_model::SbmSoilModel, dt::Float64)
    (; potential_transpiration) = soil_model.boundary_conditions
    v = soil_model.variables
    p = soil_model.parameters

    rooting_depth = get_rootingdepth(soil_model)
    n = length(rooting_depth)

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
            if k == v.n_unsatlayers[i] && v.water_table_depth[i] < rooting_depth[i]
                rootlength = min(p.actual_layer_thickness[i][k], rooting_depth[i] - p.cumulative_layer_depth[i][k])
                rootfraction_unsat =
                    p.rootfraction[i][k] * (v.unsaturated_layer_thickness[i][k] / rootlength)
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
                rooting_depth[i] > 0.0 ?
                max((1.0 / sum_rootfraction_unsat), 1.0) * rootfraction_unsat : 0.0
            volumetric_water_content = max(v.unsaturated_layer_depth[i][k] / v.unsaturated_layer_thickness[i][k], 1e-7)
            head = head_brooks_corey(volumetric_water_content, p.theta_s[i], p.theta_r[i], p.brooks_corey_exponent[i][k], p.air_entry_pressure[i])
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
                    (rooting_depth[i] - p.cumulative_layer_depth[i][k]) / v.unsaturated_layer_thickness[i][k],
                ),
            )
            maxextr = v.unsaturated_layer_depth[i][k] * availcap / dt
            actevapustore_layer =
                min(alpha * rootfraction_unsat_scaled * potential_transpiration[i], maxextr)
            unsaturated_layer_depth = v.unsaturated_layer_depth[i][k] - actevapustore_layer * dt
            actevapustore += actevapustore_layer
            v.unsaturated_layer_depth[i] = setindex(v.unsaturated_layer_depth[i], unsaturated_layer_depth, k)
        end

        # transpiration from saturated store
        wetroots = scurve(v.water_table_depth[i], rooting_depth[i], 1.0, p.wet_root_distribution_parameter[i])
        alpha = rwu_reduction_feddes(
            Float64(0.0),
            p.h1[i],
            p.h2[i],
            v.h3[i],
            p.h4[i],
            p.alpha_h1[i],
        )
        restpottrans = potential_transpiration[i] - actevapustore
        actual_evaporation_saturated_zone = min(restpottrans * wetroots * alpha, v.drainable_water_depth[i] / dt)

        v.actual_evaporation_unsaturated_store[i] = actevapustore
        v.actual_evaporation_saturated_zone[i] = actual_evaporation_saturated_zone
        v.drainable_water_depth[i] -= actual_evaporation_saturated_zone * dt
        v.transpiration[i] = actevapustore + actual_evaporation_saturated_zone
    end
    return nothing
end

"""
    actual_infiltration!(soil_model::SbmSoilModel)

Update the actual infiltration rate `actual_infiltration` of the SBM soil model for a single timestep.

A soil water balance check is performed. Unsaturated storage that exceeds the maximum
storage per unsaturated soil layer is transferred to the layer above (or surface), from the
bottom to the top unsaturated soil layer. The resulting excess water `ustoredepth_excess` is
subtracted from the infiltration rate `infiltration`.
"""
function actual_infiltration!(soil_model::SbmSoilModel, dt::Float64)
    v = soil_model.variables
    p = soil_model.parameters

    n = length(v.actual_infiltration)
    threaded_foreach(1:n; basesize = 1000) do i
        # check soil moisture balance per layer
        ustoredepth_excess = 0.0
        for k in v.n_unsatlayers[i]:-1:1
            ustoredepth_excess = max(
                0.0,
                v.unsaturated_layer_depth[i][k] -
                v.unsaturated_layer_thickness[i][k] * (p.theta_s[i] - p.theta_r[i]),
            )
            v.unsaturated_layer_depth[i] = setindex(
                v.unsaturated_layer_depth[i],
                v.unsaturated_layer_depth[i][k] - ustoredepth_excess,
                k,
            )
            if k > 1
                v.unsaturated_layer_depth[i] = setindex(
                    v.unsaturated_layer_depth[i],
                    v.unsaturated_layer_depth[i][k - 1] + ustoredepth_excess,
                    k - 1,
                )
            end
        end
        v.actual_infiltration[i] = v.infiltration[i] - ustoredepth_excess / dt
    end
    return nothing
end

"""
    actual_infiltration_soil_path!(soil_model::SbmSoilModel)

Update the actual infiltration rate for soil `actual_infiltration_soil` and paved area `actual_infiltration_compacted_soil`
of the SBM soil model for a single timestep.
"""
function actual_infiltration_soil_path!(soil_model::SbmSoilModel)
    v = soil_model.variables
    p = soil_model.parameters
    (; water_flux_surface) = soil_model.boundary_conditions

    n = length(water_flux_surface)
    threaded_foreach(1:n; basesize = 1000) do i
        v.actual_infiltration_soil[i], v.actual_infiltration_compacted_soil[i] = actual_infiltration_soil_path(
            water_flux_surface[i],
            v.actual_infiltration[i],
            p.compacted_soil_area_fraction[i],
            p.infiltration_capacity_soil[i],
            p.infiltration_capacity_compacted_soil[i],
            v.f_infiltration_reduction[i],
        )
    end
    return nothing
end

"""
    capillary_flux!(soil_model::SbmSoilModel)

Update the capillary flux `actual_capillary_flux` of the SBM soil model for a single timestep.
"""
function capillary_flux!(soil_model::SbmSoilModel, dt::Float64)
    v = soil_model.variables
    p = soil_model.parameters
    rooting_depth = get_rootingdepth(soil_model)

    n = length(rooting_depth)
    threaded_foreach(1:n; basesize = 1000) do i
        if v.n_unsatlayers[i] > 0
            ksat = hydraulic_conductivity_at_depth(
                p.kv_profile,
                p.vertical_hydraulic_conductivity_factor,
                v.water_table_depth[i],
                i,
                v.n_unsatlayers[i],
            )
            maxcapflux = max(
                0.0,
                min(
                    ksat,
                    v.actual_evaporation_unsaturated_store[i],
                    v.unsaturated_store_capacity[i] / dt,
                    v.drainable_water_depth[i] / dt,
                ),
            )

            capflux = if v.water_table_depth[i] > rooting_depth[i]
                maxcapflux *
                pow(1.0 - min(v.water_table_depth[i], p.cap_hmax[i]) / (p.cap_hmax[i]), p.cap_n[i])
            else
                0.0
            end
            netcapflux = capflux
            actual_capillary_flux = 0.0
            for k in v.n_unsatlayers[i]:-1:1
                toadd = min(
                    netcapflux,
                    max(
                        (
                            v.unsaturated_layer_thickness[i][k] * (p.theta_s[i] - p.theta_r[i]) - v.unsaturated_layer_depth[i][k]
                        ) / dt,
                        0.0,
                    ),
                )
                v.unsaturated_layer_depth[i] = setindex(
                    v.unsaturated_layer_depth[i],
                    v.unsaturated_layer_depth[i][k] + toadd * dt,
                    k,
                )
                netcapflux -= toadd
                actual_capillary_flux += toadd
            end
            v.actual_capillary_flux[i] = actual_capillary_flux
        else
            v.actual_capillary_flux[i] = 0.0
        end
    end
    return nothing
end

"""
    leakage!(soil_model::SbmSoilModel, dt::Float64)

Update the actual leakage rate `actual_leakage` of the SBM soil model for a single timestep.
"""
function leakage!(soil_model::SbmSoilModel, dt::Float64)
    v = soil_model.variables
    p = soil_model.parameters

    n = length(v.actual_leakage)
    threaded_foreach(1:n; basesize = 1000) do i
        deepksat = hydraulic_conductivity_at_depth(
            p.kv_profile,
            p.vertical_hydraulic_conductivity_factor,
            p.soil_thickness[i],
            i,
            p.number_of_layers[i],
        )

        deeptransfer = min(v.drainable_water_depth[i] / dt, deepksat)
        v.actual_leakage[i] = max(0.0, min(p.maximum_leakage[i], deeptransfer))
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
    @. v.saturation_excess_water = (water_flux_surface - v.actual_infiltration) - v.infiltration_excess
    actual_infiltration_soil_path!(soil_model)
    @. v.excess_water_soil =
        max(water_flux_surface * (1.0 - p.compacted_soil_area_fraction) - v.actual_infiltration_soil, 0.0)
    @. v.excess_water_compacted_soil = max(water_flux_surface * p.compacted_soil_area_fraction - v.actual_infiltration_compacted_soil, 0.0)
    # recompute the unsaturated store and unsaturated_store_capacity (for capillary flux)
    unsaturated_store_depth!(soil_model)
    @. v.unsaturated_store_capacity = p.soil_water_capacity - v.saturated_water_depth - v.unsaturated_store_depth
    # capillary flux and leakage
    capillary_flux!(soil_model, dt)
    leakage!(soil_model, dt)
    # recharge rate to the saturated store
    @. v.recharge =
        (v.transfer - v.actual_capillary_flux - v.actual_leakage - v.actual_evaporation_saturated_zone - v.soil_evaporation_saturated_zone)
    # total actual evapotranspiration
    v.actual_evapotranspiration .=
        v.soil_evaporation .+ v.transpiration .+ runoff.variables.actual_open_water_evaporation_river .+
        runoff.variables.actual_open_water_evaporation_land .+ get_evaporation(demand.paddy)
    return nothing
end

function update_ustorelayerdepth!(soil, zi_prev, water_table_depth, i)
    v = soil.variables
    p = soil.parameters

    n_unsatlayers_prev = v.n_unsatlayers[i]
    ustorelayerthickness_prev = v.unsaturated_layer_thickness[i]
    unsaturated_layer_thickness = set_layerthickness(water_table_depth, p.cumulative_layer_depth[i], p.actual_layer_thickness[i])
    n_unsatlayers = number_of_active_layers(unsaturated_layer_thickness)
    unsaturated_layer_depth = v.unsaturated_layer_depth[i]
    if water_table_depth < zi_prev
        for k in n_unsatlayers:n_unsatlayers_prev
            k == 0 && continue
            if isnan(unsaturated_layer_thickness[k])
                unsaturated_layer_depth = setindex(unsaturated_layer_depth, 0.0, k)
            else
                hydraulic_conductivity_scale_parameter = unsaturated_layer_thickness[k] / ustorelayerthickness_prev[k]
                unsaturated_layer_depth = setindex(unsaturated_layer_depth, hydraulic_conductivity_scale_parameter * unsaturated_layer_depth[k], k)
            end
        end
    else
        for k in n_unsatlayers_prev:n_unsatlayers
            k == 0 && continue
            thickness_prev =
                isnan(ustorelayerthickness_prev[k]) ? 0.0 : ustorelayerthickness_prev[k]
            delta_thickness = unsaturated_layer_thickness[k] - thickness_prev
            unsaturated_layer_depth = setindex(
                unsaturated_layer_depth,
                unsaturated_layer_depth[k] + delta_thickness * (p.theta_fc[i] - p.theta_r[i]),
                k,
            )
        end
    end
    v.n_unsatlayers[i] = n_unsatlayers
    v.unsaturated_layer_depth[i] = unsaturated_layer_depth
    v.unsaturated_layer_thickness[i] = unsaturated_layer_thickness
    v.water_table_depth[i] = water_table_depth
end

"""
    update_ustorelayerdepth!(soil_model::SbmSoilModel, subsurface_flow)

Update the `SbmSoilModel` variables unsaturated store depth of soil layers
`unsaturated_layer_depth`, number of unsaturated zone soil layers `n_unsatlayers`, thickness of
unsaturated zone soil layers `unsaturated_layer_thickness` and water table depth `water_table_depth`, based on the
water table change computed by a subsurface flow model.
"""
function update_ustorelayerdepth!(soil_model::SbmSoilModel, subsurface_flow)
    p = soil_model.parameters
    v = soil_model.variables

    water_table_depth = get_water_depth(subsurface_flow)

    n = length(v.water_table_depth)
    threaded_foreach(1:n; basesize = 1000) do i
        zi_prev = v.water_table_depth[i]
        update_ustorelayerdepth!(soil_model, zi_prev, water_table_depth[i], i)
    end
end

"""
    update_soil_water_storage!(soil_model::SbmSoilModel, external_models::NamedTuple)

Update the SBM soil model for a single timestep based on the update of a subsurface flow
model, resulting in a change in water table depth and an exfiltration rate `exfiltwater_average`.

The available water in unsaturated zone `unsaturated_store_depth`, unsaturated store capacity
`unsaturated_store_capacity`, `total_soil_water_storage`, land `runoff` and `net_runoff`, the saturated
store `saturated_water_depth` and the water exfiltrating during saturation excess conditions
`exfiltration_saturated_water` are updated. Additionally, volumetric water content per soil layer and for
the root zone are updated.
"""
function update_soil_water_storage!(
    soil_model::SbmSoilModel,
    external_models::NamedTuple,
    dt::Float64,
)
    (; runoff, demand, subsurface_flow) = external_models
    (; runoff_land, actual_open_water_evaporation_land) = runoff.variables
    p = soil_model.parameters
    v = soil_model.variables
    exfiltration_saturated_water = subsurface_flow.variables.exfiltwater_average

    rooting_depth = get_rootingdepth(soil_model)

    n = length(v.water_table_depth)
    threaded_foreach(1:n; basesize = 1000) do i
        unsaturated_layer_depth = v.unsaturated_layer_depth[i]
        unsaturated_layer_thickness = v.unsaturated_layer_thickness[i]
        unsaturated_store_depth = sum(@view unsaturated_layer_depth[1:(v.n_unsatlayers[i])])
        sbm_runoff = max(
            0.0,
            exfiltration_saturated_water[i] + v.saturation_excess_water[i] + runoff_land[i] + v.infiltration_excess[i],
        )

        # volumetric water content per soil layer and root zone
        volumetric_water_content = v.volumetric_water_content[i]
        relative_volumetric_water_content = v.relative_volumetric_water_content[i]
        for k in 1:p.number_of_layers[i]
            if k <= v.n_unsatlayers[i]
                volumetric_water_content = setindex(
                    volumetric_water_content,
                    (
                        unsaturated_layer_depth[k] +
                        (p.actual_layer_thickness[i][k] - unsaturated_layer_thickness[k]) *
                        (p.theta_s[i] - p.theta_r[i])
                    ) / p.actual_layer_thickness[i][k] + p.theta_r[i],
                    k,
                )
            else
                volumetric_water_content = setindex(volumetric_water_content, p.theta_s[i], k)
            end
            relative_volumetric_water_content = setindex(relative_volumetric_water_content, from_SI(volumetric_water_content[k] / p.theta_s[i], PERCENTAGE), k)
        end

        rootstore_unsat = 0
        for k in 1:(v.n_unsatlayers[i])
            rootstore_unsat +=
                min(
                    1.0,
                    (
                        max(0.0, rooting_depth[i] - p.cumulative_layer_depth[i][k]) /
                        unsaturated_layer_thickness[k]
                    ),
                ) * unsaturated_layer_depth[k]
        end

        rootstore_sat = max(0.0, rooting_depth[i] - v.water_table_depth[i]) * (p.theta_s[i] - p.theta_r[i])
        root_zone_storage = rootstore_sat + rootstore_unsat
        volumetric_water_content_root_zone = root_zone_storage / rooting_depth[i] + p.theta_r[i]
        relative_volumetric_water_content_root_zone = from_SI(volumetric_water_content_root_zone / p.theta_s[i], PERCENTAGE)
        saturated_water_depth = (p.soil_thickness[i] - v.water_table_depth[i]) * (p.theta_s[i] - p.theta_r[i])
        drainable_water_depth =
            (p.soil_thickness[i] - v.water_table_depth[i]) *
            lower_bound_drainable_porosity(p.theta_s[i], p.theta_fc[i])
        unsaturated_store_capacity = p.soil_water_capacity[i] - saturated_water_depth - unsaturated_store_depth

        # update the outputs and states
        v.unsaturated_store_capacity[i] = unsaturated_store_capacity
        v.unsaturated_store_depth[i] = unsaturated_store_depth
        v.saturated_water_depth[i] = saturated_water_depth
        v.drainable_water_depth[i] = drainable_water_depth
        v.exfiltration_saturated_water[i] = exfiltration_saturated_water[i]
        v.runoff[i] = sbm_runoff
        v.volumetric_water_content[i] = volumetric_water_content
        v.relative_volumetric_water_content[i] = relative_volumetric_water_content
        v.root_zone_storage[i] = root_zone_storage
        v.volumetric_water_content_root_zone[i] = volumetric_water_content_root_zone
        v.relative_volumetric_water_content_root_zone[i] = relative_volumetric_water_content_root_zone
        v.total_soil_water_storage[i] = saturated_water_depth + unsaturated_store_depth
    end
    # update runoff and net_runoff (the runoff rate depends on the presence of paddy fields
    # and the h_max parameter of a paddy field)
    update_runoff!(demand.paddy, v.runoff, dt)
    @. v.net_runoff = v.runoff - actual_open_water_evaporation_land
    return nothing
end

"""
    update_diagnostic_vars!(soil_model::SbmSoilModel)

Update diagnostic variables of `SbmSoilModel` that are critical for subsequent computations
and depend on state variables `saturated_water_depth` and `unsaturated_layer_depth`.
"""
function update_diagnostic_vars!(soil_model::SbmSoilModel)
    (;
        water_table_depth,
        saturated_water_depth,
        drainable_water_depth,
        unsaturated_layer_thickness,
        unsaturated_store_capacity,
        unsaturated_store_depth,
        total_soil_water_storage,
        n_unsatlayers,
    ) = soil_model.variables
    (;
        soil_thickness,
        theta_s,
        theta_r,
        theta_fc,
        soil_water_capacity,
        cumulative_layer_depth,
        actual_layer_thickness,
    ) = soil_model.parameters

    unsaturated_store_depth!(soil_model)
    @. water_table_depth = max(0.0, soil_thickness - saturated_water_depth / (theta_s - theta_r))
    @. drainable_water_depth =
        (soil_thickness - water_table_depth) * lower_bound_drainable_porosity(theta_s, theta_fc)
    @. unsaturated_store_capacity = soil_water_capacity - saturated_water_depth - unsaturated_store_depth
    @. unsaturated_layer_thickness = set_layerthickness(water_table_depth, cumulative_layer_depth, actual_layer_thickness)
    @. n_unsatlayers = number_of_active_layers(unsaturated_layer_thickness)
    @. total_soil_water_storage = saturated_water_depth + unsaturated_store_depth
end

# wrapper method
get_rootingdepth(soil_model::SbmSoilModel) =
    soil_model.parameters.vegetation_parameter_set.rooting_depth
