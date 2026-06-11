abstract type AbstractSoilModel end

"Struct for storing SBM soil model variables"
@with_kw struct SbmSoilVariables{N}
    n_cells::Int
    # Calculated soil water pressure head h3 of the root water uptake reduction function (Feddes) [m]
    h3::Vector{Float64} = fill(MISSING_VALUE, n_cells)
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
    transpiration::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Actual evaporation from unsaturated store [m s⁻¹]
    actual_evaporation_unsaturated_store::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Soil evaporation from unsaturated and saturated store [m s⁻¹]
    soil_evaporation::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Soil evaporation from saturated store [m s⁻¹]
    soil_evaporation_saturated_zone::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Actual capillary rise [m s⁻¹]
    actual_capillary_flux::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Actual transpiration from saturated store [m s⁻¹]
    actual_evaporation_saturated_zone::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Total actual evapotranspiration [m s⁻¹]
    actual_evapotranspiration::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Actual infiltration into the unsaturated zone [m s⁻¹]
    actual_infiltration::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Actual infiltration non-compacted fraction [m s⁻¹]
    actual_infiltration_soil::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Actual infiltration compacted fraction [m s⁻¹]
    actual_infiltration_compacted_soil::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Actual infiltration (compacted and the non-compacted areas) [m s⁻¹]
    infiltration::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Infiltration excess water [m s⁻¹]
    infiltration_excess::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Water that cannot infiltrate due to saturated soil (saturation excess) [m s⁻¹]
    saturation_excess_water::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Water exfiltrating during saturation excess conditions [m s⁻¹]
    exfiltration_saturated_water::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Excess water for non-compacted fraction [m s⁻¹]
    excess_water_soil::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Excess water for compacted fraction [m s⁻¹]
    excess_water_compacted_soil::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Total surface runoff from infiltration and saturation excess (excluding actual open water evaporation) [m s⁻¹]
    runoff::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Net surface runoff (surface runoff - actual open water evaporation) [m s⁻¹]
    net_runoff::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Volumetric water content [-] per soil layer (including theta_r and saturated zone)
    volumetric_water_content::Vector{SVector{N, Float64}}
    # Volumetric water content [%] per soil layer (including theta_r and saturated zone)
    relative_volumetric_water_content::Vector{SVector{N, Float64}}
    # Root water storage [m] in unsaturated and saturated zone (excluding theta_r)
    root_zone_storage::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Volumetric water content [-] in root zone (including theta_r and saturated zone)
    volumetric_water_content_root_zone::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Volumetric water content [%] in root zone (including theta_r and saturated zone)
    relative_volumetric_water_content_root_zone::Vector{Float64} =
        fill(MISSING_VALUE, n_cells)
    # Amount of available water in the unsaturated zone [m]
    unsaturated_store_depth::Vector{Float64} = zeros(n_cells)
    # Downward flux from unsaturated to saturated zone [m s⁻¹]
    transfer::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Net recharge to saturated store [m s⁻¹]
    recharge::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Actual leakage from saturated store [m s⁻¹]
    actual_leakage::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Total water storage (excluding floodplain volume and reservoirs) [m]
    total_storage::Vector{Float64} = zeros(Float64, n_cells)
    # Total soil water storage [m]
    total_soil_water_storage::Vector{Float64}
    # Top soil temperature [K]
    soil_surface_temperature::Vector{Float64} =
        to_SI.(fill(10.0, n_cells), Ref(ABSOLUTE_DEGREES))
    # Soil infiltration reduction factor (when soil is frozen) [-]
    f_infiltration_reduction::Vector{Float64} = ones(n_cells)
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
    water_table_depth =
        @. max(0.0, soil_thickness - saturated_water_depth / (theta_s - theta_r))
    drainable_water_depth = @. (soil_thickness - water_table_depth) *
       lower_bound_drainable_porosity(theta_s, theta_fc)
    unsaturated_layer_thickness = set_layerthickness.(
        water_table_depth,
        cumulative_layer_depth,
        actual_layer_thickness,
    )
    n_unsatlayers = number_of_active_layers.(unsaturated_layer_thickness)
    volumetric_water_content = fill(MISSING_VALUE, maximum_number_of_layers, n)
    relative_volumetric_water_content = fill(MISSING_VALUE, maximum_number_of_layers, n)
    total_soil_water_storage = saturated_water_depth .+ unsaturated_store_depth

    vars = SbmSoilVariables(;
        n_cells = n,
        unsaturated_layer_depth = zero(actual_layer_thickness),
        unsaturated_store_capacity = soil_water_capacity .- saturated_water_depth,
        unsaturated_layer_thickness,
        saturated_water_depth,
        drainable_water_depth,
        water_table_depth,
        n_unsatlayers,
        volumetric_water_content = svectorscopy(
            volumetric_water_content,
            Val{maximum_number_of_layers}(),
        ),
        relative_volumetric_water_content = svectorscopy(
            relative_volumetric_water_content,
            Val{maximum_number_of_layers}(),
        ),
        total_soil_water_storage,
    )
    return vars
end

"Struct for storing SBM soil model boundary conditions"
@with_kw struct SbmSoilBC
    n_cells::Int
    # Water flux at the soil surface [m s⁻¹]
    water_flux_surface::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Potential transpiration rate [m s⁻¹]
    potential_transpiration::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Potential soil evaporation rate [m s⁻¹]
    potential_soilevaporation::Vector{Float64} = fill(MISSING_VALUE, n_cells)
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
    n_cells = length(indices)
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
            error(
                "$parname needs a layer dimension of size $maximum_number_of_layers, but is $size1",
            )
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
            nlayers_kv = fill(0, n_cells)
            for cell_idx in eachindex(nlayers_kv)
                layers =
                    @view cumulative_layer_depth[cell_idx][2:number_of_layers[cell_idx]]
                k = argmin(abs.(z_layered[cell_idx] .- layers))
                nlayers_kv[cell_idx] = k
                z_layered[cell_idx] = layers[k]
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
        error(
            "$parname needs a layer dimension of size $maximum_number_of_layers, but is $size1",
        )
    end
    brooks_corey_exponent =
        svectorscopy(brooks_corey_exponent, Val{maximum_number_of_layers}())
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
        error(
            "$parname needs a layer dimension of size $maximum_number_of_layers, but is $size1",
        )
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
        theta_fc = field_capacity.(
            actual_layer_thickness,
            number_of_layers,
            theta_s,
            theta_r,
            brooks_corey_exponent,
            air_entry_pressure,
        )
    end

    # optional root fraction
    rootfraction_name = "soil_root__length_density_fraction"
    if haskey(config.input.static, rootfraction_name)
        rootfraction =
            ncread(dataset, config, rootfraction_name, LandHydrologySBM; sel = indices)
    else
        n_cells = length(indices)
        (; rooting_depth) = vegetation_parameter_set
        # default root fraction
        rootfraction = zeros(maximum_number_of_layers, n_cells)
        for cell_idx in 1:n_cells
            if rooting_depth[cell_idx] > 0.0
                for soil_layer_idx in 1:maximum_number_of_layers
                    if (
                        rooting_depth[cell_idx] -
                        cumulative_layer_depth[cell_idx][soil_layer_idx]
                    ) >= actual_layer_thickness[cell_idx][soil_layer_idx]
                        rootfraction[soil_layer_idx, cell_idx] =
                            actual_layer_thickness[cell_idx][soil_layer_idx] /
                            rooting_depth[cell_idx]
                    else
                        rootfraction[soil_layer_idx, cell_idx] = max(
                            (
                                rooting_depth[cell_idx] -
                                cumulative_layer_depth[cell_idx][soil_layer_idx]
                            ) / rooting_depth[cell_idx],
                            0.0,
                        )
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

    n_cells = length(indices)
    sbm_params = SbmSoilParameters(;
        maximum_number_of_layers,
        number_of_layers,
        soil_water_capacity,
        theta_s,
        theta_r,
        theta_fc,
        vertical_hydraulic_conductivity_factor = svectorscopy(
            vertical_hydraulic_conductivity_factor,
            Val{maximum_number_of_layers}(),
        ),
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
        soil_fraction = fill(MISSING_VALUE, n_cells),
        kv_profile,
        vegetation_parameter_set,
    )
    return sbm_params
end

"SBM soil model"
@with_kw struct SbmSoilModel{N, M, Kv} <: AbstractSoilModel
    n_cells::Int
    boundary_conditions::SbmSoilBC = SbmSoilBC(; n_cells)
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
    n_cells = length(indices)
    parameters = SbmSoilParameters(dataset, config, vegetation_parameter_set, indices, dt)
    variables = SbmSoilVariables(n_cells, parameters)
    soil_model = SbmSoilModel(; n_cells, parameters, variables)
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
    water_flux_surface .= max.(
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
    @. v.soil_surface_temperature =
        soil_temperature(v.soil_surface_temperature, p.w_soil, temperature)
    return nothing
end

soil_temperature!(::SbmSoilModel, ::NoSnowModel, ::Vector{Float64}) = nothing

"Update total available water in the unsaturated zone of the SBM soil model for a single timestep"
function unsaturated_store_depth!(soil_model::SbmSoilModel)
    v = soil_model.variables
    p = soil_model.parameters
    for cell_idx in eachindex(v.unsaturated_layer_depth)
        v.unsaturated_store_depth[cell_idx] =
            sum(@view v.unsaturated_layer_depth[cell_idx][1:p.number_of_layers[cell_idx]])
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

    n_cells = length(v.soil_surface_temperature)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        v.f_infiltration_reduction[cell_idx] = infiltration_reduction_factor(
            v.soil_surface_temperature[cell_idx],
            p.cf_soil[cell_idx];
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

    n_cells = length(v.infiltration)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        v.infiltration[cell_idx], v.infiltration_excess[cell_idx] = infiltration(
            water_flux_surface[cell_idx],
            p.compacted_soil_area_fraction[cell_idx],
            p.infiltration_capacity_soil[cell_idx],
            p.infiltration_capacity_compacted_soil[cell_idx],
            v.unsaturated_store_capacity[cell_idx],
            v.f_infiltration_reduction[cell_idx],
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

    n_cells = length(v.transfer)
    threaded_foreach(1:n_cells; basesize = 250) do cell_idx
        if v.n_unsatlayers[cell_idx] > 0
            # Brooks-Corey approach
            z = cumsum(v.unsaturated_layer_thickness[cell_idx])
            flow_rate = 0.0
            for soil_layer_idx in 1:v.n_unsatlayers[cell_idx]
                l_sat =
                    v.unsaturated_layer_thickness[cell_idx][soil_layer_idx] *
                    (p.theta_s[cell_idx] - p.theta_r[cell_idx])
                kv_z = hydraulic_conductivity_at_depth(
                    p.kv_profile,
                    p.vertical_hydraulic_conductivity_factor,
                    z[soil_layer_idx],
                    cell_idx,
                    soil_layer_idx,
                )
                unsaturated_layer_depth = if soil_layer_idx == 1
                    v.unsaturated_layer_depth[cell_idx][soil_layer_idx] +
                    v.infiltration[cell_idx] * dt
                else
                    v.unsaturated_layer_depth[cell_idx][soil_layer_idx] + flow_rate * dt
                end
                unsaturated_layer_depth, flow_rate = unsatzone_flow_layer(
                    unsaturated_layer_depth,
                    kv_z,
                    l_sat,
                    p.brooks_corey_exponent[cell_idx][soil_layer_idx],
                    dt,
                )
                v.unsaturated_layer_depth[cell_idx] = setindex(
                    v.unsaturated_layer_depth[cell_idx],
                    unsaturated_layer_depth,
                    soil_layer_idx,
                )
            end
            v.transfer[cell_idx] = flow_rate
        else
            v.transfer[cell_idx] = 0.0
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

    n_cells = length(potential_soilevaporation)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        potsoilevap = potential_soilevaporation[cell_idx]
        # First calculate the evaporation of unsaturated storage into the
        # atmosphere from the upper layer.
        soilevapunsat = soil_evaporation_unsaturated_store(
            potsoilevap,
            v.unsaturated_layer_depth[cell_idx][1],
            v.unsaturated_layer_thickness[cell_idx][1],
            v.n_unsatlayers[cell_idx],
            v.water_table_depth[cell_idx],
            p.theta_s[cell_idx] - p.theta_r[cell_idx],
        )
        # Ensure that the unsaturated evaporation rate does not exceed the
        # available unsaturated moisture
        soilevapunsat = min(soilevapunsat, v.unsaturated_layer_depth[cell_idx][1] / dt)
        # Update the additional atmospheric demand
        potsoilevap -= soilevapunsat
        v.unsaturated_layer_depth[cell_idx] = setindex(
            v.unsaturated_layer_depth[cell_idx],
            v.unsaturated_layer_depth[cell_idx][1] - soilevapunsat * dt,
            1,
        )
        theta_drainable =
            lower_bound_drainable_porosity(p.theta_s[cell_idx], p.theta_fc[cell_idx])
        soil_evaporation_saturated_zone = soil_evaporation_saturated_store(
            potsoilevap,
            v.n_unsatlayers[cell_idx],
            p.actual_layer_thickness[cell_idx][1],
            v.water_table_depth[cell_idx],
            theta_drainable,
            dt,
        )
        v.soil_evaporation_saturated_zone[cell_idx] = soil_evaporation_saturated_zone
        v.soil_evaporation[cell_idx] = soilevapunsat + soil_evaporation_saturated_zone
        v.drainable_water_depth[cell_idx] -= soil_evaporation_saturated_zone * dt
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
    n_cells = length(rooting_depth)

    threaded_foreach(1:n_cells; basesize = 250) do cell_idx
        v.h3[cell_idx] = feddes_h3(
            p.h3_high[cell_idx],
            p.h3_low[cell_idx],
            potential_transpiration[cell_idx],
        )

        # compute sum of root fraction in unsaturated soil layers and adapt root fraction
        # lowest unsaturated soil layer if water table depth intersects the unsaturated root
        # zone
        sum_rootfraction_unsat = 0.0
        rootfraction_unsat_lowest = 0.0
        for soil_layer_idx in 1:v.n_unsatlayers[cell_idx]
            # the root fraction is valid for the root length in a soil layer, if zi decreases
            # the root length the root fraction needs to be adapted
            if soil_layer_idx == v.n_unsatlayers[cell_idx] &&
               v.water_table_depth[cell_idx] < rooting_depth[cell_idx]
                rootlength = min(
                    p.actual_layer_thickness[cell_idx][soil_layer_idx],
                    rooting_depth[cell_idx] -
                    p.cumulative_layer_depth[cell_idx][soil_layer_idx],
                )
                rootfraction_unsat =
                    p.rootfraction[cell_idx][soil_layer_idx] *
                    (v.unsaturated_layer_thickness[cell_idx][soil_layer_idx] / rootlength)
            else
                rootfraction_unsat = p.rootfraction[cell_idx][soil_layer_idx]
            end
            sum_rootfraction_unsat += rootfraction_unsat

            # rootfraction lowest unsaturated layer
            rootfraction_unsat_lowest = rootfraction_unsat
        end

        actevapustore = 0.0
        for soil_layer_idx in 1:v.n_unsatlayers[cell_idx]
            # scale rootfraction soil layer unsaturated zone based on sum of rootfraction in
            # unsaturated zone
            if soil_layer_idx < v.n_unsatlayers[cell_idx]
                rootfraction_unsat = p.rootfraction[cell_idx][soil_layer_idx]
            else
                rootfraction_unsat = rootfraction_unsat_lowest
            end
            rootfraction_unsat_scaled =
                rooting_depth[cell_idx] > 0.0 ?
                max((1.0 / sum_rootfraction_unsat), 1.0) * rootfraction_unsat : 0.0
            volumetric_water_content = max(
                v.unsaturated_layer_depth[cell_idx][soil_layer_idx] /
                v.unsaturated_layer_thickness[cell_idx][soil_layer_idx],
                1e-7,
            )
            head = head_brooks_corey(
                volumetric_water_content,
                p.theta_s[cell_idx],
                p.theta_r[cell_idx],
                p.brooks_corey_exponent[cell_idx][soil_layer_idx],
                p.air_entry_pressure[cell_idx],
            )
            alpha = rwu_reduction_feddes(
                head,
                p.h1[cell_idx],
                p.h2[cell_idx],
                v.h3[cell_idx],
                p.h4[cell_idx],
                p.alpha_h1[cell_idx],
            )
            availcap = min(
                1.0,
                max(
                    0.0,
                    (
                        rooting_depth[cell_idx] -
                        p.cumulative_layer_depth[cell_idx][soil_layer_idx]
                    ) / v.unsaturated_layer_thickness[cell_idx][soil_layer_idx],
                ),
            )
            maxextr = v.unsaturated_layer_depth[cell_idx][soil_layer_idx] * availcap / dt
            actevapustore_layer = min(
                alpha * rootfraction_unsat_scaled * potential_transpiration[cell_idx],
                maxextr,
            )
            unsaturated_layer_depth =
                v.unsaturated_layer_depth[cell_idx][soil_layer_idx] -
                actevapustore_layer * dt
            actevapustore += actevapustore_layer
            v.unsaturated_layer_depth[cell_idx] = setindex(
                v.unsaturated_layer_depth[cell_idx],
                unsaturated_layer_depth,
                soil_layer_idx,
            )
        end

        # transpiration from saturated store
        wetroots = scurve(
            v.water_table_depth[cell_idx],
            rooting_depth[cell_idx],
            1.0,
            p.wet_root_distribution_parameter[cell_idx],
        )
        alpha = rwu_reduction_feddes(
            Float64(0.0),
            p.h1[cell_idx],
            p.h2[cell_idx],
            v.h3[cell_idx],
            p.h4[cell_idx],
            p.alpha_h1[cell_idx],
        )
        restpottrans = potential_transpiration[cell_idx] - actevapustore
        actual_evaporation_saturated_zone =
            min(restpottrans * wetroots * alpha, v.drainable_water_depth[cell_idx] / dt)

        v.actual_evaporation_unsaturated_store[cell_idx] = actevapustore
        v.actual_evaporation_saturated_zone[cell_idx] = actual_evaporation_saturated_zone
        v.drainable_water_depth[cell_idx] -= actual_evaporation_saturated_zone * dt
        v.transpiration[cell_idx] = actevapustore + actual_evaporation_saturated_zone
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

    n_cells = length(v.actual_infiltration)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        # check soil moisture balance per layer
        ustoredepth_excess = 0.0
        for soil_layer_idx in v.n_unsatlayers[cell_idx]:-1:1
            ustoredepth_excess = max(
                0.0,
                v.unsaturated_layer_depth[cell_idx][soil_layer_idx] -
                v.unsaturated_layer_thickness[cell_idx][soil_layer_idx] *
                (p.theta_s[cell_idx] - p.theta_r[cell_idx]),
            )
            v.unsaturated_layer_depth[cell_idx] = setindex(
                v.unsaturated_layer_depth[cell_idx],
                v.unsaturated_layer_depth[cell_idx][soil_layer_idx] - ustoredepth_excess,
                soil_layer_idx,
            )
            if soil_layer_idx > 1
                v.unsaturated_layer_depth[cell_idx] = setindex(
                    v.unsaturated_layer_depth[cell_idx],
                    v.unsaturated_layer_depth[cell_idx][soil_layer_idx - 1] +
                    ustoredepth_excess,
                    soil_layer_idx - 1,
                )
            end
        end
        v.actual_infiltration[cell_idx] = v.infiltration[cell_idx] - ustoredepth_excess / dt
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

    n_cells = length(water_flux_surface)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        v.actual_infiltration_soil[cell_idx],
        v.actual_infiltration_compacted_soil[cell_idx] = actual_infiltration_soil_path(
            water_flux_surface[cell_idx],
            v.actual_infiltration[cell_idx],
            p.compacted_soil_area_fraction[cell_idx],
            p.infiltration_capacity_soil[cell_idx],
            p.infiltration_capacity_compacted_soil[cell_idx],
            v.f_infiltration_reduction[cell_idx],
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

    n_cells = length(rooting_depth)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        if v.n_unsatlayers[cell_idx] > 0
            ksat = hydraulic_conductivity_at_depth(
                p.kv_profile,
                p.vertical_hydraulic_conductivity_factor,
                v.water_table_depth[cell_idx],
                cell_idx,
                v.n_unsatlayers[cell_idx],
            )
            maxcapflux = max(
                0.0,
                min(
                    ksat,
                    v.actual_evaporation_unsaturated_store[cell_idx],
                    v.unsaturated_store_capacity[cell_idx] / dt,
                    v.drainable_water_depth[cell_idx] / dt,
                ),
            )

            capflux = if v.water_table_depth[cell_idx] > rooting_depth[cell_idx]
                maxcapflux * pow(
                    1.0 -
                    min(v.water_table_depth[cell_idx], p.cap_hmax[cell_idx]) /
                    (p.cap_hmax[cell_idx]),
                    p.cap_n[cell_idx],
                )
            else
                0.0
            end
            netcapflux = capflux
            actual_capillary_flux = 0.0
            for soil_layer_idx in v.n_unsatlayers[cell_idx]:-1:1
                toadd = min(
                    netcapflux,
                    max(
                        (
                            v.unsaturated_layer_thickness[cell_idx][soil_layer_idx] *
                            (p.theta_s[cell_idx] - p.theta_r[cell_idx]) -
                            v.unsaturated_layer_depth[cell_idx][soil_layer_idx]
                        ) / dt,
                        0.0,
                    ),
                )
                v.unsaturated_layer_depth[cell_idx] = setindex(
                    v.unsaturated_layer_depth[cell_idx],
                    v.unsaturated_layer_depth[cell_idx][soil_layer_idx] + toadd * dt,
                    soil_layer_idx,
                )
                netcapflux -= toadd
                actual_capillary_flux += toadd
            end
            v.actual_capillary_flux[cell_idx] = actual_capillary_flux
        else
            v.actual_capillary_flux[cell_idx] = 0.0
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

    n_cells = length(v.actual_leakage)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        deepksat = hydraulic_conductivity_at_depth(
            p.kv_profile,
            p.vertical_hydraulic_conductivity_factor,
            p.soil_thickness[cell_idx],
            cell_idx,
            p.number_of_layers[cell_idx],
        )

        deeptransfer = min(v.drainable_water_depth[cell_idx] / dt, deepksat)
        v.actual_leakage[cell_idx] =
            max(0.0, min(p.maximum_leakage[cell_idx], deeptransfer))
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
    @. v.saturation_excess_water =
        (water_flux_surface - v.actual_infiltration) - v.infiltration_excess
    actual_infiltration_soil_path!(soil_model)
    @. v.excess_water_soil = max(
        water_flux_surface * (1.0 - p.compacted_soil_area_fraction) -
        v.actual_infiltration_soil,
        0.0,
    )
    @. v.excess_water_compacted_soil = max(
        water_flux_surface * p.compacted_soil_area_fraction -
        v.actual_infiltration_compacted_soil,
        0.0,
    )
    # recompute the unsaturated store and unsaturated_store_capacity (for capillary flux)
    unsaturated_store_depth!(soil_model)
    @. v.unsaturated_store_capacity =
        p.soil_water_capacity - v.saturated_water_depth - v.unsaturated_store_depth
    # capillary flux and leakage
    capillary_flux!(soil_model, dt)
    leakage!(soil_model, dt)
    # recharge rate to the saturated store
    @. v.recharge = (
        v.transfer - v.actual_capillary_flux - v.actual_leakage -
        v.actual_evaporation_saturated_zone - v.soil_evaporation_saturated_zone
    )
    # total actual evapotranspiration
    v.actual_evapotranspiration .=
        v.soil_evaporation .+ v.transpiration .+
        runoff.variables.actual_open_water_evaporation_river .+
        runoff.variables.actual_open_water_evaporation_land .+ get_evaporation(demand.paddy)
    return nothing
end

function update_ustorelayerdepth!(soil, zi_prev, water_table_depth, i)
    v = soil.variables
    p = soil.parameters

    n_unsatlayers_prev = v.n_unsatlayers[i]
    ustorelayerthickness_prev = v.unsaturated_layer_thickness[i]
    unsaturated_layer_thickness = set_layerthickness(
        water_table_depth,
        p.cumulative_layer_depth[i],
        p.actual_layer_thickness[i],
    )
    n_unsatlayers = number_of_active_layers(unsaturated_layer_thickness)
    unsaturated_layer_depth = v.unsaturated_layer_depth[i]
    if water_table_depth < zi_prev
        for soil_layer_idx in n_unsatlayers:n_unsatlayers_prev
            soil_layer_idx == 0 && continue
            if isnan(unsaturated_layer_thickness[soil_layer_idx])
                unsaturated_layer_depth =
                    setindex(unsaturated_layer_depth, 0.0, soil_layer_idx)
            else
                hydraulic_conductivity_scale_parameter =
                    unsaturated_layer_thickness[soil_layer_idx] /
                    ustorelayerthickness_prev[soil_layer_idx]
                unsaturated_layer_depth = setindex(
                    unsaturated_layer_depth,
                    hydraulic_conductivity_scale_parameter *
                    unsaturated_layer_depth[soil_layer_idx],
                    soil_layer_idx,
                )
            end
        end
    else
        for soil_layer_idx in n_unsatlayers_prev:n_unsatlayers
            soil_layer_idx == 0 && continue
            thickness_prev =
                isnan(ustorelayerthickness_prev[soil_layer_idx]) ? 0.0 :
                ustorelayerthickness_prev[soil_layer_idx]
            delta_thickness = unsaturated_layer_thickness[soil_layer_idx] - thickness_prev
            unsaturated_layer_depth = setindex(
                unsaturated_layer_depth,
                unsaturated_layer_depth[soil_layer_idx] +
                delta_thickness * (p.theta_fc[i] - p.theta_r[i]),
                soil_layer_idx,
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

    n_cells = length(v.water_table_depth)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        zi_prev = v.water_table_depth[cell_idx]
        update_ustorelayerdepth!(soil_model, zi_prev, water_table_depth[cell_idx], cell_idx)
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

    n_cells = length(v.water_table_depth)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        unsaturated_layer_depth = v.unsaturated_layer_depth[cell_idx]
        unsaturated_layer_thickness = v.unsaturated_layer_thickness[cell_idx]
        unsaturated_store_depth =
            sum(@view unsaturated_layer_depth[1:(v.n_unsatlayers[cell_idx])])
        sbm_runoff = max(
            0.0,
            exfiltration_saturated_water[cell_idx] +
            v.saturation_excess_water[cell_idx] +
            runoff_land[cell_idx] +
            v.infiltration_excess[cell_idx],
        )

        # volumetric water content per soil layer and root zone
        volumetric_water_content = v.volumetric_water_content[cell_idx]
        relative_volumetric_water_content = v.relative_volumetric_water_content[cell_idx]
        for soil_layer_idx in 1:p.number_of_layers[cell_idx]
            if soil_layer_idx <= v.n_unsatlayers[cell_idx]
                volumetric_water_content = setindex(
                    volumetric_water_content,
                    (
                        unsaturated_layer_depth[soil_layer_idx] +
                        (
                            p.actual_layer_thickness[cell_idx][soil_layer_idx] -
                            unsaturated_layer_thickness[soil_layer_idx]
                        ) * (p.theta_s[cell_idx] - p.theta_r[cell_idx])
                    ) / p.actual_layer_thickness[cell_idx][soil_layer_idx] +
                    p.theta_r[cell_idx],
                    soil_layer_idx,
                )
            else
                volumetric_water_content =
                    setindex(volumetric_water_content, p.theta_s[cell_idx], soil_layer_idx)
            end
            relative_volumetric_water_content = setindex(
                relative_volumetric_water_content,
                from_SI(
                    volumetric_water_content[soil_layer_idx] / p.theta_s[cell_idx],
                    PERCENTAGE,
                ),
                soil_layer_idx,
            )
        end

        rootstore_unsat = 0
        for soil_layer_idx in 1:(v.n_unsatlayers[cell_idx])
            rootstore_unsat +=
                min(
                    1.0,
                    (
                        max(
                            0.0,
                            rooting_depth[cell_idx] -
                            p.cumulative_layer_depth[cell_idx][soil_layer_idx],
                        ) / unsaturated_layer_thickness[soil_layer_idx]
                    ),
                ) * unsaturated_layer_depth[soil_layer_idx]
        end

        rootstore_sat =
            max(0.0, rooting_depth[cell_idx] - v.water_table_depth[cell_idx]) *
            (p.theta_s[cell_idx] - p.theta_r[cell_idx])
        root_zone_storage = rootstore_sat + rootstore_unsat
        volumetric_water_content_root_zone =
            root_zone_storage / rooting_depth[cell_idx] + p.theta_r[cell_idx]
        relative_volumetric_water_content_root_zone =
            from_SI(volumetric_water_content_root_zone / p.theta_s[cell_idx], PERCENTAGE)
        saturated_water_depth =
            (p.soil_thickness[cell_idx] - v.water_table_depth[cell_idx]) *
            (p.theta_s[cell_idx] - p.theta_r[cell_idx])
        drainable_water_depth =
            (p.soil_thickness[cell_idx] - v.water_table_depth[cell_idx]) *
            lower_bound_drainable_porosity(p.theta_s[cell_idx], p.theta_fc[cell_idx])
        unsaturated_store_capacity =
            p.soil_water_capacity[cell_idx] - saturated_water_depth -
            unsaturated_store_depth

        # update the outputs and states
        v.unsaturated_store_capacity[cell_idx] = unsaturated_store_capacity
        v.unsaturated_store_depth[cell_idx] = unsaturated_store_depth
        v.saturated_water_depth[cell_idx] = saturated_water_depth
        v.drainable_water_depth[cell_idx] = drainable_water_depth
        v.exfiltration_saturated_water[cell_idx] = exfiltration_saturated_water[cell_idx]
        v.runoff[cell_idx] = sbm_runoff
        v.volumetric_water_content[cell_idx] = volumetric_water_content
        v.relative_volumetric_water_content[cell_idx] = relative_volumetric_water_content
        v.root_zone_storage[cell_idx] = root_zone_storage
        v.volumetric_water_content_root_zone[cell_idx] = volumetric_water_content_root_zone
        v.relative_volumetric_water_content_root_zone[cell_idx] =
            relative_volumetric_water_content_root_zone
        v.total_soil_water_storage[cell_idx] =
            saturated_water_depth + unsaturated_store_depth
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
    @. water_table_depth =
        max(0.0, soil_thickness - saturated_water_depth / (theta_s - theta_r))
    @. drainable_water_depth =
        (soil_thickness - water_table_depth) *
        lower_bound_drainable_porosity(theta_s, theta_fc)
    @. unsaturated_store_capacity =
        soil_water_capacity - saturated_water_depth - unsaturated_store_depth
    @. unsaturated_layer_thickness = set_layerthickness(
        water_table_depth,
        cumulative_layer_depth,
        actual_layer_thickness,
    )
    @. n_unsatlayers = number_of_active_layers(unsaturated_layer_thickness)
    @. total_soil_water_storage = saturated_water_depth + unsaturated_store_depth
end

# wrapper method
get_rootingdepth(soil_model::SbmSoilModel) =
    soil_model.parameters.vegetation_parameter_set.rooting_depth
