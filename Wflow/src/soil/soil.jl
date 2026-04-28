abstract type AbstractSoilModel end

"Struct for storing SBM soil model variables"
@with_kw struct SbmSoilVariables{N}
    n_land_cells::Int
    # Calculated soil water pressure head h3 of the root water uptake reduction function (Feddes) [cm]
    h3::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Unsaturated store capacity [mm]
    ustorecapacity::Vector{Float64}
    # Amount of water in the unsaturated store, per layer [mm]
    ustorelayerdepth::Vector{SVector{N,Float64}}
    # Thickness of unsaturated zone, per layer [mm]
    ustorelayerthickness::Vector{SVector{N,Float64}}
    # Saturated store [mm]
    satwaterdepth::Vector{Float64}
    # Drainable water store [mm]
    drainable_waterdepth::Vector{Float64}
    # Pseudo-water table depth [mm] (top of the saturated zone)
    zi::Vector{Float64}
    # Number of unsaturated soil layers
    n_unsatlayers::Vector{Int}
    # Transpiration [mm Δt⁻¹]
    transpiration::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Actual evaporation from unsaturated store [mm Δt⁻¹]
    ae_ustore::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Soil evaporation from unsaturated and saturated store [mm Δt⁻¹]
    soilevap::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Soil evaporation from saturated store [mm Δt⁻¹]
    soilevapsat::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Actual capillary rise [mm Δt⁻¹]
    actcapflux::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Actual transpiration from saturated store [mm Δt⁻¹]
    actevapsat::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Total actual evapotranspiration [mm Δt⁻¹]
    actevap::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Actual infiltration into the unsaturated zone [mm Δt⁻¹]
    actinfilt::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Actual infiltration non-compacted fraction [mm Δt⁻¹]
    actinfiltsoil::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Actual infiltration compacted fraction [mm Δt⁻¹]
    actinfiltpath::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Actual infiltration (compacted and the non-compacted areas) [mm Δt⁻¹]
    infiltsoilpath::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Infiltration excess water [mm Δt⁻¹]
    infiltexcess::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Water that cannot infiltrate due to saturated soil (saturation excess) [mm Δt⁻¹]
    excesswater::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Water exfiltrating during saturation excess conditions [mm Δt⁻¹]
    exfiltsatwater::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Excess water for non-compacted fraction [mm Δt⁻¹]
    excesswatersoil::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Excess water for compacted fraction [mm Δt⁻¹]
    excesswaterpath::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Total surface runoff from infiltration and saturation excess (excluding actual open water evaporation) [mm Δt⁻¹]
    runoff::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Net surface runoff (surface runoff - actual open water evaporation) [mm Δt⁻¹]
    net_runoff::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Volumetric water content [-] per soil layer (including theta_r and saturated zone)
    vwc::Vector{SVector{N,Float64}}
    # Volumetric water content [%] per soil layer (including theta_r and saturated zone)
    vwc_perc::Vector{SVector{N,Float64}}
    # Root water storage [mm] in unsaturated and saturated zone (excluding theta_r)
    rootstore::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Volumetric water content [-] in root zone (including theta_r and saturated zone)
    vwc_root::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Volumetric water content [%] in root zone (including theta_r and saturated zone)
    vwc_percroot::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Amount of available water in the unsaturated zone [mm]
    ustoredepth::Vector{Float64} = zeros(n_land_cells)
    # Downward flux from unsaturated to saturated zone [mm Δt⁻¹]
    transfer::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Net recharge to saturated store [mm Δt⁻¹]
    recharge::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Actual leakage from saturated store [mm Δt⁻¹]
    actleakage::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Total water storage (excluding floodplain volume and reservoirs) [mm]
    total_storage::Vector{Float64} = zeros(n_land_cells)
    # Total soil water storage [mm]
    total_soilwater_storage::Vector{Float64}
    # Top soil temperature [ᵒC]
    tsoil::Vector{Float64} = fill(10.0, n_land_cells)
    # Soil infiltration reduction factor (when soil is frozen) [-]
    f_infiltration_reduction::Vector{Float64} = ones(n_land_cells)
end

"Struct for storing SBM soil model parameters"
@with_kw struct SbmSoilParameters{N,M,Kv}
    # Maximum number of soil layers [-]
    maxlayers::Int
    # Number of soil layers [-]
    nlayers::Vector{Int}
    # Saturated water content (porosity) [-]
    theta_s::Vector{Float64}
    # Residual water content [-]
    theta_r::Vector{Float64}
    # Field capacity water content [-]
    theta_fc::Vector{Float64}
    # Soilwater capacity [mm]
    soilwatercapacity::Vector{Float64}
    # Multiplication factor [-] applied to kv_z (vertical flow)
    kvfrac::Vector{SVector{N,Float64}}
    # Air entry pressure [cm] of soil (Brooks-Corey)
    hb::Vector{Float64}
    # Soil thickness [mm]
    soilthickness::Vector{Float64}
    # Thickness of soil layers [mm]
    act_thickl::Vector{SVector{N,Float64}}
    # Cumulative sum of soil layers [mm], starting at soil surface (0)
    sumlayers::Vector{SVector{M,Float64}}
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
    c::Vector{SVector{N,Float64}}
    # Soil temperature smooth factor [-]
    w_soil::Vector{Float64}
    # Controls soil infiltration reduction factor when soil is frozen [-]
    cf_soil::Vector{Float64}
    # Fraction of compacted area  [-]
    pathfrac::Vector{Float64}
    # Controls how roots are linked to water table [-]
    rootdistpar::Vector{Float64}
    # Fraction of the root length density in each soil layer [-]
    rootfraction::Vector{SVector{N,Float64}}
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
function SbmSoilVariables(n_land_cells::Int, parameters::SbmSoilParameters)
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
    ustoredepth = zeros(n_land_cells)
    zi = @. max(0.0, soilthickness - satwaterdepth / (theta_s - theta_r))
    drainable_waterdepth =
        @. (soilthickness - zi) * lower_bound_drainable_porosity(theta_s, theta_fc)
    ustorelayerthickness = set_layerthickness.(zi, sumlayers, act_thickl)
    n_unsatlayers = number_of_active_layers.(ustorelayerthickness)

    vwc = fill(MISSING_VALUE, maxlayers, n_land_cells)
    vwc_perc = fill(MISSING_VALUE, maxlayers, n_land_cells)
    total_soilwater_storage = satwaterdepth .+ ustoredepth

    vars = SbmSoilVariables(;
        n_land_cells,
        ustorelayerdepth=zero(act_thickl),
        ustorecapacity=soilwatercapacity .- satwaterdepth,
        ustorelayerthickness,
        satwaterdepth,
        drainable_waterdepth,
        zi,
        n_unsatlayers,
        vwc=svectorscopy(vwc, Val{maxlayers}()),
        vwc_perc=svectorscopy(vwc_perc, Val{maxlayers}()),
        total_soilwater_storage,
    )
    return vars
end

"Struct for storing SBM soil model boundary conditions"
@with_kw struct SbmSoilBC
    n_land_cells::Int
    # Water flux at the soil surface [mm Δt⁻¹]
    water_flux_surface::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Potential transpiration rate [mm Δt⁻¹]
    potential_transpiration::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
    # Potential soil evaporation rate [mm Δt⁻¹]
    potential_soilevaporation::Vector{Float64} = fill(MISSING_VALUE, n_land_cells)
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
    kv::Vector{SVector{N,Float64}}
end

"Layered exponential depth profile of vertical hydraulic conductivity"
struct KvLayeredExponential{N}
    # A scaling parameter [mm⁻¹] (controls exponential decline of kv_0)
    f::Vector{Float64}
    # Vertical hydraulic conductivity [mm Δt⁻¹] per soil layer
    kv::Vector{SVector{N,Float64}}
    # Number of soil layers [-] with vertical hydraulic conductivity value `kv`
    nlayers_kv::Vector{Int}
    # Depth [mm] from soil surface for which layered profile is valid
    z_layered::Vector{Float64}
end

"Initialize SBM soil model hydraulic conductivity depth profile"
function sbm_kv_profiles(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
    kv_0::Vector{Float64},
    f::Vector{Float64},
    maxlayers::Int,
    nlayers::Vector{Int},
    sumlayers::Vector,
    dt::Second,
)
    kv_profile_type = config.model.saturated_hydraulic_conductivity_profile
    n_land_cells = length(land_indices_2d)
    if kv_profile_type == VerticalConductivityProfile.exponential
        kv_profile = KvExponential(kv_0, f)
    elseif kv_profile_type == VerticalConductivityProfile.exponential_constant
        z_exp = ncread(
            dataset,
            config,
            "soil_exponential_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth",
            LandHydrologySBM;
            sel=land_indices_2d,
        )
        exp_profile = KvExponential(kv_0, f)
        kv_profile = KvExponentialConstant(exp_profile, z_exp)
    elseif kv_profile_type == VerticalConductivityProfile.layered ||
           kv_profile_type == VerticalConductivityProfile.layered_exponential
        kv =
            ncread(
                dataset,
                config,
                "soil_layer_water__vertical_saturated_hydraulic_conductivity",
                LandHydrologySBM;
                sel=land_indices_2d,
            ) .* (dt / BASETIMESTEP)
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
                sel=land_indices_2d,
            )
            nlayers_kv = zeros(Int, n_land_cells)
            for land_cell_idx in 1:n_land_cells
                layers = @view sumlayers[land_cell_idx][2:nlayers[land_cell_idx]]
                _, n_soil_layers = findmin(abs.(z_layered[land_cell_idx] .- layers))
                nlayers_kv[land_cell_idx] = n_soil_layers
                z_layered[land_cell_idx] = layers[n_soil_layers]
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
    land_indices_2d::Vector{CartesianIndex{2}},
    dt::Second,
)
    config_soil_layer_thickness = config.model.soil_layer__thickness

    soil_layer_thickness =
        SVector(Tuple(push!(Float64.(config_soil_layer_thickness), MISSING_VALUE)))
    cum_depth_layers = pushfirst(cumsum(soil_layer_thickness), 0.0)
    maxlayers = length(soil_layer_thickness) # max number of soil layers

    @info "Using `$(maxlayers - 1)` soil layers with the following thickness: `$config_soil_layer_thickness`"

    w_soil =
        ncread(
            dataset,
            config,
            "soil_surface_temperature__weight_coefficient",
            LandHydrologySBM;
            sel=land_indices_2d,
        ) .* (dt / BASETIMESTEP)
    cf_soil = ncread(
        dataset,
        config,
        "soil_surface_water__infiltration_reduction_parameter",
        LandHydrologySBM;
        sel=land_indices_2d,
    )

    # soil parameters
    theta_s = ncread(
        dataset,
        config,
        "soil_water__saturated_volume_fraction",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    theta_r = ncread(
        dataset,
        config,
        "soil_water__residual_volume_fraction",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    kv_0 =
        ncread(
            dataset,
            config,
            "soil_surface_water__vertical_saturated_hydraulic_conductivity",
            LandHydrologySBM;
            sel=land_indices_2d,
        ) .* (dt / BASETIMESTEP)
    f = ncread(
        dataset,
        config,
        "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    hb = ncread(
        dataset,
        config,
        "soil_water__air_entry_pressure_head",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    h1 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h1",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    h2 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h2",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    h3_high = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h3_high",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    h3_low = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h3_low",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    h4 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h4",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    alpha_h1 = ncread(
        dataset,
        config,
        "vegetation_root__feddes_critical_pressure_head_h1_reduction_coefficient",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    soilthickness =
        ncread(dataset, config, "soil__thickness", LandHydrologySBM; sel=land_indices_2d)
    infiltcappath =
        ncread(
            dataset,
            config,
            "compacted_soil_surface_water__infiltration_capacity",
            LandHydrologySBM;
            sel=land_indices_2d,
        ) .* (dt / BASETIMESTEP)
    maxleakage =
        ncread(
            dataset,
            config,
            "soil_water_saturated_zone_bottom__max_leakage_volume_flux",
            LandHydrologySBM;
            sel=land_indices_2d,
        ) .* (dt / BASETIMESTEP)
    c = ncread(
        dataset,
        config,
        "soil_layer_water__brooks_corey_exponent",
        LandHydrologySBM;
        sel=land_indices_2d,
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
        sel=land_indices_2d,
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
        sel=land_indices_2d,
    )

    # vegetation parameters
    rootdistpar = ncread(
        dataset,
        config,
        "soil_wet_root__sigmoid_function_shape_parameter",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    cap_hmax = ncread(
        dataset,
        config,
        "soil_water_saturated_zone_top__capillary_rise_max_water_table_depth",
        LandHydrologySBM;
        sel=land_indices_2d,
    )
    cap_n = ncread(
        dataset,
        config,
        "soil_water_saturated_zone_top__capillary_rise_averianov_exponent",
        LandHydrologySBM;
        sel=land_indices_2d,
    )

    act_thickl =
        set_layerthickness.(soilthickness, (cum_depth_layers,), (soil_layer_thickness,))
    sumlayers = @. pushfirst(cumsum(act_thickl), 0.0)
    nlayers = number_of_active_layers.(act_thickl)

    if haskey(config.input.static, "soil_water__field_capacity_volume_fraction")
        theta_fc = ncread(
            dataset,
            config,
            "soil_water__field_capacity_volume_fraction";
            optional=false,
            sel=land_indices_2d,
            type=Float64,
        )
    else
        theta_fc = field_capacity.(act_thickl, nlayers, theta_s, theta_r, c, hb)
    end

    n_land_cells = length(land_indices_2d)

    # optional root fraction
    rootfraction_name = "soil_root__length_density_fraction"
    if haskey(config.input.static, rootfraction_name)
        rootfraction = ncread(
            dataset,
            config,
            rootfraction_name,
            LandHydrologySBM;
            sel=land_indices_2d,
        )
    else
        (; rootingdepth) = vegetation_parameter_set
        # default root fraction
        rootfraction = zeros(maxlayers, n_land_cells)
        for land_cell_idx in 1:n_land_cells
            if rootingdepth[land_cell_idx] > 0.0
                for soil_layer_idx in 1:maxlayers
                    if (
                        rootingdepth[land_cell_idx] -
                        sumlayers[land_cell_idx][soil_layer_idx]
                    ) >= act_thickl[land_cell_idx][soil_layer_idx]
                        rootfraction[soil_layer_idx, land_cell_idx] =
                            act_thickl[land_cell_idx][soil_layer_idx] /
                            rootingdepth[land_cell_idx]
                    else
                        rootfraction[soil_layer_idx, land_cell_idx] =
                            max(
                                rootingdepth[land_cell_idx] -
                                sumlayers[land_cell_idx][soil_layer_idx],
                                0.0,
                            ) / rootingdepth[land_cell_idx]
                    end
                end
            end
        end
    end

    kv_profile = sbm_kv_profiles(
        dataset,
        config,
        land_indices_2d,
        kv_0,
        f,
        maxlayers,
        nlayers,
        sumlayers,
        dt,
    )

    soilwatercapacity = @. soilthickness * (theta_s - theta_r)

    sbm_params = SbmSoilParameters(;
        maxlayers,
        nlayers,
        soilwatercapacity,
        theta_s,
        theta_r,
        theta_fc,
        kvfrac=svectorscopy(kvfrac, Val{maxlayers}()),
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
        rootfraction=svectorscopy(rootfraction, Val{maxlayers}()),
        cap_hmax,
        cap_n,
        c,
        w_soil,
        cf_soil,
        soil_fraction=fill(MISSING_VALUE, n_land_cells),
        kv_profile,
        vegetation_parameter_set,
    )
    return sbm_params
end

"SBM soil model"
@with_kw struct SbmSoilModel{N,M,Kv} <: AbstractSoilModel
    n_land_cells::Int
    boundary_conditions::SbmSoilBC = SbmSoilBC(; n_land_cells)
    parameters::SbmSoilParameters{N,M,Kv}
    variables::SbmSoilVariables{N}
end

"Initialize SBM soil model"
function SbmSoilModel(
    dataset::NCDataset,
    config::Config,
    vegetation_parameter_set::VegetationParameters,
    land_indices_2d::Vector{CartesianIndex{2}},
    dt::Second,
)
    n_land_cells = length(land_indices_2d)
    parameters =
        SbmSoilParameters(dataset, config, vegetation_parameter_set, land_indices_2d, dt)
    variables = SbmSoilVariables(n_land_cells, parameters)
    soil_model = SbmSoilModel(; n_land_cells, parameters, variables)
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
)
    (; interception, runoff, demand, allocation) = external_models
    (; potential_transpiration, water_flux_surface, potential_soilevaporation) =
        soil_model.boundary_conditions

    potential_transpiration .= get_potential_transpiration(interception)

    @. potential_soilevaporation =
        soil_model.parameters.soil_fraction * atmospheric_forcing.potential_evaporation
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
    (; n_land_cells) = soil_model
    v = soil_model.variables
    p = soil_model.parameters
    for land_cell_idx in 1:n_land_cells
        v.ustoredepth[land_cell_idx] =
            sum(@view v.ustorelayerdepth[land_cell_idx][1:p.nlayers[land_cell_idx]])
    end
    return nothing
end

"Update the infiltration reduction factor of the SBM soil model for a single timestep"
function infiltration_reduction_factor!(
    soil_model::SbmSoilModel;
    modelsnow=false,
    soil_infiltration_reduction=false,
)
    (; n_land_cells) = soil_model
    v = soil_model.variables
    p = soil_model.parameters

    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        v.f_infiltration_reduction[land_cell_idx] = infiltration_reduction_factor(
            v.tsoil[land_cell_idx],
            p.cf_soil[land_cell_idx];
            modelsnow,
            soil_infiltration_reduction,
        )
    end
    return nothing
end

"""
    infiltration!(soil_model::SbmSoilMsoil

Update the infiltration rate `infiltsoilpath` and infiltration excess water rate
`infiltexcess` of the SBM soil model for a single timestep.
"""
function infiltration!(soil_model::SbmSoilModel)
    (; n_land_cells) = soil_model

    v = soil_model.variables
    p = soil_model.parameters
    (; water_flux_surface) = soil_model.boundary_conditions
    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        v.infiltsoilpath[land_cell_idx], v.infiltexcess[land_cell_idx] = infiltration(
            water_flux_surface[land_cell_idx],
            p.pathfrac[land_cell_idx],
            p.infiltcapsoil[land_cell_idx],
            p.infiltcappath[land_cell_idx],
            v.ustorecapacity[land_cell_idx],
            v.f_infiltration_reduction[land_cell_idx],
        )
    end
    return nothing
end

"""
    unsaturated_zone_flow!(soil_model::SbmSoilModel)

Update unsaturated storage `ustorelayerdepth` and the `transfer` of water from the unsaturated
to the saturated store of the SBM soil model for a single timestep, based on the Brooks-Corey
approach.
"""
function unsaturated_zone_flow!(soil_model::SbmSoilModel)
    (; n_land_cells) = soil_model
    v = soil_model.variables
    p = soil_model.parameters

    threaded_foreach(1:n_land_cells; basesize=250) do land_cell_idx
        if v.n_unsatlayers[land_cell_idx] > 0
            # Brooks-Corey approach
            z = cumsum(v.ustorelayerthickness[land_cell_idx])
            flow_rate = 0.0
            for soil_layer_idx in 1:v.n_unsatlayers[land_cell_idx]
                l_sat =
                    v.ustorelayerthickness[land_cell_idx][soil_layer_idx] *
                    (p.theta_s[land_cell_idx] - p.theta_r[land_cell_idx])
                kv_z = hydraulic_conductivity_at_depth(
                    p.kv_profile,
                    p.kvfrac,
                    z[soil_layer_idx],
                    land_cell_idx,
                    soil_layer_idx,
                )
                ustorelayerdepth = if isone(soil_layer_idx)
                    v.ustorelayerdepth[land_cell_idx][soil_layer_idx] +
                    v.infiltsoilpath[land_cell_idx]
                else
                    v.ustorelayerdepth[land_cell_idx][soil_layer_idx] + flow_rate
                end
                ustorelayerdepth, flow_rate = unsatzone_flow_layer(
                    ustorelayerdepth,
                    kv_z,
                    l_sat,
                    p.c[land_cell_idx][soil_layer_idx],
                )
                v.ustorelayerdepth[land_cell_idx] = setindex(
                    v.ustorelayerdepth[land_cell_idx],
                    ustorelayerdepth,
                    soil_layer_idx,
                )
            end
            v.transfer[land_cell_idx] = flow_rate
        else
            v.transfer[land_cell_idx] = 0.0
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
function soil_evaporation!(soil_model::SbmSoilModel)
    (; n_land_cells) = soil_model
    (; potential_soilevaporation) = soil_model.boundary_conditions
    v = soil_model.variables
    p = soil_model.parameters

    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        potsoilevap = potential_soilevaporation[land_cell_idx]
        # First calculate the evaporation of unsaturated storage into the
        # atmosphere from the upper layer.
        soilevapunsat = soil_evaporation_unsatured_store(
            potsoilevap,
            v.ustorelayerdepth[land_cell_idx][1],
            v.ustorelayerthickness[land_cell_idx][1],
            v.n_unsatlayers[land_cell_idx],
            v.zi[land_cell_idx],
            p.theta_s[land_cell_idx] - p.theta_r[land_cell_idx],
        )
        # Ensure that the unsaturated evaporation rate does not exceed the
        # available unsaturated moisture
        soilevapunsat = min(soilevapunsat, v.ustorelayerdepth[land_cell_idx][1])
        # Update the additional atmospheric demand
        potsoilevap -= soilevapunsat
        v.ustorelayerdepth[land_cell_idx] = setindex(
            v.ustorelayerdepth[land_cell_idx],
            v.ustorelayerdepth[land_cell_idx][1] - soilevapunsat,
            1,
        )

        theta_drainable = lower_bound_drainable_porosity(
            p.theta_s[land_cell_idx],
            p.theta_fc[land_cell_idx],
        )
        soilevapsat = soil_evaporation_satured_store(
            potsoilevap,
            v.n_unsatlayers[land_cell_idx],
            p.act_thickl[land_cell_idx][1],
            v.zi[land_cell_idx],
            theta_drainable,
        )

        v.soilevapsat[land_cell_idx] = soilevapsat
        v.soilevap[land_cell_idx] = soilevapunsat + soilevapsat
        v.drainable_waterdepth[land_cell_idx] =
            v.drainable_waterdepth[land_cell_idx] - soilevapsat
    end
    return nothing
end

"""
    transpiration!(soil_model::SbmSoilModel, dt)

Update total `transpiration`, transpiration from the unsaturated store `ae_ustore` and
saturated store `actevapsat` of the SBM soil model for a single timestep. Also unsaturated
storage `ustorelayerdepth` and the saturated store `satwaterdepth` are updated.
"""
function transpiration!(soil_model::SbmSoilModel, dt::Float64)
    (; n_land_cells) = soil_model
    (; potential_transpiration) = soil_model.boundary_conditions
    v = soil_model.variables
    p = soil_model.parameters

    rootingdepth = get_rootingdepth(soil_model)

    threaded_foreach(1:n_land_cells; basesize=250) do land_cell_idx
        v.h3[land_cell_idx] = feddes_h3(
            p.h3_high[land_cell_idx],
            p.h3_low[land_cell_idx],
            potential_transpiration[land_cell_idx],
            dt,
        )

        # compute sum of root fraction in unsaturated soil layers and adapt root fraction
        # lowest unsaturated soil layer if water table depth intersects the unsaturated root
        # zone
        sum_rootfraction_unsat = 0.0
        rootfraction_unsat_lowest = 0.0
        for soil_layer_idx in 1:v.n_unsatlayers[land_cell_idx]
            # the root fraction is valid for the root length in a soil layer, if zi decreases
            # the root length the root fraction needs to be adapted
            if soil_layer_idx == v.n_unsatlayers[land_cell_idx] &&
               v.zi[land_cell_idx] < rootingdepth[land_cell_idx]
                rootlength = min(
                    p.act_thickl[land_cell_idx][soil_layer_idx],
                    rootingdepth[land_cell_idx] -
                    p.sumlayers[land_cell_idx][soil_layer_idx],
                )
                rootfraction_unsat =
                    p.rootfraction[land_cell_idx][soil_layer_idx] *
                    (v.ustorelayerthickness[land_cell_idx][soil_layer_idx] / rootlength)
                sum_rootfraction_unsat += rootfraction_unsat
            else
                rootfraction_unsat = p.rootfraction[land_cell_idx][soil_layer_idx]
                sum_rootfraction_unsat += rootfraction_unsat
            end
            # rootfraction lowest unsaturated layer
            rootfraction_unsat_lowest = rootfraction_unsat
        end

        actevapustore = 0.0
        for soil_layer_idx in 1:v.n_unsatlayers[land_cell_idx]
            # scale rootfraction soil layer unsaturated zone based on sum of rootfraction in
            # unsaturated zone
            if soil_layer_idx < v.n_unsatlayers[land_cell_idx]
                rootfraction_unsat = p.rootfraction[land_cell_idx][soil_layer_idx]
            else
                rootfraction_unsat = rootfraction_unsat_lowest
            end
            rootfraction_unsat_scaled =
                rootingdepth[land_cell_idx] > 0.0 ?
                max((1.0 / sum_rootfraction_unsat), 1.0) * rootfraction_unsat : 0.0

            vwc = max(
                v.ustorelayerdepth[land_cell_idx][soil_layer_idx] /
                v.ustorelayerthickness[land_cell_idx][soil_layer_idx],
                Float64(0.0000001),
            )
            head = head_brooks_corey(
                vwc,
                p.theta_s[land_cell_idx],
                p.theta_r[land_cell_idx],
                p.c[land_cell_idx][soil_layer_idx],
                p.hb[land_cell_idx],
            )
            alpha = rwu_reduction_feddes(
                head,
                p.h1[land_cell_idx],
                p.h2[land_cell_idx],
                v.h3[land_cell_idx],
                p.h4[land_cell_idx],
                p.alpha_h1[land_cell_idx],
            )

            availcap = min(
                1.0,
                max(
                    0.0,
                    (
                        rootingdepth[land_cell_idx] -
                        p.sumlayers[land_cell_idx][soil_layer_idx]
                    ) / v.ustorelayerthickness[land_cell_idx][soil_layer_idx],
                ),
            )
            maxextr = v.ustorelayerdepth[land_cell_idx][soil_layer_idx] * availcap
            actevapustore_layer = min(
                alpha * rootfraction_unsat_scaled * potential_transpiration[land_cell_idx],
                maxextr,
            )
            ustorelayerdepth =
                v.ustorelayerdepth[land_cell_idx][soil_layer_idx] - actevapustore_layer
            actevapustore += actevapustore_layer
            v.ustorelayerdepth[land_cell_idx] = setindex(
                v.ustorelayerdepth[land_cell_idx],
                ustorelayerdepth,
                soil_layer_idx,
            )
        end

        # transpiration from saturated store
        wetroots = scurve(
            v.zi[land_cell_idx],
            rootingdepth[land_cell_idx],
            Float64(1.0),
            p.rootdistpar[land_cell_idx],
        )
        alpha = rwu_reduction_feddes(
            Float64(0.0),
            p.h1[land_cell_idx],
            p.h2[land_cell_idx],
            v.h3[land_cell_idx],
            p.h4[land_cell_idx],
            p.alpha_h1[land_cell_idx],
        )
        restpottrans = potential_transpiration[land_cell_idx] - actevapustore
        actevapsat =
            min(restpottrans * wetroots * alpha, v.drainable_waterdepth[land_cell_idx])

        v.ae_ustore[land_cell_idx] = actevapustore
        v.actevapsat[land_cell_idx] = actevapsat
        v.drainable_waterdepth[land_cell_idx] =
            v.drainable_waterdepth[land_cell_idx] - actevapsat
        v.transpiration[land_cell_idx] = actevapustore + actevapsat
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
function actual_infiltration!(soil_model::SbmSoilModel)
    (; n_land_cells) = soil_model
    v = soil_model.variables
    p = soil_model.parameters

    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        # check soil moisture balance per layer
        ustoredepth_excess = 0.0
        for soil_layer_idx in v.n_unsatlayers[land_cell_idx]:-1:1
            ustoredepth_excess = max(
                0.0,
                v.ustorelayerdepth[land_cell_idx][soil_layer_idx] -
                v.ustorelayerthickness[land_cell_idx][soil_layer_idx] *
                (p.theta_s[land_cell_idx] - p.theta_r[land_cell_idx]),
            )
            v.ustorelayerdepth[land_cell_idx] = setindex(
                v.ustorelayerdepth[land_cell_idx],
                v.ustorelayerdepth[land_cell_idx][soil_layer_idx] - ustoredepth_excess,
                soil_layer_idx,
            )
            if soil_layer_idx > 1
                v.ustorelayerdepth[land_cell_idx] = setindex(
                    v.ustorelayerdepth[land_cell_idx],
                    v.ustorelayerdepth[land_cell_idx][soil_layer_idx-1] +
                    ustoredepth_excess,
                    soil_layer_idx - 1,
                )
            end
        end

        v.actinfilt[land_cell_idx] = v.infiltsoilpath[land_cell_idx] - ustoredepth_excess
    end
    return nothing
end

"""
    actual_infiltration_soil_path!(soil_model::SbmSoilModel)

Update the actual infiltration rate for soil `actinfiltsoil` and paved area `actinfiltpath`
of the SBM soil model for a single timestep.
"""
function actual_infiltration_soil_path!(soil_model::SbmSoilModel)
    (; n_land_cells) = soil_model
    v = soil_model.variables
    p = soil_model.parameters
    (; water_flux_surface) = soil_model.boundary_conditions

    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        v.actinfiltsoil[land_cell_idx], v.actinfiltpath[land_cell_idx] =
            actual_infiltration_soil_path(
                water_flux_surface[land_cell_idx],
                v.actinfilt[land_cell_idx],
                p.pathfrac[land_cell_idx],
                p.infiltcapsoil[land_cell_idx],
                p.infiltcappath[land_cell_idx],
                v.f_infiltration_reduction[land_cell_idx],
            )
    end
    return nothing
end

"""
    capillary_flux!(soil_model::SbmSoilModel)

Update the capillary flux `actcapflux` of the SBM soil model for a single timestep.
"""
function capillary_flux!(soil_model::SbmSoilModel)
    (; n_land_cells) = soil_model
    v = soil_model.variables
    p = soil_model.parameters
    rootingdepth = get_rootingdepth(soil_model)

    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        if v.n_unsatlayers[land_cell_idx] > 0
            ksat = hydraulic_conductivity_at_depth(
                p.kv_profile,
                p.kvfrac,
                v.zi[land_cell_idx],
                land_cell_idx,
                v.n_unsatlayers[land_cell_idx],
            )
            maxcapflux = max(
                0.0,
                min(
                    ksat,
                    v.ae_ustore[land_cell_idx],
                    v.ustorecapacity[land_cell_idx],
                    v.drainable_waterdepth[land_cell_idx],
                ),
            )

            if v.zi[land_cell_idx] > rootingdepth[land_cell_idx]
                capflux =
                    maxcapflux * pow(
                        1.0 -
                        min(v.zi[land_cell_idx], p.cap_hmax[land_cell_idx]) /
                        (p.cap_hmax[land_cell_idx]),
                        p.cap_n[land_cell_idx],
                    )
            else
                capflux = 0.0
            end

            netcapflux = capflux
            actcapflux = 0.0
            for soil_layer_idx in v.n_unsatlayers[land_cell_idx]:-1:1
                toadd = min(
                    netcapflux,
                    max(
                        v.ustorelayerthickness[land_cell_idx][soil_layer_idx] *
                        (p.theta_s[land_cell_idx] - p.theta_r[land_cell_idx]) -
                        v.ustorelayerdepth[land_cell_idx][soil_layer_idx],
                        0.0,
                    ),
                )
                v.ustorelayerdepth[land_cell_idx] = setindex(
                    v.ustorelayerdepth[land_cell_idx],
                    v.ustorelayerdepth[land_cell_idx][soil_layer_idx] + toadd,
                    soil_layer_idx,
                )
                netcapflux -= toadd
                actcapflux += toadd
            end
            v.actcapflux[land_cell_idx] = actcapflux
        else
            v.actcapflux[land_cell_idx] = 0.0
        end
    end
    return nothing
end

"""
    leakage!(soil_model::SbmSoilModel)

Update the actual leakage rate `actleakage` of the SBM soil model for a single timestep.
"""
function leakage!(soil_model::SbmSoilModel)
    (; n_land_cells) = soil_model
    v = soil_model.variables
    p = soil_model.parameters

    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        deepksat = hydraulic_conductivity_at_depth(
            p.kv_profile,
            p.kvfrac,
            p.soilthickness[land_cell_idx],
            land_cell_idx,
            p.nlayers[land_cell_idx],
        )
        deeptransfer = min(v.drainable_waterdepth[land_cell_idx], deepksat)
        v.actleakage[land_cell_idx] =
            max(0.0, min(p.maxleakage[land_cell_idx], deeptransfer))
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
        modelsnow=config.model.snow__flag,
        soil_infiltration_reduction=config.model.soil_infiltration_reduction__flag,
    )
    infiltration!(soil_model)
    # unsaturated zone flow
    unsaturated_zone_flow!(soil_model)
    # soil evaporation and transpiration
    soil_evaporation!(soil_model)
    transpiration!(soil_model, dt)
    # actual infiltration and excess water
    actual_infiltration!(soil_model)
    @. v.excesswater = water_flux_surface - v.actinfilt - v.infiltexcess
    actual_infiltration_soil_path!(soil_model)
    @. v.excesswatersoil =
        max(water_flux_surface * (1.0 - p.pathfrac) - v.actinfiltsoil, 0.0)
    @. v.excesswaterpath = max(water_flux_surface * p.pathfrac - v.actinfiltpath, 0.0)
    # recompute the unsaturated store and ustorecapacity (for capillary flux)
    ustoredepth!(soil_model)
    @. v.ustorecapacity = p.soilwatercapacity - v.satwaterdepth - v.ustoredepth
    # capillary flux and leakage
    capillary_flux!(soil_model)
    leakage!(soil_model)
    # recharge rate to the saturated store
    @. v.recharge =
        (v.transfer - v.actcapflux - v.actleakage - v.actevapsat - v.soilevapsat)
    # total actual evapotranspiration
    v.actevap .=
        v.soilevap .+ v.transpiration .+ runoff.variables.ae_openw_r .+
        runoff.variables.ae_openw_l .+ get_evaporation(demand.paddy)
    return nothing
end

function update_ustorelayerdepth!(soil, zi_prev, zi, land_cell_idx)
    v = soil.variables
    p = soil.parameters

    n_unsatlayers_prev = v.n_unsatlayers[land_cell_idx]
    ustorelayerthickness_prev = v.ustorelayerthickness[land_cell_idx]
    ustorelayerthickness =
        set_layerthickness(zi, p.sumlayers[land_cell_idx], p.act_thickl[land_cell_idx])
    n_unsatlayers = number_of_active_layers(ustorelayerthickness)
    ustorelayerdepth = v.ustorelayerdepth[land_cell_idx]
    if zi < zi_prev
        for soil_layer_idx in n_unsatlayers:n_unsatlayers_prev
            iszero(soil_layer_idx) && continue
            if isnan(ustorelayerthickness[soil_layer_idx])
                ustorelayerdepth = setindex(ustorelayerdepth, 0.0, soil_layer_idx)
            else
                f =
                    ustorelayerthickness[soil_layer_idx] /
                    ustorelayerthickness_prev[soil_layer_idx]
                ustorelayerdepth = setindex(
                    ustorelayerdepth,
                    f * ustorelayerdepth[soil_layer_idx],
                    soil_layer_idx,
                )
            end
        end
    else
        for soil_layer_idx in n_unsatlayers_prev:n_unsatlayers
            iszero(soil_layer_idx) && continue
            thickness_prev =
                isnan(ustorelayerthickness_prev[soil_layer_idx]) ? 0.0 :
                ustorelayerthickness_prev[soil_layer_idx]
            delta_thickness = ustorelayerthickness[soil_layer_idx] - thickness_prev
            ustorelayerdepth = setindex(
                ustorelayerdepth,
                ustorelayerdepth[soil_layer_idx] +
                delta_thickness * (p.theta_fc[land_cell_idx] - p.theta_r[land_cell_idx]),
                soil_layer_idx,
            )
        end
    end
    v.n_unsatlayers[land_cell_idx] = n_unsatlayers
    v.ustorelayerdepth[land_cell_idx] = ustorelayerdepth
    v.ustorelayerthickness[land_cell_idx] = ustorelayerthickness
    v.zi[land_cell_idx] = zi
end

"""
    update_ustorelayerdepth!(soil_model::SbmSoilModel, subsurface_flow)

Update the `SbmSoilModel` variables unsaturated store depth of soil layers
`ustorelayerdepth`, number of unsaturated zone soil layers `n_unsatlayers`, thickness of
unsaturated zone soil layers `ustorelayerthickness` and water table depth `zi`, based on the
water table change computed by a subsurface flow model.
"""
function update_ustorelayerdepth!(soil_model::SbmSoilModel, subsurface_flow)
    (; n_land_cells) = soil_model
    p = soil_model.parameters
    v = soil_model.variables

    zi = get_water_depth(subsurface_flow) * 1000.0 # convert from [m] to [mm]

    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        zi_prev = v.zi[land_cell_idx]
        update_ustorelayerdepth!(soil_model, zi_prev, zi[land_cell_idx], land_cell_idx)
    end
end

"""
    update_soil_water_storage!(soil_model::SbmSoilModel, external_models::NamedTuple)

Update the SBM soil model for a single timestep based on the update of a subsurface flow
model, resulting in a change in water table depth and an exfiltration rate `exfiltwater`.

The available water in unsaturated zone `ustoredepth`, unsaturated store capacity
`ustorecapacity`, `total_soilwater_storage`, land `runoff` and `net_runoff`, the saturated
store `satwaterdepth` and the water exfiltrating during saturation excess conditions
`exfiltsatwater` are updated. Additionally, volumetric water content per soil layer and for
the root zone are updated.
"""
function update_soil_water_storage!(soil_model::SbmSoilModel, external_models::NamedTuple)
    (; n_land_cells) = soil_model
    (; runoff, demand, subsurface_flow) = external_models
    (; runoff_land, ae_openw_l) = runoff.variables
    p = soil_model.parameters
    v = soil_model.variables

    exfiltsatwater = subsurface_flow.variables.exfiltwater * 1000.0 # convert from [m] to [mm]
    rootingdepth = get_rootingdepth(soil_model)

    threaded_foreach(1:n_land_cells; basesize=1000) do land_cell_idx
        ustorelayerdepth = v.ustorelayerdepth[land_cell_idx]
        ustorelayerthickness = v.ustorelayerthickness[land_cell_idx]
        ustoredepth = sum(@view ustorelayerdepth[1:(v.n_unsatlayers[land_cell_idx])])
        sbm_runoff = max(
            0.0,
            exfiltsatwater[land_cell_idx] +
            v.excesswater[land_cell_idx] +
            runoff_land[land_cell_idx] +
            v.infiltexcess[land_cell_idx],
        )

        # volumetric water content per soil layer and root zone
        vwc = v.vwc[land_cell_idx]
        vwc_perc = v.vwc_perc[land_cell_idx]
        for soil_layer_idx in 1:p.nlayers[land_cell_idx]
            if soil_layer_idx <= v.n_unsatlayers[land_cell_idx]
                vwc = setindex(
                    vwc,
                    (
                        ustorelayerdepth[soil_layer_idx] +
                        (
                            p.act_thickl[land_cell_idx][soil_layer_idx] -
                            ustorelayerthickness[soil_layer_idx]
                        ) * (p.theta_s[land_cell_idx] - p.theta_r[land_cell_idx])
                    ) / p.act_thickl[land_cell_idx][soil_layer_idx] +
                    p.theta_r[land_cell_idx],
                    soil_layer_idx,
                )
            else
                vwc = setindex(vwc, p.theta_s[land_cell_idx], soil_layer_idx)
            end
            vwc_perc = setindex(
                vwc_perc,
                (vwc[soil_layer_idx] / p.theta_s[land_cell_idx]) * 100.0,
                soil_layer_idx,
            )
        end

        rootstore_unsat = 0
        for soil_layer_idx in 1:(v.n_unsatlayers[land_cell_idx])
            rootstore_unsat =
                rootstore_unsat +
                min(
                    1.0,
                    (
                        max(
                            0.0,
                            rootingdepth[land_cell_idx] -
                            p.sumlayers[land_cell_idx][soil_layer_idx],
                        ) / ustorelayerthickness[soil_layer_idx]
                    ),
                ) * ustorelayerdepth[soil_layer_idx]
        end

        rootstore_sat =
            max(0.0, rootingdepth[land_cell_idx] - v.zi[land_cell_idx]) *
            (p.theta_s[land_cell_idx] - p.theta_r[land_cell_idx])
        rootstore = rootstore_sat + rootstore_unsat
        vwc_root = rootstore / rootingdepth[land_cell_idx] + p.theta_r[land_cell_idx]
        vwc_percroot = (vwc_root / p.theta_s[land_cell_idx]) * 100.0

        satwaterdepth =
            (p.soilthickness[land_cell_idx] - v.zi[land_cell_idx]) *
            (p.theta_s[land_cell_idx] - p.theta_r[land_cell_idx])
        drainable_waterdepth =
            (p.soilthickness[land_cell_idx] - v.zi[land_cell_idx]) *
            lower_bound_drainable_porosity(
                p.theta_s[land_cell_idx],
                p.theta_fc[land_cell_idx],
            )
        ustorecapacity = p.soilwatercapacity[land_cell_idx] - satwaterdepth - ustoredepth

        # update the outputs and states
        v.ustorecapacity[land_cell_idx] = ustorecapacity
        v.ustoredepth[land_cell_idx] = ustoredepth
        v.satwaterdepth[land_cell_idx] = satwaterdepth
        v.drainable_waterdepth[land_cell_idx] = drainable_waterdepth
        v.exfiltsatwater[land_cell_idx] = exfiltsatwater[land_cell_idx]
        v.runoff[land_cell_idx] = sbm_runoff
        v.vwc[land_cell_idx] = vwc
        v.vwc_perc[land_cell_idx] = vwc_perc
        v.rootstore[land_cell_idx] = rootstore
        v.vwc_root[land_cell_idx] = vwc_root
        v.vwc_percroot[land_cell_idx] = vwc_percroot
        v.total_soilwater_storage[land_cell_idx] = satwaterdepth + ustoredepth
    end
    # update runoff and net_runoff (the runoff rate depends on the presence of paddy fields
    # and the h_max parameter of a paddy field)
    update_runoff!(demand.paddy, v.runoff)
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
