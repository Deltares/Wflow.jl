abstract type AbstractSoilModel end

"Struct for storing SBM soil model variables"
@with_kw struct SbmSoilVariables{N}
    # Calculated soil water pressure head h3 of the root water uptake reduction function (Feddes) [cm]
    h3::Vector{Float}
    # Unsaturated store capacity [mm]
    ustorecapacity::Vector{Float}
    # Amount of water in the unsaturated store, per layer [mm]
    ustorelayerdepth::Vector{SVector{N, Float}}
    # Thickness of unsaturated zone, per layer [mm]
    ustorelayerthickness::Vector{SVector{N, Float}}
    # Saturated store [mm]
    satwaterdepth::Vector{Float}
    # Pseudo-water table depth [mm] (top of the saturated zone)
    zi::Vector{Float}
    # Number of unsaturated soil layers
    n_unsatlayers::Vector{Int}
    # Transpiration [mm Δt⁻¹]
    transpiration::Vector{Float}
    # Actual evaporation from unsaturated store [mm Δt⁻¹]
    ae_ustore::Vector{Float}
    # Soil evaporation from unsaturated and saturated store [mm Δt⁻¹]
    soilevap::Vector{Float}
    # Soil evaporation from saturated store [mm Δt⁻¹]
    soilevapsat::Vector{Float}
    # Actual capillary rise [mm Δt⁻¹]
    actcapflux::Vector{Float}
    # Actual transpiration from saturated store [mm Δt⁻¹]
    actevapsat::Vector{Float}
    # Total actual evapotranspiration [mm Δt⁻¹]
    actevap::Vector{Float}
    # Actual infiltration into the unsaturated zone [mm Δt⁻¹]
    actinfilt::Vector{Float}
    # Actual infiltration non-compacted fraction [mm Δt⁻¹]
    actinfiltsoil::Vector{Float}
    # Actual infiltration compacted fraction [mm Δt⁻¹]
    actinfiltpath::Vector{Float}
    # Actual infiltration (compacted and the non-compacted areas) [mm Δt⁻¹]
    infiltsoilpath::Vector{Float}
    # Infiltration excess water [mm Δt⁻¹]
    infiltexcess::Vector{Float}
    # Water that cannot infiltrate due to saturated soil (saturation excess) [mm Δt⁻¹]
    excesswater::Vector{Float}
    # Water exfiltrating during saturation excess conditions [mm Δt⁻¹]
    exfiltsatwater::Vector{Float}
    # Water exfiltrating from unsaturated store because of change in water table [mm Δt⁻¹]
    exfiltustore::Vector{Float}
    # Excess water for non-compacted fraction [mm Δt⁻¹]
    excesswatersoil::Vector{Float}
    # Excess water for compacted fraction [mm Δt⁻¹]
    excesswaterpath::Vector{Float}
    # Total surface runoff from infiltration and saturation excess (excluding actual open water evaporation) [mm Δt⁻¹]
    runoff::Vector{Float}
    # Net surface runoff (surface runoff - actual open water evaporation) [mm Δt⁻¹]
    net_runoff::Vector{Float}
    # Volumetric water content [-] per soil layer (including theta_r and saturated zone)
    vwc::Vector{SVector{N, Float}}
    # Volumetric water content [%] per soil layer (including theta_r and saturated zone)
    vwc_perc::Vector{SVector{N, Float}}
    # Root water storage [mm] in unsaturated and saturated zone (excluding theta_r)
    rootstore::Vector{Float}
    # Volumetric water content [-] in root zone (including theta_r and saturated zone)
    vwc_root::Vector{Float}
    # Volumetric water content [%] in root zone (including theta_r and saturated zone)
    vwc_percroot::Vector{Float}
    # Amount of available water in the unsaturated zone [mm]
    ustoredepth::Vector{Float}
    # Downward flux from unsaturated to saturated zone [mm Δt⁻¹]
    transfer::Vector{Float}
    # Net recharge to saturated store [mm Δt⁻¹]
    recharge::Vector{Float}
    # Actual leakage from saturated store [mm Δt⁻¹]
    actleakage::Vector{Float}
    # Total water storage (excluding floodplain volume, lakes and reservoirs) [mm]
    total_storage::Vector{Float}
    # Top soil temperature [ᵒC]
    tsoil::Vector{Float}
    # Soil infiltration reduction factor (when soil is frozen) [-]
    f_infiltration_reduction::Vector{Float}
end

"Struct for storing SBM soil model parameters"
@with_kw struct SbmSoilParameters{N, M, Kv}
    # Maximum number of soil layers [-]
    maxlayers::Int
    # Number of soil layers [-]
    nlayers::Vector{Int}
    # Saturated water content (porosity) [-]
    theta_s::Vector{Float}
    # Residual water content [-]
    theta_r::Vector{Float}
    # Soilwater capacity [mm]
    soilwatercapacity::Vector{Float}
    # Muliplication factor [-] applied to kv_z (vertical flow)
    kvfrac::Vector{SVector{N, Float}}
    # Air entry pressure [cm] of soil (Brooks-Corey)
    hb::Vector{Float}
    # Soil thickness [mm]
    soilthickness::Vector{Float}
    # Thickness of soil layers [mm]
    act_thickl::Vector{SVector{N, Float}}
    # Cumulative sum of soil layers [mm], starting at soil surface (0)
    sumlayers::Vector{SVector{M, Float}}
    # Infiltration capacity of the compacted areas [mm Δt⁻¹]
    infiltcappath::Vector{Float}
    # Soil infiltration capacity [mm Δt⁻¹]
    infiltcapsoil::Vector{Float}
    # Maximum leakage [mm Δt⁻¹] from saturated zone
    maxleakage::Vector{Float}
    # Parameter [mm] controlling capillary rise
    cap_hmax::Vector{Float}
    # Coefficient [-] controlling capillary rise
    cap_n::Vector{Float}
    # Brooks-Corey power coefﬁcient [-] for each soil layer
    c::Vector{SVector{N, Float}}
    # Soil temperature smooth factor [-]
    w_soil::Vector{Float}
    # Controls soil infiltration reduction factor when soil is frozen [-]
    cf_soil::Vector{Float}
    # Fraction of compacted area  [-]
    pathfrac::Vector{Float}
    # Controls how roots are linked to water table [-]
    rootdistpar::Vector{Float}
    # Fraction of the root length density in each soil layer [-]
    rootfraction::Vector{SVector{N, Float}}
    # Soil water pressure head h1 of the root water uptake reduction function (Feddes) [cm]
    h1::Vector{Float}
    # Soil water pressure head h2 of the root water uptake reduction function (Feddes) [cm]
    h2::Vector{Float}
    # Soil water pressure head h3_high of the root water uptake reduction function (Feddes) [cm]
    h3_high::Vector{Float}
    # Soil water pressure head h3_low of the root water uptake reduction function (Feddes) [cm]
    h3_low::Vector{Float}
    # Soil water pressure head h4 of the root water uptake reduction function (Feddes) [cm]
    h4::Vector{Float}
    # Root water uptake reduction at soil water pressure head h1 (0.0 or 1.0) [-]
    alpha_h1::Vector{Float}
    # Soil fraction [-]
    soil_fraction::Vector{Float}
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
    zi = @. max(0.0, soilthickness - satwaterdepth / (theta_s - theta_r))
    ustorelayerthickness = set_layerthickness.(zi, sumlayers, act_thickl)
    n_unsatlayers = number_of_active_layers.(ustorelayerthickness)

    vwc = fill(MISSING_VALUE, maxlayers, n)
    vwc_perc = fill(MISSING_VALUE, maxlayers, n)

    vars = SbmSoilVariables(;
        ustorelayerdepth = zero(act_thickl),
        ustorecapacity = soilwatercapacity .- satwaterdepth,
        ustorelayerthickness,
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
        vwc = svectorscopy(vwc, Val{Int64(maxlayers)}()),
        vwc_perc = svectorscopy(vwc_perc, Val{Int64(maxlayers)}()),
        rootstore = fill(MISSING_VALUE, n),
        vwc_root = fill(MISSING_VALUE, n),
        vwc_percroot = fill(MISSING_VALUE, n),
        ustoredepth = fill(MISSING_VALUE, n),
        transfer = fill(MISSING_VALUE, n),
        recharge = fill(MISSING_VALUE, n),
        actleakage = fill(MISSING_VALUE, n),
        tsoil = fill(10.0, n),
        total_storage = zeros(Float, n), # Set the total water storage from initialized values
    )
    return vars
end

"Struct for storing SBM soil model boundary conditions"
@with_kw struct SbmSoilBC
    # Water flux at the soil surface [mm Δt⁻¹]
    water_flux_surface::Vector{Float}
    # Potential transpiration rate [mm Δt⁻¹]
    potential_transpiration::Vector{Float}
    # Potential soil evaporation rate [mm Δt⁻¹]
    potential_soilevaporation::Vector{Float}
end

"Initialize SBM soil model boundary conditions"
function SbmSoilBC(
    n::Int;
    water_flux_surface::Vector{Float} = fill(MISSING_VALUE, n),
    potential_transpiration::Vector{Float} = fill(MISSING_VALUE, n),
    potential_soilevaporation::Vector{Float} = fill(MISSING_VALUE, n),
)
    return SbmSoilBC(;
        water_flux_surface = water_flux_surface,
        potential_transpiration = potential_transpiration,
        potential_soilevaporation = potential_soilevaporation,
    )
end

"Exponential depth profile of vertical hydraulic conductivity at the soil surface"
struct KvExponential
    # Vertical hydraulic conductivity [mm Δt⁻¹] at soil surface
    kv_0::Vector{Float}
    # A scaling parameter [mm⁻¹] (controls exponential decline of kv_0)
    f::Vector{Float}
end

"Exponential constant depth profile of vertical hydraulic conductivity"
struct KvExponentialConstant
    exponential::KvExponential
    # Depth [mm] from soil surface for which exponential decline of kv_0 is valid
    z_exp::Vector{Float}
end

"Layered depth profile of vertical hydraulic conductivity"
struct KvLayered{N}
    # Vertical hydraulic conductivity [mm Δt⁻¹] per soil layer
    kv::Vector{SVector{N, Float}}
end

"Layered exponential depth profile of vertical hydraulic conductivity"
struct KvLayeredExponential{N}
    # A scaling parameter [mm⁻¹] (controls exponential decline of kv_0)
    f::Vector{Float}
    # Vertical hydraulic conductivity [mm Δt⁻¹] per soil layer
    kv::Vector{SVector{N, Float}}
    # Number of soil layers [-] with vertical hydraulic conductivity value `kv`
    nlayers_kv::Vector{Int}
    # Depth [mm] from soil surface for which layered profile is valid
    z_layered::Vector{Float}
end

"Initialize SBM soil model hydraulic conductivity depth profile"
function sbm_kv_profiles(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    kv_0::Vector{Float},
    f::Vector{Float},
    maxlayers::Int,
    nlayers::Vector{Int},
    sumlayers::Vector,
    dt::Second,
)
    kv_profile_type =
        get(config.model, "saturated_hydraulic_conductivity_profile", "exponential")::String
    n = length(indices)
    if kv_profile_type == "exponential"
        kv_profile = KvExponential(kv_0, f)
    elseif kv_profile_type == "exponential_constant"
        lens = lens_input_parameter(
            config,
            "soil_vertical_saturated_hydraulic_conductivity_profile~exponential_below-surface__depth";
            optional = false,
        )
        z_exp = ncread(dataset, config, lens; sel = indices, type = Float)
        exp_profile = KvExponential(kv_0, f)
        kv_profile = KvExponentialConstant(exp_profile, z_exp)
    elseif kv_profile_type == "layered" || kv_profile_type == "layered_exponential"
        lens = lens_input_parameter(
            config,
            "soil_layer_water__vertical_saturated_hydraulic_conductivity",
        )
        kv =
            ncread(dataset, config, lens; sel = indices, type = Float, dimname = :layer) .*
            (dt / BASETIMESTEP)
        if size(kv, 1) != maxlayers
            parname = lens(config)
            size1 = size(kv, 1)
            error("$parname needs a layer dimension of size $maxlayers, but is $size1")
        end
        if kv_profile_type == "layered"
            kv_profile = KvLayered(svectorscopy(kv, Val{Int64(maxlayers)}()))
        else
            lens = lens_input_parameter(
                config,
                "soil_vertical_saturated_hydraulic_conductivity_profile~layered_below-surface__depth";
                optional = false,
            )
            z_layered = ncread(dataset, config, lens; sel = indices, type = Float)
            nlayers_kv = fill(0, n)
            for i in eachindex(nlayers_kv)
                layers = @view sumlayers[i][2:nlayers[i]]
                _, k = findmin(abs.(z_layered[i] .- layers))
                nlayers_kv[i] = k
                z_layered[i] = layers[k]
            end
            kv_profile = KvLayeredExponential(
                f,
                svectorscopy(kv, Val{Int64(maxlayers)}()),
                nlayers_kv,
                z_layered,
            )
        end
    else
        error("""An unknown "saturated_hydraulic_conductivity_profile" is specified in the
              TOML file ($ksat_profile). This should be "exponential",
              "exponential_constant", "layered" or "layered_exponential".
              """)
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
    config_soil_layer_thickness = get(config.model, "soil_layer__thickness", Float[])

    if length(config_soil_layer_thickness) > 0
        soil_layer_thickness =
            SVector(Tuple(push!(Float.(config_soil_layer_thickness), MISSING_VALUE)))
        cum_depth_layers = pushfirst(cumsum(soil_layer_thickness), 0.0)
        maxlayers = Int(length(soil_layer_thickness)) # max number of soil layers
    else
        soil_layer_thickness = SVector.(soilthickness)
        cum_depth_layers = pushfirst(cumsum(soil_layer_thickness), 0.0)
        maxlayers = Int(1)
    end

    lens = lens_input_parameter(config, "soil_surface_temperature__weight_coefficient")
    w_soil =
        ncread(dataset, config, lens; sel = indices, defaults = 0.1125, type = Float) .*
        (dt / BASETIMESTEP)

    lens =
        lens_input_parameter(config, "soil_surface_water__infiltration_reduction_parameter")
    cf_soil = ncread(dataset, config, lens; sel = indices, defaults = 0.038, type = Float)

    # soil parameters
    lens = lens_input_parameter(config, "soil_water__saturated_volume_fraction")
    theta_s = ncread(dataset, config, lens; sel = indices, type = Float)

    lens = lens_input_parameter(config, "soil_water__residual_volume_fraction")
    theta_r = ncread(dataset, config, lens; sel = indices, type = Float)

    lens = lens_input_parameter(
        config,
        "soil_surface_water__vertical_saturated_hydraulic_conductivity",
    )
    kv_0 =
        ncread(dataset, config, lens; sel = indices, type = Float) .*
        Float(dt / BASETIMESTEP)

    lens = lens_input_parameter(
        config,
        "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter",
    )
    f = ncread(dataset, config, lens; sel = indices, type = Float)

    lens = lens_input_parameter(config, "soil_water__air_entry_pressure_head")
    hb = ncread(dataset, config, lens; sel = indices, defaults = -10.0, type = Float)

    lens = lens_input_parameter(config, "vegetation_root__feddes_critial_pressure_head_h~1")
    h1 = ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float)

    lens = lens_input_parameter(config, "vegetation_root__feddes_critial_pressure_head_h~2")
    h2 = ncread(dataset, config, lens; sel = indices, defaults = -100.0, type = Float)

    lens = lens_input_parameter(
        config,
        "vegetation_root__feddes_critial_pressure_head_h~3~high",
    )
    h3_high = ncread(dataset, config, lens; sel = indices, defaults = -400.0, type = Float)

    lens = lens_input_parameter(
        config,
        "vegetation_root__feddes_critial_pressure_head_h~3~low",
    )
    h3_low = ncread(dataset, config, lens; sel = indices, defaults = -1000.0, type = Float)

    lens = lens_input_parameter(config, "vegetation_root__feddes_critial_pressure_head_h~4")
    h4 = ncread(dataset, config, lens; sel = indices, defaults = -16000.0, type = Float)

    lens = lens_input_parameter(
        config,
        "vegetation_root__feddes_critial_pressure_head_h~1_reduction_coefficient",
    )
    alpha_h1 = ncread(dataset, config, lens; sel = indices, defaults = 1.0, type = Float)

    lens = lens_input_parameter(config, "soil__thickness")
    soilthickness = ncread(dataset, config, lens; sel = indices, type = Float)

    lens =
        lens_input_parameter(config, "soil~compacted_surface_water__infiltration_capacity")
    infiltcappath =
        ncread(dataset, config, lens; sel = indices, defaults = 10.0, type = Float) .*
        (dt / BASETIMESTEP)
    lens =
        lens_input_parameter(config, "soil_water_sat-zone_bottom__max_leakage_volume_flux")
    maxleakage =
        ncread(dataset, config, lens; sel = indices, defaults = 0.0, type = Float) .*
        (dt / BASETIMESTEP)

    lens = lens_input_parameter(config, "soil_layer_water__brooks-corey_exponent")
    c = ncread(dataset, config, lens; sel = indices, type = Float, dimname = :layer)
    if size(c, 1) != maxlayers
        parname = lens(config)
        size1 = size(c, 1)
        error("$parname needs a layer dimension of size $maxlayers, but is $size1")
    end

    lens = lens_input_parameter(
        config,
        "soil_layer_water__vertical_saturated_hydraulic_conductivity_factor",
    )
    kvfrac = ncread(
        dataset,
        config,
        lens;
        sel = indices,
        defaults = 1.0,
        type = Float,
        dimname = :layer,
    )
    if size(kvfrac, 1) != maxlayers
        parname = lens(config)
        size1 = size(kvfrac, 1)
        error("$parname needs a layer dimension of size $maxlayers, but is $size1")
    end

    # soil infiltration capacity based on kv_0 and kvfrac upper soil layer
    infiltcapsoil = kv_0 .* @view kvfrac[1, :]
    # fraction compacted area
    lens = lens_input_parameter(config, "soil~compacted__area_fraction")
    pathfrac = ncread(dataset, config, lens; sel = indices, type = Float)

    # vegetation parameters
    lens = lens_input_parameter(config, "soil_root~wet__sigmoid_function_shape_parameter")
    rootdistpar =
        ncread(dataset, config, lens; sel = indices, defaults = -500.0, type = Float)
    lens = lens_input_parameter(
        config,
        "soil_water_sat-zone_top_capillary-rise__max_water-table_depth",
    )
    cap_hmax = ncread(dataset, config, lens; sel = indices, defaults = 2000.0, type = Float)

    lens = lens_input_parameter(
        config,
        "soil_water_sat-zone_top_capillary-rise__averianov_exponent",
    )
    cap_n = ncread(dataset, config, lens; sel = indices, defaults = 2.0, type = Float)

    act_thickl =
        set_layerthickness.(soilthickness, (cum_depth_layers,), (soil_layer_thickness,))
    sumlayers = @. pushfirst(cumsum(act_thickl), 0.0)
    nlayers = Vector{Int}(number_of_active_layers.(act_thickl))

    if length(config_soil_layer_thickness) > 0
        # root fraction read from dataset file, in case of multiple soil layers and TOML file
        # includes "soil_root__length_density_fraction"
        par_name = "soil_root__length_density_fraction"
        do_cyclic = haskey(config.input, "cyclic")
        do_root_fraction =
            do_cyclic ? haskey(config.input.cyclic, par_name) :
            haskey(config.input.static, par_name)
        if do_root_fraction
            lens = lens_input_parameter(config, par_name; optional = false)
            rootfraction =
                ncread(dataset, config, lens; sel = indices, type = Float, dimname = :layer)
        else
            n = length(indices)
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
        kvfrac = svectorscopy(kvfrac, Val{Int64(maxlayers)}()),  # casting N to 64-bit int is required because of bug; https://github.com/JuliaArrays/StaticArrays.jl/issues/1233
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
        rootfraction = svectorscopy(rootfraction, Val{Int64(maxlayers)}()),
        cap_hmax,
        cap_n,
        c = svectorscopy(c, Val{Int64(maxlayers)}()),
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
    n = Int(length(indices))
    params = SbmSoilParameters(dataset, config, vegetation_parameter_set, indices, dt)
    vars = SbmSoilVariables(n, params)
    bc = SbmSoilBC(n)
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
    temperature::Vector{Float},
)
    v = model.variables
    p = model.parameters
    @. v.tsoil = soil_temperature(v.tsoil, p.w_soil, temperature)
    return nothing
end

soil_temperature!(model::SbmSoilModel, snow::NoSnowModel, temperature::Vector{Float}) =
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
    unsaturated_zone_flow!(model::SbmSoilModel; topog_sbm_transfer = false)

Update unsaturated storage `ustorelayerdepth` and the `transfer` of water from the
unsaturated to the saturated store of the SBM soil model for a single timestep, based on the
original Topog_SBM formulation (`topog_sbm_transfer = true` and one soil layer) or the
Brooks-Corey approach (`topog_sbm_transfer = false` and one or multiple soil layers).
"""
function unsaturated_zone_flow!(model::SbmSoilModel; topog_sbm_transfer = false)
    v = model.variables
    p = model.parameters

    n = length(v.transfer)
    threaded_foreach(1:n; basesize = 250) do i
        if v.n_unsatlayers[i] > 0
            if topog_sbm_transfer && p.maxlayers == 1
                # original Topog_SBM formulation
                ustorelayerdepth = v.ustorelayerdepth[i][1] + v.infiltsoilpath[i]
                kv_z =
                    hydraulic_conductivity_at_depth(p.kv_profile, p.kvfrac, v.zi[i], i, 1)
                ustorelayerdepth, v.transfer[i] = unsatzone_flow_sbm(
                    ustorelayerdepth,
                    p.soilwatercapacity[i],
                    v.satwaterdepth[i],
                    kv_z,
                    v.ustorelayerthickness[1],
                    p.theta_s[i],
                    p.theta_r[i],
                )
                v.ustorelayerdepth[i] = setindex(v.ustorelayerdepth[i], ustorelayerdepth, 1)
            else
                # Brooks-Corey approach
                z = cumsum(v.ustorelayerthickness[i])
                flow_rate = 0.0
                for m in 1:v.n_unsatlayers[i]
                    l_sat = v.ustorelayerthickness[i][m] * (p.theta_s[i] - p.theta_r[i])
                    kv_z =
                        hydraulic_conductivity_at_depth(p.kv_profile, p.kvfrac, z[m], i, m)
                    ustorelayerdepth = if m == 1
                        v.ustorelayerdepth[i][m] + v.infiltsoilpath[i]
                    else
                        v.ustorelayerdepth[i][m] + flow_rate
                    end
                    ustorelayerdepth, flow_rate =
                        unsatzone_flow_layer(ustorelayerdepth, kv_z, l_sat, p.c[i][m])
                    v.ustorelayerdepth[i] =
                        setindex(v.ustorelayerdepth[i], ustorelayerdepth, m)
                end
                v.transfer[i] = flow_rate
            end
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
        if p.maxlayers == 1
            saturationdeficit = p.soilwatercapacity[i] - v.satwaterdepth[i]
            soilevapunsat =
                potsoilevap * min(1.0, saturationdeficit / p.soilwatercapacity[i])
        else
            soilevapunsat = soil_evaporation_unsatured_store(
                potsoilevap,
                v.ustorelayerdepth[i][1],
                v.ustorelayerthickness[i][1],
                v.n_unsatlayers[i],
                v.zi[i],
                p.theta_s[i] - p.theta_r[i],
            )
        end
        # Ensure that the unsaturated evaporation rate does not exceed the
        # available unsaturated moisture
        soilevapunsat = min(soilevapunsat, v.ustorelayerdepth[i][1])
        # Update the additional atmospheric demand
        potsoilevap = potsoilevap - soilevapunsat
        v.ustorelayerdepth[i] =
            setindex(v.ustorelayerdepth[i], v.ustorelayerdepth[i][1] - soilevapunsat, 1)

        if p.maxlayers == 1
            soilevapsat = 0.0
        else
            soilevapsat = soil_evaporation_satured_store(
                potsoilevap,
                v.n_unsatlayers[i],
                p.act_thickl[i][1],
                v.zi[i],
                p.theta_s[i] - p.theta_r[i],
            )
        end
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
function transpiration!(model::SbmSoilModel, dt::Float)
    (; potential_transpiration) = model.boundary_conditions
    v = model.variables
    p = model.parameters

    rootingdepth = get_rootingdepth(model)
    n = length(rootingdepth)
    threaded_foreach(1:n; basesize = 250) do i
        actevapustore = 0.0
        rootfraction_unsat = 0.0
        v.h3[i] = feddes_h3(p.h3_high[i], p.h3_low[i], potential_transpiration[i], dt)
        for k in 1:v.n_unsatlayers[i]
            vwc = max(
                v.ustorelayerdepth[i][k] / v.ustorelayerthickness[i][k],
                Float(0.0000001),
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
            # the rootfraction is valid for the root length in a soil layer, if zi decreases the root length
            # the rootfraction needs to be adapted
            if k == v.n_unsatlayers[i] && v.zi[i] < rootingdepth[i]
                rootlength = min(p.act_thickl[i][k], rootingdepth[i] - p.sumlayers[i][k])
                rootfraction_act =
                    p.rootfraction[i][k] * (v.ustorelayerthickness[i][k] / rootlength)
            else
                rootfraction_act = p.rootfraction[i][k]
            end
            actevapustore_layer =
                min(alpha * rootfraction_act * potential_transpiration[i], maxextr)
            rootfraction_unsat = rootfraction_unsat + rootfraction_act
            ustorelayerdepth = v.ustorelayerdepth[i][k] - actevapustore_layer
            actevapustore = actevapustore + actevapustore_layer
            v.ustorelayerdepth[i] = setindex(v.ustorelayerdepth[i], ustorelayerdepth, k)
        end

        # transpiration from saturated store
        wetroots = scurve(v.zi[i], rootingdepth[i], Float(1.0), p.rootdistpar[i])
        alpha = rwu_reduction_feddes(
            Float(0.0),
            p.h1[i],
            p.h2[i],
            v.h3[i],
            p.h4[i],
            p.alpha_h1[i],
        )
        # include remaining root fraction if rooting depth is below water table zi
        if v.zi[i] >= rootingdepth[i]
            f_roots = wetroots
            restevap = potential_transpiration[i] - actevapustore
        else
            f_roots = wetroots * (1.0 - rootfraction_unsat)
            restevap = potential_transpiration[i]
        end
        actevapsat = min(restevap * f_roots * alpha, v.satwaterdepth[i])

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
                netcapflux = netcapflux - toadd
                actcapflux = actcapflux + toadd
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
        dt::Float,
    )

Update the SBM soil model (infiltration, unsaturated zone flow, soil evaporation and
transpiration, capillary flux and leakage) for a single timestep.
"""
function update!(
    model::SbmSoilModel,
    atmospheric_forcing::AtmosphericForcing,
    external_models::NamedTuple,
    config::Config,
    dt::Float,
)
    soil_infiltration_reduction =
        get(config.model, "soil_infiltration_reduction__flag", false)::Bool
    modelsnow = get(config.model, "snow__flag", false)::Bool
    topog_sbm_transfer = get(config.model, "topog_sbm_transfer__flag", false)::Bool

    (; snow, runoff, demand) = external_models
    (; temperature) = atmospheric_forcing
    (; water_flux_surface) = model.boundary_conditions
    v = model.variables
    p = model.parameters

    ustoredepth!(model)
    @. v.ustorecapacity = p.soilwatercapacity - v.satwaterdepth - v.ustoredepth
    @. v.ustorelayerthickness = set_layerthickness(v.zi, p.sumlayers, p.act_thickl)
    @. v.n_unsatlayers = number_of_active_layers(v.ustorelayerthickness)

    # infiltration
    soil_temperature!(model, snow, temperature)
    infiltration_reduction_factor!(model; modelsnow, soil_infiltration_reduction)
    infiltration!(model)
    # unsaturated zone flow
    unsaturated_zone_flow!(model; topog_sbm_transfer)
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

    zi = get_water_depth(subsurface_flow) * 1000.0
    exfiltsatwater = get_exfiltwater(subsurface_flow) * 1000.0
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
        sbm_runoff =
            exfiltustore +
            exfiltsatwater[i] +
            v.excesswater[i] +
            runoff_land[i] +
            v.infiltexcess[i]

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

        # update the outputs and states
        v.n_unsatlayers[i] = n_unsatlayers
        v.ustorelayerdepth[i] = ustorelayerdepth
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
    end
    # update runoff and net_runoff (the runoff rate depends on the presence of paddy fields
    # and the h_max parameter of a paddy field)
    update_runoff!(demand.paddy, v.runoff)
    @. v.net_runoff = v.runoff - ae_openw_l
    return nothing
end

# wrapper method
get_rootingdepth(model::SbmSoilModel) =
    model.parameters.vegetation_parameter_set.rootingdepth