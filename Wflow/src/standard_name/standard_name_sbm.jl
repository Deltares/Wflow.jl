"""
Mapping of (CSDMS) standard names to model variables and units for models with a land model
of type `LandHydrologySBM`. The `lens` allows access to a nested model variable.
"""
const sbm_standard_name_map = Dict{String, VariableMetadata}(
    "atmosphere_water__precipitation_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.atmospheric_forcing.precipitation),
        unit = Unit(; mm = 1, dt = -1),
        description = "Precipitation",
    ),
    "land_surface_water__potential_evaporation_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.atmospheric_forcing.potential_evaporation),
        unit = Unit(; mm = 1, dt = -1),
        description = "Potential evaporation",
    ),
    "land_surface__evapotranspiration_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.actevap),
        unit = Unit(; mm = 1, dt = -1),
        description = "Actual evapotranspiration",
    ),
    "land_water_storage__total_depth" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.total_storage),
        unit = Unit(; mm = 1),
        description = "Total land water storage depth",
    ),
    "land_water_mass_balance_error__volume_flux" => VariableMetadata(;
        lens = @optic(_.mass_balance.land_water_balance.error),
        unit = Unit(; mm = 1, dt = -1),
        description = "Land water mass balance error",
    ),
    "land_water_mass_balance_relative_error__number" => VariableMetadata(;
        lens = @optic(_.mass_balance.land_water_balance.relative_error),
        description = "Land water mass balance relative error",
    ),
    "atmosphere_air__temperature" => VariableMetadata(;
        lens = @optic(_.land.atmospheric_forcing.temperature),
        unit = Unit(; degC = 1, absolute_temperature = true),
        description = "Air temperature",
    ),
    "vegetation_root__depth" => VariableMetadata(;
        lens = @optic(_.land.vegetation_parameters.rootingdepth),
        unit = Unit(; mm = 1),
        default = 750.0,
        description = "Rooting depth",
    ),
    "vegetation__crop_factor" => VariableMetadata(;
        lens = @optic(_.land.vegetation_parameters.kc),
        default = 1.0,
        description = "Crop factor",
    ),
    "vegetation__specific_leaf_storage" => VariableMetadata(;
        lens = @optic(_.land.vegetation_parameters.storage_specific_leaf),
        unit = Unit(; mm = 1),
        default = 0.5,
        description = "Specific leaf storage capacity",
    ),
    "vegetation_wood_water__storage_capacity" => VariableMetadata(;
        lens = @optic(_.land.vegetation_parameters.storage_wood),
        unit = Unit(; mm = 1),
        default = 0.0,
        description = "Wood water storage capacity",
    ),
    "vegetation_canopy__light_extinction_coefficient" => VariableMetadata(;
        lens = @optic(_.land.vegetation_parameters.kext),
        default = 0.5,
        description = "Light extinction coefficient",
    ),
    "vegetation_canopy__gap_fraction" => VariableMetadata(;
        lens = @optic(_.land.vegetation_parameters.canopygapfraction),
        default = 0.1,
        description = "Canopy gap fraction",
    ),
    "vegetation_canopy_water__mean_evaporation_to_mean_precipitation_ratio" =>
        VariableMetadata(;
            lens = @optic(_.land.interception.parameters.e_r),
            default = 0.25,
            description = "Evaporation to precipitation ratio",
        ),
    "vegetation_water__storage_capacity" => VariableMetadata(;
        lens = @optic(_.land.vegetation_parameters.cmax),
        unit = Unit(; mm = 1),
        default = 1.0,
        description = "Maximum canopy storage",
    ),
    "vegetation__leaf_area_index" => VariableMetadata(;
        lens = @optic(_.land.vegetation_parameters.leaf_area_index),
        unit = Unit(; m = (2, 2)),
        default = 3.0,
        description = "Leaf area index",
    ),
    "vegetation_canopy_water__depth" => VariableMetadata(;
        lens = @optic(_.land.interception.variables.canopy_storage),
        unit = Unit(; mm = 1),
        description = "Canopy water storage depth",
    ),
    "vegetation_canopy_water__stemflow_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.interception.variables.stemflow),
        unit = Unit(; mm = 1, dt = -1),
        description = "Stemflow flux",
    ),
    "vegetation_canopy_water__throughfall_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.interception.variables.throughfall),
        unit = Unit(; mm = 1, dt = -1),
        description = "Throughfall flux",
    ),
    "vegetation_canopy_water__interception_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.interception.variables.interception_rate),
        unit = Unit(; mm = 1, dt = -1),
        description = "Interception rate",
    ),
    "atmosphere_air__snowfall_temperature_threshold" => VariableMetadata(;
        lens = @optic(_.land.snow.parameters.tt),
        unit = Unit(; degC = 1, absolute_temperature = true),
        default = 0.0,
        description = "Temperature threshold for snowfall",
    ),
    "atmosphere_air__snowfall_temperature_interval" => VariableMetadata(;
        lens = @optic(_.land.snow.parameters.tti),
        unit = Unit(; degC = 1),
        default = 1.0,
        description = "Temperature interval for mixed precipitation",
    ),
    "snowpack__melting_temperature_threshold" => VariableMetadata(;
        lens = @optic(_.land.snow.parameters.ttm),
        unit = Unit(; degC = 1, absolute_temperature = true),
        default = 0.0,
        description = "Temperature threshold for snowmelt",
    ),
    "snowpack__liquid_water_holding_capacity" => VariableMetadata(;
        lens = @optic(_.land.snow.parameters.whc),
        default = 0.1,
        description = "Liquid water holding capacity of snowpack",
    ),
    "snowpack__degree_day_coefficient" => VariableMetadata(;
        lens = @optic(_.land.snow.parameters.cfmax),
        unit = Unit(; mm = 1, degC = -1, d = -1),
        default = 3.75,
        description = "Degree-day coefficient for snowmelt",
    ),
    "snowpack__leq_depth" => VariableMetadata(;
        lens = @optic(_.land.snow.variables.swe),
        unit = Unit(; mm = 1),
        description = "Snow water equivalent depth",
    ),
    "snowpack_meltwater__volume_flux" => VariableMetadata(;
        lens = @optic(_.land.snow.variables.snow_melt),
        unit = Unit(; mm = 1, dt = -1),
        description = "Snowmelt flux",
    ),
    "snowpack_water__runoff_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.snow.variables.runoff),
        unit = Unit(; mm = 1, dt = -1),
        description = "Snow runoff flux",
    ),
    "river_water__external_inflow_volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.external_inflow),
        unit = Unit(; m = 3, s = -1),
        description = "External inflow to river",
    ),
    "river_water__external_abstraction_volume_flow_rate" => VariableMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.actual_external_abstraction_av
        ),
        unit = Unit(; m = 3, s = -1),
        description = "External abstraction from river",
    ),
    "river_water__lateral_inflow_volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.inwater),
        unit = Unit(; m = 3, s = -1),
        description = "Lateral inflow to river",
    ),
    "river_water__instantaneous_volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.variables.q),
        unit = Unit(; m = 3, s = -1),
        description = "Instantaneous river discharge",
    ),
    "river_water__volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.variables.q_av),
        unit = Unit(; m = 3, s = -1),
        description = "Average river discharge",
    ),
    "river_water__depth" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.variables.h),
        unit = Unit(; m = 1),
        description = "River water depth",
    ),
    "river_water__volume" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.variables.storage),
        unit = Unit(; m = 3),
        description = "River water volume",
    ),
    "river_water_mass_balance_error__volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.mass_balance.routing.river_water_balance.error),
        unit = Unit(; m = 3, s = -1),
        description = "River water mass balance error",
    ),
    "river_water_mass_balance_relative_error__number" => VariableMetadata(;
        lens = @optic(_.mass_balance.routing.river_water_balance.relative_error),
        description = "River water mass balance relative error",
    ),
    "land_surface_water__abstraction_volume_flux" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.allocation.variables.act_surfacewater_abst),
        unit = Unit(; mm = 1, dt = -1),
        description = "Surface water abstraction",
    ),
    "floodplain_water__volume" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.floodplain.variables.storage),
        unit = Unit(; m = 3),
        description = "Floodplain water volume",
    ),
    "floodplain_water__depth" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.floodplain.variables.h),
        unit = Unit(; m = 1),
        description = "Floodplain water depth",
    ),
    "floodplain_water__instantaneous_volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.floodplain.variables.q),
        unit = Unit(; m = 3, s = -1),
        description = "Instantaneous floodplain discharge",
    ),
    "floodplain_water__volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.floodplain.variables.q_av),
        unit = Unit(; m = 3, s = -1),
        description = "Average floodplain discharge",
    ),
    "reservoir_water__target_min_volume_fraction" => VariableMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.targetminfrac
        ),
        default = 0.1,
        description = "Target minimum reservoir volume fraction",
    ),
    "reservoir_water__target_full_volume_fraction" => VariableMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.targetfullfrac
        ),
        default = 0.8,
        description = "Target full reservoir volume fraction",
    ),
    "reservoir_water_demand__required_downstream_volume_flow_rate" =>
        VariableMetadata(;
            lens = @optic(
                _.routing.river_flow.boundary_conditions.reservoir.parameters.demand
            ),
            unit = Unit(; m = 3, s = -1),
            default = 0.0,
            description = "Required downstream flow rate",
        ),
    "reservoir_water_release_below_spillway__max_volume_flow_rate" =>
        VariableMetadata(;
            lens = @optic(
                _.routing.river_flow.boundary_conditions.reservoir.parameters.maxrelease
            ),
            unit = Unit(; m = 3, s = -1),
            default = 0.0,
            description = "Maximum release below spillway",
        ),
    "reservoir_water__volume" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.storage),
        unit = Unit(; m = 3),
        description = "Reservoir water volume",
    ),
    "reservoir_water__outgoing_volume_flow_rate" => VariableMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.outflow_av
        ),
        unit = Unit(; m = 3, s = -1),
        description = "Reservoir outgoing flow rate",
    ),
    "reservoir_water__outgoing_observed_volume_flow_rate" => VariableMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.outflow_obs
        ),
        unit = Unit(; m = 3, s = -1),
        description = "Observed reservoir outgoing flow rate",
    ),
    "reservoir_water__incoming_volume_flow_rate" => VariableMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.inflow
        ),
        unit = Unit(; m = 3, s = -1),
        description = "Reservoir incoming flow rate",
    ),
    "reservoir_water__evaporation_volume_flux" => VariableMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.actevap),
        unit = Unit(; mm = 1, dt = -1),
        description = "Reservoir evaporation flux",
    ),
    "reservoir_water__precipitation_volume_flux" => VariableMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation
        ),
        unit = Unit(; mm = 1, dt = -1),
        description = "Reservoir precipitation flux",
    ),
    "reservoir_water__potential_evaporation_volume_flux" => VariableMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.evaporation
        ),
        unit = Unit(; mm = 1, dt = -1),
        description = "Reservoir potential evaporation flux",
    ),
    "reservoir_water_surface__elevation" => VariableMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.waterlevel
        ),
        unit = Unit(; m = 1),
        description = "Reservoir water surface elevation",
    ),
    "reservoir_water__external_inflow_volume_flow_rate" => VariableMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.external_inflow
        ),
        unit = Unit(; m = 3, s = -1),
        description = "External inflow to reservoir",
    ),
    "reservoir_water_mass_balance_error__volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.mass_balance.routing.reservoir_water_balance.error),
        unit = Unit(; m = 3, s = -1),
        description = "Reservoir water mass balance error",
    ),
    "reservoir_water_mass_balance_relative_error__number" => VariableMetadata(;
        lens = @optic(_.mass_balance.routing.reservoir_water_balance.relative_error),
        description = "Reservoir water mass balance relative error",
    ),
    "soil_water__infiltration_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.actinfilt),
        unit = Unit(; mm = 1, dt = -1),
        description = "Soil water infiltration flux",
    ),
    "soil_water__transpiration_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.transpiration),
        unit = Unit(; mm = 1, dt = -1),
        description = "Soil water transpiration flux",
    ),
    "soil_surface_temperature__weight_coefficient" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.w_soil),
        default = 0.1125,
        description = "Soil surface temperature weight coefficient",
    ),
    "soil_surface_water__infiltration_reduction_parameter" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.cf_soil),
        default = 0.038,
        description = "Soil surface infiltration reduction parameter",
    ),
    "soil_water__saturated_volume_fraction" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.theta_s),
        default = 0.6,
        description = "Soil saturated water content",
    ),
    "soil_water__residual_volume_fraction" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.theta_r),
        default = 0.01,
        description = "Soil residual water content",
    ),
    "soil_surface_water__vertical_saturated_hydraulic_conductivity" =>
        VariableMetadata(;
            lens = @optic(_.land.soil.parameters.kv_profile.kv_0),
            unit = Unit(; mm = 1, dt = -1),
            default = 3000.0,
            description = "Vertical saturated hydraulic conductivity",
        ),
    "soil_water__air_entry_pressure_head" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.hb),
        unit = Unit(; cm = 1),
        default = -10.0,
        description = "Air entry pressure head",
    ),
    "vegetation_root__feddes_critical_pressure_head_h1" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.h1),
        unit = Unit(; cm = 1),
        default = -10.0,
        description = "Feddes critical pressure head h1",
    ),
    "vegetation_root__feddes_critical_pressure_head_h2" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.h2),
        unit = Unit(; cm = 1),
        default = -25.0,
        description = "Feddes critical pressure head h2",
    ),
    "vegetation_root__feddes_critical_pressure_head_h3_high" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.h3_high),
        unit = Unit(; cm = 1),
        default = -400.0,
        description = "Feddes critical pressure head h3 high",
    ),
    "vegetation_root__feddes_critical_pressure_head_h3_low" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.h3_low),
        unit = Unit(; cm = 1),
        default = -1000.0,
        description = "Feddes critical pressure head h3 low",
    ),
    "vegetation_root__feddes_critical_pressure_head_h4" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.h4),
        unit = Unit(; cm = 1),
        default = -15849.0,
        description = "Feddes critical pressure head h4",
    ),
    "vegetation_root__feddes_critical_pressure_head_h1_reduction_coefficient" =>
        VariableMetadata(;
            lens = @optic(_.land.soil.parameters.alpha_h1),
            default = 1.0,
            description = "Feddes h1 reduction coefficient",
        ),
    "soil__thickness" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.soilthickness),
        unit = Unit(; mm = 1),
        default = 2000.0,
        description = "Soil thickness",
    ),
    "compacted_soil_surface_water__infiltration_capacity" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.infiltcappath),
        unit = Unit(; mm = 1, dt = -1),
        default = 10.0,
        description = "Compacted soil infiltration capacity",
    ),
    "soil_water_saturated_zone_bottom__max_leakage_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.maxleakage),
        unit = Unit(; mm = 1, dt = -1),
        default = 0.0,
        description = "Maximum leakage from saturated zone bottom",
    ),
    "soil_layer_water__brooks_corey_exponent" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.c),
        default = 10.0,
        description = "Brooks-Corey exponent",
    ),
    "soil_layer_water__vertical_saturated_hydraulic_conductivity_factor" =>
        VariableMetadata(;
            lens = @optic(_.land.soil.parameters.kvfrac),
            default = 1.0,
            description = "Vertical saturated hydraulic conductivity factor",
        ),
    "compacted_soil__area_fraction" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.pathfrac),
        default = 0.01,
        description = "Compacted soil area fraction",
    ),
    "soil_wet_root__sigmoid_function_shape_parameter" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.rootdistpar),
        default = -500.0,
        description = "Root distribution shape parameter",
    ),
    "soil_water_saturated_zone_top__capillary_rise_max_water_table_depth" =>
        VariableMetadata(;
            lens = @optic(_.land.soil.parameters.cap_hmax),
            unit = Unit(; mm = 1),
            default = 2000.0,
            description = "Maximum water table depth for capillary rise",
        ),
    "soil_water_saturated_zone_top__capillary_rise_averianov_exponent" =>
        VariableMetadata(;
            lens = @optic(_.land.soil.parameters.cap_n),
            default = 2.0,
            description = "Averianov exponent for capillary rise",
        ),
    "soil_root__length_density_fraction" => VariableMetadata(;
        lens = @optic(_.land.soil.parameters.rootfraction),
        default = 1.0,
        description = "Root length density fraction",
    ),
    "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter" =>
        VariableMetadata(;
            lens = nothing,
            unit = Unit(; mm = -1),
            description = "Vertical hydraulic conductivity scale parameter",
        ),
    "soil_layer_water__vertical_saturated_hydraulic_conductivity" => VariableMetadata(;
        lens = nothing,
        unit = Unit(; mm = 1, dt = -1),
        description = "Layer vertical saturated hydraulic conductivity",
    ),
    "soil_exponential_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth" =>
        VariableMetadata(;
            lens = nothing,
            unit = Unit(; mm = 1),
            description = "Exponential profile depth",
        ),
    "soil_layered_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth" =>
        VariableMetadata(;
            lens = nothing,
            unit = Unit(; mm = 1),
            description = "Layered profile depth",
        ),
    "soil_surface_water__runoff_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.runoff),
        unit = Unit(; mm = 1, dt = -1),
        description = "Soil surface runoff flux",
    ),
    "soil_surface_water__net_runoff_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.net_runoff),
        unit = Unit(; mm = 1, dt = -1),
        description = "Net soil surface runoff flux",
    ),
    "soil_surface_water_unsaturated_zone__exfiltration_volume_flux" =>
        VariableMetadata(;
            lens = @optic(_.land.soil.variables.exfiltustore),
            unit = Unit(; mm = 1, dt = -1),
            description = "Exfiltration from unsaturated zone",
        ),
    "soil_surface_water_saturated_zone__exfiltration_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.exfiltustore),
        unit = Unit(; mm = 1, dt = -1),
        description = "Exfiltration from saturated zone",
    ),
    "compacted_soil_surface_water__excess_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.excesswatersoil),
        unit = Unit(; mm = 1, dt = -1),
        description = "Excess water flux from compacted soil",
    ),
    "non_compacted_soil_surface_water__excess_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.excesswaterpath),
        unit = Unit(; mm = 1, dt = -1),
        description = "Excess water flux from non-compacted soil",
    ),
    "soil_layer_water__volume_fraction" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.vwc),
        unit = Unit(; m = (3, 3)),
        description = "Soil layer water volume fraction",
    ),
    "soil_layer_water__volume_percentage" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.vwc_perc),
        unit = Unit(; percentage = 1),
        description = "Soil layer water volume percentage",
    ),
    "soil_water_root_zone__volume_fraction" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.vwc_root),
        unit = Unit(; m = (3, 3)),
        description = "Root zone water volume fraction",
    ),
    "soil_water_root_zone__volume_percentage" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.vwc_percroot),
        unit = Unit(; percentage = 1),
        description = "Root zone water volume percentage",
    ),
    "soil_water_root_zone__depth" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.rootstore),
        unit = Unit(; mm = 1),
        description = "Root zone water depth",
    ),
    "soil_layer_water_unsaturated_zone__depth" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.ustorelayerdepth),
        unit = Unit(; mm = 1),
        description = "Unsaturated zone layer water depth",
    ),
    "soil_water_unsaturated_zone__depth" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.ustoredepth),
        unit = Unit(; mm = 1),
        description = "Unsaturated zone water depth",
    ),
    "soil_water_saturated_zone_top__capillary_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.actcapflux),
        unit = Unit(; mm = 1, dt = -1),
        description = "Capillary flux at saturated zone top",
    ),
    "soil_water_saturated_zone_top__recharge_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.transfer),
        unit = Unit(; mm = 1, dt = -1),
        description = "Recharge flux to saturated zone",
    ),
    "soil_water_saturated_zone_top__net_recharge_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.recharge),
        unit = Unit(; mm = 1, dt = -1),
        description = "Net recharge flux to saturated zone",
    ),
    "soil_water_saturated_zone_bottom__leakage_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.actleakage),
        unit = Unit(; mm = 1, dt = -1),
        description = "Leakage flux from saturated zone bottom",
    ),
    "soil_water_saturated_zone__depth" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.satwaterdepth),
        unit = Unit(; mm = 1),
        description = "Saturated zone water depth",
    ),
    "soil_water_saturated_zone_top__depth" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.zi),
        unit = Unit(; mm = 1),
        description = "Water table depth",
    ),
    "soil_surface__temperature" => VariableMetadata(;
        lens = @optic(_.land.soil.variables.tsoil),
        unit = Unit(; degC = 1, absolute_temperature = true),
        description = "Soil surface temperature",
    ),
    "subsurface_water_saturated_zone_top__depth" => VariableMetadata(;
        lens = @optic(_.routing.subsurface_flow.variables.zi),
        unit = Unit(; m = 1),
        description = "Subsurface water table depth",
    ),
    "subsurface_water__exfiltration_volume_flux" => VariableMetadata(;
        lens = @optic(_.routing.subsurface_flow.variables.exfiltwater),
        unit = Unit(; m = 1, dt = -1),
        description = "Subsurface water exfiltration flux",
    ),
    "subsurface_water__volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.subsurface_flow.variables.ssf),
        unit = Unit(; m = 3, d = -1),
        description = "Subsurface water flow rate",
    ),
    "subsurface_water__to_river_volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.subsurface_flow.variables.to_river),
        unit = Unit(; m = 3, d = -1),
        description = "Subsurface water flow to river",
    ),
    "subsurface_water_mass_balance_error__volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.mass_balance.routing.subsurface_water_balance.error),
        unit = Unit(; m = 3, d = -1),
        description = "Subsurface water mass balance error",
    ),
    "subsurface_water_mass_balance_relative_error__number" => VariableMetadata(;
        lens = @optic(_.mass_balance.routing.subsurface_water_balance.relative_error),
        description = "Subsurface water mass balance relative error",
    ),
    "snowpack_liquid_water__depth" => VariableMetadata(;
        lens = @optic(_.land.snow.variables.snow_water),
        unit = Unit(; mm = 1),
        description = "Liquid water depth in snowpack",
    ),
    "snowpack_dry_snow__leq_depth" => VariableMetadata(;
        lens = @optic(_.land.snow.variables.snow_storage),
        unit = Unit(; mm = 1),
        description = "Dry snow storage depth",
    ),
    "glacier_ice__leq_depth" => VariableMetadata(;
        lens = @optic(_.land.glacier.variables.glacier_store),
        unit = Unit(; mm = 1),
        description = "Glacier ice storage depth",
    ),
    "glacier_ice__melt_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.glacier.variables.glacier_melt),
        unit = Unit(; mm = 1, dt = -1),
        description = "Glacier ice melt flux",
    ),
    "glacier_ice__initial_leq_depth" => VariableMetadata(;
        lens = @optic(_.land.glacier.variables.glacier_store),
        unit = Unit(; mm = 1),
        description = "Initial glacier ice storage depth",
    ),
    "glacier_ice__melting_temperature_threshold" => VariableMetadata(;
        lens = @optic(_.land.glacier.parameters.g_ttm),
        unit = Unit(; degC = 1, absolute_temperature = true),
        default = 0.0,
        description = "Glacier ice melting temperature threshold",
    ),
    "glacier_ice__degree_day_coefficient" => VariableMetadata(;
        lens = @optic(_.land.glacier.parameters.g_cfmax),
        unit = Unit(; mm = 1, degC = -1, dt = -1),
        default = 3.0,
        description = "Glacier ice degree-day coefficient",
    ),
    "glacier_firn_accumulation__snowpack_dry_snow_leq_depth_fraction" =>
        VariableMetadata(;
            lens = @optic(_.land.glacier.parameters.g_sifrac),
            unit = Unit(; dt = -1),
            default = 0.001,
            description = "Firn accumulation fraction",
        ),
    "glacier_surface__area_fraction" => VariableMetadata(;
        lens = @optic(_.land.glacier.parameters.glacier_frac),
        default = 0.0,
        description = "Glacier surface area fraction",
    ),
    "land_surface_water__instantaneous_volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.overland_flow.variables.q),
        unit = Unit(; m = 3, s = -1),
        description = "Instantaneous overland flow rate",
    ),
    "land_surface_water__volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.overland_flow.variables.q_av),
        unit = Unit(; m = 3, s = -1),
        description = "Average overland flow rate",
    ),
    "land_surface_water__to_river_volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.overland_flow.variables.to_river),
        unit = Unit(; m = 3, s = -1),
        description = "Overland flow to river",
    ),
    "land_surface_water__depth" => VariableMetadata(;
        lens = @optic(_.routing.overland_flow.variables.h),
        unit = Unit(; m = 1),
        description = "Overland flow water depth",
    ),
    "land_surface_water__volume" => VariableMetadata(;
        lens = @optic(_.routing.overland_flow.variables.storage),
        unit = Unit(; m = 3),
        description = "Overland flow water volume",
    ),
    "land_surface_water__x_component_of_instantaneous_volume_flow_rate" =>
        VariableMetadata(;
            lens = @optic(_.routing.overland_flow.variables.qx),
            unit = Unit(; m = 3, s = -1),
            description = "X-component of instantaneous overland flow",
        ),
    "land_surface_water__y_component_of_instantaneous_volume_flow_rate" =>
        VariableMetadata(;
            lens = @optic(_.routing.overland_flow.variables.qy),
            unit = Unit(; m = 3, s = -1),
            description = "Y-component of instantaneous overland flow",
        ),
    "land_surface_water_mass_balance_error__volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.mass_balance.routing.overland_water_balance.error),
        unit = Unit(; m = 3, s = -1),
        description = "Overland flow mass balance error",
    ),
    "land_surface_water_mass_balance_relative_error__number" => VariableMetadata(;
        lens = @optic(_.mass_balance.routing.overland_water_balance.relative_error),
        description = "Overland flow mass balance relative error",
    ),
    "paddy_surface_water__depth" => VariableMetadata(;
        lens = @optic(_.land.demand.paddy.variables.h),
        unit = Unit(; mm = 1),
        description = "Paddy surface water depth",
    ),
    "domestic__gross_water_demand_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.demand.domestic.demand.demand_gross),
        unit = Unit(; mm = 1, dt = -1),
        description = "Domestic gross water demand flux",
    ),
    "domestic__net_water_demand_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.demand.domestic.demand.demand_net),
        unit = Unit(; mm = 1, dt = -1),
        description = "Domestic net water demand flux",
    ),
    "industry__gross_water_demand_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.demand.industry.demand.demand_gross),
        unit = Unit(; mm = 1, dt = -1),
        description = "Industry gross water demand flux",
    ),
    "industry__net_water_demand_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.demand.industry.demand.demand_net),
        unit = Unit(; mm = 1, dt = -1),
        description = "Industry net water demand flux",
    ),
    "livestock__gross_water_demand_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.demand.livestock.demand.demand_gross),
        unit = Unit(; mm = 1, dt = -1),
        description = "Livestock gross water demand flux",
    ),
    "livestock__net_water_demand_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.demand.livestock.demand.demand_net),
        unit = Unit(; mm = 1, dt = -1),
        description = "Livestock net water demand flux",
    ),
    "irrigated_paddy__min_depth" => VariableMetadata(;
        lens = @optic(_.land.demand.paddy.parameters.h_min),
        unit = Unit(; mm = 1),
        description = "Minimum paddy water depth",
    ),
    "irrigated_paddy__optimal_depth" => VariableMetadata(;
        lens = @optic(_.land.demand.paddy.parameters.h_opt),
        unit = Unit(; mm = 1),
        description = "Optimal paddy water depth",
    ),
    "irrigated_paddy__max_depth" => VariableMetadata(;
        lens = @optic(_.land.demand.paddy.parameters.h_max),
        unit = Unit(; mm = 1),
        description = "Maximum paddy water depth",
    ),
    "irrigated_paddy__irrigation_efficiency" => VariableMetadata(;
        lens = @optic(_.land.demand.paddy.parameters.irrigation_efficiency),
        description = "Paddy irrigation efficiency",
    ),
    "irrigated_paddy_area__count" => VariableMetadata(;
        lens = @optic(_.land.demand.paddy.parameters.irrigation_areas),
        description = "Number of paddy irrigation areas",
    ),
    "irrigated_paddy__irrigation_trigger_flag" => VariableMetadata(;
        lens = @optic(_.land.demand.paddy.parameters.irrigation_trigger),
        description = "Paddy irrigation trigger flag",
    ),
    "irrigated_paddy__max_irrigation_rate" => VariableMetadata(;
        lens = @optic(_.land.demand.paddy.parameters.maximum_irrigation_rate),
        unit = Unit(; mm = 1, dt = -1),
        description = "Maximum paddy irrigation rate",
    ),
    "irrigated_paddy__gross_water_demand_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.demand.paddy.variables.demand_gross),
        unit = Unit(; mm = 1, dt = -1),
        description = "Paddy gross water demand flux",
    ),
    "irrigated_non_paddy__irrigation_efficiency" => VariableMetadata(;
        lens = @optic(_.land.demand.nonpaddy.parameters.irrigation_efficiency),
        description = "Non-paddy irrigation efficiency",
    ),
    "irrigated_non_paddy_area__count" => VariableMetadata(;
        lens = @optic(_.land.demand.nonpaddy.parameters.irrigation_areas),
        description = "Number of non-paddy irrigation areas",
    ),
    "irrigated_non_paddy__irrigation_trigger_flag" => VariableMetadata(;
        lens = @optic(_.land.demand.nonpaddy.parameters.irrigation_trigger),
        description = "Non-paddy irrigation trigger flag",
    ),
    "irrigated_non_paddy__max_irrigation_rate" => VariableMetadata(;
        lens = @optic(_.land.demand.nonpaddy.parameters.maximum_irrigation_rate),
        unit = Unit(; mm = 1, dt = -1),
        description = "Maximum non-paddy irrigation rate",
    ),
    "irrigated_non_paddy__gross_water_demand_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.demand.nonpaddy.variables.demand_gross),
        unit = Unit(; mm = 1, dt = -1),
        description = "Non-paddy gross water demand flux",
    ),
    "land__allocated_irrigation_water_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.allocation.variables.irri_alloc),
        unit = Unit(; mm = 1, dt = -1),
        description = "Allocated irrigation water flux",
    ),
    "land__allocated_non_irrigation_water_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.allocation.variables.nonirri_alloc),
        unit = Unit(; mm = 1, dt = -1),
        description = "Allocated non-irrigation water flux",
    ),
    "subsurface_water__abstraction_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.allocation.variables.act_groundwater_abst),
        unit = Unit(; mm = 1, dt = -1),
        description = "Subsurface water abstraction flux",
    ),
    "land__non_irrigation_return_flow_volume_flux" => VariableMetadata(;
        lens = @optic(_.land.allocation.variables.nonirri_returnflow),
        unit = Unit(; mm = 1, dt = -1),
        description = "Non-irrigation return flow flux",
    ),
    "subsurface_water__hydraulic_head" => VariableMetadata(;
        lens = @optic(_.routing.subsurface_flow.aquifer.variables.head),
        unit = Unit(; m = 1),
        description = "Subsurface water hydraulic head",
    ),
    "subsurface_water_saturated_zone_top__net_recharge_volume_flow_rate" =>
        VariableMetadata(;
            lens = @optic(_.routing.subsurface_flow.boundaries.recharge.variables.flux_av),
            unit = Unit(; m = 3, d = -1),
            description = "Net recharge volume flow rate to saturated zone",
        ),
    "land_drain_water__to_subsurface_volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.variables.flux_av),
        unit = Unit(; m = 3, d = -1),
        description = "Land drain water to subsurface volume flow rate",
    ),
    "river_water__to_subsurface_volume_flow_rate" => VariableMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.river.variables.flux_av),
        unit = Unit(; m = 3, d = -1),
        description = "River water to subsurface volume flow rate",
    ),
    "land_surface_water__withdrawal_fraction" => VariableMetadata(;
        lens = @optic(_.land.allocation.parameters.frac_sw_used),
        description = "Surface water withdrawal fraction",
    ),
    "land_water_allocation_area__count" => VariableMetadata(;
        lens = @optic(_.land.allocation.parameters.areas),
        description = "Number of water allocation areas",
    ),
)
