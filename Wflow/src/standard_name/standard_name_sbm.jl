"""
Mapping of (CSDMS) standard names to model variables and units for models with a land model
of type `LandHydrologySBM`. The `lens` of the NamedTuple allows access to a nested model
variable.
"""
const sbm_standard_name_map = Dict{String, NamedTuple}(
    "atmosphere_water__precipitation_volume_flux" => (
        lens = @optic(_.land.atmospheric_forcing.precipitation),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "land_surface_water__potential_evaporation_volume_flux" => (
        lens = @optic(_.land.atmospheric_forcing.potential_evaporation),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "land_surface__evapotranspiration_volume_flux" =>
        (lens = @optic(_.land.soil.variables.actevap), unit = Unit(; mm = 1, dt = -1)),
    "land_water_storage__total_depth" =>
        (lens = @optic(_.land.soil.variables.total_storage), unit = Unit(; mm = 1)),
    "land_water_mass_balance_error__volume_flux" => (
        lens = @optic(_.mass_balance.land_water_balance.error),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "land_water_mass_balance_relative_error__number" => (
        lens = @optic(_.mass_balance.land_water_balance.relative_error),
        unit = Unit(),
    ),
    "atmosphere_air__temperature" => (
        lens = @optic(_.land.atmospheric_forcing.temperature),
        unit = Unit(; degC = 1, absolute_temperature = true),
    ),
    "vegetation_root__depth" => (
        lens = @optic(_.land.vegetation_parameters.rootingdepth),
        unit = Unit(; mm = 1),
    ),
    "vegetation__crop_factor" =>
        (lens = @optic(_.land.vegetation_parameters.kc), unit = Unit()),
    "vegetation__specific_leaf_storage" => (
        lens = @optic(_.land.vegetation_parameters.storage_specific_leaf),
        unit = Unit(; mm = 1),
    ),
    "vegetation_wood_water__storage_capacity" => (
        lens = @optic(_.land.vegetation_parameters.storage_wood),
        unit = Unit(; mm = 1),
    ),
    "vegetation_canopy__light_extinction_coefficient" =>
        (lens = @optic(_.land.vegetation_parameters.kext), unit = Unit()),
    "vegetation_canopy__gap_fraction" =>
        (lens = @optic(_.land.vegetation_parameters.canopygapfraction), unit = Unit()),
    "vegetation_canopy_water__mean_evaporation_to_mean_precipitation_ratio" =>
        (lens = @optic(_.land.interception.parameters.e_r), unit = Unit()),
    "vegetation_water__storage_capacity" =>
        (lens = @optic(_.land.vegetation_parameters.cmax), unit = Unit(; mm = 1)),
    "vegetation__leaf_area_index" => (
        lens = @optic(_.land.vegetation_parameters.leaf_area_index),
        unit = Unit(; m = (2, 2)),
    ),
    "vegetation_canopy_water__depth" => (
        lens = @optic(_.land.interception.variables.canopy_storage),
        unit = Unit(; mm = 1),
    ),
    "vegetation_canopy_water__stemflow_volume_flux" => (
        lens = @optic(_.land.interception.variables.stemflow),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "vegetation_canopy_water__throughfall_volume_flux" => (
        lens = @optic(_.land.interception.variables.throughfall),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "vegetation_canopy_water__interception_volume_flux" => (
        lens = @optic(_.land.interception.variables.interception_rate),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "atmosphere_air__snowfall_temperature_threshold" => (
        lens = @optic(_.land.snow.parameters.tt),
        unit = Unit(; degC = 1, absolute_temperature = true),
    ),
    "atmosphere_air__snowfall_temperature_interval" =>
        (lens = @optic(_.land.snow.parameters.tti), unit = Unit(; degC = 1)),
    "snowpack__melting_temperature_threshold" => (
        lens = @optic(_.land.snow.parameters.ttm),
        unit = Unit(; degC = 1, absolute_temperature = true),
    ),
    "snowpack__liquid_water_holding_capacity" =>
        (lens = @optic(_.land.snow.parameters.whc), unit = Unit()),
    "snowpack__degree_day_coefficient" => (
        lens = @optic(_.land.snow.parameters.cfmax),
        unit = Unit(; mm = 1, degC = -1, d = -1),
    ),
    "snowpack__leq_depth" =>
        (lens = @optic(_.land.snow.variables.swe), unit = Unit(; mm = 1)),
    "snowpack_meltwater__volume_flux" => (
        lens = @optic(_.land.snow.variables.snow_melt),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "snowpack_water__runoff_volume_flux" =>
        (lens = @optic(_.land.snow.variables.runoff), unit = Unit(; mm = 1, dt = -1)),
    "river_water__external_inflow_volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.external_inflow),
        unit = Unit(; m = 3, s = -1),
    ),
    "river_water__external_abstraction_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.actual_external_abstraction_av
        ),
        unit = Unit(; m = 3, s = -1),
    ),
    "river_water__lateral_inflow_volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.inwater),
        unit = Unit(; m = 3, s = -1),
    ),
    "river_water__instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.river_flow.variables.q), unit = Unit(; m = 3, s = -1)),
    "river_water__volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.variables.q_av),
        unit = Unit(; m = 3, s = -1),
    ),
    "river_water__depth" =>
        (lens = @optic(_.routing.river_flow.variables.h), unit = Unit(; m = 1)),
    "river_water__volume" =>
        (lens = @optic(_.routing.river_flow.variables.storage), unit = Unit(; m = 3)),
    "river_water_mass_balance_error__volume_flow_rate" => (
        lens = @optic(_.mass_balance.routing.river_water_balance.error),
        unit = Unit(; m = 3, s = -1),
    ),
    "river_water_mass_balance_relative_error__number" => (
        lens = @optic(_.mass_balance.routing.river_water_balance.relative_error),
        unit = Unit(),
    ),
    "land_surface_water__abstraction_volume_flux" => (
        lens = @optic(_.routing.river_flow.allocation.variables.act_surfacewater_abst),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "floodplain_water__volume" => (
        lens = @optic(_.routing.river_flow.floodplain.variables.storage),
        unit = Unit(; m = 3),
    ),
    "floodplain_water__depth" => (
        lens = @optic(_.routing.river_flow.floodplain.variables.h),
        unit = Unit(; m = 1),
    ),
    "floodplain_water__instantaneous_volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.floodplain.variables.q),
        unit = Unit(; m = 3, s = -1),
    ),
    "floodplain_water__volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.floodplain.variables.q_av),
        unit = Unit(; m = 3, s = -1),
    ),
    "reservoir_water__target_min_volume_fraction" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.targetminfrac
        ),
        unit = Unit(),
    ),
    "reservoir_water__target_full_volume_fraction" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.targetfullfrac
        ),
        unit = Unit(),
    ),
    "reservoir_water_demand__required_downstream_volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.parameters.demand),
        unit = Unit(; m = 3, s = -1),
    ),
    "reservoir_water_release_below_spillway__max_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.maxrelease
        ),
        unit = Unit(; m = 3, s = -1),
    ),
    "reservoir_water__volume" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.storage),
        unit = Unit(; m = 3),
    ),
    "reservoir_water__outgoing_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.outflow_av
        ),
        unit = Unit(; m = 3, s = -1),
    ),
    "reservoir_water__outgoing_observed_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.outflow_obs
        ),
        unit = Unit(; m = 3, s = -1),
    ),
    "reservoir_water__incoming_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.inflow
        ),
        unit = Unit(; m = 3, s = -1),
    ),
    "reservoir_water__evaporation_volume_flux" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.actevap),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "reservoir_water__precipitation_volume_flux" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation
        ),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "reservoir_water__potential_evaporation_volume_flux" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.evaporation
        ),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "reservoir_water_surface__elevation" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.waterlevel
        ),
        unit = Unit(; m = 1),
    ),
    "reservoir_water__external_inflow_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.external_inflow
        ),
        unit = Unit(; m = 3, s = -1),
    ),
    "reservoir_water_mass_balance_error__volume_flow_rate" => (
        lens = @optic(_.mass_balance.routing.reservoir_water_balance.error),
        unit = Unit(; m = 3, s = -1),
    ),
    "reservoir_water_mass_balance_relative_error__number" => (
        lens = @optic(_.mass_balance.routing.reservoir_water_balance.relative_error),
        unit = Unit(),
    ),
    "soil_water__infiltration_volume_flux" => (
        lens = @optic(_.land.soil.variables.actinfilt),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_water__transpiration_volume_flux" => (
        lens = @optic(_.land.soil.variables.transpiration),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_surface_temperature__weight_coefficient" =>
        (lens = @optic(_.land.soil.parameters.w_soil), unit = Unit()),
    "soil_surface_water__infiltration_reduction_parameter" =>
        (lens = @optic(_.land.soil.parameters.cf_soil), unit = Unit()),
    "soil_water__saturated_volume_fraction" =>
        (lens = @optic(_.land.soil.parameters.theta_s), unit = Unit()),
    "soil_water__residual_volume_fraction" =>
        (lens = @optic(_.land.soil.parameters.theta_r), unit = Unit()),
    "soil_surface_water__vertical_saturated_hydraulic_conductivity" => (
        lens = @optic(_.land.soil.parameters.kv_profile.kv_0),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_water__air_entry_pressure_head" =>
        (lens = @optic(_.land.soil.parameters.hb), unit = Unit(; cm = 1)),
    "vegetation_root__feddes_critical_pressure_head_h1" =>
        (lens = @optic(_.land.soil.parameters.h1), unit = Unit(; cm = 1)),
    "vegetation_root__feddes_critical_pressure_head_h2" =>
        (lens = @optic(_.land.soil.parameters.h2), unit = Unit(; cm = 1)),
    "vegetation_root__feddes_critical_pressure_head_h3_high" =>
        (lens = @optic(_.land.soil.parameters.h3_high), unit = Unit(; cm = 1)),
    "vegetation_root__feddes_critical_pressure_head_h3_low" =>
        (lens = @optic(_.land.soil.parameters.h3_low), unit = Unit(; cm = 1)),
    "vegetation_root__feddes_critical_pressure_head_h4" =>
        (lens = @optic(_.land.soil.parameters.h4), unit = Unit(; cm = 1)),
    "vegetation_root__feddes_critical_pressure_head_h1_reduction_coefficient" =>
        (lens = @optic(_.land.soil.parameters.alpha_h1), unit = Unit()),
    "soil__thickness" =>
        (lens = @optic(_.land.soil.parameters.soilthickness), unit = Unit(; mm = 1)),
    "compacted_soil_surface_water__infiltration_capacity" => (
        lens = @optic(_.land.soil.parameters.infiltcappath),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_water_saturated_zone_bottom__max_leakage_volume_flux" => (
        lens = @optic(_.land.soil.parameters.maxleakage),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_layer_water__brooks_corey_exponent" =>
        (lens = @optic(_.land.soil.parameters.c), unit = Unit()),
    "soil_layer_water__vertical_saturated_hydraulic_conductivity_factor" =>
        (lens = @optic(_.land.soil.parameters.kvfrac), unit = Unit()),
    "compacted_soil__area_fraction" =>
        (lens = @optic(_.land.soil.parameters.pathfrac), unit = Unit()),
    "soil_wet_root__sigmoid_function_shape_parameter" =>
        (lens = @optic(_.land.soil.parameters.rootdistpar), unit = Unit()),
    "soil_water_saturated_zone_top__capillary_rise_max_water_table_depth" =>
        (lens = @optic(_.land.soil.parameters.cap_hmax), unit = Unit(; mm = 1)),
    "soil_water_saturated_zone_top__capillary_rise_averianov_exponent" =>
        (lens = @optic(_.land.soil.parameters.cap_n), unit = Unit()),
    "soil_root__length_density_fraction" =>
        (lens = @optic(_.land.soil.parameters.rootfraction), unit = Unit()),
    "soil_water__vertical_saturated_hydraulic_conductivity_scale_parameter" =>
        (lens = nothing, unit = Unit(; mm = -1)),
    "soil_layer_water__vertical_saturated_hydraulic_conductivity" =>
        (lens = nothing, unit = Unit(; mm = 1, dt = -1)),
    "soil_exponential_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth" =>
        (lens = nothing, unit = Unit(; mm = 1)),
    "soil_layered_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth" =>
        (lens = nothing, unit = Unit(; mm = 1)),
    "soil_surface_water__runoff_volume_flux" =>
        (lens = @optic(_.land.soil.variables.runoff), unit = Unit(; mm = 1, dt = -1)),
    "soil_surface_water__net_runoff_volume_flux" => (
        lens = @optic(_.land.soil.variables.net_runoff),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_surface_water_unsaturated_zone__exfiltration_volume_flux" => (
        lens = @optic(_.land.soil.variables.exfiltustore),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_surface_water_saturated_zone__exfiltration_volume_flux" => (
        lens = @optic(_.land.soil.variables.exfiltustore),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "compacted_soil_surface_water__excess_volume_flux" => (
        lens = @optic(_.land.soil.variables.excesswatersoil),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "non_compacted_soil_surface_water__excess_volume_flux" => (
        lens = @optic(_.land.soil.variables.excesswaterpath),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_layer_water__volume_fraction" =>
        (lens = @optic(_.land.soil.variables.vwc), unit = Unit(; m = (3, 3))),
    "soil_layer_water__volume_percentage" =>
        (lens = @optic(_.land.soil.variables.vwc_perc), unit = Unit(; percentage = 1)),
    "soil_water_root_zone__volume_fraction" =>
        (lens = @optic(_.land.soil.variables.vwc_root), unit = Unit(; m = (3, 3))),
    "soil_water_root_zone__volume_percentage" => (
        lens = @optic(_.land.soil.variables.vwc_percroot),
        unit = Unit(; percentage = 1),
    ),
    "soil_water_root_zone__depth" =>
        (lens = @optic(_.land.soil.variables.rootstore), unit = Unit(; mm = 1)),
    "soil_layer_water_unsaturated_zone__depth" =>
        (lens = @optic(_.land.soil.variables.ustorelayerdepth), unit = Unit(; mm = 1)),
    "soil_water_unsaturated_zone__depth" =>
        (lens = @optic(_.land.soil.variables.ustoredepth), unit = Unit(; mm = 1)),
    "soil_water_saturated_zone_top__capillary_volume_flux" => (
        lens = @optic(_.land.soil.variables.actcapflux),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_water_saturated_zone_top__recharge_volume_flux" =>
        (lens = @optic(_.land.soil.variables.transfer), unit = Unit(; mm = 1, dt = -1)),
    "soil_water_saturated_zone_top__net_recharge_volume_flux" =>
        (lens = @optic(_.land.soil.variables.recharge), unit = Unit(; mm = 1, dt = -1)),
    "soil_water_saturated_zone_bottom__leakage_volume_flux" => (
        lens = @optic(_.land.soil.variables.actleakage),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_water_saturated_zone__depth" =>
        (lens = @optic(_.land.soil.variables.satwaterdepth), unit = Unit(; mm = 1)),
    "soil_water_saturated_zone_top__depth" =>
        (lens = @optic(_.land.soil.variables.zi), unit = Unit(; mm = 1)),
    "soil_surface__temperature" => (
        lens = @optic(_.land.soil.variables.tsoil),
        unit = Unit(; degC = 1, absolute_temperature = true),
    ),
    "subsurface_water_saturated_zone_top__depth" =>
        (lens = @optic(_.routing.subsurface_flow.variables.zi), unit = Unit(; m = 1)),
    "subsurface_water__exfiltration_volume_flux" => (
        lens = @optic(_.routing.subsurface_flow.variables.exfiltwater),
        unit = Unit(; m = 1, dt = -1),
    ),
    "subsurface_water__volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.variables.ssf),
        unit = Unit(; m = 3, d = -1),
    ),
    "subsurface_water__to_river_volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.variables.to_river),
        unit = Unit(; m = 3, d = -1),
    ),
    "subsurface_water_mass_balance_error__volume_flow_rate" => (
        lens = @optic(_.mass_balance.routing.subsurface_water_balance.error),
        unit = Unit(; m = 3, d = -1),
    ),
    "subsurface_water_mass_balance_relative_error__number" => (
        lens = @optic(_.mass_balance.routing.subsurface_water_balance.relative_error),
        unit = Unit(),
    ),
    "snowpack_liquid_water__depth" =>
        (lens = @optic(_.land.snow.variables.snow_water), unit = Unit(; mm = 1)),
    "snowpack_dry_snow__leq_depth" =>
        (lens = @optic(_.land.snow.variables.snow_storage), unit = Unit(; mm = 1)),
    "glacier_ice__leq_depth" =>
        (lens = @optic(_.land.glacier.variables.glacier_store), unit = Unit(; mm = 1)),
    "glacier_ice__melt_volume_flux" => (
        lens = @optic(_.land.glacier.variables.glacier_melt),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "glacier_ice__initial_leq_depth" =>
        (lens = @optic(_.land.glacier.variables.glacier_store), unit = Unit(; mm = 1)),
    "glacier_ice__melting_temperature_threshold" => (
        lens = @optic(_.land.glacier.parameters.g_ttm),
        unit = Unit(; degC = 1, absolute_temperature = true),
    ),
    "glacier_ice__degree_day_coefficient" => (
        lens = @optic(_.land.glacier.parameters.g_cfmax),
        unit = Unit(; mm = 1, degC = -1, dt = -1),
    ),
    "glacier_firn_accumulation__snowpack_dry_snow_leq_depth_fraction" =>
        (lens = @optic(_.land.glacier.parameters.g_sifrac), unit = Unit(; dt = -1)),
    "glacier_surface__area_fraction" =>
        (lens = @optic(_.land.glacier.parameters.glacier_frac), unit = Unit()),
    "land_surface_water__instantaneous_volume_flow_rate" => (
        lens = @optic(_.routing.overland_flow.variables.q),
        unit = Unit(; m = 3, s = -1),
    ),
    "land_surface_water__volume_flow_rate" => (
        lens = @optic(_.routing.overland_flow.variables.q_av),
        unit = Unit(; m = 3, s = -1),
    ),
    "land_surface_water__to_river_volume_flow_rate" => (
        lens = @optic(_.routing.overland_flow.variables.to_river),
        unit = Unit(; m = 3, s = -1),
    ),
    "land_surface_water__depth" =>
        (lens = @optic(_.routing.overland_flow.variables.h), unit = Unit(; m = 1)),
    "land_surface_water__volume" => (
        lens = @optic(_.routing.overland_flow.variables.storage),
        unit = Unit(; m = 3),
    ),
    "land_surface_water__x_component_of_instantaneous_volume_flow_rate" => (
        lens = @optic(_.routing.overland_flow.variables.qx),
        unit = Unit(; m = 3, s = -1),
    ),
    "land_surface_water__y_component_of_instantaneous_volume_flow_rate" => (
        lens = @optic(_.routing.overland_flow.variables.qy),
        unit = Unit(; m = 3, s = -1),
    ),
    "land_surface_water_mass_balance_error__volume_flow_rate" => (
        lens = @optic(_.mass_balance.routing.overland_water_balance.error),
        unit = Unit(; m = 3, s = -1),
    ),
    "land_surface_water_mass_balance_relative_error__number" => (
        lens = @optic(_.mass_balance.routing.overland_water_balance.relative_error),
        unit = Unit(),
    ),
    "paddy_surface_water__depth" =>
        (lens = @optic(_.land.demand.paddy.variables.h), unit = Unit(; mm = 1)),
    "domestic__gross_water_demand_volume_flux" => (
        lens = @optic(_.land.demand.domestic.demand.demand_gross),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "domestic__net_water_demand_volume_flux" => (
        lens = @optic(_.land.demand.domestic.demand.demand_net),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "industry__gross_water_demand_volume_flux" => (
        lens = @optic(_.land.demand.industry.demand.demand_gross),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "industry__net_water_demand_volume_flux" => (
        lens = @optic(_.land.demand.industry.demand.demand_net),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "livestock__gross_water_demand_volume_flux" => (
        lens = @optic(_.land.demand.livestock.demand.demand_gross),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "livestock__net_water_demand_volume_flux" => (
        lens = @optic(_.land.demand.livestock.demand.demand_net),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "irrigated_paddy__min_depth" =>
        (lens = @optic(_.land.demand.paddy.parameters.h_min), unit = Unit(; mm = 1)),
    "irrigated_paddy__optimal_depth" =>
        (lens = @optic(_.land.demand.paddy.parameters.h_opt), unit = Unit(; mm = 1)),
    "irrigated_paddy__max_depth" =>
        (lens = @optic(_.land.demand.paddy.parameters.h_max), unit = Unit(; mm = 1)),
    "irrigated_paddy__irrigation_efficiency" => (
        lens = @optic(_.land.demand.paddy.parameters.irrigation_efficiency),
        unit = Unit(),
    ),
    "irrigated_paddy_area__count" =>
        (lens = @optic(_.land.demand.paddy.parameters.irrigation_areas), unit = Unit()),
    "irrigated_paddy__irrigation_trigger_flag" => (
        lens = @optic(_.land.demand.paddy.parameters.irrigation_trigger),
        unit = Unit(),
    ),
    "irrigated_paddy__max_irrigation_rate" => (
        lens = @optic(_.land.demand.paddy.parameters.maximum_irrigation_rate),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "irrigated_paddy__gross_water_demand_volume_flux" => (
        lens = @optic(_.land.demand.paddy.variables.demand_gross),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "irrigated_non_paddy__irrigation_efficiency" => (
        lens = @optic(_.land.demand.nonpaddy.parameters.irrigation_efficiency),
        unit = Unit(),
    ),
    "irrigated_non_paddy_area__count" => (
        lens = @optic(_.land.demand.nonpaddy.parameters.irrigation_areas),
        unit = Unit(),
    ),
    "irrigated_non_paddy__irrigation_trigger_flag" => (
        lens = @optic(_.land.demand.nonpaddy.parameters.irrigation_trigger),
        unit = Unit(),
    ),
    "irrigated_non_paddy__max_irrigation_rate" => (
        lens = @optic(_.land.demand.nonpaddy.parameters.maximum_irrigation_rate),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "irrigated_non_paddy__gross_water_demand_volume_flux" => (
        lens = @optic(_.land.demand.nonpaddy.variables.demand_gross),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "land__allocated_irrigation_water_volume_flux" => (
        lens = @optic(_.land.allocation.variables.irri_alloc),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "land__allocated_non_irrigation_water_volume_flux" => (
        lens = @optic(_.land.allocation.variables.nonirri_alloc),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "subsurface_water__abstraction_volume_flux" => (
        lens = @optic(_.land.allocation.variables.act_groundwater_abst),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "land__non_irrigation_return_flow_volume_flux" => (
        lens = @optic(_.land.allocation.variables.nonirri_returnflow),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "subsurface_water__hydraulic_head" => (
        lens = @optic(_.routing.subsurface_flow.aquifer.variables.head),
        unit = Unit(; m = 1),
    ),
    "subsurface_water_saturated_zone_top__net_recharge_volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.recharge.variables.flux_av),
        unit = Unit(; m = 3, d = -1),
    ),
    "land_drain_water__to_subsurface_volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.variables.flux_av),
        unit = Unit(; m = 3, d = -1),
    ),
    "river_water__to_subsurface_volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.river.variables.flux_av),
        unit = Unit(; m = 3, d = -1),
    ),
    "land_surface_water__withdrawal_fraction" =>
        (lens = @optic(_.land.allocation.parameters.frac_sw_used), unit = Unit()),
    "land_water_allocation_area__count" =>
        (lens = @optic(_.land.allocation.parameters.areas), unit = Unit()),
)
