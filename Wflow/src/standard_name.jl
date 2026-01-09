"""
Mapping of (CSDMS) standard names to model variables and units for models with a land model
of type `LandHydrologySBM`. The `lens` of the NamedTuple allows access to a nested model
variable.
"""
const sbm_standard_name_map = Dict{String, NamedTuple}(
    "atmosphere_water__precipitation_volume_flux" =>
        (lens = @optic(_.land.atmospheric_forcing.precipitation), unit = "mm dt-1"),
    "land_surface_water__potential_evaporation_volume_flux" => (
        lens = @optic(_.land.atmospheric_forcing.potential_evaporation),
        unit = "mm dt-1",
    ),
    "land_surface__evapotranspiration_volume_flux" =>
        (lens = @optic(_.land.soil.variables.actevap), unit = "mm dt-1"),
    "land_water_storage__total_depth" =>
        (lens = @optic(_.land.soil.variables.total_storage), unit = "mm"),
    "land_water_mass_balance_error__volume_flux" =>
        (lens = @optic(_.mass_balance.land_water_balance.error), unit = "mm dt-1"),
    "land_water_mass_balance_relative_error__number" =>
        (lens = @optic(_.mass_balance.land_water_balance.relative_error), unit = "-"),
    "atmosphere_air__temperature" =>
        (lens = @optic(_.land.atmospheric_forcing.temperature), unit = "°C"),
    "vegetation__leaf_area_index" =>
        (lens = @optic(_.land.vegetation_parameters.leaf_area_index), unit = "m2 m-2"),
    "vegetation_canopy_water__depth" =>
        (lens = @optic(_.land.interception.variables.canopy_storage), unit = "mm"),
    "vegetation_canopy_water__stemflow_volume_flux" =>
        (lens = @optic(_.land.interception.variables.stemflow), unit = "mm dt-1"),
    "vegetation_canopy_water__throughfall_volume_flux" =>
        (lens = @optic(_.land.interception.variables.throughfall), unit = "mm dt-1"),
    "vegetation_canopy_water__interception_volume_flux" => (
        lens = @optic(_.land.interception.variables.interception_rate),
        unit = "mm dt-1",
    ),
    "snowpack__leq_depth" => (lens = @optic(_.land.snow.variables.swe), unit = "mm"),
    "snowpack_meltwater__volume_flux" =>
        (lens = @optic(_.land.snow.variables.snow_melt), unit = "mm dt-1"),
    "snowpack_water__runoff_volume_flux" =>
        (lens = @optic(_.land.snow.variables.runoff), unit = "mm dt-1"),
    "river_water__external_inflow_volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.external_inflow),
        unit = "m3 s-1",
    ),
    "river_water__external_abstraction_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.actual_external_abstraction_av
        ),
        unit = "m3 s-1",
    ),
    "river_water__lateral_inflow_volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.inwater),
        unit = "m3 s-1",
    ),
    "river_water__instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.river_flow.variables.q), unit = "m3 s-1"),
    "river_water__volume_flow_rate" =>
        (lens = @optic(_.routing.river_flow.variables.q_av), unit = "m3 s-1"),
    "river_water__depth" =>
        (lens = @optic(_.routing.river_flow.variables.h), unit = "m"),
    "river_water__volume" =>
        (lens = @optic(_.routing.river_flow.variables.storage), unit = "m3"),
    "river_water_mass_balance_error__volume_flow_rate" => (
        lens = @optic(_.mass_balance.routing.river_water_balance.error),
        unit = "m3 s-1",
    ),
    "river_water_mass_balance_relative_error__number" => (
        lens = @optic(_.mass_balance.routing.river_water_balance.relative_error),
        unit = "-",
    ),
    "land_surface_water__abstraction_volume_flux" => (
        lens = @optic(_.routing.river_flow.allocation.variables.act_surfacewater_abst),
        unit = "mm dt-1",
    ),
    "floodplain_water__volume" =>
        (lens = @optic(_.routing.river_flow.floodplain.variables.storage), unit = "m3"),
    "floodplain_water__depth" =>
        (lens = @optic(_.routing.river_flow.floodplain.variables.h), unit = "m"),
    "floodplain_water__instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.river_flow.floodplain.variables.q), unit = "m3 s-1"),
    "floodplain_water__volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.floodplain.variables.q_av),
        unit = "m3 s-1",
    ),
    "reservoir_water__target_min_volume_fraction" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.targetminfrac
        ),
        unit = "-",
    ),
    "reservoir_water__target_full_volume_fraction" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.targetfullfrac
        ),
        unit = "-",
    ),
    "reservoir_water_demand__required_downstream_volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.parameters.demand),
        unit = "m3 s-1",
    ),
    "reservoir_water_release_below_spillway__max_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.maxrelease
        ),
        unit = "m3 s-1",
    ),
    "reservoir_water__volume" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.storage),
        unit = "m3",
    ),
    "reservoir_water__outgoing_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.outflow_av
        ),
        unit = "m3 s-1",
    ),
    "reservoir_water__outgoing_observed_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.outflow_obs
        ),
        unit = "m3 s-1",
    ),
    "reservoir_water__incoming_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.inflow
        ),
        unit = "m3 s-1",
    ),
    "reservoir_water__evaporation_volume_flux" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.actevap),
        unit = "mm dt-1",
    ),
    "reservoir_water__precipitation_volume_flux" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation
        ),
        unit = "mm dt-1",
    ),
    "reservoir_water__potential_evaporation_volume_flux" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.evaporation
        ),
        unit = "mm dt-1",
    ),
    "reservoir_water_surface__elevation" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.waterlevel
        ),
        unit = "m",
    ),
    "reservoir_water__external_inflow_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.external_inflow
        ),
        unit = "m3 s-1",
    ),
    "reservoir_water_mass_balance_error__volume_flow_rate" => (
        lens = @optic(_.mass_balance.routing.reservoir_water_balance.error),
        unit = "m3 s-1",
    ),
    "reservoir_water_mass_balance_relative_error__number" => (
        lens = @optic(_.mass_balance.routing.reservoir_water_balance.relative_error),
        unit = "-",
    ),
    "soil_water__infiltration_volume_flux" =>
        (lens = @optic(_.land.soil.variables.actinfilt), unit = "mm dt-1"),
    "soil_water__transpiration_volume_flux" =>
        (lens = @optic(_.land.soil.variables.transpiration), unit = "mm dt-1"),
    "soil_surface_water__runoff_volume_flux" =>
        (lens = @optic(_.land.soil.variables.runoff), unit = "mm dt-1"),
    "soil_surface_water__net_runoff_volume_flux" =>
        (lens = @optic(_.land.soil.variables.net_runoff), unit = "mm dt-1"),
    "soil_surface_water_unsaturated_zone__exfiltration_volume_flux" =>
        (lens = @optic(_.land.soil.variables.exfiltustore), unit = "mm dt-1"),
    "soil_surface_water_saturated_zone__exfiltration_volume_flux" =>
        (lens = @optic(_.land.soil.variables.exfiltsatwater), unit = "mm dt-1"),
    "compacted_soil_surface_water__excess_volume_flux" =>
        (lens = @optic(_.land.soil.variables.excesswatersoil), unit = "mm dt-1"),
    "non_compacted_soil_surface_water__excess_volume_flux" =>
        (lens = @optic(_.land.soil.variables.excesswaterpath), unit = "mm dt-1"),
    "soil_layer_water__volume_fraction" =>
        (lens = @optic(_.land.soil.variables.vwc), unit = "m3 m-3"),
    "soil_layer_water__volume_percentage" =>
        (lens = @optic(_.land.soil.variables.vwc_perc), unit = "%"),
    "soil_water_root_zone__volume_fraction" =>
        (lens = @optic(_.land.soil.variables.vwc_root), unit = "m3 m-3"),
    "soil_water_root_zone__volume_percentage" =>
        (lens = @optic(_.land.soil.variables.vwc_percroot), unit = "%"),
    "soil_water_root_zone__depth" =>
        (lens = @optic(_.land.soil.variables.rootstore), unit = "mm"),
    "soil_layer_water_unsaturated_zone__depth" =>
        (lens = @optic(_.land.soil.variables.ustorelayerdepth), unit = "mm"),
    "soil_water_unsaturated_zone__depth" =>
        (lens = @optic(_.land.soil.variables.ustoredepth), unit = "mm"),
    "soil_water_saturated_zone_top__capillary_volume_flux" =>
        (lens = @optic(_.land.soil.variables.actcapflux), unit = "mm dt-1"),
    "soil_water_saturated_zone_top__recharge_volume_flux" =>
        (lens = @optic(_.land.soil.variables.transfer), unit = "mm dt-1"),
    "soil_water_saturated_zone_top__net_recharge_volume_flux" =>
        (lens = @optic(_.land.soil.variables.recharge), unit = "mm dt-1"),
    "soil_water_saturated_zone_bottom__leakage_volume_flux" =>
        (lens = @optic(_.land.soil.variables.actleakage), unit = "mm dt-1"),
    "soil_water_saturated_zone__depth" =>
        (lens = @optic(_.land.soil.variables.satwaterdepth), unit = "mm"),
    "soil_water_saturated_zone_top__depth" =>
        (lens = @optic(_.land.soil.variables.zi), unit = "mm"),
    "soil_surface__temperature" =>
        (lens = @optic(_.land.soil.variables.tsoil), unit = "ᵒC"),
    "subsurface_water_saturated_zone_top__depth" =>
        (lens = @optic(_.routing.subsurface_flow.variables.zi), unit = "m"),
    "subsurface_water__exfiltration_volume_flux" => (
        lens = @optic(_.routing.subsurface_flow.variables.exfiltwater),
        unit = "m dt-1",
    ),
    "subsurface_water__volume_flow_rate" =>
        (lens = @optic(_.routing.subsurface_flow.variables.ssf), unit = "m3 d-1"),
    "subsurface_water__to_river_volume_flow_rate" =>
        (lens = @optic(_.routing.subsurface_flow.variables.to_river), unit = "m3 d-1"),
    "subsurface_water_mass_balance_error__volume_flow_rate" => (
        lens = @optic(_.mass_balance.routing.subsurface_water_balance.error),
        unit = "m3 d-1",
    ),
    "subsurface_water_mass_balance_relative_error__number" => (
        lens = @optic(_.mass_balance.routing.subsurface_water_balance.relative_error),
        unit = "-",
    ),
    "snowpack_liquid_water__depth" =>
        (lens = @optic(_.land.snow.variables.snow_water), unit = "mm"),
    "snowpack_dry_snow__leq_depth" =>
        (lens = @optic(_.land.snow.variables.snow_storage), unit = "mm"),
    "glacier_ice__leq_depth" =>
        (lens = @optic(_.land.glacier.variables.glacier_store), unit = "mm"),
    "glacier_ice__melt_volume_flux" =>
        (lens = @optic(_.land.glacier.variables.glacier_melt), unit = "mm dt-1"),
    "land_surface_water__instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.overland_flow.variables.q), unit = "m3 s-1"),
    "land_surface_water__volume_flow_rate" =>
        (lens = @optic(_.routing.overland_flow.variables.q_av), unit = "m3 s-1"),
    "land_surface_water__to_river_volume_flow_rate" =>
        (lens = @optic(_.routing.overland_flow.variables.to_river), unit = "m3 s-1"),
    "land_surface_water__depth" =>
        (lens = @optic(_.routing.overland_flow.variables.h), unit = "m"),
    "land_surface_water__volume" =>
        (lens = @optic(_.routing.overland_flow.variables.storage), unit = "m3"),
    "land_surface_water__x_component_of_instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.overland_flow.variables.qx), unit = "m3 s-1"),
    "land_surface_water__y_component_of_instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.overland_flow.variables.qy), unit = "m3 s-1"),
    "land_surface_water_mass_balance_error__volume_flow_rate" => (
        lens = @optic(_.mass_balance.routing.overland_water_balance.error),
        unit = "m3 s-1",
    ),
    "land_surface_water_mass_balance_relative_error__number" => (
        lens = @optic(_.mass_balance.routing.overland_water_balance.relative_error),
        unit = "-",
    ),
    "paddy_surface_water__depth" =>
        (lens = @optic(_.land.demand.paddy.variables.h), unit = "mm"),
    "domestic__gross_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.domestic.demand.demand_gross), unit = "mm dt-1"),
    "domestic__net_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.domestic.demand.demand_net), unit = "mm dt-1"),
    "industry__gross_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.industry.demand.demand_gross), unit = "mm dt-1"),
    "industry__net_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.industry.demand.demand_net), unit = "mm dt-1"),
    "livestock__gross_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.livestock.demand.demand_gross), unit = "mm dt-1"),
    "livestock__net_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.livestock.demand.demand_net), unit = "mm dt-1"),
    "irrigated_paddy__irrigation_trigger_flag" => (
        lens = @optic(_.land.demand.paddy.parameters.irrigation_trigger),
        unit = "mm dt-1",
    ),
    "irrigated_paddy__gross_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.paddy.variables.demand_gross), unit = "mm dt-1"),
    "irrigated_non_paddy__irrigation_trigger_flag" => (
        lens = @optic(_.land.demand.nonpaddy.parameters.irrigation_trigger),
        unit = "mm dt-1",
    ),
    "irrigated_non_paddy__gross_water_demand_volume_flux" => (
        lens = @optic(_.land.demand.nonpaddy.variables.demand_gross),
        unit = "mm dt-1",
    ),
    "land__allocated_irrigation_water_volume_flux" =>
        (lens = @optic(_.land.allocation.variables.irri_alloc), unit = "mm dt-1"),
    "land__allocated_non_irrigation_water_volume_flux" =>
        (lens = @optic(_.land.allocation.variables.nonirri_alloc), unit = "mm dt-1"),
    "subsurface_water__abstraction_volume_flux" => (
        lens = @optic(_.land.allocation.variables.act_groundwater_abst),
        unit = "mm dt-1",
    ),
    "land__non_irrigation_return_flow_volume_flux" => (
        lens = @optic(_.land.allocation.variables.nonirri_returnflow),
        unit = "mm dt-1",
    ),
    "subsurface_water__hydraulic_head" =>
        (lens = @optic(_.routing.subsurface_flow.aquifer.variables.head), unit = "m"),
    "subsurface_water_saturated_zone_top__net_recharge_volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.recharge.variables.flux_av),
        unit = "m3 d-1",
    ),
    "land_drain_water__to_subsurface_volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.variables.flux_av),
        unit = "m3 d-1",
    ),
    "river_water__to_subsurface_volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.river.variables.flux_av),
        unit = "m3 d-1",
    ),
)

"""
Mapping of (CSDMS) standard names to model variables and units for models with a land model
of type `SoilLoss`. The `lens` of the NamedTuple allows access to a nested model variable.
"""
const sediment_standard_name_map = Dict{String, NamedTuple}(
    "atmosphere_water__precipitation_volume_flux" =>
        (lens = @optic(_.land.atmospheric_forcing.precipitation), unit = "mm dt-1"),
    "river_water__volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.hydrological_forcing.q_river),
        unit = "m3 s-1",
    ),
    "river_water__depth" => (
        lens = @optic(_.routing.river_flow.hydrological_forcing.waterlevel_river),
        unit = "m",
    ),
    "land_surface_water__volume_flow_rate" =>
        (lens = @optic(_.land.hydrological_forcing.q_land), unit = "m3 s-1"),
    "land_surface_water__depth" =>
        (lens = @optic(_.land.hydrological_forcing.waterlevel_land), unit = "m"),
    "vegetation_canopy_water__interception_volume_flux" =>
        (lens = @optic(_.land.hydrological_forcing.interception), unit = "mm dt-1"),
    "rainfall_soil_erosion__mass_flow_rate" => (
        lens = @optic(_.land.rainfall_erosion.variables.soil_erosion_rate),
        unit = "t dt-1",
    ),
    "overland_flow_soil_erosion__mass_flow_rate" => (
        lens = @optic(_.land.overland_flow_erosion.variables.soil_erosion_rate),
        unit = "t dt-1",
    ),
    "soil_erosion__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.soil_erosion_rate),
        unit = "t dt-1",
    ),
    "soil_erosion_clay__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.clay_erosion_rate),
        unit = "t dt-1",
    ),
    "soil_erosion_silt__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.silt_erosion_rate),
        unit = "t dt-1",
    ),
    "soil_erosion_sand__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.sand_erosion_rate),
        unit = "t dt-1",
    ),
    "soil_erosion_small_aggregates__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.sagg_erosion_rate),
        unit = "t dt-1",
    ),
    "soil_erosion_large_aggregates__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.lagg_erosion_rate),
        unit = "t dt-1",
    ),
    "land_surface_water_sediment_transport_capacity__mass_flow_rate" => (
        lens = @optic(
            _.routing.overland_flow.transport_capacity.variables.sediment_transport_capacity
        ),
        unit = "t dt-1",
    ),
    "land_surface_water_sediment__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.sediment_rate),
        unit = "t dt-1",
    ),
    "land_surface_water_clay__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.clay_rate),
        unit = "t dt-1",
    ),
    "land_surface_water_silt__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.silt_rate),
        unit = "t dt-1",
    ),
    "land_surface_water_sand__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.sand_rate),
        unit = "t dt-1",
    ),
    "land_surface_water_small_aggregates__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.sagg_rate),
        unit = "t dt-1",
    ),
    "land_surface_water_large_aggregates__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.lagg_rate),
        unit = "t dt-1",
    ),
    "land_surface_water_sediment__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.sediment_flux.variables.sediment_rate),
        unit = "t dt-1",
    ),
    "river_water_sediment__bedload_mass_concentration" => (
        lens = @optic(_.routing.river_flow.concentrations.variables.bed),
        unit = "g m-3",
    ),
    "river_water_sediment__suspended_mass_concentration" => (
        lens = @optic(_.routing.river_flow.concentrations.variables.suspended),
        unit = "g m-3",
    ),
    "river_water_sediment__mass_concentration" => (
        lens = @optic(_.routing.river_flow.concentrations.variables.total),
        unit = "g m-3",
    ),
    "river_water_clay__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_clay),
        unit = "t",
    ),
    "river_bed_clay__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_clay),
        unit = "t",
    ),
    "river_water_gravel__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_gravel),
        unit = "t",
    ),
    "river_bed_gravel__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_gravel),
        unit = "t",
    ),
    "river_water_large_aggregates__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_lagg),
        unit = "t",
    ),
    "river_bed_large_aggregates__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_lagg),
        unit = "t",
    ),
    "river_water_clay__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.clay_rate),
        unit = "t dt-1",
    ),
    "river_water_gravel__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.gravel_rate),
        unit = "t dt-1",
    ),
    "river_water_large_aggregates__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.lagg_rate),
        unit = "t dt-1",
    ),
    "river_water_small_aggregates__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.sagg_rate),
        unit = "t dt-1",
    ),
    "river_water_sand__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.sand_rate),
        unit = "t dt-1",
    ),
    "river_water_silt__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.silt_rate),
        unit = "t dt-1",
    ),
    "river_water_small_aggregates__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_sagg),
        unit = "t",
    ),
    "river_bed_small_aggregates__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_sagg),
        unit = "t",
    ),
    "river_water_sand__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_sand),
        unit = "t",
    ),
    "river_bed_sand__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_sand),
        unit = "t",
    ),
    "river_water_silt__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_silt),
        unit = "t",
    ),
    "river_bed_silt__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_silt),
        unit = "t",
    ),
    "river_water_sediment_erosion__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.erosion),
        unit = "t",
    ),
    "river_water_sediment_deposition__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.deposition),
        unit = "t",
    ),
)

# wrapper methods for standard name mapping
standard_name_map(::LandHydrologySBM) = sbm_standard_name_map
standard_name_map(::SoilLoss) = sediment_standard_name_map
get_lens(name::AbstractString, T::AbstractLandModel) = standard_name_map(T)[name].lens
