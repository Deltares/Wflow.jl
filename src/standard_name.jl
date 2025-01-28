const sbm_standard_name_map = Dict{String, ComposedFunction}(
    "atmosphere_water__precipitation_volume_flux" =>
        @optic(_.land.atmospheric_forcing.precipitation),
    "land_surface_water__potential_evaporation_volume_flux" =>
        @optic(_.land.atmospheric_forcing.potential_evaporation),
    "land_surface__evapotranspiration_volume_flux" =>
        @optic(_.land.soil.variables.actevap),
    "land_water__total_storage" => @optic(_.land.soil.variables.total_storage),
    "atmosphere_air__temperature" => @optic(_.land.atmospheric_forcing.temperature),
    "vegetation__leaf-area_index" =>
        @optic(_.land.vegetation_parameter_set.leaf_area_index),
    "vegetation_canopy_water__storage" =>
        @optic(_.land.interception.variables.canopy_storage),
    "vegetation_canopy_water__stemflow_volume_flux" =>
        @optic(_.land.interception.variables.stemflow),
    "vegetation_canopy_water__throughfall_volume_flux" =>
        @optic(_.land.interception.variables.throughfall),
    "snowpack__liquid-equivalent_depth" => @optic(_.land.snow.variables.swe),
    "snowpack_meltwater__volume_flux" => @optic(_.land.snow.variables.snow_melt),
    "snowpack_water__runoff_volume_flux" => @optic(_.land.snow.variables.runoff),
    "river_water_inflow~external__volume_flow_rate" =>
        @optic(_.routing.river_flow.boundary_conditions.inflow),
    "river_water_inflow~lateral__volume_flow_rate" =>
        @optic(_.routing.river_flow.boundary_conditions.inwater),
    "river_water__instantaneous_volume_flow_rate" =>
        @optic(_.routing.river_flow.variables.q),
    "river_water__volume_flow_rate" => @optic(_.routing.river_flow.variables.q_av),
    "river_water__instantaneous_depth" => @optic(_.routing.river_flow.variables.h),
    "river_water__depth" => @optic(_.routing.river_flow.variables.h_av),
    "river_water__volume" => @optic(_.routing.river_flow.variables.volume),
    "land_surface_water_abstraction__volume_flux" =>
        @optic(_.routing.river_flow.allocation.variables.act_surfacewater_abst),
    "floodplain_water__volume" =>
        @optic(_.routing.river_flow.floodplain.variables.volume),
    "floodplain_water__instantaneous_depth" =>
        @optic(_.routing.river_flow.floodplain.variables.h),
    "floodplain_water__instantaneous_volume_flow_rate" =>
        @optic(_.routing.river_flow.floodplain.variables.q),
    "floodplain_water__volume_flow_rate" =>
        @optic(_.routing.river_flow.floodplain.variables.q_av),
    "reservoir_water__instantaneous_volume" =>
        @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.storage),
    "reservoir_water__volume" =>
        @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.storage_av),
    "reservoir_water~outgoing__volume_flow_rate" =>
        @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.outflow_av),
    "lake_water__volume" =>
        @optic(_.routing.river_flow.boundary_conditions.lake.variables.storage_av),
    "lake_water_level__instantaneous_elevation" =>
        @optic(_.routing.river_flow.boundary_conditions.lake.variables.waterlevel),
    "lake_water_level__elevation" =>
        @optic(_.routing.river_flow.boundary_conditions.lake.variables.waterlevel_av),
    "lake_water~outgoing__volume_flow_rate" =>
        @optic(_.routing.river_flow.boundary_conditions.lake.variables.outflow_av),
    "soil_water__infiltration_volume_flux" => @optic(_.land.soil.variables.actinfilt),
    "soil_water__transpiration_volume_flux" =>
        @optic(_.land.soil.variables.transpiration),
    "soil_surface_water__runoff_volume_flux" => @optic(_.land.soil.variables.runoff),
    "soil_surface_water__net_runoff_volume_flux" =>
        @optic(_.land.soil.variables.net_runoff),
    "soil_layer_water__volume_fraction" => @optic(_.land.soil.variables.vwc),
    "soil_layer_water__volume_percentage" => @optic(_.land.soil.variables.vwc_perc),
    "soil_water_root-zone__volume_fraction" => @optic(_.land.soil.variables.vwc_root),
    "soil_water_root-zone__volume_percentage" =>
        @optic(_.land.soil.variables.vwc_percroot),
    "soil_water_sat-zone_top__recharge_volume_flux" =>
        @optic(_.land.soil.variables.recharge),
    "soil_layer_water_unsat-zone__depth" =>
        @optic(_.land.soil.variables.ustorelayerdepth),
    "soil_water_sat-zone_top__capillary_volume_flux" =>
        @optic(_.land.soil.variables.actcapflux),
    "soil_water_sat-zone_top__recharge_volume_flux" =>
        @optic(_.land.soil.variables.transfer),
    "soil_water_sat-zone_top__net_recharge_volume_flux" =>
        @optic(_.land.soil.variables.recharge),
    "soil_water_sat-zone_bottom__leakage_volume_flux" =>
        @optic(_.land.soil.variables.actleakage),
    "soil_water_sat-zone__depth" => @optic(_.land.soil.variables.satwaterdepth),
    "soil_water_sat-zone_top__depth" => @optic(_.land.soil.variables.zi),
    "soil_surface__temperature" => @optic(_.land.soil.variables.tsoil),
    "subsurface_water__volume_flow_rate" =>
        @optic(_.routing.subsurface_flow.variables.ssf),
    "snowpack~liquid__depth" => @optic(_.land.snow.variables.snow_water),
    "snowpack~dry__leq-depth" => @optic(_.land.snow.variables.snow_storage),
    "glacier_ice__leq-volume" => @optic(_.land.glacier.variables.glacier_store),
    "glacier_ice__melt_volume_flux" => @optic(_.land.glacier.variables.glacier_melt),
    "land_surface_water__instantaneous_volume_flow_rate" =>
        @optic(_.routing.overland_flow.variables.q),
    "land_surface_water__volume_flow_rate" =>
        @optic(_.routing.overland_flow.variables.q_av),
    "land_surface_water__instantaneous_depth" =>
        @optic(_.routing.overland_flow.variables.h),
    "land_surface_water__depth" => @optic(_.routing.overland_flow.variables.h_av),
    "land_surface_water__volume" => @optic(_.routing.overland_flow.variables.volume),
    "land_surface_water__x_component_of_instantaneous_volume_flow_rate" =>
        @optic(_.routing.overland_flow.variables.qx),
    "land_surface_water__y_component_of_instantaneous_volume_flow_rate" =>
        @optic(_.routing.overland_flow.variables.qy),
    "land_surface_water~paddy__depth" => @optic(_.land.demand.paddy.variables.h),
    "land~domestic__gross_water_demand_volume_flux" =>
        @optic(_.land.demand.domestic.demand.demand_gross),
    "land~domestic__net_water_demand_volume_flux" =>
        @optic(_.land.demand.domestic.demand.demand_net),
    "land~industry__gross_water_demand_volume_flux" =>
        @optic(_.land.demand.industry.demand.demand_gross),
    "land~industry__net_water_demand_volume_flux" =>
        @optic(_.land.demand.industry.demand.demand_net),
    "land~livestock__gross_water_demand_volume_flux" =>
        @optic(_.land.demand.livestock.demand.demand_gross),
    "land~livestock__net_water_demand_volume_flux" =>
        @optic(_.land.demand.livestock.demand.demand_net),
    "land~irrigated-paddy__irrigation_trigger_flag" =>
        @optic(_.land.demand.paddy.parameters.irrigation_trigger),
    "land~irrigated-paddy__gross_water_demand_volume_flux" =>
        @optic(_.land.demand.paddy.variables.demand_gross),
    "land~irrigated-non-paddy__irrigation_trigger_flag" =>
        @optic(_.land.demand.nonpaddy.parameters.irrigation_trigger),
    "land~irrigated-non-paddy__gross_water_demand_volume_flux" =>
        @optic(_.land.demand.nonpaddy.variables.demand_gross),
    "land~irrigated__allocated_water_volume_flux" =>
        @optic(_.land.allocation.variables.irri_alloc),
    "land~non-irrigated__allocated_water_volume_flux" =>
        @optic(_.land.allocation.variables.nonirri_alloc),
    "subsurface_water_abstraction__volume_flux" =>
        @optic(_.land.allocation.variables.act_groundwater_abst),
    "subsurface_water__hydraulic_head" =>
        @optic(_.routing.subsurface_flow.aquifer.variables.head),
    "subsurface_water_sat-zone_top__net_recharge_volume_flux" =>
        @optic(_.routing.subsurface_flow.boundaries.recharge.variables.flux),
    "land_drain_water~to-subsurface__volume_flow_rate" =>
        @optic(_.routing.subsurface_flow.boundaries.drain.variables.flux),
    "river_water~to-subsurface__volume_flow_rate" =>
        @optic(_.routing.subsurface_flow.boundaries.river.variables.flux),
)

const sediment_standard_name_map = Dict{String, ComposedFunction}(
    "atmosphere_water__precipitation_volume_flux" =>
        @optic(_.land.atmospheric_forcing.precipitation),
    "river_water__volume_flow_rate" =>
        @optic(_.routing.river_flow.hydrological_forcing.q_river),
    "river_water__depth" =>
        @optic(_.routing.river_flow.hydrological_forcing.waterlevel_river),
    "land_surface_water__volume_flow_rate" =>
        @optic(_.land.hydrological_forcing.q_land),
    "land_surface_water__depth" => @optic(_.land.hydrological_forcing.waterlevel_land),
    "vegetation_canopy_water__interception_volume_flux" =>
        @optic(_.land.hydrological_forcing.interception),
    "soil_erosion~rainfall__mass_flow_rate" =>
        @optic(_.land.rainfall_erosion.variables.amount),
    "soil_erosion~overland_flow__mass_flow_rate" =>
        @optic(_.land.overland_flow_erosion.variables.amount),
    "soil_erosion__mass_flow_rate" => @optic(_.land.soil_erosion.variables.amount),
    "soil_erosion_clay__mass_flow_rate" => @optic(_.land.soil_erosion.variables.clay),
    "soil_erosion_silt__mass_flow_rate" => @optic(_.land.soil_erosion.variables.silt),
    "soil_erosion_sand__mass_flow_rate" => @optic(_.land.soil_erosion.variables.sand),
    "soil_erosion_aggregates~small__mass_flow_rate" =>
        @optic(_.land.soil_erosion.variables.sagg),
    "soil_erosion_aggregates~large__mass_flow_rate" =>
        @optic(_.land.soil_erosion.variables.lagg),
    "land_surface_water_sediment_transport_capacity__mass_flow_rate" =>
        @optic(_.routing.overland_flow.transport_capacity.variables.amount),
    "land_surface_water_clay_transport_capacity__mass_flow_rate" =>
        @optic(_.routing.overland_flow.transport_capacity.variables.clay),
    "land_surface_water_sediment~to-river__mass_flow_rate" =>
        @optic(_.routing.overland_flow.to_river.variables.amount),
    "land_surface_water_clay~to-river__mass_flow_rate" =>
        @optic(_.routing.overland_flow.to_river.variables.clay),
    "land_surface_water_silt~to-river__mass_flow_rate" =>
        @optic(_.routing.overland_flow.to_river.variables.silt),
    "land_surface_water_sand~to-river__mass_flow_rate" =>
        @optic(_.routing.overland_flow.to_river.variables.sand),
    "land_surface_water_aggregates~small~to-river__mass_flow_rate" =>
        @optic(_.routing.overland_flow.to_river.variables.sagg),
    "land_surface_water_aggregates~large~to-river__mass_flow_rate" =>
        @optic(_.routing.overland_flow.to_river.variables.lagg),
    "land_surface_water_sediment__mass_flow_rate" =>
        @optic(_.routing.overland_flow.sediment_flux.variables.amount),
    "river_water_sediment~bedload__mass_concentration" =>
        @optic(_.routing.river_flow.concentrations.variables.bed),
    "river_water_sediment~suspended__mass_concentration" =>
        @optic(_.routing.river_flow.concentrations.variables.suspended),
    "river_water_sediment__mass_concentration" =>
        @optic(_.routing.river_flow.concentrations.variables.total),
    "river_water_clay__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.leftover_clay),
    "river_bed_clay__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.store_clay),
    "river_water_gravel__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.leftover_gravel),
    "river_bed_gravel__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.store_gravel),
    "river_water_aggregates~large__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.leftover_lagg),
    "river_bed_aggregates~large__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.store_lagg),
    "river_water_clay__mass_flow_rate" =>
        @optic(_.routing.river_flow.sediment_flux.variables.clay),
    "river_water_gravel__mass_flow_rate" =>
        @optic(_.routing.river_flow.sediment_flux.variables.gravel),
    "river_water_aggregates~large__mass_flow_rate" =>
        @optic(_.routing.river_flow.sediment_flux.variables.lagg),
    "river_water_aggregates~small__mass_flow_rate" =>
        @optic(_.routing.river_flow.sediment_flux.variables.sagg),
    "river_water_sand__mass_flow_rate" =>
        @optic(_.routing.river_flow.sediment_flux.variables.sand),
    "river_water_silt__mass_flow_rate" =>
        @optic(_.routing.river_flow.sediment_flux.variables.silt),
    "river_water_aggregates~small__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.leftover_sagg),
    "river_bed_aggregates~small__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.store_sagg),
    "river_water_sand__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.leftover_sand),
    "river_bed_sand__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.store_sand),
    "river_water_silt__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.leftover_silt),
    "river_bed_silt__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.store_silt),
    "river_water_sediment_erosion__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.erosion),
    "river_water_sediment_deposition__mass" =>
        @optic(_.routing.river_flow.sediment_flux.variables.deposition),
)

standard_name_map(model::LandHydrologySBM) = sbm_standard_name_map
standard_name_map(model::SoilLoss) = sediment_standard_name_map
