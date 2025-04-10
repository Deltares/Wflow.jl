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
    "land_water~storage~total__depth" =>
        (lens = @optic(_.land.soil.variables.total_storage), unit = "mm"),
    "atmosphere_air__temperature" =>
        (lens = @optic(_.land.atmospheric_forcing.temperature), unit = "°C"),
    "vegetation__leaf-area_index" =>
        (lens = @optic(_.land.vegetation_parameters.leaf_area_index), unit = "m2 m-2"),
    "vegetation_canopy_water__depth" =>
        (lens = @optic(_.land.interception.variables.canopy_storage), unit = "mm"),
    "vegetation_canopy_water__stemflow_volume_flux" =>
        (lens = @optic(_.land.interception.variables.stemflow), unit = "mm dt-1"),
    "vegetation_canopy_water__throughfall_volume_flux" =>
        (lens = @optic(_.land.interception.variables.throughfall), unit = "mm dt-1"),
    "snowpack__leq-depth" => (lens = @optic(_.land.snow.variables.swe), unit = "mm"),
    "snowpack_meltwater__volume_flux" =>
        (lens = @optic(_.land.snow.variables.snow_melt), unit = "mm dt-1"),
    "snowpack_water__runoff_volume_flux" =>
        (lens = @optic(_.land.snow.variables.runoff), unit = "mm dt-1"),
    "river_water_inflow~external__volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.inflow),
        unit = "m3 s-1",
    ),
    "river_water_inflow~lateral__volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.inwater),
        unit = "m3 s-1",
    ),
    "river_water__instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.river_flow.variables.q), unit = "m3 s-1"),
    "river_water__volume_flow_rate" =>
        (lens = @optic(_.routing.river_flow.variables.q_av), unit = "m3 s-1"),
    "river_water__instantaneous_depth" =>
        (lens = @optic(_.routing.river_flow.variables.h), unit = "m"),
    "river_water__depth" =>
        (lens = @optic(_.routing.river_flow.variables.h_av), unit = "m"),
    "river_water__volume" =>
        (lens = @optic(_.routing.river_flow.variables.storage_av), unit = "m3"),
    "land_surface_water_abstraction__volume_flux" => (
        lens = @optic(_.routing.river_flow.allocation.variables.act_surfacewater_abst),
        unit = "mm dt-1",
    ),
    "floodplain_water__volume" => (
        lens = @optic(_.routing.river_flow.floodplain.variables.storage_av),
        unit = "m3",
    ),
    "floodplain_water__instantaneous_depth" =>
        (lens = @optic(_.routing.river_flow.floodplain.variables.h), unit = "m"),
    "floodplain_water__instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.river_flow.floodplain.variables.q), unit = "m3 s-1"),
    "floodplain_water__volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.floodplain.variables.q_av),
        unit = "m3 s-1",
    ),
    "reservoir_water__instantaneous_volume" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.storage),
        unit = "m3",
    ),
    "reservoir_water__volume" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.storage_av
        ),
        unit = "m3",
    ),
    "reservoir_water~outgoing__volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.outflow_av
        ),
        unit = "m3 s-1",
    ),
    "reservoir_water~incoming__volume_flow_rate" => (
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
    "lake_water__volume" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.lake.variables.storage_av),
        unit = "m3",
    ),
    "lake_water_surface__instantaneous_elevation" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.lake.variables.waterlevel),
        unit = "m",
    ),
    "lake_water_surface__elevation" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.lake.variables.waterlevel_av
        ),
        unit = "m",
    ),
    "lake_water~outgoing__volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.lake.variables.outflow_av),
        unit = "m3 s-1",
    ),
    "lake_water~incoming__volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.lake.boundary_conditions.inflow
        ),
        unit = "m3 s-1",
    ),
    "lake_water__evaporation_volume_flux" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.lake.variables.actevap),
        unit = "mm dt-1",
    ),
    "lake_water__precipitation_volume_flux" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.lake.boundary_conditions.precipitation
        ),
        unit = "mm dt-1",
    ),
    "lake_water__potential_evaporation_volume_flux" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.lake.boundary_conditions.evaporation
        ),
        unit = "mm dt-1",
    ),
    "soil_water__infiltration_volume_flux" =>
        (lens = @optic(_.land.soil.variables.actinfilt), unit = "mm dt-1"),
    "soil_water__transpiration_volume_flux" =>
        (lens = @optic(_.land.soil.variables.transpiration), unit = "mm dt-1"),
    "soil_surface_water__runoff_volume_flux" =>
        (lens = @optic(_.land.soil.variables.runoff), unit = "mm dt-1"),
    "soil_surface_water__net_runoff_volume_flux" =>
        (lens = @optic(_.land.soil.variables.net_runoff), unit = "mm dt-1"),
    "soil_layer_water__volume_fraction" =>
        (lens = @optic(_.land.soil.variables.vwc), unit = "m3 m-3"),
    "soil_layer_water__volume_percentage" =>
        (lens = @optic(_.land.soil.variables.vwc_perc), unit = "%"),
    "soil_water_root-zone__volume_fraction" =>
        (lens = @optic(_.land.soil.variables.vwc_root), unit = "m3 m-3"),
    "soil_water_root-zone__volume_percentage" =>
        (lens = @optic(_.land.soil.variables.vwc_percroot), unit = "%"),
    "soil_water_root-zone__depth" =>
        (lens = @optic(_.land.soil.variables.rootstore), unit = "mm"),
    "soil_layer_water_unsat-zone__depth" =>
        (lens = @optic(_.land.soil.variables.ustorelayerdepth), unit = "mm"),
    "soil_water_unsat-zone__depth" =>
        (lens = @optic(_.land.soil.variables.ustoredepth), unit = "mm"),
    "soil_water_sat-zone_top__capillary_volume_flux" =>
        (lens = @optic(_.land.soil.variables.actcapflux), unit = "mm dt-1"),
    "soil_water_sat-zone_top__recharge_volume_flux" =>
        (lens = @optic(_.land.soil.variables.transfer), unit = "mm dt-1"),
    "soil_water_sat-zone_top__net_recharge_volume_flux" =>
        (lens = @optic(_.land.soil.variables.recharge), unit = "mm dt-1"),
    "soil_water_sat-zone_bottom__leakage_volume_flux" =>
        (lens = @optic(_.land.soil.variables.actleakage), unit = "mm dt-1"),
    "soil_water_sat-zone__depth" =>
        (lens = @optic(_.land.soil.variables.satwaterdepth), unit = "mm"),
    "soil_water_sat-zone_top__depth" =>
        (lens = @optic(_.land.soil.variables.zi), unit = "mm"),
    "soil_surface__temperature" =>
        (lens = @optic(_.land.soil.variables.tsoil), unit = "ᵒC"),
    "subsurface_water_sat-zone_top__depth" =>
        (lens = @optic(_.routing.subsurface_flow.variables.zi), unit = "m"),
    "subsurface_water__exfiltration_volume_flux" => (
        lens = @optic(_.routing.subsurface_flow.variables.exfiltwater),
        unit = "m dt-1",
    ),
    "subsurface_water__volume_flow_rate" =>
        (lens = @optic(_.routing.subsurface_flow.variables.ssf), unit = "m3 d-1"),
    "subsurface_water~to-river__volume_flow_rate" =>
        (lens = @optic(_.routing.subsurface_flow.variables.to_river), unit = "m3 d-1"),
    "snowpack~liquid__depth" =>
        (lens = @optic(_.land.snow.variables.snow_water), unit = "mm"),
    "snowpack~dry__leq-depth" =>
        (lens = @optic(_.land.snow.variables.snow_storage), unit = "mm"),
    "glacier_ice__leq-depth" =>
        (lens = @optic(_.land.glacier.variables.glacier_store), unit = "mm"),
    "glacier_ice__melt_volume_flux" =>
        (lens = @optic(_.land.glacier.variables.glacier_melt), unit = "mm dt-1"),
    "land_surface_water__instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.overland_flow.variables.q), unit = "m3 s-1"),
    "land_surface_water__volume_flow_rate" =>
        (lens = @optic(_.routing.overland_flow.variables.q_av), unit = "m3 s-1"),
    "land_surface_water__instantaneous_depth" =>
        (lens = @optic(_.routing.overland_flow.variables.h), unit = "m"),
    "land_surface_water__depth" =>
        (lens = @optic(_.routing.overland_flow.variables.h_av), unit = "m"),
    "land_surface_water__volume" =>
        (lens = @optic(_.routing.overland_flow.variables.storage_av), unit = "m3"),
    "land_surface_water__x_component_of_instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.overland_flow.variables.qx), unit = "m3 s-1"),
    "land_surface_water__y_component_of_instantaneous_volume_flow_rate" =>
        (lens = @optic(_.routing.overland_flow.variables.qy), unit = "m3 s-1"),
    "land_surface_water~paddy__depth" =>
        (lens = @optic(_.land.demand.paddy.variables.h), unit = "mm"),
    "land~domestic__gross_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.domestic.demand.demand_gross), unit = "mm dt-1"),
    "land~domestic__net_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.domestic.demand.demand_net), unit = "mm dt-1"),
    "land~industry__gross_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.industry.demand.demand_gross), unit = "mm dt-1"),
    "land~industry__net_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.industry.demand.demand_net), unit = "mm dt-1"),
    "land~livestock__gross_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.livestock.demand.demand_gross), unit = "mm dt-1"),
    "land~livestock__net_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.livestock.demand.demand_net), unit = "mm dt-1"),
    "land~irrigated-paddy__irrigation_trigger_flag" => (
        lens = @optic(_.land.demand.paddy.parameters.irrigation_trigger),
        unit = "mm dt-1",
    ),
    "land~irrigated-paddy__gross_water_demand_volume_flux" =>
        (lens = @optic(_.land.demand.paddy.variables.demand_gross), unit = "mm dt-1"),
    "land~irrigated-non-paddy__irrigation_trigger_flag" => (
        lens = @optic(_.land.demand.nonpaddy.parameters.irrigation_trigger),
        unit = "mm dt-1",
    ),
    "land~irrigated-non-paddy__gross_water_demand_volume_flux" => (
        lens = @optic(_.land.demand.nonpaddy.variables.demand_gross),
        unit = "mm dt-1",
    ),
    "land~irrigated__allocated_water_volume_flux" =>
        (lens = @optic(_.land.allocation.variables.irri_alloc), unit = "mm dt-1"),
    "land~non-irrigated__allocated_water_volume_flux" =>
        (lens = @optic(_.land.allocation.variables.nonirri_alloc), unit = "mm dt-1"),
    "subsurface_water_abstraction__volume_flux" => (
        lens = @optic(_.land.allocation.variables.act_groundwater_abst),
        unit = "mm dt-1",
    ),
    "subsurface_water__hydraulic_head" =>
        (lens = @optic(_.routing.subsurface_flow.aquifer.variables.head), unit = "m"),
    "subsurface_water_sat-zone_top__net_recharge_volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.recharge.variables.flux),
        unit = "m3 d-1",
    ),
    "land_drain_water~to-subsurface__volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.variables.flux),
        unit = "mm dt-1",
    ),
    "river_water~to-subsurface__volume_flow_rate" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.river.variables.flux),
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
    "soil_erosion~rainfall__mass_flow_rate" =>
        (lens = @optic(_.land.rainfall_erosion.variables.amount), unit = "t dt-1"),
    "soil_erosion~overland_flow__mass_flow_rate" =>
        (lens = @optic(_.land.overland_flow_erosion.variables.amount), unit = "t dt-1"),
    "soil_erosion__mass_flow_rate" =>
        (lens = @optic(_.land.soil_erosion.variables.amount), unit = "t dt-1"),
    "soil_erosion_clay__mass_flow_rate" =>
        (lens = @optic(_.land.soil_erosion.variables.clay), unit = "t dt-1"),
    "soil_erosion_silt__mass_flow_rate" =>
        (lens = @optic(_.land.soil_erosion.variables.silt), unit = "t dt-1"),
    "soil_erosion_sand__mass_flow_rate" =>
        (lens = @optic(_.land.soil_erosion.variables.sand), unit = "t dt-1"),
    "soil_erosion_aggregates~small__mass_flow_rate" =>
        (lens = @optic(_.land.soil_erosion.variables.sagg), unit = "t dt-1"),
    "soil_erosion_aggregates~large__mass_flow_rate" =>
        (lens = @optic(_.land.soil_erosion.variables.lagg), unit = "t dt-1"),
    "land_surface_water_sediment_transport_capacity__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.transport_capacity.variables.amount),
        unit = "t dt-1",
    ),
    "land_surface_water_clay_transport_capacity__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.transport_capacity.variables.clay),
        unit = "t dt-1",
    ),
    "land_surface_water_sediment~to-river__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.amount),
        unit = "t dt-1",
    ),
    "land_surface_water_clay~to-river__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.clay),
        unit = "t dt-1",
    ),
    "land_surface_water_silt~to-river__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.silt),
        unit = "t dt-1",
    ),
    "land_surface_water_sand~to-river__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.sand),
        unit = "t dt-1",
    ),
    "land_surface_water_aggregates~small~to-river__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.sagg),
        unit = "t dt-1",
    ),
    "land_surface_water_aggregates~large~to-river__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.lagg),
        unit = "t dt-1",
    ),
    "land_surface_water_sediment__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.sediment_flux.variables.amount),
        unit = "t dt-1",
    ),
    "river_water_sediment~bedload__mass_concentration" => (
        lens = @optic(_.routing.river_flow.concentrations.variables.bed),
        unit = "g m-3",
    ),
    "river_water_sediment~suspended__mass_concentration" => (
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
    "river_water_aggregates~large__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_lagg),
        unit = "t",
    ),
    "river_bed_aggregates~large__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_lagg),
        unit = "t",
    ),
    "river_water_clay__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.clay),
        unit = "t dt-1",
    ),
    "river_water_gravel__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.gravel),
        unit = "t dt-1",
    ),
    "river_water_aggregates~large__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.lagg),
        unit = "t dt-1",
    ),
    "river_water_aggregates~small__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.sagg),
        unit = "t dt-1",
    ),
    "river_water_sand__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.sand),
        unit = "t dt-1",
    ),
    "river_water_silt__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.silt),
        unit = "t dt-1",
    ),
    "river_water_aggregates~small__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_sagg),
        unit = "t",
    ),
    "river_bed_aggregates~small__mass" => (
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
standard_name_map(model::LandHydrologySBM) = sbm_standard_name_map
standard_name_map(model::SoilLoss) = sediment_standard_name_map
