
const sbm_standard_name_map = Dict{String, ComposedFunction}(
    "atmosphere_water__precipitation_volume_flux" =>
        @optic(_.vertical.atmospheric_forcing.precipitation),
    "land_surface_water__potential_evaporation_volume_flux" =>
        @optic(_.vertical.atmospheric_forcing.potential_evaporation),
    "atmosphere_air__temperature" => @optic(_.vertical.atmospheric_forcing.temperature),
    "vegetation__leaf-area_index" =>
        @optic(_.vertical.vegetation_parameter_set.leaf_area_index),
    "vegetation_canopy_water__storage" =>
        @optic(_.vertical.interception.variables.canopy_storage),
    "river_water__volume_inflow_rate" =>
        @optic(_.lateral.river.boundary_conditions.inflow),
    "river_water__volume_flow_rate" => @optic(_.lateral.river.variables.q),
    "river_water__time_average_of_volume_flow_rate" =>
        @optic(_.lateral.river.variables.q_av),
    "river_water__depth" => @optic(_.lateral.river.variables.h),
    "river_water__time_average_of_depth" => @optic(_.lateral.river.variables.h_av),
    "river_water__volume" => @optic(_.lateral.river.variables.volume),
    "floodplain_water__volume" => @optic(_.lateral.river.floodplain.variables.volume),
    "floodplain_water__depth" => @optic(_.lateral.river.floodplain.variables.h),
    "floodplain_water__volume_flow_rate" =>
        @optic(_.lateral.river.floodplain.variables.q),
    "reservoir_water__volume" =>
        @optic(_.lateral.river.boundary_conditions.reservoir.variables.volume),
    "soil_water_sat-zone_top__recharge_volume_flux" =>
        @optic(_.vertical.soil.variables.recharge),
    "soil_water_unsat-zone__depth-per-soil_layer" =>
        @optic(_.vertical.soil.variables.ustorelayerdepth),
    "soil_water_sat-zone__depth" => @optic(_.vertical.soil.variables.satwaterdepth),
    "soil_water_sat-zone_top__depth" => @optic(_.vertical.soil.variables.zi),
    "soil_surface__temperature" => @optic(_.vertical.soil.variables.tsoil),
    "subsurface_water__volume_flow_rate" => @optic(_.lateral.subsurface.variables.ssf),
    "snowpack~liquid__depth" => @optic(_.vertical.snow.variables.snow_water),
    "snowpack~dry__leq-depth" => @optic(_.vertical.snow.variables.snow_storage),
    "glacier_ice__leq-volume" => @optic(_.vertical.glacier.variables.glacier_store),
    "land_surface_water__volume_flow_rate" => @optic(_.lateral.land.variables.q),
    "land_surface_water__depth" => @optic(_.lateral.land.variables.h),
    "land_surface_water__time_average_of_depth" =>
        @optic(_.lateral.land.variables.h_av),
    "land_surface_water__volume" => @optic(_.lateral.land.variables.volume),
    "land_surface_water__x_component_of_volume_flow_rate" =>
        @optic(_.lateral.land.variables.qx),
    "land_surface_water__y_component_of_volume_flow_rate" =>
        @optic(_.lateral.land.variables.qy),
    "land_surface_water~paddy__depth" => @optic(_.vertical.demand.paddy.variables.h),
    "land~domestic__gross_water_demand_flux" =>
        @optic(_.vertical.demand.domestic.demand.demand_gross),
    "land~domestic__net_water_demand_flux" =>
        @optic(_.vertical.demand.domestic.demand.demand_net),
    "land~industry__gross_water_demand_flux" =>
        @optic(_.vertical.demand.industry.demand.demand_gross),
    "land~industry__net_water_demand_flux" =>
        @optic(_.vertical.demand.industry.demand.demand_net),
    "land~livestock__gross_water_demand_flux" =>
        @optic(_.vertical.demand.livestock.demand.demand_gross),
    "land~livestock__net_water_demand_flux" =>
        @optic(_.vertical.demand.livestock.demand.demand_net),
    "land~irrigated-paddy__irrigation_trigger_flag" =>
        @optic(_.vertical.demand.paddy.parameters.irrigation_trigger),
    "land~irrigated-non-paddy__irrigation_trigger_flag" =>
        @optic(_.vertical.demand.nonpaddy.parameters.irrigation_trigger),
    "land~irrigated__allocated_water_volume_flux" =>
        @optic(_.vertical.allocation.variables.irri_alloc),
    "subsurface_water__hydraulic_head" =>
        @optic(_.lateral.subsurface.flow.aquifer.variables.head),
    "subsurface_water_sat-zone_top__net-recharge_volume_flux" =>
        @optic(_.lateral.subsurface.recharge.variables.flux),
    "land_drain_water~to-subsurface__volume_flow_rate" =>
        @optic(_.lateral.subsurface.drain.variables.flux),
    "river_water~to-subsurface__volume_flow_rate" =>
        @optic(_.lateral.subsurface.river.variables.flux),
)

const sediment_standard_name_map = Dict{String, ComposedFunction}(
    "atmosphere_water__precipitation_volume_flux" =>
        @optic(_.vertical.atmospheric_forcing.precipitation),
    "river_water__volume_flow_rate" =>
        @optic(_.lateral.river.hydrological_forcing.q_river),
    "river_water__depth" =>
        @optic(_.lateral.river.hydrological_forcing.waterlevel_river),
    "land_surface_water__volume_flow_rate" =>
        @optic(_.vertical.hydrological_forcing.q_land),
    "land_surface_water__depth" =>
        @optic(_.vertical.hydrological_forcing.waterlevel_land),
    "vegetation_canopy_water__interception_volume_flux" =>
        @optic(_.vertical.hydrological_forcing.interception),
    "soil_erosion~rainfall__mass_flow_rate" =>
        @optic(_.vertical.rainfall_erosion.variables.amount),
    "soil_erosion~overland_flow__mass_flow_rate" =>
        @optic(_.vertical.overland_flow_erosion.variables.amount),
    "soil_erosion__mass_flow_rate" => @optic(_.vertical.soil_erosion.variables.amount),
    "soil_erosion_clay__mass_flow_rate" =>
        @optic(_.vertical.soil_erosion.variables.clay),
    "land_surface_water_sediment_transport_capacity__mass_flow_rate" =>
        @optic(_.lateral.land.transport_capacity.variables.amount),
    "land_surface_water_clay_transport_capacity__mass_flow_rate" =>
        @optic(_.lateral.land.transport_capacity.variables.clay),
    "land_surface_water_sediment~to-river__mass_flow_rate" =>
        @optic(_.lateral.land.to_river.variables.amount),
    "land_surface_water_clay~to-river__mass_flow_rate" =>
        @optic(_.lateral.land.to_river.variables.clay),
    "land_surface_water_sediment__mass_flow_rate" =>
        @optic(_.lateral.land.sediment_flux.variables.amount),
    "land_surface_water_clay__mass_flow_rate" =>
        @optic(_.lateral.land.sediment_flux.variables.clay),
    "river_water_sediment~bedload__mass_concentration" =>
        @optic(_.lateral.river.concentrations.variables.bed),
    "river_water_sediment~suspended__mass_concentration" =>
        @optic(_.lateral.river.concentrations.variables.suspended),
    "river_water_sediment__mass_concentration" =>
        @optic(_.lateral.river.concentrations.variables.total),
    "river_water_clay__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.leftover_clay),
    "river_bed_clay__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.store_clay),
    "river_water_gravel__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.leftover_gravel),
    "river_bed_gravel__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.store_gravel),
    "river_water_aggregates~large__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.leftover_lagg),
    "river_bed_aggregates~large__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.store_lagg),
    "river_water_clay__mass_flow_rate" =>
        @optic(_.lateral.river.sediment_flux.variables.clay),
    "river_water_gravel__mass_flow_rate" =>
        @optic(_.lateral.river.sediment_flux.variables.gravel),
    "river_water_aggregates~large__mass_flow_rate" =>
        @optic(_.lateral.river.sediment_flux.variables.lagg),
    "river_water_aggregates~small__mass_flow_rate" =>
        @optic(_.lateral.river.sediment_flux.variables.sagg),
    "river_water_sand__mass_flow_rate" =>
        @optic(_.lateral.river.sediment_flux.variables.sand),
    "river_water_silt__mass_flow_rate" =>
        @optic(_.lateral.river.sediment_flux.variables.silt),
    "river_water_aggregates~small__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.leftover_sagg),
    "river_bed_aggregates~small__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.store_sagg),
    "river_water_sand__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.leftover_sand),
    "river_bed_sand__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.store_sand),
    "river_water_silt__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.leftover_silt),
    "river_bed_silt__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.store_silt),
    "river_water_sediment_erosion__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.erosion),
    "river_water_sediment_deposition__mass" =>
        @optic(_.lateral.river.sediment_flux.variables.deposition),
)

standard_name_map(model::LandHydrologySBM) = sbm_standard_name_map
standard_name_map(model::SoilLoss) = sediment_standard_name_map
