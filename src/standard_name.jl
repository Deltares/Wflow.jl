
const standard_name_map = Dict{String, ComposedFunction}(
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