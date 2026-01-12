# NOTE: The order of the entries determines the order in the docs tables
"""
Mapping of (CSDMS) standard names to model variables and units for models with a land model
of type `SoilLoss`. The `lens` allows access to a nested model variable.
"""
const sediment_standard_name_map = OrderedDict{String, ParameterMetadata}(
    "atmosphere_water__precipitation_volume_flux" => ParameterMetadata(;
        lens = @optic(_.land.atmospheric_forcing.precipitation),
        unit = Unit(; mm = 1, dt = -1),
        description = "Precipitation",
        tags = [:atmospheric_forcing],
    ),
    "soil_clay__mass_fraction" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.parameters.clay_fraction),
        default = 0.4,
        description = "Soil content clay",
        tags = [:soil_erosion_input],
    ),
    "soil_silt__mass_fraction" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.parameters.silt_fraction),
        default = 0.3,
        description = "Soil content silt",
        tags = [:soil_erosion_input],
    ),
    "soil_sand__mass_fraction" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.parameters.sand_fraction),
        default = 0.3,
        description = "Soil content sand",
        tags = [:soil_erosion_input],
    ),
    "soil_small_aggregates__mass_fraction" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.parameters.sagg_fraction),
        default = 0.0,
        description = "Soil content small aggregates",
        tags = [:soil_erosion_input],
    ),
    "soil_large_aggregates__mass_fraction" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.parameters.lagg_fraction),
        default = 0.0,
        description = "Soil content large aggregates",
        tags = [:soil_erosion_input],
    ),
    "river_bottom_and_bank_sediment__median_diameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.potential_erosion.parameters.d50),
        unit = Unit(; mm = 1),
        default = 0.1,
        description = "Median diameter in the river bed/bank",
    ),
    "soil_erosion__answers_overland_flow_factor" => ParameterMetadata(;
        lens = @optic(_.land.overland_flow_erosion.parameters.answers_overland_flow_factor),
        default = 0.9,
        description = "Answers overland flow erosion factor",
        tags = [:overland_flow_erosion_input],
    ),
    "soil_erosion__rainfall_soil_detachability_factor" => ParameterMetadata(;
        lens = @optic(_.land.rainfall_erosion.parameters.soil_detachability),
        unit = Unit(; g = 1, J = -1),
        default = 0.6,
        description = "Soil detachability factor",
        tags = [:rainfall_erosion_input],
    ),
    "soil_erosion__eurosem_exponent" => ParameterMetadata(;
        lens = @optic(_.land.rainfall_erosion.parameters.eurosem_exponent),
        unit = Unit(; m = -1),
        default = 2.0,
        description = "Exponent EUROSEM",
        tags = [:rainfall_erosion_input],
    ),
    "vegetation_canopy__height" => ParameterMetadata(;
        lens = @optic(_.land.rainfall_erosion.parameters.canopyheight),
        unit = Unit(; m = 1),
        default = 0.5,
        description = "Canopy height",
        tags = [:rainfall_erosion_input],
    ),
    "vegetation_canopy__gap_fraction" => ParameterMetadata(;
        lens = @optic(_.land.rainfall_erosion.parameters.canopygapfraction),
        default = 0.1,
        description = "Canopy gap fraction",
        tags = [:rainfall_erosion_input],
    ),
    "compacted_soil__area_fraction" => ParameterMetadata(;
        lens = @optic(_.land.rainfall_erosion.parameters.soilcover_fraction),
        default = 0.01,
        description = "Fraction of the soil that is covered (eg paved, snow, etc)",
        tags = [:rainfall_erosion_input],
    ),
    "soil_erosion__usle_k_factor" => ParameterMetadata(;
        lens = @optic(_.land.overland_flow_erosion.parameters.usle_k),
        default = 0.1,
        description = "USLE soil erodibility factor",
        tags = [:overland_flow_erosion_input],
    ),
    "soil_erosion__usle_c_factor" => ParameterMetadata(;
        lens = @optic(_.land.overland_flow_erosion.parameters.usle_c),
        default = 0.01,
        description = "USLE crop management factor",
        tags = [:overland_flow_erosion_input],
    ),
    "soil_erosion__answers_rainfall_factor" => ParameterMetadata(;
        lens = @optic(_.land.rainfall_erosion.parameters.answers_rainfall_factor),
        default = 0.108,
        description = "Answers rainfall erosion factor",
        tags = [:rainfall_erosion_input],
    ),
    "river_bottom_and_bank_clay__mass_fraction" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.clay_fraction),
        default = 0.15,
        description = "River bed/bank content clay",
    ),
    "river_bottom_and_bank_silt__mass_fraction" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.silt_fraction),
        default = 0.65,
        description = "River bed/bank content silt",
    ),
    "river_bottom_and_bank_sand__mass_fraction" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.sand_fraction),
        default = 0.15,
        description = "River bed/bank content sand",
    ),
    "river_bottom_and_bank_gravel__mass_fraction" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.gravel_fraction),
        default = 0.05,
        description = "River bed/bank content gravel",
    ),
    "clay__mean_diameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_clay),
        unit = Unit(; μm = 1),
        default = 2.0,
        description = "Clay mean diameter",
    ),
    "silt__mean_diameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_silt),
        unit = Unit(; μm = 1),
        default = 10.0,
        description = "Silt mean diameter",
    ),
    "sand__mean_diameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_sand),
        unit = Unit(; μm = 1),
        default = 200.0,
        description = "Sand mean diameter",
    ),
    "sediment_small_aggregates__mean_diameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_sagg),
        unit = Unit(; μm = 1),
        default = 30.0,
        description = "Small aggregates mean diameter",
    ),
    "sediment_large_aggregates__mean_diameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_lagg),
        unit = Unit(; μm = 1),
        default = 500.0,
        description = "Large aggregates mean diameter",
    ),
    "gravel__mean_diameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_gravel),
        unit = Unit(; μm = 1),
        default = 2000.0,
        description = "Gravel mean diameter",
    ),
    "reservoir_location__count" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.reservoir_outlet),
        fill = 0,
        description = "Reservoir location ids",
    ),
    "reservoir_surface__area" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.reservoir_area),
        unit = Unit(; m = 2),
        fill = 0.0,
        description = "Reservoir surface area",
    ),
    "reservoir_water_sediment__bedload_trapping_efficiency" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.sediment_flux.parameters.reservoir_trapping_efficiency
        ),
        default = 1.0,
        fill = 0.0,
        description = "Reservoir sediment bedload trapping efficiency",
    ),
    "river_water__volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.hydrological_forcing.q_river),
        unit = Unit(; m = 3, s = -1),
        description = "River discharge",
        tags = [:hydrological_forcing],
    ),
    "river_water__depth" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.hydrological_forcing.waterlevel_river),
        unit = Unit(; m = 1),
        description = "River water depth",
        tags = [:hydrological_forcing],
    ),
    "land_surface_water__volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.land.hydrological_forcing.q_land),
        unit = Unit(; m = 3, s = -1),
        description = "Overland flow discharge",
        tags = [:hydrological_forcing],
    ),
    "land_surface_water__depth" => ParameterMetadata(;
        lens = @optic(_.land.hydrological_forcing.waterlevel_land),
        unit = Unit(; m = 1),
        description = "Overland flow water depth",
        tags = [:hydrological_forcing],
    ),
    "vegetation_canopy_water__interception_volume_flux" => ParameterMetadata(;
        lens = @optic(_.land.hydrological_forcing.interception),
        unit = Unit(; mm = 1, dt = -1),
        description = "Rainfall interception by the vegetation",
        tags = [:hydrological_forcing],
    ),
    "rainfall_soil_erosion__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.land.rainfall_erosion.variables.soil_erosion_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Total soil erosion from rainfall (splash)",
        tags = [:rainfall_erosion_output],
    ),
    "overland_flow_soil_erosion__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.land.overland_flow_erosion.variables.soil_erosion_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Total soil erosion from overland flow",
        tags = [:overland_flow_erosion_output],
    ),
    "soil_erosion__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.variables.soil_erosion_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Total soil erosion",
        tags = [:soil_erosion_output],
    ),
    "soil_erosion_clay__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.variables.clay_erosion_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Total clay erosion",
        tags = [:soil_erosion_output],
    ),
    "soil_erosion_silt__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.variables.silt_erosion_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Total silt erosion",
        tags = [:soil_erosion_output],
    ),
    "soil_erosion_sand__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.variables.sand_erosion_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Total sand erosion",
        tags = [:soil_erosion_output],
    ),
    "soil_erosion_small_aggregates__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.variables.sagg_erosion_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Total small aggregates erosion",
        tags = [:soil_erosion_output],
    ),
    "soil_erosion_large_aggregates__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.land.soil_erosion.variables.lagg_erosion_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Total large aggregates erosion",
        tags = [:soil_erosion_output],
    ),
    "land_surface_water_sediment_transport_capacity__mass_flow_rate" =>
        ParameterMetadata(;
            lens = @optic(
                _.routing.overland_flow.transport_capacity.variables.sediment_transport_capacity
            ),
            unit = Unit(; t = 1, dt = -1),
            description = "Total sediment transport capacity",
        ),
    "land_surface_water_sediment__to_river_mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.to_river.variables.sediment_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Total sediment flux flowing into the river",
    ),
    "land_surface_water_clay__to_river_mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.to_river.variables.clay_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Clay flux flowing into the river",
    ),
    "land_surface_water_silt__to_river_mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.to_river.variables.silt_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Silt flux flowing into the river",
    ),
    "land_surface_water_sand__to_river_mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.to_river.variables.sand_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Sand flux flowing into the river",
    ),
    "land_surface_water_small_aggregates__to_river_mass_flow_rate" =>
        ParameterMetadata(;
            lens = @optic(_.routing.overland_flow.to_river.variables.sagg_rate),
            unit = Unit(; t = 1, dt = -1),
            description = "Small aggregates flux flowing into the river",
        ),
    "land_surface_water_large_aggregates__to_river_mass_flow_rate" =>
        ParameterMetadata(;
            lens = @optic(_.routing.overland_flow.to_river.variables.lagg_rate),
            unit = Unit(; t = 1, dt = -1),
            description = "Large aggregates flux flowing into the river",
        ),
    "land_surface_water_sediment__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.sediment_flux.variables.sediment_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Overland flow sediment flux",
    ),
    "river_water_sediment__bedload_mass_concentration" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.concentrations.variables.bed),
        unit = Unit(; g = 1, m = -3),
        description = "Bedload sediment concentration in river",
    ),
    "river_water_sediment__suspended_mass_concentration" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.concentrations.variables.suspended),
        unit = Unit(; g = 1, m = -3),
        description = "Suspended sediment concentration in river",
    ),
    "river_water_sediment__mass_concentration" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.concentrations.variables.total),
        unit = Unit(; g = 1, m = -3),
        description = "Total sediment concentration in river",
    ),
    "river_water_clay__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_clay),
        unit = Unit(; t = 1),
        description = "Clay mass in river water",
    ),
    "river_bed_clay__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_clay),
        unit = Unit(; t = 1),
        description = "Clay mass in river bed",
    ),
    "river_water_gravel__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_gravel),
        unit = Unit(; t = 1),
        description = "Gravel mass in river water",
    ),
    "river_bed_gravel__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_gravel),
        unit = Unit(; t = 1),
        description = "Gravel mass in river bed",
    ),
    "river_water_large_aggregates__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_lagg),
        unit = Unit(; t = 1),
        description = "Large aggregates mass in river water",
    ),
    "river_bed_large_aggregates__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_lagg),
        unit = Unit(; t = 1),
        description = "Large aggregates mass in river bed",
    ),
    "river_water_clay__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.clay_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Clay mass flow rate in river",
    ),
    "river_water_gravel__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.gravel_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Gravel mass flow rate in river",
    ),
    "river_water_large_aggregates__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.lagg_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Large aggregates mass flow rate in river",
    ),
    "river_water_small_aggregates__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.sagg_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Small aggregates mass flow rate in river",
    ),
    "river_water_sand__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.sand_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Sand mass flow rate in river",
    ),
    "river_water_silt__mass_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.silt_rate),
        unit = Unit(; t = 1, dt = -1),
        description = "Silt mass flow rate in river",
    ),
    "river_water_small_aggregates__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_sagg),
        unit = Unit(; t = 1),
        description = "Small aggregates mass in river water",
    ),
    "river_bed_small_aggregates__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_sagg),
        unit = Unit(; t = 1),
        description = "Small aggregates mass in river bed",
    ),
    "river_water_sand__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_sand),
        unit = Unit(; t = 1),
        description = "Sand mass in river water",
    ),
    "river_bed_sand__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_sand),
        unit = Unit(; t = 1),
        description = "Sand mass in river bed",
    ),
    "river_water_silt__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_silt),
        unit = Unit(; t = 1),
        description = "Silt mass in river water",
    ),
    "river_bed_silt__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_silt),
        unit = Unit(; t = 1),
        description = "Silt mass in river bed",
    ),
    "river_water_sediment_erosion__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.erosion),
        unit = Unit(; t = 1),
        description = "Sediment erosion mass in river",
    ),
    "river_water_sediment_deposition__mass" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.sediment_flux.variables.deposition),
        unit = Unit(; t = 1),
        description = "Sediment deposition mass in river",
    ),
    "sediment__particle_density" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.transport_capacity.parameters.density),
        unit = Unit(; kg = 1, m = -3),
        default = 2650.0,
        description = "Particle density of sediment",
    ),
    "land_surface_sediment__median_diameter" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.transport_capacity.parameters.d50),
        unit = Unit(; mm = 1),
        default = 0.1,
        description = "Median diameter of land surface sediment",
    ),
    "land_surface_water_sediment__govers_transport_capacity_coefficient" =>
        ParameterMetadata(;
            lens = @optic(_.routing.overland_flow.transport_capacity.parameters.c_govers),
            default = 0.0000005,
            description = "Govers transport capacity coefficient for overland flow",
        ),
    "land_surface_water_sediment__govers_transport_capacity_exponent" =>
        ParameterMetadata(;
            lens = @optic(_.routing.overland_flow.transport_capacity.parameters.n_govers),
            default = 1.5,
            description = "Govers transport capacity exponent for overland flow",
        ),
    "river_sediment__median_diameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.transport_capacity.parameters.d50),
        unit = Unit(; mm = 1),
        default = 1.0,
        description = "Median diameter of river sediment",
    ),
    "river_water_sediment__bagnold_transport_capacity_coefficient" =>
        ParameterMetadata(;
            lens = @optic(_.routing.river_flow.transport_capacity.parameters.c_bagnold),
            default = 0.0017,
            description = "Bagnold transport capacity coefficient for river flow",
        ),
    "river_water_sediment__bagnold_transport_capacity_exponent" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.transport_capacity.parameters.e_bagnold),
        default = 1.8,
        description = "Bagnold transport capacity exponent for river flow",
    ),
    "river_water_sediment__kodatie_transport_capacity_a_coefficient" =>
        ParameterMetadata(;
            lens = @optic(_.routing.river_flow.transport_capacity.parameters.a_kodatie),
            default = 5.0,
            description = "Kodatie transport capacity coefficient A for river flow",
        ),
    "river_water_sediment__kodatie_transport_capacity_b_coefficient" =>
        ParameterMetadata(;
            lens = @optic(_.routing.river_flow.transport_capacity.parameters.b_kodatie),
            default = 0.6,
            description = "Kodatie transport capacity coefficient B for river flow",
        ),
    "river_water_sediment__kodatie_transport_capacity_c_coefficient" =>
        ParameterMetadata(;
            lens = @optic(_.routing.river_flow.transport_capacity.parameters.c_kodatie),
            default = 0.05,
            description = "Kodatie transport capacity coefficient C for river flow",
        ),
    "river_water_sediment__kodatie_transport_capacity_d_coefficient" =>
        ParameterMetadata(;
            lens = @optic(_.routing.river_flow.transport_capacity.parameters.d_kodatie),
            default = 2.0,
            description = "Kodatie transport capacity coefficient D for river flow",
        ),
)
