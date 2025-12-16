"""
Mapping of (CSDMS) standard names to model variables and units for models with a land model
of type `SoilLoss`. The `lens` of the NamedTuple allows access to a nested model variable.
"""
const sediment_standard_name_map = Dict{String, NamedTuple}(
    "atmosphere_water__precipitation_volume_flux" => (
        lens = @optic(_.land.atmospheric_forcing.precipitation),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "soil_clay__mass_fraction" =>
        (lens = @optic(_.land.soil_erosion.parameters.clay_fraction), unit = Unit()),
    "soil_silt__mass_fraction" =>
        (lens = @optic(_.land.soil_erosion.parameters.silt_fraction), unit = Unit()),
    "soil_sand__mass_fraction" =>
        (lens = @optic(_.land.soil_erosion.parameters.sand_fraction), unit = Unit()),
    "soil_small_aggregates__mass_fraction" =>
        (lens = @optic(_.land.soil_erosion.parameters.sagg_fraction), unit = Unit()),
    "soil_large_aggregates__mass_fraction" =>
        (lens = @optic(_.land.soil_erosion.parameters.lagg_fraction), unit = Unit()),
    "river_bottom_and_bank_sediment__median_diameter" => (
        lens = @optic(_.routing.potential_erosion.parameters.d50),
        unit = Unit(; mm = 1),
    ),
    "soil_erosion__usle_k_factor" =>
        (lens = @optic(_.land.overland_flow_erosion.parameters.usle_k), unit = Unit()),
    "soil_erosion__usle_c_factor" =>
        (lens = @optic(_.land.overland_flow_erosion.parameters.usle_c), unit = Unit()),
    "soil_erosion__answers_overland_flow_factor" => (
        lens = @optic(_.land.overland_flow_erosion.parameters.answers_overland_flow_factor),
        unit = Unit(),
    ),
    "soil_erosion__rainfall_soil_detachability_factor" => (
        lens = @optic(_.land.rainfall_erosion.parameters.soil_detachability),
        unit = Unit(; g = 1, J = -1),
    ),
    "soil_erosion__eurosem_exponent" => (
        lens = @optic(_.land.rainfall_erosion.parameters.eurosem_exponent),
        unit = Unit(; m = -1),
    ),
    "vegetation_canopy__height" => (
        lens = @optic(_.land.rainfall_erosion.parameters.canopyheight),
        unit = Unit(; m = 1),
    ),
    "vegetation_canopy__gap_fraction" => (
        lens = @optic(_.land.rainfall_erosion.parameters.canopygapfraction),
        unit = Unit(),
    ),
    "compacted_soil__area_fraction" => (
        lens = @optic(_.land.rainfall_erosion.parameters.soilcover_fraction),
        unit = Unit(),
    ),
    "soil_erosion__answers_rainfall_factor" => (
        lens = @optic(_.land.rainfall_erosion.parameters.answers_rainfall_factor),
        unit = Unit(),
    ),
    "river_bottom_and_bank_clay__mass_fraction" => (
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.clay_fraction),
        unit = Unit(),
    ),
    "river_bottom_and_bank_silt__mass_fraction" => (
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.silt_fraction),
        unit = Unit(),
    ),
    "river_bottom_and_bank_sand__mass_fraction" => (
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.sand_fraction),
        unit = Unit(),
    ),
    "river_bottom_and_bank_gravel__mass_fraction" => (
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.gravel_fraction),
        unit = Unit(),
    ),
    "clay__mean_diameter" => (
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_clay),
        unit = Unit(; μm = 1),
    ),
    "silt__mean_diameter" => (
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_silt),
        unit = Unit(; μm = 1),
    ),
    "sand__mean_diameter" => (
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_sand),
        unit = Unit(; μm = 1),
    ),
    "sediment_small_aggregates__mean_diameter" => (
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_sagg),
        unit = Unit(; μm = 1),
    ),
    "sediment_large_aggregates__mean_diameter" => (
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_lagg),
        unit = Unit(; μm = 1),
    ),
    "gravel__mean_diameter" => (
        lens = @optic(_.routing.river_flow.concentrations.parameters.dm_gravel),
        unit = Unit(; μm = 1),
    ),
    "reservoir_location__count" => (
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.reservoir_outlet),
        unit = Unit(),
    ),
    "reservoir_surface__area" => (
        lens = @optic(_.routing.river_flow.sediment_flux.parameters.reservoir_area),
        unit = Unit(; m = 2),
    ),
    "reservoir_water_sediment__bedload_trapping_efficiency" => (
        lens = @optic(
            _.routing.river_flow.sediment_flux.parameters.reservoir_trapping_efficiency
        ),
        unit = Unit(),
    ),
    "river_water__volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.hydrological_forcing.q_river),
        unit = Unit(; m = 3, s = -1),
    ),
    "river_water__depth" => (
        lens = @optic(_.routing.river_flow.hydrological_forcing.waterlevel_river),
        unit = Unit(; m = 1),
    ),
    "land_surface_water__volume_flow_rate" => (
        lens = @optic(_.land.hydrological_forcing.q_land),
        unit = Unit(; m = 3, s = -1),
    ),
    "land_surface_water__depth" => (
        lens = @optic(_.land.hydrological_forcing.waterlevel_land),
        unit = Unit(; m = 1),
    ),
    "vegetation_canopy_water__interception_volume_flux" => (
        lens = @optic(_.land.hydrological_forcing.interception),
        unit = Unit(; mm = 1, dt = -1),
    ),
    "rainfall_soil_erosion__mass_flow_rate" => (
        lens = @optic(_.land.rainfall_erosion.variables.amount),
        unit = Unit(; t = 1, dt = -1),
    ),
    "overland_flow_soil_erosion__mass_flow_rate" => (
        lens = @optic(_.land.overland_flow_erosion.variables.amount),
        unit = Unit(; t = 1, dt = -1),
    ),
    "soil_erosion__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.amount),
        unit = Unit(; t = 1, dt = -1),
    ),
    "soil_erosion_clay__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.clay),
        unit = Unit(; t = 1, dt = -1),
    ),
    "soil_erosion_silt__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.silt),
        unit = Unit(; t = 1, dt = -1),
    ),
    "soil_erosion_sand__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.sand),
        unit = Unit(; t = 1, dt = -1),
    ),
    "soil_erosion_small_aggregates__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.sagg),
        unit = Unit(; t = 1, dt = -1),
    ),
    "soil_erosion_large_aggregates__mass_flow_rate" => (
        lens = @optic(_.land.soil_erosion.variables.lagg),
        unit = Unit(; t = 1, dt = -1),
    ),
    "land_surface_water_sediment_transport_capacity__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.transport_capacity.variables.amount),
        unit = Unit(; t = 1, dt = -1),
    ),
    "land_surface_water_sediment__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.amount),
        unit = Unit(; t = 1, dt = -1),
    ),
    "land_surface_water_clay__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.clay),
        unit = Unit(; t = 1, dt = -1),
    ),
    "land_surface_water_silt__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.silt),
        unit = Unit(; t = 1, dt = -1),
    ),
    "land_surface_water_sand__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.sand),
        unit = Unit(; t = 1, dt = -1),
    ),
    "land_surface_water_small_aggregates__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.sagg),
        unit = Unit(; t = 1, dt = -1),
    ),
    "land_surface_water_large_aggregates__to_river_mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.to_river.variables.lagg),
        unit = Unit(; t = 1, dt = -1),
    ),
    "land_surface_water_sediment__mass_flow_rate" => (
        lens = @optic(_.routing.overland_flow.sediment_flux.variables.amount),
        unit = Unit(; t = 1, dt = -1),
    ),
    "river_water_sediment__bedload_mass_concentration" => (
        lens = @optic(_.routing.river_flow.concentrations.variables.bed),
        unit = Unit(; g = 1, m = -3),
    ),
    "river_water_sediment__suspended_mass_concentration" => (
        lens = @optic(_.routing.river_flow.concentrations.variables.suspended),
        unit = Unit(; g = 1, m = -3),
    ),
    "river_water_sediment__mass_concentration" => (
        lens = @optic(_.routing.river_flow.concentrations.variables.total),
        unit = Unit(; g = 1, m = -3),
    ),
    "river_water_clay__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_clay),
        unit = Unit(; t = 1),
    ),
    "river_bed_clay__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_clay),
        unit = Unit(; t = 1),
    ),
    "river_water_gravel__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_gravel),
        unit = Unit(; t = 1),
    ),
    "river_bed_gravel__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_gravel),
        unit = Unit(; t = 1),
    ),
    "river_water_large_aggregates__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_lagg),
        unit = Unit(; t = 1),
    ),
    "river_bed_large_aggregates__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_lagg),
        unit = Unit(; t = 1),
    ),
    "river_water_clay__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.clay),
        unit = Unit(; t = 1, dt = -1),
    ),
    "river_water_gravel__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.gravel),
        unit = Unit(; t = 1, dt = -1),
    ),
    "river_water_large_aggregates__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.lagg),
        unit = Unit(; t = 1, dt = -1),
    ),
    "river_water_small_aggregates__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.sagg),
        unit = Unit(; t = 1, dt = -1),
    ),
    "river_water_sand__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.sand),
        unit = Unit(; t = 1, dt = -1),
    ),
    "river_water_silt__mass_flow_rate" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.silt),
        unit = Unit(; t = 1, dt = -1),
    ),
    "river_water_small_aggregates__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_sagg),
        unit = Unit(; t = 1),
    ),
    "river_bed_small_aggregates__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_sagg),
        unit = Unit(; t = 1),
    ),
    "river_water_sand__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_sand),
        unit = Unit(; t = 1),
    ),
    "river_bed_sand__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_sand),
        unit = Unit(; t = 1),
    ),
    "river_water_silt__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.leftover_silt),
        unit = Unit(; t = 1),
    ),
    "river_bed_silt__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.store_silt),
        unit = Unit(; t = 1),
    ),
    "river_water_sediment_erosion__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.erosion),
        unit = Unit(; t = 1),
    ),
    "river_water_sediment_deposition__mass" => (
        lens = @optic(_.routing.river_flow.sediment_flux.variables.deposition),
        unit = Unit(; t = 1),
    ),
    "sediment__particle_density" => (
        lens = @optic(_.routing.overland_flow.transport_capacity.parameters.density),
        unit = Unit(; kg = 1, m = -3),
    ),
    "land_surface_sediment__median_diameter" => (
        lens = @optic(_.routing.overland_flow.transport_capacity.parameters.d50),
        unit = Unit(; mm = 1),
    ),
    "land_surface_water_sediment__govers_transport_capacity_coefficient" => (
        lens = @optic(_.routing.overland_flow.transport_capacity.parameters.c_govers),
        unit = Unit(),
    ),
    "land_surface_water_sediment__govers_transport_capacity_exponent" => (
        lens = @optic(_.routing.overland_flow.transport_capacity.parameters.n_govers),
        unit = Unit(),
    ),
    "river_sediment__median_diameter" => (
        lens = @optic(_.routing.river_flow.transport_capacity.parameters.d50),
        unit = Unit(; mm = 1),
    ),
    "river_water_sediment__bagnold_transport_capacity_coefficient" => (
        lens = @optic(_.routing.river_flow.transport_capacity.parameters.c_bagnold),
        unit = Unit(),
    ),
    "river_water_sediment__bagnold_transport_capacity_exponent" => (
        lens = @optic(_.routing.river_flow.transport_capacity.parameters.e_bagnold),
        unit = Unit(),
    ),
    "river_water_sediment__kodatie_transport_capacity_a_coefficient" => (
        lens = @optic(_.routing.river_flow.transport_capacity.parameters.a_kodatie),
        unit = Unit(),
    ),
    "river_water_sediment__kodatie_transport_capacity_b_coefficient" => (
        lens = @optic(_.routing.river_flow.transport_capacity.parameters.b_kodatie),
        unit = Unit(),
    ),
    "river_water_sediment__kodatie_transport_capacity_c_coefficient" => (
        lens = @optic(_.routing.river_flow.transport_capacity.parameters.c_kodatie),
        unit = Unit(),
    ),
    "river_water_sediment__kodatie_transport_capacity_d_coefficient" => (
        lens = @optic(_.routing.river_flow.transport_capacity.parameters.d_kodatie),
        unit = Unit(),
    ),
)
