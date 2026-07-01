# NOTE: The order of the entries determines the order in the docs tables
"""
Mapping of (CSDMS) standard names to the metadata associated with the corresponding
domain parameter. For more details and default values see `ParameterMetadata`.
"""
const domain_standard_name_map = OrderedDict{String, ParameterMetadata}(
    "subbasin_location__count" => ParameterMetadata(;
        type = Int,
        description = "Subbasin ids",
        allow_missing = true,
        tags = [:generic_input_map],
    ),
    "basin__local_drain_direction" => ParameterMetadata(;
        type = Int,
        description = "Local drain direction (1-9)",
        allow_missing = true,
        tags = [:generic_input_map],
    ),
    "basin_pit_location__mask" => ParameterMetadata(;
        fill = false,
        description = "Pit location mask",
        tags = [:generic_input_map],
    ),
    "river_location__mask" => ParameterMetadata(;
        fill = false,
        description = "River mask (0-1)",
        tags = [:generic_input_map],
    ),
    "land_water_allocation_area__count" => ParameterMetadata(;
        default = 1,
        description = "Water allocation area ids",
        tags = [:generic_input_map],
    ),
    "land_surface__albedo" => ParameterMetadata(;
        lens = @optic(_.domain.land.parameters.albedo),
        description = "Land surface albedo",
        default = 0.2,
        allow_dynamic_input = true,
        tags = [:generic_input_map],
    ),
    "atmosphere_bottom_air_flow__roughness_length" => ParameterMetadata(;
        lens = @optic(_.domain.land.parameters.z0m),
        unit = Unit(; m = 1),
        description = "Land surface aerodynamic roughness length for momentum transfer",
        allow_dynamic_input = true,
        tags = [:generic_input_map],
    ),
    "atmosphere_bottom_air_heat_flow__roughness_length" => ParameterMetadata(;
        lens = @optic(_.domain.land.parameters.z0h),
        unit = Unit(; m = 1),
        description = "Land surface aerodynamic roughness length for heat transfer",
        allow_dynamic_input = true,
        tags = [:generic_input_map],
    ),
    "atmosphere_bottom_air_flow__zero_plane_displacement_length" => ParameterMetadata(;
        lens = @optic(_.domain.land.parameters.d0),
        unit = Unit(; m = 1),
        description = "Zero-place displacement height",
        allow_dynamic_input = true,
        tags = [:generic_input_map],
    ),
)
