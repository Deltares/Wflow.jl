# NOTE: The order of the entries determines the order in the docs tables
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
    "basin_pit_location__mask" =>
        ParameterMetadata(; fill = false, description = "Pit location mask"),
    "river_location__mask" => ParameterMetadata(;
        fill = false,
        description = "River mask (0-1)",
        tags = [:generic_input_map],
    ),
    "land_water_allocation_area__count" =>
        ParameterMetadata(; default = 1, description = "Water allocation area ids"),
)
