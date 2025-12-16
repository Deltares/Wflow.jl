const domain_standard_name_map = Dict{String, VariableMetadata}(
    "subbasin_location__count" => VariableMetadata(; description = "Subbasin ids"),
    "basin__local_drain_direction" =>
        VariableMetadata(; description = "Local drain direction (1-9)"),
    "basin_pit_location__mask" => VariableMetadata(; description = "Pit location mask"),
    "river_location__mask" => VariableMetadata(; description = "River mask (0-1)"),
    "reservoir_location__count" =>
        VariableMetadata(; description = "Reservoir location ids"),
    "reservoir_area__count" => VariableMetadata(; description = "Reservoir area ids"),
    "land_water_allocation_area__count" =>
        VariableMetadata(; description = "Water allocation area ids"),
    "land_surface__slope" => VariableMetadata(;
        lens = @optic(_.domain.land.parameters.slope),
        unit = Unit(; m = (1, 1)),
        description = "Land slope",
    ),
    "river__width" =>
        VariableMetadata(; unit = Unit(; m = 1), description = "River width"),
    "river__length" =>
        VariableMetadata(; unit = Unit(; m = 1), description = "River length"),
    "river__slope" =>
        VariableMetadata(; unit = Unit(; m = (1, 1)), description = "River slope"),
    "land_water_covered__area_fraction" =>
        VariableMetadata(; description = "Land water covered area fraction"),
    "land_drain_location__mask" =>
        VariableMetadata(; description = "Drain location mask"),
)
