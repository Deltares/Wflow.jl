const domain_standard_name_map = Dict{String, ParameterMetadata}(
    "subbasin_location__count" => ParameterMetadata(;
        type = Int,
        description = "Subbasin ids",
        allow_missing = true,
    ),
    "basin__local_drain_direction" => ParameterMetadata(;
        type = Int,
        description = "Local drain direction (1-9)",
        allow_missing = true,
    ),
    "basin_pit_location__mask" =>
        ParameterMetadata(; fill = false, description = "Pit location mask"),
    "river_location__mask" =>
        ParameterMetadata(; fill = false, description = "River mask (0-1)"),
    "reservoir_location__count" =>
        ParameterMetadata(; fill = 0, description = "Reservoir location ids"),
    "reservoir_area__count" => ParameterMetadata(;
        type = Int,
        allow_missing = true,
        description = "Reservoir area ids",
    ),
    "land_water_allocation_area__count" =>
        ParameterMetadata(; default = 1, description = "Water allocation area ids"),
    "land_surface__slope" => ParameterMetadata(;
        lens = @optic(_.domain.land.parameters.slope),
        unit = Unit(; m = (1, 1)),
        description = "Land slope",
    ),
    "river__width" => ParameterMetadata(;
        unit = Unit(; m = 1),
        fill = 0.0,
        description = "River width",
    ),
    "river__length" => ParameterMetadata(;
        unit = Unit(; m = 1),
        fill = 0.0,
        description = "River length",
    ),
    "river__slope" =>
        ParameterMetadata(; unit = Unit(; m = (1, 1)), description = "River slope"),
    "land_water_covered__area_fraction" => ParameterMetadata(;
        default = 0.0,
        description = "Land water covered area fraction",
    ),
    "land_drain_location__mask" =>
        ParameterMetadata(; fill = false, description = "Drain location mask"),
)
