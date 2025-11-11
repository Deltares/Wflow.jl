const domain_standard_name_map = Dict{String, NamedTuple}(
    "subbasin_location__count" => (lens = nothing, unit = Unit()),
    "basin__local_drain_direction" => (lens = nothing, unit = Unit()),
    "basin_pit_location__mask" => (lens = nothing, unit = Unit()),
    "river_location__mask" => (lens = nothing, unit = Unit()),
    "reservoir_location__count" => (lens = nothing, unit = Unit()),
    "reservoir_area__count" => (lens = nothing, unit = Unit()),
    "land_water_allocation_area__count" => (lens = nothing, unit = Unit()),
    "land_surface__slope" =>
        (lens = @optic(_.domain.land.parameters.slope), unit = Unit(; m = (1, 1))),
    "river_location__mask" => (lens = nothing, unit = Unit()),
    "river__width" => (lens = nothing, unit = Unit(; m = 1)),
    "river__length" => (lens = nothing, unit = Unit(; m = 1)),
    "river__slope" => (lens = nothing, unit = Unit(; m = (1, 1))),
    "land_water_covered__area_fraction" => (lens = nothing, unit = Unit()),
    "land_drain_location__mask" => (lens = nothing, unit = Unit()),
)
