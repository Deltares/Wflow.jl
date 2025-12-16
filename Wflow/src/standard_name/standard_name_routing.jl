
const routing_standard_name_map = Dict{String, NamedTuple}(
    # Subsurface flow parameters
    "subsurface_water__horizontal_to_vertical_saturated_hydraulic_conductivity_ratio" =>
        (lens = @optic(_.routing.subsurface_flow.parameters.khfrac), unit = Unit()),

    # Groundwater flow parameters
    "subsurface_surface_water__horizontal_saturated_hydraulic_conductivity" => (
        lens = @optic(_.routing.subsurface_flow.aquifer.parameters.k),
        unit = Unit(; m = 1, d = -1),
    ),
    "subsurface_water__specific_yield" => (
        lens = @optic(_.routing.subsurface_flow.aquifer.parameters.specific_yield),
        unit = Unit(),
    ),
    "subsurface__horizontal_saturated_hydraulic_conductivity_scale_parameter" => (
        lens = @optic(_.routing.subsurface_flow.aquifer.parameters.f),
        unit = Unit(; m = -1),
    ),
    "model_constant_boundary_condition__hydraulic_head" => (
        lens = @optic(_.routing.subsurface_flow.constanthead.variables.head),
        unit = Unit(; m = 1),
    ),

    # Groundwater boundary condition parameters
    "river_water__infiltration_conductance" => (
        lens = @optic(
            _.routing.subsurface_flow.boundaries.river.parameters.infiltration_conductance
        ),
        unit = Unit(; m = 2, d = -1),
    ),
    "river_water__exfiltration_conductance" => (
        lens = @optic(
            _.routing.subsurface_flow.boundaries.river.parameters.exfiltration_conductance
        ),
        unit = Unit(; m = 2, d = -1),
    ),
    "river_bottom__elevation" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.river.parameters.bottom),
        unit = Unit(; m = 1),
    ),
    "land_drain__elevation" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.parameters.elevation),
        unit = Unit(; m = 1),
    ),
    "land_drain__conductance" => (
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.parameters.conductance),
        unit = Unit(; m = 2, d = -1),
    ),

    # Reservoir parameters
    "reservoir_surface__area" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.parameters.area),
        unit = Unit(; m = 2),
    ),
    "reservoir_water_surface__initial_elevation" =>
        (lens = nothing, unit = Unit(; m = 1)),
    "reservoir_water__storage_curve_type_count" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.storfunc
        ),
        unit = Unit(),
    ),
    "reservoir_water__rating_curve_type_count" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.outflowfunc
        ),
        unit = Unit(),
    ),
    "reservoir_lower_location__count" => (lens = nothing, unit = Unit()),
    "reservoir_location__count" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.parameters.id),
        unit = Unit(),
    ),
    "reservoir_water_flow_threshold_level__elevation" =>
        (lens = nothing, unit = Unit(; m = 1)),
    "reservoir_water__rating_curve_coefficient" => (lens = nothing, unit = Unit()),
    "reservoir_water__rating_curve_exponent" => (lens = nothing, unit = Unit()),
    "reservoir_water_demand__required_downstream_volume_flow_rate" =>
        (lens = nothing, unit = Unit(; m = 3, s = -1)),
    "reservoir_water_release_below_spillway__max_volume_flow_rate" =>
        (lens = nothing, unit = Unit(; m = 3, s = -1)),
    "reservoir_water__max_volume" => (lens = nothing, unit = Unit(; m = 3)),
    "reservoir_water__target_full_volume_fraction" => (lens = nothing, unit = Unit()),
    "reservoir_water__target_min_volume_fraction" => (lens = nothing, unit = Unit()),
    "reservoir_water__outgoing_observed_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.outflow_obs
        ),
        unit = Unit(; m = 3, s = -1),
    ),
    "reservoir_water__external_inflow_volume_flow_rate" => (
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.external_inflow
        ),
        unit = Unit(; m = 3, s = -1),
    ),

    # River flow parameters
    "model_boundary_condition_river__length" => (lens = nothing, unit = Unit(; m = 1)),
    "river_bank_water__elevation" => (lens = nothing, unit = Unit(; m = 1)),
    "river_bank_water__depth" => (
        lens = @optic(_.routing.river_flow.parameters.bankfull_depth),
        unit = Unit(; m = 1),
    ),
    "model_boundary_condition_river_bank_water__depth" =>
        (lens = nothing, unit = Unit(; m = 1)),
    "river_water_flow__manning_n_parameter" => (
        lens = @optic(_.routing.river_flow.parameters.flow.mannings_n),
        unit = Unit(; s = 1, m = -1 // 3),
    ),
    "river_water__external_inflow_volume_flow_rate" => (
        lens = @optic(_.routing.river_flow.boundary_conditions.external_inflow),
        unit = Unit(; m = 3, s = -1),
    ),

    # Land/overland flow parameters
    "land_surface_water_flow__manning_n_parameter" => (
        lens = @optic(_.routing.overland_flow.parameters.mannings_n),
        unit = Unit(; s = 1, m = -1 // 3),
    ),
    "land_surface_water_flow__ground_elevation" =>
        (lens = nothing, unit = Unit(; m = 1)),
    "land_surface__elevation" => (lens = nothing, unit = Unit(; m = 1)),
    "floodplain_water__sum_of_volume_per_depth" =>
        (lens = nothing, unit = Unit(; m = 3)),
    "floodplain_water_flow__manning_n_parameter" =>
        (lens = nothing, unit = Unit(; s = 1, m = -1 // 3)),
)
