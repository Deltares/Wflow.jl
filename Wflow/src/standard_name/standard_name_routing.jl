
const routing_standard_name_map = Dict{String, ParameterMetadata}(
    # Subsurface flow parameters
    "subsurface_water__horizontal_to_vertical_saturated_hydraulic_conductivity_ratio" =>
        ParameterMetadata(;
            lens = @optic(_.routing.subsurface_flow.parameters.khfrac),
            description = "A muliplication factor applied to vertical hydraulic conductivity",
        ),

    # Groundwater flow parameters
    "subsurface_surface_water__horizontal_saturated_hydraulic_conductivity" =>
        ParameterMetadata(;
            lens = @optic(_.routing.subsurface_flow.aquifer.parameters.k),
            unit = Unit(; m = 1, d = -1),
            description = "Horizontal conductivity",
        ),
    "subsurface_water__specific_yield" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.aquifer.parameters.specific_yield),
        description = "Specific yield",
    ),
    "subsurface__horizontal_saturated_hydraulic_conductivity_scale_parameter" =>
        ParameterMetadata(;
            lens = @optic(_.routing.subsurface_flow.aquifer.parameters.f),
            unit = Unit(; m = -1),
            description = "Factor controlling the reduction of horizontal conductivity with depth",
        ),
    "model_constant_boundary_condition__hydraulic_head" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.constanthead.variables.head),
        unit = Unit(; m = 1),
        description = "Head of the boundary",
    ),

    # Groundwater boundary condition parameters
    "river_water__infiltration_conductance" => ParameterMetadata(;
        lens = @optic(
            _.routing.subsurface_flow.boundaries.river.parameters.infiltration_conductance
        ),
        unit = Unit(; m = 2, d = -1),
        description = "River bed infiltration conductance",
    ),
    "river_water__exfiltration_conductance" => ParameterMetadata(;
        lens = @optic(
            _.routing.subsurface_flow.boundaries.river.parameters.exfiltration_conductance
        ),
        unit = Unit(; m = 2, d = -1),
        description = "River bed exfiltration conductance",
    ),
    "river_bottom__elevation" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.river.parameters.bottom),
        unit = Unit(; m = 1),
        description = "River bottom elevation",
    ),
    "land_drain__elevation" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.parameters.elevation),
        unit = Unit(; m = 1),
        description = "Drain elevation",
    ),
    "land_drain__conductance" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.parameters.conductance),
        unit = Unit(; m = 2, d = -1),
        description = "Drain conductance",
    ),

    # Reservoir parameters
    "reservoir_surface__area" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.parameters.area),
        unit = Unit(; m = 2),
        description = "Area of the reservoir",
    ),
    "reservoir_water_surface__initial_elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Water level of reservoir (used for initialization)",
    ),
    "reservoir_water__storage_curve_type_count" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.storfunc
        ),
        description = "Type of reservoir storage curve",
    ),
    "reservoir_water__rating_curve_type_count" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.outflowfunc
        ),
        description = "Type of reservoir rating curve",
    ),
    "reservoir_lower_location__count" => ParameterMetadata(;
        default = 0,
        description = "Index of lower reservoir (linked reservoir)",
    ),
    "reservoir_location__count" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.parameters.id),
        description = "Outlet of the reservoirs in which each reservoir has a unique id",
    ),
    "reservoir_water_flow_threshold_level__elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Water level threshold, below this level outflow is zero",
    ),
    "reservoir_water__rating_curve_coefficient" =>
        ParameterMetadata(; description = "Rating curve coefficient"),
    "reservoir_water__rating_curve_exponent" =>
        ParameterMetadata(; description = "Rating curve exponent"),
    "reservoir_water_demand__required_downstream_volume_flow_rate" =>
        ParameterMetadata(;
            unit = Unit(; m = 3, s = -1),
            description = "Minimum (environmental) flow released from reservoir",
        ),
    "reservoir_water_release_below_spillway__max_volume_flow_rate" =>
        ParameterMetadata(;
            unit = Unit(; m = 3, s = -1),
            description = "Maximum amount that can be released if below spillway",
        ),
    "reservoir_water__max_volume" => ParameterMetadata(;
        unit = Unit(; m = 3),
        description = "Maximum volume (above which water is spilled)",
    ),
    "reservoir_water__target_full_volume_fraction" =>
        ParameterMetadata(; description = "Target fraction full (of max storage)"),
    "reservoir_water__target_min_volume_fraction" => ParameterMetadata(;
        description = "Target minimum full fraction (of max storage)",
    ),
    "reservoir_water__outgoing_observed_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.outflow_obs
        ),
        unit = Unit(; m = 3, s = -1),
        description = "Observed outflow reservoir",
    ),
    "reservoir_water__external_inflow_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.external_inflow
        ),
        unit = Unit(; m = 3, s = -1),
        description = "External inflow reservoir (negative for abstractions)",
    ),

    # River flow parameters
    "model_boundary_condition_river__length" => ParameterMetadata(;
        unit = Unit(; m = 1),
        default = 1.0e4,
        description = "Boundary condition river length downstream river outlets",
    ),
    "river_bank_water__elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Bankfull elevation of the river",
    ),
    "river_bank_water__depth" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.parameters.bankfull_depth),
        unit = Unit(; m = 1),
        default = 1.0,
        description = "Bankfull river depth",
    ),
    "model_boundary_condition_river_bank_water__depth" => ParameterMetadata(;
        unit = Unit(; m = 1),
        default = 0.0,
        description = "Boundary condition bankfull depth downstream river outlets",
    ),
    "river_water_flow__manning_n_parameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.parameters.flow.mannings_n),
        unit = Unit(; s = 1, m = -1 // 3),
        default = 0.036,
        description = "Manning's roughness",
    ),
    "river_water__external_inflow_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.external_inflow),
        unit = Unit(; m = 3, s = -1),
        default = 0.0,
        description = "External inflow into the river (negative for abstractions)",
    ),

    # Land/overland flow parameters
    "land_surface_water_flow__manning_n_parameter" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.parameters.mannings_n),
        unit = Unit(; s = 1, m = -1 // 3),
        default = 0.072,
        description = "Manning's roughness",
    ),
    "land_surface_water_flow__ground_elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Elevation of each cell",
    ),
    "land_surface__elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Land surface elevation",
    ),
    "floodplain_water__sum_of_volume_per_depth" => ParameterMetadata(;
        unit = Unit(; m = 3),
        description = "Floodplain profile (cumulative volume per flood depth)",
    ),
    "floodplain_water_flow__manning_n_parameter" => ParameterMetadata(;
        unit = Unit(; s = 1, m = -1 // 3),
        default = 0.072,
        description = "Manning's roughness",
    ),
)
