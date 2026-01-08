# NOTE: The order of the entries determines the order in the docs tables
const routing_standard_name_map = OrderedDict{String, ParameterMetadata}(
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
        fill = MISSING_VALUE,
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
        fill = MISSING_VALUE,
        description = "Drain elevation",
    ),
    "land_drain__conductance" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.parameters.conductance),
        unit = Unit(; m = 2, d = -1),
        fill = MISSING_VALUE,
        description = "Drain conductance",
    ),

    ## Reservoir parameters
    #### Generic input
    "reservoir_location__count" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.parameters.id),
        fill = 0,
        description = "Outlet of the reservoirs in which each reservoir has a unique id",
        flags = [:reservoir_generic_input],
    ),
    "reservoir_area__count" => ParameterMetadata(;
        type = Int,
        allow_missing = true,
        description = "Reservoir coverage",
        flags = [:reservoir_generic_input],
    ),
    #### Input
    "reservoir_surface__area" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.parameters.area),
        unit = Unit(; m = 2),
        description = "Area of the reservoir",
        flags = [:reservoir_input],
    ),
    "reservoir_water__max_volume" => ParameterMetadata(;
        unit = Unit(; m = 3),
        description = "Maximum volume (above which water is spilled)",
        flags = [:reservoir_input],
    ),
    "reservoir_water__rating_curve_coefficient" => ParameterMetadata(;
        description = "Rating curve coefficient",
        flags = [:reservoir_input],
    ),
    "reservoir_water__rating_curve_exponent" => ParameterMetadata(;
        description = "Rating curve exponent",
        flags = [:reservoir_input],
    ),
    "reservoir_water_flow_threshold_level__elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Water level threshold, below this level outflow is zero",
        flags = [:reservoir_input],
    ),
    "reservoir_lower_location__count" => ParameterMetadata(;
        default = 0,
        fill = 0,
        description = "Index of lower reservoir (linked reservoir)",
        flags = [:reservoir_input],
    ),
    "reservoir_water__storage_curve_type_count" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.storfunc
        ),
        type = Int,
        description = "Type of reservoir storage curve",
        flags = [:reservoir_input],
    ),
    "reservoir_water__rating_curve_type_count" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.outflowfunc
        ),
        type = Int,
        description = "Type of reservoir rating curve",
        flags = [:reservoir_input],
    ),
    "reservoir_water_surface__initial_elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Water level of reservoir (used for initialization)",
        flags = [:reservoir_input],
    ),
    #### Static or cyclic/forcing input
    "reservoir_water_demand__required_downstream_volume_flow_rate" =>
        ParameterMetadata(;
            unit = Unit(; m = 3, s = -1),
            description = "Minimum (environmental) flow released from reservoir",
            flags = [:reservoir_static_cyclic_forcing_input],
        ),
    "reservoir_water_release_below_spillway__max_volume_flow_rate" =>
        ParameterMetadata(;
            unit = Unit(; m = 3, s = -1),
            description = "Maximum amount that can be released if below spillway",
            flags = [:reservoir_static_cyclic_forcing_input],
        ),
    "reservoir_water__target_full_volume_fraction" => ParameterMetadata(;
        description = "Target fraction full (of max storage)",
        flags = [:reservoir_static_cyclic_forcing_input],
    ),
    "reservoir_water__target_min_volume_fraction" => ParameterMetadata(;
        description = "Target minimum full fraction (of max storage)",
        flags = [:reservoir_static_cyclic_forcing_input],
    ),
    #### States and Output
    "reservoir_water__volume" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.storage),
        unit = Unit(; m = 3),
        description = "Reservoir water volume",
        flags = [:reservoir_output],
    ),
    "reservoir_water_surface__elevation" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.waterlevel
        ),
        unit = Unit(; m = 1),
        description = "Reservoir water level",
        flags = [:reservoir_state, :reservoir_output],
    ),
    "reservoir_water__outgoing_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.outflow),
        unit = Unit(; m = 3, s = -1),
        description = "Outflow of the reservoir",
        flags = [:reservoir_output],
    ),
    "reservoir_water__incoming_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.inflow
        ),
        unit = Unit(; m = 3, s = -1),
        description = "Reservoir incoming flow rate",
        flags = [:reservoir_output],
    ),
    "reservoir_water__evaporation_volume_flux" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.actevap),
        unit = Unit(; mm = 1, dt = -1),
        description = "Average actual evaporation over the reservoir area",
        flags = [:reservoir_output],
    ),
    "reservoir_water__precipitation_volume_flux" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation
        ),
        unit = Unit(; mm = 1, dt = -1),
        description = "Reservoir precipitation flux",
        flags = [:reservoir_output],
    ),
    "reservoir_water__potential_evaporation_volume_flux" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.evaporation
        ),
        unit = Unit(; mm = 1, dt = -1),
        description = "Reservoir potential evaporation flux",
        flags = [:reservoir_output],
    ),

    ## Kinematic wave
    ### River flow
    #### Input
    "river__length" => ParameterMetadata(;
        unit = Unit(; m = 1),
        fill = 0.0,
        description = "River length",
        flags = [:kinematic_wave_river_flow_input],
    ),
    "river__width" => ParameterMetadata(;
        unit = Unit(; m = 1),
        fill = 0.0,
        description = "River width",
        flags = [:kinematic_wave_river_flow_input],
    ),
    "river_water_flow__manning_n_parameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.parameters.flow.mannings_n),
        unit = Unit(; s = 1, m = -1 // 3),
        default = 0.036,
        description = "Manning's roughness",
        flags = [:kinematic_wave_river_flow_input],
    ),
    "river_bank_water__depth" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.parameters.bankfull_depth),
        unit = Unit(; m = 1),
        default = 1.0,
        fill = 0.0,
        description = "Bankfull river depth",
        flags = [:kinematic_wave_river_flow_input],
    ),
    "river__slope" => ParameterMetadata(;
        unit = Unit(; m = (1, 1)),
        description = "River slope",
        flags = [:kinematic_wave_river_flow_input],
    ),
    #### States
    "river_water__instantaneous_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.variables.q),
        unit = Unit(; m = 3, s = -1),
        description = "River discharge",
        flags = [:kinematic_wave_river_state],
    ),
    "river_water__depth" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.variables.h),
        unit = Unit(; m = 1),
        description = "River water depth",
    ),
    #### Output

    # River flow parameters
    "model_boundary_condition_river__length" => ParameterMetadata(;
        unit = Unit(; m = 1),
        default = 1.0e4,
        description = "Boundary condition river length downstream river outlets",
    ),
    "river_bank_water__elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        fill = 0.0,
        description = "Bankfull elevation of the river",
    ),
    "model_boundary_condition_river_bank_water__depth" =>
        ParameterMetadata(; default = 0, unit = Unit(; m = 1)),
    # Land/overland flow parameters
    "land_surface_water_flow__manning_n_parameter" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.parameters.mannings_n),
        unit = Unit(; s = 1, m = -1 // 3),
        default = 0.072,
        description = "Manning's roughness",
    ),
    "land_surface_water_flow__ground_elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        fill = 0.0,
        description = "Elevation of each cell",
    ),
    "land_surface__elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Land surface elevation",
    ),
    "floodplain_water__sum_of_volume_per_depth" => ParameterMetadata(;
        unit = Unit(; m = 3),
        dimname = :flood_depth,
        description = "Floodplain profile (cumulative volume per flood depth)",
    ),
    "floodplain_water_flow__manning_n_parameter" => ParameterMetadata(;
        unit = Unit(; s = 1, m = -1 // 3),
        default = 0.072,
        description = "Manning's roughness",
    ),
)
