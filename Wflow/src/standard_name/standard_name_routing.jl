# NOTE: The order of the entries determines the order in the docs tables
const routing_standard_name_map = OrderedDict{String, ParameterMetadata}(
    ## Reservoir parameters
    #### Generic input
    "reservoir_location__count" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.parameters.id),
        fill = 0,
        description = "Outlet of the reservoirs in which each reservoir has a unique id",
        tags = [:reservoir_generic_input],
    ),
    "reservoir_area__count" => ParameterMetadata(;
        type = Int,
        allow_missing = true,
        description = "Reservoir coverage",
        tags = [:reservoir_generic_input],
    ),
    #### Input
    "reservoir_surface__area" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.parameters.area),
        unit = Unit(; m = 2),
        description = "Area of the reservoir",
        tags = [:reservoir_input],
    ),
    "reservoir_water__max_volume" => ParameterMetadata(;
        unit = Unit(; m = 3),
        description = "Maximum volume (above which water is spilled)",
        tags = [:reservoir_input],
    ),
    "reservoir_water__rating_curve_coefficient" => ParameterMetadata(;
        description = "Rating curve coefficient",
        tags = [:reservoir_input],
    ),
    "reservoir_water__rating_curve_exponent" => ParameterMetadata(;
        description = "Rating curve exponent",
        tags = [:reservoir_input],
    ),
    "reservoir_water_flow_threshold_level__elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Water level threshold, below this level outflow is zero",
        tags = [:reservoir_input],
    ),
    "reservoir_lower_location__count" => ParameterMetadata(;
        default = 0,
        fill = 0,
        description = "Index of lower reservoir (linked reservoir)",
        tags = [:reservoir_input],
    ),
    "reservoir_water__storage_curve_type_count" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.storfunc
        ),
        type = Int,
        description = "Type of reservoir storage curve",
        tags = [:reservoir_input],
    ),
    "reservoir_water__rating_curve_type_count" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.parameters.outflowfunc
        ),
        type = Int,
        description = "Type of reservoir rating curve",
        tags = [:reservoir_input],
    ),
    "reservoir_water_surface__initial_elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Water level of reservoir (used for initialization)",
        tags = [:reservoir_input],
    ),
    #### Static or cyclic/forcing input
    "reservoir_water_demand__required_downstream_volume_flow_rate" =>
        ParameterMetadata(;
            unit = Unit(; m = 3, s = -1),
            description = "Minimum (environmental) flow released from reservoir",
            tags = [:reservoir_static_cyclic_forcing_input],
        ),
    "reservoir_water_release_below_spillway__max_volume_flow_rate" =>
        ParameterMetadata(;
            unit = Unit(; m = 3, s = -1),
            description = "Maximum amount that can be released if below spillway",
            tags = [:reservoir_static_cyclic_forcing_input],
        ),
    "reservoir_water__target_full_volume_fraction" => ParameterMetadata(;
        description = "Target fraction full (of max storage)",
        tags = [:reservoir_static_cyclic_forcing_input],
    ),
    "reservoir_water__target_min_volume_fraction" => ParameterMetadata(;
        description = "Target minimum full fraction (of max storage)",
        tags = [:reservoir_static_cyclic_forcing_input],
    ),
    #### States and Output
    "reservoir_water__volume" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.storage),
        unit = Unit(; m = 3),
        description = "Reservoir water volume",
        tags = [:reservoir_output],
    ),
    "reservoir_water_surface__elevation" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.variables.waterlevel
        ),
        unit = Unit(; m = 1),
        description = "Reservoir water level",
        tags = [
            :reservoir_state,
            :reservoir_output,
            :kinematic_wave_river_static_cyclic_forcing_input,
        ],
    ),
    "reservoir_water__outgoing_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.outflow),
        unit = Unit(; m = 3, s = -1),
        description = "Outflow of the reservoir",
        tags = [:reservoir_output],
    ),
    "reservoir_water__incoming_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.inflow
        ),
        unit = Unit(; m = 3, s = -1),
        description = "Inflow into the reservoir",
        tags = [:reservoir_output],
    ),
    "reservoir_water__evaporation_volume_flux" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.reservoir.variables.actevap),
        unit = Unit(; mm = 1, dt = -1),
        description = "Average actual evaporation over the reservoir area",
        tags = [:reservoir_output],
    ),
    "reservoir_water__precipitation_volume_flux" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation
        ),
        unit = Unit(; mm = 1, dt = -1),
        description = "Average precipitation over the reservoir area",
        tags = [:reservoir_output],
    ),
    "reservoir_water__potential_evaporation_volume_flux" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.evaporation
        ),
        unit = Unit(; mm = 1, dt = -1),
        description = "Average potential evaporation over the reservoir area",
        tags = [:reservoir_output],
    ),

    ## Kinematic wave
    ### River flow
    #### Input
    "river__length" => ParameterMetadata(;
        unit = Unit(; m = 1),
        fill = 0.0,
        description = "River length",
        tags = [:kinematic_wave_river_flow_input, :local_inertial_river_input],
    ),
    "river__width" => ParameterMetadata(;
        unit = Unit(; m = 1),
        fill = 0.0,
        description = "River width",
        tags = [:kinematic_wave_river_flow_input, :local_inertial_river_input],
    ),
    "river_water_flow__manning_n_parameter" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.parameters.flow.mannings_n),
        unit = Unit(; s = 1, m = -1 // 3),
        default = 0.036,
        description = "Manning's roughness",
        tags = [:kinematic_wave_river_flow_input, :local_inertial_river_input],
    ),
    "river_bank_water__depth" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.parameters.bankfull_depth),
        unit = Unit(; m = 1),
        default = 1.0,
        fill = 0.0,
        description = "Bankfull river depth",
        tags = [:kinematic_wave_river_flow_input, :local_inertial_river_input],
    ),
    "river_water__external_inflow_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.external_inflow),
        unit = Unit(; m = 3, s = -1),
        default = 0.0,
        description = "External inflow into the river (negative for abstractions)",
        tags = [
            :kinematic_wave_river_static_cyclic_forcing_input,
            :local_inertial_river_input,
        ],
    ),
    "river__slope" => ParameterMetadata(;
        unit = Unit(; m = 1 // 1),
        description = "River slope",
        tags = [:kinematic_wave_river_flow_input],
    ),
    #### States
    "river_water__instantaneous_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.variables.q),
        unit = Unit(; m = 3, s = -1),
        description = "River discharge",
        tags = [
            :kinematic_wave_river_state,
            :kinematic_wave_river_output,
            :local_inertial_river_state,
            :local_inertial_river_output,
        ],
    ),
    "river_water__depth" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.variables.h),
        unit = Unit(; m = 1),
        description = "River water depth",
        tags = [
            :kinematic_wave_river_state,
            :kinematic_wave_river_output,
            :local_inertial_river_state,
            :local_inertial_river_output,
        ],
    ),
    #### Output
    "river_water__volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.variables.q_av),
        unit = Unit(; m = 3, s = -1),
        description = "River (+ floodplain) discharge",
        tags = [:kinematic_wave_river_output, :local_inertial_river_output],
    ),
    "river_water__volume" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.variables.storage),
        unit = Unit(; m = 3),
        description = "River water volume",
        tags = [:kinematic_wave_river_output, :local_inertial_river_output],
    ),
    "river_water__lateral_inflow_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.boundary_conditions.inwater),
        unit = Unit(; m = 3, s = -1),
        description = "Lateral inflow to river",
        tags = [:kinematic_wave_river_output],
    ),
    "river_water__external_abstraction_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(
            _.routing.river_flow.boundary_conditions.actual_external_abstraction_av
        ),
        unit = Unit(; m = 3, s = -1),
        description = "Actual abstraction based on external negative inflow",
        tags = [:kinematic_wave_river_output, :local_inertial_river_output],
    ),
    ### Land/overland flow parameters
    #### input
    "land_surface_water_flow__manning_n_parameter" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.parameters.mannings_n),
        unit = Unit(; s = 1, m = -1 // 3),
        default = 0.072,
        description = "Manning's roughness",
        tags = [:kinematic_wave_overland_input, :local_inertial_overland_input],
    ),
    #### States
    "land_surface_water__instantaneous_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.variables.q),
        unit = Unit(; m = 3, s = -1),
        description = "Discharge overland flow",
        tags = [:kinematic_wave_overland_state, :kinematic_wave_overland_output],
    ),
    "land_surface_water__depth" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.variables.h),
        unit = Unit(; m = 1),
        description = "Water depth",
        tags = [
            :kinematic_wave_overland_state,
            :kinematic_wave_overland_output,
            :local_inertial_overland_state,
            :local_inertial_overland_output,
        ],
    ),
    #### Output
    "land_surface_water__volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.variables.q_av),
        unit = Unit(; m = 3, s = -1),
        description = "Discharge overland flow",
        tags = [:kinematic_wave_overland_output],
    ),
    "land_surface_water__to_river_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.variables.to_river),
        unit = Unit(; m = 3, s = -1),
        description = "Discharge overland flow that flows to the river",
        tags = [:kinematic_wave_overland_output],
    ),
    "land_surface_water__volume" => ParameterMetadata(;
        lens = @optic(_.routing.overland_flow.variables.storage),
        unit = Unit(; m = 3),
        description = "Total surface water storage of cell (including river storage for river cells)",
        tags = [:kinematic_wave_overland_output, :local_inertial_overland_output],
    ),
    ### Lateral subsurface_flow
    #### Input
    "subsurface_water__horizontal_to_vertical_saturated_hydraulic_conductivity_ratio" =>
        ParameterMetadata(;
            lens = @optic(_.routing.subsurface_flow.parameters.khfrac),
            description = "A muliplication factor applied to vertical hydraulic conductivity",
            tags = [:kinematic_lateral_subsurface_input],
        ),
    #### States
    "subsurface_water__volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.variables.ssf),
        unit = Unit(; m = 3, d = -1),
        description = "Subsurface flow",
        tags = [:kinematic_lateral_subsurface_state, :kinematic_lateral_subsurface_output],
    ),
    #### Output
    "subsurface_water_saturated_zone_top__depth" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.variables.zi),
        unit = Unit(; m = 1),
        description = "Pseudo-water table depth (top of the saturated zone)",
        tags = [:kinematic_lateral_subsurface_output],
    ),
    "subsurface_water__exfiltration_volume_flux" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.variables.exfiltwater),
        unit = Unit(; m = 1, dt = -1),
        description = "Exfiltration (groundwater above surface level, saturated excess conditions)",
        tags = [:kinematic_lateral_subsurface_output],
    ),
    "subsurface_water__to_river_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.variables.to_river),
        unit = Unit(; m = 3, d = -1),
        description = "Part of subsurface flow that flows to the river",
        tags = [:kinematic_lateral_subsurface_output],
    ),
    ## Local inertial
    ### River flow
    #### Input
    "model_boundary_condition_river__length" => ParameterMetadata(;
        unit = Unit(; m = 1),
        default = 1.0e4,
        description = "Boundary condition river length downstream river outlets",
        tags = [:local_inertial_river_input],
    ),
    "model_boundary_condition_river_bank_water__depth" => ParameterMetadata(;
        default = 0,
        unit = Unit(; m = 1),
        description = "Boundary condition bankfull depth downstream river outlets",
        tags = [:local_inertial_river_input],
    ),
    "river_bank_water__elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        fill = 0.0,
        description = "Bankfull elevation of the river",
        tags = [:local_inertial_river_input],
    ),
    ### 1D floodplain flow
    #### Input
    "floodplain_water__sum_of_volume_per_depth" => ParameterMetadata(;
        unit = Unit(; m = 3),
        dimname = :flood_depth,
        description = "Floodplain profile (cumulative volume per flood depth)",
        tags = [:local_inertial_floodplain_1D_flow_input],
    ),
    "floodplain_water_flow__manning_n_parameter" => ParameterMetadata(;
        unit = Unit(; s = 1, m = -1 // 3),
        default = 0.072,
        description = "Manning's roughness",
        tags = [:local_inertial_floodplain_1D_flow_input],
    ),
    #### States
    "floodplain_water__instantaneous_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.floodplain.variables.q),
        unit = Unit(; m = 3, s = -1),
        description = "Floodplain discharge",
        tags = [
            :local_inertial_floodplain_1D_flow_state,
            :local_inertial_floodplain_1D_flow_output,
        ],
    ),
    "floodplain_water__depth" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.floodplain.variables.h),
        unit = Unit(; m = 1),
        description = "Floodplain water depth",
        tags = [
            :local_inertial_floodplain_1D_flow_state,
            :local_inertial_floodplain_1D_flow_output,
        ],
    ),
    #### Output
    "floodplain_water__volume" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.floodplain.variables.storage),
        unit = Unit(; m = 3),
        description = "Floodplain water volume",
        tags = [:local_inertial_floodplain_1D_flow_output],
    ),
    "floodplain_water__volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.river_flow.floodplain.variables.q_av),
        unit = Unit(; m = 3, s = -1),
        description = "Floodplain discharge",
        tags = [:local_inertial_floodplain_1D_flow_output],
    ),
    ### Overland flow
    #### Input
    "land_surface_water_flow__ground_elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        fill = 0.0,
        description = "Elevation of each cell",
        tags = [:local_inertial_overland_input],
    ),
    "land_surface__slope" => ParameterMetadata(;
        lens = @optic(_.domain.land.parameters.slope),
        unit = Unit(; m = (1, 1)),
        description = "Land surface slope",
        tags = [:kinematic_wave_overland_input, :kinematic_lateral_subsurface_input],
    ),
    #### States
    "land_surface_water__x_component_of_instantaneous_volume_flow_rate" =>
        ParameterMetadata(;
            lens = @optic(_.routing.overland_flow.variables.qx),
            unit = Unit(; m = 3, s = -1),
            description = "Flow in x direction",
            tags = [:local_inertial_overland_state, :local_inertial_overland_output],
        ),
    "land_surface_water__y_component_of_instantaneous_volume_flow_rate" =>
        ParameterMetadata(;
            lens = @optic(_.routing.overland_flow.variables.qy),
            unit = Unit(; m = 3, s = -1),
            description = "Flow in y direction",
            tags = [:local_inertial_overland_state, :local_inertial_overland_output],
        ),
    ## Groundwater flow
    ### Unconfined aquifer
    #### Input
    "land_surface__elevation" => ParameterMetadata(;
        unit = Unit(; m = 1),
        description = "Land surface elevation",
        tags = [:groundwater_unconfined_aquifer_input],
    ),
    "subsurface_surface_water__horizontal_saturated_hydraulic_conductivity" =>
        ParameterMetadata(;
            lens = @optic(_.routing.subsurface_flow.aquifer.parameters.k),
            unit = Unit(; m = 1, d = -1),
            description = "Horizontal conductivity",
            tags = [:groundwater_unconfined_aquifer_input],
        ),
    "subsurface_water__specific_yield" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.aquifer.parameters.specific_yield),
        description = "Specific yield",
        tags = [:groundwater_unconfined_aquifer_input],
    ),
    "subsurface__horizontal_saturated_hydraulic_conductivity_scale_parameter" =>
        ParameterMetadata(;
            lens = @optic(_.routing.subsurface_flow.aquifer.parameters.f),
            unit = Unit(; m = -1),
            description = "Factor controlling the reduction of horizontal conductivity with depth",
            tags = [:groundwater_unconfined_aquifer_input],
        ),
    #### States
    "subsurface_water__hydraulic_head" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.aquifer.variables.head),
        unit = Unit(; m = 1),
        description = "Groundwater head",
        tags = [
            :groundwater_unconfined_aquifer_state,
            :groundwater_unconfined_aquifer_output,
        ],
    ),
    ### River boundary
    #### Input
    "river_water__infiltration_conductance" => ParameterMetadata(;
        lens = @optic(
            _.routing.subsurface_flow.boundaries.river.parameters.infiltration_conductance
        ),
        unit = Unit(; m = 2, d = -1),
        description = "River bed infiltration conductance",
        tags = [:groundwater_river_boundary_input],
    ),
    "river_water__exfiltration_conductance" => ParameterMetadata(;
        lens = @optic(
            _.routing.subsurface_flow.boundaries.river.parameters.exfiltration_conductance
        ),
        unit = Unit(; m = 2, d = -1),
        description = "River bed exfiltration conductance",
        tags = [:groundwater_river_boundary_input],
    ),
    "river_bottom__elevation" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.river.parameters.bottom),
        unit = Unit(; m = 1),
        description = "River bottom elevation",
        tags = [:groundwater_river_boundary_input],
    ),
    #### Output
    "river_water__to_subsurface_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.river.variables.flux_av),
        unit = Unit(; m = 3, d = -1),
        description = "Exchange flux (river to aquifer)",
        tags = [:groundwater_river_boundary_output],
    ),
    ### Drainage boundary
    #### Input
    "land_drain__elevation" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.parameters.elevation),
        unit = Unit(; m = 1),
        fill = MISSING_VALUE,
        description = "Drain elevation",
        tags = [:groundwater_drainage_boundary_input],
    ),
    "land_drain__conductance" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.parameters.conductance),
        unit = Unit(; m = 2, d = -1),
        fill = MISSING_VALUE,
        description = "Drain conductance",
        tags = [:groundwater_drainage_boundary_input],
    ),
    "land_drain_location__mask" => ParameterMetadata(;
        fill = false,
        description = "Drain location",
        tags = [:groundwater_drainage_boundary_input],
    ),
    #### Output
    "land_drain_water__to_subsurface_volume_flow_rate" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.boundaries.drain.variables.flux_av),
        unit = Unit(; m = 3, d = -1),
        description = "Exchange flux (drain to aquifer)",
        tags = [:groundwater_drainage_boundary_output],
    ),
    ### Recharge boundary
    #### Output
    "subsurface_water_saturated_zone_top__net_recharge_volume_flow_rate" =>
        ParameterMetadata(;
            lens = @optic(_.routing.subsurface_flow.boundaries.recharge.variables.flux_av),
            unit = Unit(; m = 3, d = -1),
            description = "Net recharge flux",
            tags = [:groundwater_recharge_boundary_output],
        ),
    ### Constant head boundary
    #### Input
    "model_constant_boundary_condition__hydraulic_head" => ParameterMetadata(;
        lens = @optic(_.routing.subsurface_flow.constanthead.variables.head),
        unit = Unit(; m = 1),
        fill = MISSING_VALUE,
        description = "Head of the boundary",
        tags = [:groundwater_constant_head_boundary_input],
    ),
)
