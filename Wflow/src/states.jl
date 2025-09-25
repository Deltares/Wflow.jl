"""
    get_snow_states(model_type::AbstractString)

Extract required snow model states, given a certain `model_type`. Returns a tuple with the
required states (standard name).
"""
function get_snow_states(model_type::AbstractString)
    if model_type == "sbm" || model_type == "sbm_gwf"
        states = ("snowpack_dry_snow__leq_depth", "snowpack_liquid_water__depth")
    elseif model_type == "sediment"
        states = ()
    else
        throw(ArgumentError("Unknown model_type provided (`$model_type`)"))
    end
    return states
end

"""
    get_glacier_states(model_type::AbstractString)

Extract required glacier model states, given a certain `model_type`. Returns a tuple with
the required states (standard name).
"""
function get_glacier_states(model_type::AbstractString)
    if model_type == "sbm" || model_type == "sbm_gwf"
        states = ("glacier_ice__leq_depth",)
    elseif model_type == "sediment"
        states = ()
    else
        throw(ArgumentError("Unknown model_type provided (`$model_type`)"))
    end
    return states
end

"""
    get_interception_states(model_type::AbstractString)

Extract required interception model states, given a certain `model_type`. Returns a tuple
with the required states (standard name).
"""
function get_interception_states(model_type::AbstractString)
    if model_type == "sbm" || model_type == "sbm_gwf"
        states = ("vegetation_canopy_water__depth",)
    elseif model_type == "sediment"
        states = ()
    else
        throw(ArgumentError("Unknown model_type provided (`$model_type`)"))
    end
    return states
end

"""
    get_soil_states(model_type::AbstractString; snow = false)

Extract required soil model states, given a certain `model_type` and whether `snow` is
modelled. Returns a tuple with the required states (internal names as symbols).
"""
function get_soil_states(model_type::AbstractString; snow = false)
    if model_type == "sbm" || model_type == "sbm_gwf"
        if snow
            states = (
                "soil_water_saturated_zone__depth",
                "soil_surface__temperature",
                "soil_layer_water_unsaturated_zone__depth",
            )
        else
            states = (
                "soil_water_saturated_zone__depth",
                "soil_layer_water_unsaturated_zone__depth",
            )
        end
    elseif model_type == "sediment"
        states = ()
    else
        throw(ArgumentError("Unknown model_type provided (`$model_type`)"))
    end
    return states
end

function get_sediment_states()
    states = (
        "river_water_clay__mass",
        "river_bed_clay__mass",
        "river_water_gravel__mass",
        "river_bed_gravel__mass",
        "river_water_large_aggregates__mass",
        "river_bed_large_aggregates__mass",
        "river_water_clay__mass_flow_rate",
        "river_water_gravel__mass_flow_rate",
        "river_water_large_aggregates__mass_flow_rate",
        "river_water_small_aggregates__mass_flow_rate",
        "river_water_sand__mass_flow_rate",
        "river_water_silt__mass_flow_rate",
        "river_water_small_aggregates__mass",
        "river_bed_small_aggregates__mass",
        "river_water_sand__mass",
        "river_bed_sand__mass",
        "river_water_silt__mass",
        "river_bed_silt__mass",
    )
    return states
end

"""
    extract_required_states(config::Config)

Function to retrieve the required states, given the model configuration. The required states
are inferred from the model settings (mainly model_type, model options and routing types).
Returns as list of required states, in the same formats as the keys that are returned from the
`ncnames` function.
"""
function extract_required_states(config::Config)
    # Extract model type
    model_type = config.model.type::String

    # Extract model settings
    do_snow = get(config.model, "snow__flag", false)::Bool
    do_glaciers = get(config.model, "glacier__flag", false)::Bool
    do_reservoirs = get(config.model, "reservoir__flag", false)::Bool
    do_floodplains = get(config.model, "floodplain_1d__flag", false)::Bool
    do_paddy = false
    if haskey(config.model, "water_demand")
        do_paddy = get(config.model.water_demand, "paddy__flag", false)::Bool
    end

    # Extract required stated based on model configuration file
    if do_snow
        snow_states = get_snow_states(model_type)
    else
        snow_states = ()
    end
    if do_snow && do_glaciers
        glacier_states = get_glacier_states(model_type)
    else
        glacier_states = ()
    end
    interception_states = get_interception_states(model_type)
    soil_states = get_soil_states(model_type; snow = do_snow)

    # Subsurface states
    if model_type == "sbm_gwf"
        ssf_states = ("subsurface_water__instantaneous_hydraulic_head",)
    elseif model_type == "sbm"
        ssf_states = ("subsurface_water__volume_flow_rate",)
    else
        ssf_states = ()
    end

    # Land states
    if model_type == "sediment"
        land_states = ()
    else
        routing_options = ("kinematic_wave", "local_inertial")
        land_routing = get_options(
            config.model,
            "land_routing",
            routing_options,
            "kinematic_wave",
        )::String
        if land_routing == "local_inertial"
            land_states = (
                "land_surface_water__x_component_of_instantaneous_volume_flow_rate",
                "land_surface_water__y_component_of_instantaneous_volume_flow_rate",
                "land_surface_water__instantaneous_depth",
            )
        else
            land_states = (
                "land_surface_water__instantaneous_volume_flow_rate",
                "land_surface_water__instantaneous_depth",
            )
        end
    end

    # River states
    if model_type == "sediment"
        river_states = get_sediment_states()
    else
        river_states = (
            "river_water__instantaneous_volume_flow_rate",
            "river_water__instantaneous_depth",
        )
    end

    # Floodplain states
    floodplain_states =
        do_floodplains ?
        (
            "floodplain_water__instantaneous_volume_flow_rate",
            "floodplain_water__instantaneous_depth",
        ) : ()

    # Reservoir states
    reservoir_states =
        do_reservoirs && model_type !== "sediment" ?
        ("reservoir_water_surface__instantaneous_elevation",) : ()

    # Paddy states
    paddy_states = do_paddy ? ("paddy_surface_water__depth",) : ()

    # Add required states to a tuple, similar to the keys in the output of
    # `ncnames(config.state.variables)`
    required_states = snow_states...,
    glacier_states...,
    interception_states...,
    soil_states...,
    ssf_states...,
    land_states...,
    river_states...,
    floodplain_states...,
    reservoir_states...,
    paddy_states...

    return required_states
end

"""
    check_states(config::Config)

This function checks whether the states provided in the model configuration are covering all
required states. If not all states are covered, an error is thrown to warn the user
(including logging messages to explain which variables are missing). If more states are
provided than required, a warning is logged, and the unnecessary state is removed.
"""
function check_states(config::Config)
    state_ncnames = ncnames(config.state.variables)

    # Get required states
    required_states = extract_required_states(config)

    # Flag to keep track of whether to throw an error
    error = false
    missing_states = ()
    # Check if all states are covered
    for state in required_states
        if !haskey(state_ncnames, state)
            @error string(
                "State variable `$state` not provided, please ensure that it is set ",
                "correctly in the model configuration",
            )
            error = true
            missing_states = (missing_states..., state)
        end
    end
    # Throw error when not all states are covered
    if error
        throw(
            ArgumentError(
                "Not all required states are provided, these states are missing: `$missing_states`",
            ),
        )
    end

    # Check if more states are passed than required, and warn user if this is the case
    for (state, _) in state_ncnames
        if !(state in required_states)
            @warn "State variable `$state` provided, but is not used in model setup, skipping."
            delete!(state_ncnames, state)
        end
    end
    return state_ncnames
end
