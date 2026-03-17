get_states_by_tag(model_type::Type, tag::Symbol) =
    [name for (name, metadata) in get_standard_name_map(model_type) if tag in metadata.tags]

"""
    get_snow_states(model_type::ModelType.T)

Extract required snow model states, given a certain `model_type`. Returns a tuple with the
required states (standard name).
"""
function get_snow_states(model_type::ModelType.T)
    return if model_type == ModelType.sbm || model_type == ModelType.sbm_gwf
        get_states_by_tag(LandHydrologySBM, :snow_state)
    elseif model_type == ModelType.sediment
        String[]
    end
end

"""
    get_glacier_states(model_type::ModelType.T)

Extract required glacier model states, given a certain `model_type`. Returns a tuple with
the required states (standard name).
"""
function get_glacier_states(model_type::ModelType.T)
    return if model_type == ModelType.sbm || model_type == ModelType.sbm_gwf
        get_states_by_tag(LandHydrologySBM, :glacier_state)
    elseif model_type == ModelType.sediment
        String[]
    end
end

"""
    get_interception_states(model_type::ModelType.T)

Extract required interception model states, given a certain `model_type`. Returns a tuple
with the required states (standard name).
"""
function get_interception_states(model_type::ModelType.T)
    return if model_type == ModelType.sbm || model_type == ModelType.sbm_gwf
        get_states_by_tag(LandHydrologySBM, :vegetation_state)
    elseif model_type == ModelType.sediment
        String[]
    end
end

"""
    get_soil_states(model_type::ModelType.T; snow = false)

Extract required soil model states, given a certain `model_type` and whether `snow` is
modelled. Returns a tuple with the required states (internal names as symbols).
"""
function get_soil_states(model_type::ModelType.T; snow::Bool=false)
    return if model_type == ModelType.sbm || model_type == ModelType.sbm_gwf
        states = get_states_by_tag(LandHydrologySBM, :soil_state)
        snow ? states : filter(!=("soil_surface__temperature"), states)
    elseif model_type == ModelType.sediment
        String[]
    end
end

function get_sediment_states()
    return get_states_by_tag(SoilLoss, :sediment_river_transport_state)
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
    model_type = config.model.type

    # Extract model settings
    do_snow = config.model.snow__flag
    do_glaciers = config.model.glacier__flag
    do_reservoirs = config.model.reservoir__flag
    do_floodplains = config.model.floodplain_1d__flag
    do_paddy = config.model.water_demand.paddy__flag

    # Extract required stated based on model configuration file
    if do_snow
        snow_states = get_snow_states(model_type)
    else
        snow_states = String[]
    end
    if do_snow && do_glaciers
        glacier_states = get_glacier_states(model_type)
    else
        glacier_states = String[]
    end
    interception_states = get_interception_states(model_type)
    soil_states = get_soil_states(model_type; snow=do_snow)

    # Subsurface states
    ssf_states = if model_type == ModelType.sbm_gwf
        ["subsurface_water__hydraulic_head"]
    elseif model_type == ModelType.sbm
        ["subsurface_water__volume_flow_rate"]
    else # model_type == ModelType.sediment
        String[]
    end

    # Land states
    land_states = if model_type == ModelType.sediment
        String[]
    else
        if config.model.land_routing == RoutingType.local_inertial
            get_states_by_tag(Routing, :local_inertial_overland_state)
        else
            get_states_by_tag(Routing, :kinematic_wave_overland_state)
        end
    end

    # River states
    river_states = if model_type == ModelType.sediment
        river_states = get_sediment_states()
    else
        # :local_inertial_river_state yields the same results
        get_states_by_tag(Routing, :kinematic_wave_river_state)
    end

    # Floodplain states
    floodplain_states =
        do_floodplains ?
        get_states_by_tag(Routing, :local_inertial_floodplain_1D_flow_state) : String[]

    # Reservoir states
    reservoir_states =
        do_reservoirs && model_type !== ModelType.sediment ?
        get_states_by_tag(Routing, :reservoir_state) : String[]

    # Paddy states
    paddy_states =
        do_paddy ? get_states_by_tag(LandHydrologySBM, :demand_paddy_irrigation_state) :
        String[]

    required_states = vcat(
        snow_states,
        glacier_states,
        interception_states,
        soil_states,
        ssf_states,
        land_states,
        river_states,
        floodplain_states,
        reservoir_states,
        paddy_states,
    )

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
    missing_states = String[]
    # Check if all states are covered
    for state in required_states
        if !haskey(state_ncnames, state)
            @error string(
                "State variable `$state` not provided, please ensure that it is set ",
                "correctly in the model configuration",
            )
            error = true
            missing_states = push!(missing_states, state)
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
