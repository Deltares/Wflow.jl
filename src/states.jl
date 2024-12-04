"""
    get_snow_states(model_type::AbstractString)

Extract required snow model states, given a certain `model_type`. Returns a tuple with the
required states (internal names as symbols).
"""
function get_snow_states(model_type::AbstractString)
    if model_type == "sbm" || model_type == "sbm_gwf"
        states = (:snow_storage, :snow_water)
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
the required states (internal names as symbols).
"""
function get_glacier_states(model_type::AbstractString)
    if model_type == "sbm" || model_type == "sbm_gwf"
        states = (:glacier_store,)
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
with the required states (internal names as symbols).
"""
function get_interception_states(model_type::AbstractString)
    if model_type == "sbm" || model_type == "sbm_gwf"
        states = (:canopy_storage,)
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
            states = (:satwaterdepth, :tsoil, :ustorelayerdepth)
        else
            states = (:satwaterdepth, :ustorelayerdepth)
        end
    elseif model_type == "sediment"
        states = ()
    else
        throw(ArgumentError("Unknown model_type provided (`$model_type`)"))
    end
    return states
end

"""
    add_to_required_states(required_states::Tuple, key_entry::Tuple, states::Tuple)
    add_to_required_states(required_states::Tuple, key_entry::Tuple, states::Nothing)

Function to iterate through the list of `states` and add these to the `required_states`
variable. This is a tuple of tuples, that contains the correct "prefix" to the state
category (provided by the `key_entry` variable). This is added in front of each element in
the list of `states`, and the `required_states` tuple is returned.

When `nothing` is passed as the `states`, the `required_states` variable is simply returned.
"""
function add_to_required_states(required_states::Tuple, key_entry::Tuple, states::Tuple)
    for state in states
        required_states = (required_states..., (key_entry..., state))
    end

    return required_states
end

function add_to_required_states(required_states::Tuple, key_entry::Tuple, states::Nothing)
    return required_states
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
    do_snow = get(config.model, "snow", false)::Bool
    do_glaciers = get(config.model, "glacier", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_floodplains = get(config.model, "floodplain_1d", false)::Bool
    do_paddy = false
    if haskey(config.model, "water_demand")
        do_paddy = get(config.model.water_demand, "paddy", false)::Bool
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
        ssf_states = (:head,)
    else
        ssf_states = haskey(config.input.lateral, "subsurface") ? (:ssf,) : nothing
    end

    # Land states
    if model_type == "sediment"
        land_states = ()
    else
        routing_options = ("kinematic-wave", "local-inertial")
        land_routing = get_options(
            config.model,
            "land_routing",
            routing_options,
            "kinematic-wave",
        )::String
        if land_routing == "local-inertial"
            land_states = (:qx, :qy, :h, :h_av)
        elseif land_routing == "kinematic-wave"
            if model_type == "sbm" || model_type == "sbm_gwf"
                land_states = (:q, :h, :h_av)
            else
                land_states = (:q, :h)
            end
        end
    end

    # River states
    if model_type == "sediment"
        river_states = (
            :clayload,
            :siltload,
            :sandload,
            :saggload,
            :laggload,
            :gravload,
            :claystore,
            :siltstore,
            :sandstore,
            :saggstore,
            :laggstore,
            :gravstore,
            :outclay,
            :outsilt,
            :outsand,
            :outsagg,
            :outlagg,
            :outgrav,
        )
    elseif model_type == "sbm" || model_type == "sbm_gwf"
        river_states = (:q, :h, :h_av)
    else
        river_states = (:q, :h)
    end

    # Floodplain states
    floodplain_states = do_floodplains ? (:q, :h) : nothing

    # Lake and reservoir states
    lake_states = do_lakes ? (:waterlevel,) : nothing
    reservoir_states = do_reservoirs ? (:volume,) : nothing

    # Paddy states
    paddy_states = do_paddy ? (:h,) : nothing

    # Build required states in a tuple, similar to the keys in the output of
    # `ncnames(config.state)`
    required_states = ()
    # Add snow, glacier, interception and sbm soil states to dict
    required_states =
        add_to_required_states(required_states, (:vertical, :snow, :variables), snow_states)
    required_states = add_to_required_states(
        required_states,
        (:vertical, :glacier, :variables),
        glacier_states,
    )
    required_states = add_to_required_states(
        required_states,
        (:vertical, :interception, :variables),
        interception_states,
    )
    required_states =
        add_to_required_states(required_states, (:vertical, :soil, :variables), soil_states)
    # Add subsurface states to dict
    if model_type == "sbm_gwf"
        key_entry = (:lateral, :subsurface, :flow, :aquifer, :variables)
    else
        key_entry = (:lateral, :subsurface, :variables)
    end
    required_states = add_to_required_states(required_states, key_entry, ssf_states)
    # Add land states to dict
    required_states =
        add_to_required_states(required_states, (:lateral, :land, :variables), land_states)
    # Add river states to dict
    if model_type == "sediment"
        key_entry = (:lateral, :river)
    else
        key_entry = (:lateral, :river, :variables)
    end
    required_states = add_to_required_states(required_states, key_entry, river_states)
    # Add floodplain states to dict
    required_states = add_to_required_states(
        required_states,
        (:lateral, :river, :floodplain, :variables),
        floodplain_states,
    )
    # Add lake states to dict
    required_states = add_to_required_states(
        required_states,
        (:lateral, :river, :boundary_conditions, :lake, :variables),
        lake_states,
    )
    # Add reservoir states to dict
    required_states = add_to_required_states(
        required_states,
        (:lateral, :river, :boundary_conditions, :reservoir, :variables),
        reservoir_states,
    )
    # Add paddy states to dict
    required_states = add_to_required_states(
        required_states,
        (:vertical, :demand, :paddy, :variables),
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
    state_ncnames = ncnames(config.state)

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
