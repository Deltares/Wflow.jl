"""
    get_vertical_states(model_type; snow = false, glacier = false)

Function to extract all required vertical states, given a certain model type. Passes the snow
and glacier options only for the `sbm` model_type. Returns a tuple with the required states
(internal names as symbols)
"""
function get_vertical_states(model_type::AbstractString; snow = false, glacier = false)
    if model_type == "sbm" || model_type == "sbm_gwf"
        if snow && glacier
            vertical_states = (
                :satwaterdepth,
                :snow,
                :tsoil,
                :ustorelayerdepth,
                :snowwater,
                :canopystorage,
                :glacierstore,
            )
        elseif snow
            vertical_states = (
                :satwaterdepth,
                :snow,
                :tsoil,
                :ustorelayerdepth,
                :snowwater,
                :canopystorage,
            )
        else
            vertical_states = (:satwaterdepth, :ustorelayerdepth, :canopystorage)
        end
    elseif model_type == "sediment"
        vertical_states = ()
    else
        throw(ArgumentError("Unknown model_type provided (`$model_type`)"))
    end
    return vertical_states
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
    vertical_states = get_vertical_states(model_type; snow = do_snow, glacier = do_glaciers)

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
            :leftover_clay,
            :leftover_silt,
            :leftover_sand,
            :leftover_sagg,
            :leftover_lagg,
            :leftover_gravel,
            :store_clay,
            :store_silt,
            :store_sand,
            :store_sagg,
            :store_lagg,
            :store_gravel,
            :clay,
            :silt,
            :sand,
            :sagg,
            :lagg,
            :gravel,
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
    # Add vertical states to dict
    required_states = add_to_required_states(required_states, (:vertical,), vertical_states)
    # Add subsurface states to dict
    if model_type == "sbm_gwf"
        key_entry = (:lateral, :subsurface, :flow, :aquifer)
    else
        key_entry = (:lateral, :subsurface)
    end
    required_states = add_to_required_states(required_states, key_entry, ssf_states)
    # Add land states to dict
    required_states =
        add_to_required_states(required_states, (:lateral, :land), land_states)
    # Add river states to dict
    required_states =
        add_to_required_states(required_states, (:lateral, :river), river_states)
    # Add floodplain states to dict
    required_states = add_to_required_states(
        required_states,
        (:lateral, :river, :floodplain),
        floodplain_states,
    )
    # Add lake states to dict
    required_states =
        add_to_required_states(required_states, (:lateral, :river, :lake), lake_states)
    # Add reservoir states to dict
    required_states = add_to_required_states(
        required_states,
        (:lateral, :river, :reservoir),
        reservoir_states,
    )
    # Add paddy states to dict
    required_states =
        add_to_required_states(required_states, (:vertical, :paddy), paddy_states)
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
