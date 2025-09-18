argument_error(msg::String) = throw(ArgumentError(msg))

# Use convert_value wrapper to avoid type piracy
convert_value(::Type{<:PropertyDict}, dict::AbstractDict{String}) = PropertyDict(dict)
convert_value(::Type{LogLevel}, level::Int) = LogLevel(level)
convert_value(::Type{Vector{T}}, v::Vector) where {T <: AbstractConfigSection} =
    init_config_section.(T, v)

function convert_value(::Type{T}, value::Any) where {T}
    return if value isa T
        value
    elseif T isa Union
        T_ = get_something_type(T)
        if T_ <: AbstractConfigSection
            init_config_section(T_, value)
        else
            convert(T_, value)
        end
    else
        convert(T, value)
    end
end

convert_value(::Type{Union{Nothing, DateTime}}, s::AbstractString) = DateTime(s)

function convert_value(::Type{LogLevel}, level::String)
    level = get(log_level_map, lowercase(level), missing)
    ismissing(level) && argument_error(
        "Invalid loglevel = $input_string in the [logging] section of the TOML, must be one of $(keys(log_level_map)).",
    )
    return level
end

const no_override = Dict{String, Any}()

"""
Initialize the configuration section of type T<:AbstractConfigSection. If T itself contains fields of 
type <:AbstractConfigSection, this function works recursively.
"""
function init_config_section(
    ::Type{T},
    dict::AbstractDict{String};
    override::AbstractDict{String} = no_override,
)::T where {T <: AbstractConfigSection}
    args = Dict{Symbol, Any}()
    all_target_fields = fieldnames(T)

    # Check validity of override
    for key in keys(override)
        if Symbol(key) ∉ all_target_fields
            error("Found invalid override key $key for section of type $T.")
        end
    end

    for (option_string, value) in dict
        option = Symbol(option_string)
        if option in all_target_fields
            target_arg_type = fieldtype(T, option)
            # Recursively call this function if the target type is a AbstractConfigSection
            args[option] = if target_arg_type <: AbstractConfigSection
                init_config_section(
                    target_arg_type,
                    value;
                    override = get(override, option_string, no_override),
                )
            else
                # Take override value if available
                value = get(override, option_string, value)
                try
                    convert_value(target_arg_type, value)
                catch
                    argument_error(
                        "Couldn't convert the value $option_string = $value in the TOML section $(section_name(T)) " *
                        "to the target type $target_arg_type.",
                    )
                end
            end
        else
            @warn "'$option_string' is not recognized as a valid field of the $(section_name(T)) section in the TOML, " *
                  "this will be ignored."
        end
    end

    try
        return T(; args...)
    catch
        argument_error(
            "Couldn't parse the $(section_name(T)) section in the TOML with data $dict, check for missing mandatory fields.",
        )
    end
end

init_config_section_default(
    type::Type{<:AbstractConfigSection},
    dict::AbstractDict{String};
    override::AbstractDict{String} = no_override,
) = invoke(
    init_config_section,
    Tuple{Type{<:AbstractConfigSection}, AbstractDict{String}},
    type,
    dict;
    override,
)

init_config_section(
    ::Type{T},
    data::Any;
    override::AbstractDict{String} = no_override,
) where {T <: AbstractConfigSection} = argument_error(
    "Couldn't parse the $(section_name(T)) section in the TOML with data $data, must be a section.",
)

# Option 1 and 2 (see struct InputEntry)
function init_config_section(
    ::Type{InputEntry},
    dict::AbstractDict{String};
    override::AbstractDict{String} = no_override,
)
    netcdf_variable_name = get_nested(dict, ["netcdf", "variable", "name"])
    value = get(dict, "value", nothing)

    if !isnothing(netcdf_variable_name)
        dict["netcdf_variable_name"] = netcdf_variable_name
        pop!(dict, "netcdf")
        # Option 1
        len = missing

        for field_name in ("scale", "offset", "layer")
            if haskey(dict, field_name)
                value = dict[field_name]
                # Convert scalar to vector
                value = (value isa Number) ? [value] : value
                dict[field_name] = value
                if !ismissing(len) && length(value) != len
                    throw(
                        ArgumentError(
                            "In the TOML [input] section with netcdf variable name $netcdf_variable_name," *
                            "the 'scale' 'offset' and 'layer' values must all have the same length if provided.",
                        ),
                    )
                else
                    len = length(value)
                end
            end
        end

        # Add a proper amount of the default values
        len = coalesce(len, 0)
        !haskey(dict, "scale") && (dict["scale"] = ones(len))
        !haskey(dict, "offset") && (dict["offset"] = zeros(len))
        has_layer = haskey(dict, "layer")

        for (i, (scale, offset)) in enumerate(zip(dict["scale"], dict["offset"]))
            # Modifier info
            msg = "NetCDF parameter $netcdf_variable_name is modified with scale $scale and offset $offset"
            has_layer && (msg *= "at index $(dict["layer"][i])")
            @info "$msg."
        end

        # Invoke default method
        return init_config_section_default(InputEntry, dict; override)
    elseif !isnothing(value)
        # Option 2
        return InputEntry(; value)
    end
    argument_error("Couldn't parse the section under [input] in the TOML with data $dict.")
end

# Option 3 (see struct InputEntry)
init_config_section(
    ::Type{InputEntry},
    standard_name::String;
    override::AbstractDict{String} = no_override,
) = InputEntry(; standard_name)
Base.convert(::Type{InputEntry}, standard_name::String) = InputEntry(; standard_name)

convert_value(::Type{IndexSection}, i::Int) = IndexSection(; i)
convert_value(::Type{IndexSection}, dict::AbstractDict{String}) =
    init_config_section_default(IndexSection, dict)
init_config_section(::Type{IndexSection}, i::Int; override = no_override) =
    IndexSection(; i)

function init_config_section(
    ::Type{InputEntries},
    dict::AbstractDict{String};
    override::AbstractDict{String} = no_override,
)
    InputEntries(
        Dict(key => init_config_section(InputEntry, value) for (key, value) in dict),
    )
end

function init_config_section(
    ::Type{InputSection},
    dict::AbstractDict{String};
    override::AbstractDict{String} = no_override,
)
    # Move flexible part of the input section into input.flexible
    flexible = Dict{String, Any}()
    for key in collect(keys(dict))
        if key ∉ input_field_names
            flexible[key] = pop!(dict, key)
        end
    end
    dict["flexible"] = flexible
    # Invoke default method
    input = init_config_section_default(InputSection, dict; override)

    # check if there is overlap in forcing and cyclic parameters
    if !haskey(input.cyclic, "_was_not_specified") # do_cyclic
        overlap = intersect(keys(input.forcing), keys(input.cyclic))
        if !isempty(overlap)
            argument_error(
                "These parameters were specified in both the forcing and cyclic" *
                "sections of the [input] in the TOML: $overlap.",
            )
        end
    end

    return input
end

function init_config_section(
    ::Type{NetCDFScalarVariable},
    dict::AbstractDict{String};
    override::AbstractDict{String} = no_override,
)
    for key in ("name", "parameter")
        !haskey(dict, key) && throw(
            ArgumentError(
                "A variable in a [[output.netcdf_scalar.variable]] section in the TOML is missing the mandatory field $key.",
            ),
        )
    end

    if haskey(dict, "map")
        dict["_location_dim"] = string(dict["name"], '_', dict["map"])
    elseif haskey(dict, "location")
        dict["_location_dim"] = dict["location"]
    else
        throw(
            ArgumentError(
                "The variable in the [[output.netcdf_scalar.variable]] section with name = '$name'" *
                " is missing either a 'location' or a 'map' field.",
            ),
        )
    end
    # Invoke default method
    return init_config_section_default(NetCDFScalarVariable, dict; override)
end

"""
    Config(path::AbstractString
    Config(dict::AbstractDict; path::Union{Nothing, String} = nothing)

Struct that contains the parsed TOML configurations, as well as a reference to the TOML path
if it exists. The object behaves largely like an immutable nested struct, apart from the flexible
sections (see the fields with `PropertyDictType` type in `config.jl`) which are nested property
dicts (see PropertyDicts.jl) and thus mutable.
An override of the data in the TOML or dict can be passed as a dictionary,
e.g. override = Dict("time" => Dict("timestepsecs" => 42))
Fields whose name start with an underscore should not be specified in the TOML.
"""
function Config(path::AbstractString; override::AbstractDict{String} = no_override)
    dict = TOML.parsefile(path)
    Config(dict; path, override)
end

function Config(
    dict::AbstractDict;
    path::Union{Nothing, String} = nothing,
    override::AbstractDict{String} = no_override,
)
    # Add path to config
    dict["path"] = path
    config = init_config_section(Config, dict; override)
    return config
end

## Convert Config back to dict
function to_dict(
    config_section::T;
    dict = Dict{String, Any}(),
) where {T <: AbstractConfigSection}
    for field_name in fieldnames(T)
        value = getfield(config_section, field_name)
        if (field_name ∈ (:flexible, :_was_specified)) || isnothing(value)
            continue
        end
        dict[String(field_name)] = to_dict(value)
    end
    return dict
end

to_dict(input::InputSection) =
    invoke(to_dict, Tuple{AbstractConfigSection}, input; dict = deepcopy(input.flexible))
to_dict(input_entries::InputEntries) =
    Dict{String, Any}(name => to_dict(entry) for (name, entry) in input_entries.dict)

function to_dict(input_entry::InputEntry)
    (; netcdf_variable_name, standard_name) = input_entry
    dict = invoke(to_dict, Tuple{AbstractConfigSection}, input_entry)
    if !isnothing(input_entry.standard_name)
        return standard_name
    end
    if !isnothing(netcdf_variable_name)
        dict["netcdf"] = Dict{String, Any}(
            "variable" => Dict{String, Any}("name" => netcdf_variable_name),
        )
        pop!(dict, "netcdf_variable_name")
    end
    return dict
end

to_dict(data::Vector) = to_dict.(data)

to_dict(x::Any) = x
to_dict(x::EnumX.Enum) = remove_leading_underscore(String(Symbol(x)))
to_dict(x::RoutingType.T) =
    (x == RoutingType.kinematic_wave) ? "kinematic-wave" : "local-inertial"
to_dict(x::LogLevel) = lowercase(string(x))