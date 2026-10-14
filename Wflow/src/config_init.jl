argument_error(msg::String) = throw(ArgumentError(msg))

"""
    toml_string(value)::String

String representation of a value parsed from TOML, for error messages. Uses explicit type
checks to avoid dynamic dispatch on values of unknown type.
"""
function toml_string(@nospecialize(value))::String
    return if value isa String
        value
    elseif value isa Int
        string(value)
    elseif value isa Float64
        string(value)
    elseif value isa Bool
        string(value)
    elseif value isa DateTime
        format_datetime(value)
    elseif value isa AbstractDict
        "<table>"
    elseif value isa AbstractVector
        "<array>"
    else
        "<value>"
    end
end

"Length of a TOML value, 1 for scalars."
function toml_length(@nospecialize(value))::Int
    return if value isa Vector{Any}
        length(value)
    elseif value isa Vector{Float64}
        length(value)
    elseif value isa Vector{Int}
        length(value)
    else
        1
    end
end

conversion_error(::Type{T}, @nospecialize(value)) where {T} =
    argument_error("Cannot convert `$(toml_string(value))` to $T.")

missing_field_error(::Type{T}, key::String) where {T} = argument_error(
    "Couldn't parse the $(section_name(T)) section in the TOML, the mandatory field `$key` is missing.",
)

function warn_unknown_keys(::Type{T}, dict::AbstractDict{String})::Nothing where {T}
    for key in keys(dict)
        if Symbol(key) ∉ fieldnames(T)
            @warn "'$key' is not recognized as a valid field of the $(section_name(T)) section in the TOML, " *
                "this will be ignored."
        end
    end
    return nothing
end

function convert_field(
        ::Type{S},
        key::String,
        ::Type{T},
        @nospecialize(value),
    )::T where {S, T}
    try
        return convert_value(T, value)
    catch e
        e isa ArgumentError || rethrow()
        argument_error(
            "Couldn't convert the value $key = $(toml_string(value)) in the TOML section $(section_name(S)) " *
                "to the target type $T.",
        )
    end
end

#=
Conversion of values parsed from TOML to the types of the configuration fields. Methods
only dispatch on the target type, and check the value type explicitly. Values from a parsed
TOML have type `Any`, so dispatching on them would require dynamic dispatch.
=#

convert_value(::Type{Any}, @nospecialize(value)) = value
convert_value(::Type{Union{Nothing, T}}, @nospecialize(value)) where {T} = convert_value(T, value)

convert_value(::Type{String}, @nospecialize(value))::String =
    value isa String ? value : conversion_error(String, value)

convert_value(::Type{Bool}, @nospecialize(value))::Bool =
    value isa Bool ? value : conversion_error(Bool, value)

function convert_value(::Type{Int}, @nospecialize(value))::Int
    return if value isa Int
        value
    elseif value isa Float64 && isinteger(value)
        Int(value)
    else
        conversion_error(Int, value)
    end
end

function convert_value(::Type{Float64}, @nospecialize(value))::Float64
    return if value isa Float64
        value
    elseif value isa Int
        Float64(value)
    else
        conversion_error(Float64, value)
    end
end

function convert_value(::Type{DateTime}, @nospecialize(value))::DateTime
    return if value isa DateTime
        value
    elseif value isa Date
        DateTime(value)
    elseif value isa String
        DateTime(value)
    else
        conversion_error(DateTime, value)
    end
end

convert_value(::Type{VersionNumber}, @nospecialize(value))::VersionNumber =
    value isa String ? VersionNumber(value) : conversion_error(VersionNumber, value)

function convert_value(::Type{LogLevel}, @nospecialize(value))::LogLevel
    if value isa Int
        return LogLevel(value)
    elseif value isa String
        for (name, level) in log_level_map
            name == lowercase(value) && return level
        end
        argument_error(
            "Invalid loglevel = $value in the [logging] section of the TOML, must be one of debug, info, warn, error.",
        )
    else
        conversion_error(LogLevel, value)
    end
end

convert_value(::Type{PropertyDictType}, @nospecialize(value))::PropertyDictType =
    value isa Dict{String, Any} ? PropertyDict(value) :
    conversion_error(PropertyDictType, value)

function convert_value(::Type{Vector{T}}, @nospecialize(value))::Vector{T} where {T}
    value isa AbstractVector || conversion_error(Vector{T}, value)
    return if value isa Vector{T}
        value
    elseif value isa Vector{Any}
        T[convert_value(T, x) for x in value]
    elseif value isa Vector{Int}
        T[convert_value(T, x) for x in value]
    elseif value isa Vector{Float64}
        T[convert_value(T, x) for x in value]
    elseif value isa Vector{String}
        T[convert_value(T, x) for x in value]
    elseif value isa Vector{Dict{String, Any}}
        T[convert_value(T, x) for x in value]
    else
        conversion_error(Vector{T}, value)
    end
end

# Automatically convert strings into enumerator instances
function convert_value(::Type{T}, @nospecialize(value))::T where {T <: EnumX.Enum}
    value isa String || conversion_error(T, value)
    option = Symbol(add_leading_underscore(value))
    for instance in instances(T)
        if option == Symbol(instance)
            return instance
        end
    end
    throw(
        ArgumentError(
            "Invalid value $value in the TOML, must be one of $(join(option_names(T), ", ")).",
        ),
    )
end

function convert_value(::Type{T}, @nospecialize(value))::T where {T <: AbstractConfigSection}
    return if value isa Dict{String, Any}
        init_config_section(T, value)
    else
        init_config_section_from_value(T, value)
    end
end

"""
Initialize the configuration section of type T<:AbstractConfigSection. If T itself contains fields of
type <:AbstractConfigSection, this function works recursively.
"""
init_config_section(::Type{T}, dict::AbstractDict{String}) where {T <: AbstractConfigSection} =
    init_config_section_default(T, dict)

init_config_section_from_value(
    ::Type{T},
    @nospecialize(data),
) where {T <: AbstractConfigSection} =
    argument_error(
    "Couldn't parse the $(section_name(T)) section in the TOML with data $(toml_string(data)), must be a section.",
)

# Option 1 and 2 (see struct InputEntry)
function init_config_section(::Type{InputEntry}, dict::AbstractDict{String})
    netcdf_variable_name = get(dict, "netcdf_variable_name", nothing)
    value = get(dict, "value", nothing)

    if !isnothing(netcdf_variable_name)
        # Option 1
        len = -1

        for field_name in ("scale", "offset", "layer")
            if haskey(dict, field_name)
                field_value = dict[field_name]
                # Convert scalar to vector
                field_value = (field_value isa Number) ? Any[field_value] : field_value
                dict[field_name] = field_value
                field_length = toml_length(field_value)
                if len >= 0 && field_length != len
                    throw(
                        ArgumentError(
                            "In the TOML [input] section with netcdf variable name $(toml_string(netcdf_variable_name))," *
                                "the 'scale' 'offset' and 'layer' values must all have the same length if provided.",
                        ),
                    )
                else
                    len = field_length
                end
            end
        end

        # Add a proper amount of the default values
        len = max(len, 0)
        !haskey(dict, "scale") && (dict["scale"] = ones(len))
        !haskey(dict, "offset") && (dict["offset"] = zeros(len))

        # Invoke default method
        return init_config_section_default(InputEntry, dict)
    elseif value isa String
        # Option 3
        return InputEntry(; external_name = value)
    elseif value isa Number
        # Option 2
        return InputEntry(; value = [convert_value(Float64, value)])
    elseif !isnothing(value)
        # Option 2
        return InputEntry(; value = convert_value(Vector{Float64}, value))
    end
    argument_error("Couldn't parse the section under [input] in the TOML.")
end

# Option 3 (see struct InputEntry), or a uniform value
function init_config_section_from_value(
        ::Type{InputEntry},
        @nospecialize(value),
    )::InputEntry
    return if value isa String
        InputEntry(; external_name = value)
    elseif value isa Number
        InputEntry(; value = [convert_value(Float64, value)])
    else
        argument_error("Couldn't parse the section under [input] in the TOML.")
    end
end
Base.convert(::Type{InputEntry}, external_name::String) = InputEntry(; external_name)
Base.convert(::Type{InputEntry}, value::Number) = InputEntry(; value = [Float64(value)])

init_config_section_from_value(::Type{IndexSection}, @nospecialize(i))::IndexSection =
    i isa Int ? IndexSection(; i) : conversion_error(IndexSection, i)

function init_config_section(::Type{InputEntries}, dict::AbstractDict{String})
    entries = Dict{String, InputEntry}()
    for (key, value) in dict
        entries[key] = convert_value(InputEntry, value)
    end
    return InputEntries(entries)
end

function init_config_section(::Type{InputSection}, dict::AbstractDict{String})
    # Move flexible part of the input section into input.location_maps
    location_maps = Dict{String, Any}()
    for key in collect(keys(dict))
        if key ∉ input_field_names
            location_maps[key] = pop!(dict, key)
        end
    end
    dict["_location_maps"] = location_maps
    # Invoke default method
    input = init_config_section_default(InputSection, dict)

    # Check whether there is overlap in specified parameters
    overlap = intersect(keys(input.static), keys(input.forcing), keys(input.cyclic))
    if !isempty(overlap)
        argument_error(
            "These parameters were specified in more than one of the static, forcing and cyclic" *
                "sections of the [input] in the TOML: $(join(overlap, ", ")).",
        )
    end

    return input
end

function init_config_section(::Type{NetCDFScalarVariable}, dict::AbstractDict{String})
    for key in ("name", "parameter")
        !haskey(dict, key) && throw(
            ArgumentError(
                "A variable in a [[output.netcdf_scalar.variable]] section in the TOML is missing the mandatory field $key.",
            ),
        )
    end

    if haskey(dict, "map")
        dict["_location_dim"] = dict["map"]
    elseif haskey(dict, "location")
        dict["_location_dim"] = dict["location"]
    else
        throw(
            ArgumentError(
                "The variable in the [[output.netcdf_scalar.variable]] section with name = '$(toml_string(dict["name"]))'" *
                    " is missing either a 'location' or a 'map' field.",
            ),
        )
    end
    # Invoke default method
    return init_config_section_default(NetCDFScalarVariable, dict)
end

"""
    Config(path::AbstractString
    Config(dict::AbstractDict; path::AbstractString)

Struct that contains the parsed TOML configurations, as well as a reference to the TOML
path. The object behaves largely like an mutable nested struct. Fields whose name start with
an underscore should not be specified in the TOML.
"""
function Config(path::AbstractString)
    dict = TOML.parsefile(path)
    config = Config(dict; path)
    return config
end

function Config(dict::AbstractDict; path::AbstractString)
    # Add path to config
    dict["path"] = String(path)
    config = init_config_section(Config, dict)
    return config
end

## Convert Config back to dict
function to_dict(
        config_section::T;
        dict = Dict{String, Any}(),
    ) where {T <: AbstractConfigSection}
    for field_name in fieldnames(T)
        value = getfield(config_section, field_name)
        # Skip internals and unspecified fields
        if startswith(string(field_name), "_") || isnothing(value)
            continue
        end
        # Skip subsections that were not specified in the original TOML
        if value isa AbstractConfigSection &&
                hasproperty(value, :_was_specified) &&
                !getfield(value, :_was_specified)
            continue
        end
        if value isa VersionNumber
            value = string(value)
        end
        dict[String(field_name)] = to_dict(value)
    end
    return dict
end

to_dict(input::InputSection) = invoke(
    to_dict,
    Tuple{AbstractConfigSection},
    input;
    dict = deepcopy(input._location_maps),
)
to_dict(input_entries::InputEntries) =
    Dict{String, Any}(name => to_dict(entry) for (name, entry) in input_entries.dict)

function to_dict(input_entry::InputEntry)
    (; external_name) = input_entry
    dict = invoke(to_dict, Tuple{AbstractConfigSection}, input_entry)
    if !isnothing(external_name)
        return external_name
    end
    return dict
end

function to_dict(index::IndexSection)
    # Return compact integer form when only `i` is set
    if !isnothing(index.i) && isnothing(index.x) && isnothing(index.y)
        return index.i
    end
    return invoke(to_dict, Tuple{AbstractConfigSection}, index)
end

to_dict(data::Vector) = to_dict.(data)

to_dict(x::Any) = x
to_dict(x::EnumX.Enum) = remove_leading_underscore(String(Symbol(x)))
to_dict(x::LogLevel) = lowercase(string(x))

Base.setproperty!(
    config_section::T,
    field::Symbol,
    value,
) where {T <: AbstractConfigSection} =
    setfield!(config_section, field, convert_value(fieldtype(T, field), value))
