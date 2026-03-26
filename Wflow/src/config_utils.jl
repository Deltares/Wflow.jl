Base.convert(::Type{String}, calendar_type::CalendarType.T) =
    remove_leading_underscore(String(Symbol(calendar_type)))

# e.g. section_name(InputSection) = "[input]"
function section_name(::Type{<:AbstractConfigSection})
    name = lowercase(replace(String(nameof(InputSection)), "Section" => ""))
    return "[$name]"
end

# Create a nested PropertyDict to retrieve data from a nested dict as
# if it were a nested struct
function nested_property_dict(dict::AbstractDict)
    new_dict = Dict{String, Any}()
    for (key, value) in dict
        if value isa AbstractDict
            new_dict[key] = nested_property_dict(value)
        else
            new_dict[key] = value
        end
    end

    return PropertyDict(new_dict)
end

# Handle options starting with a number
add_leading_underscore(s::String) = isdigit(first(s)) ? "_$s" : s
remove_leading_underscore(s::String) = (first(s) == '_') ? s[2:end] : s

# Get printable tuple of possible options
option_names(::Type{T}) where {T <: EnumX.Enum} =
    remove_leading_underscore.(String.(Symbol.(instances(T))))

get_something_type(::Type{T}) where {T} = last(Base.uniontypes(T))

# Pretty printing AbstractConfigurationSection instances
function Base.show(io::IO, config_section::AbstractConfigSection)
    first = true
    for field_name in fieldnames(typeof(config_section))
        (field_name == :_was_specified) && continue
        value = getfield(config_section, field_name)
        if !isnothing(value)
            first && (first = false; println(io))
            println(io, "\t\t$field_name\t= $value")
        end
    end
end

# Specialized printing for InputEntry
function Base.show(io::IO, input_entry::InputEntry)
    (; netcdf_variable_name, scale, offset, value, external_name) = input_entry
    if !isnothing(netcdf_variable_name)
        if scale == [1.0] && offset == [0.0]
            print(io, netcdf_variable_name)
        else
            print(io, netcdf_variable_name * " (scale = $scale, offset = $offset)")
        end
    elseif !isnothing(value)
        print(io, value)
    else
        print(io, external_name)
    end
end

# Specialized printing for CoordinateSection
Base.show(io::IO, cs::CoordinateSection) = print(io, "(x = $(cs.x), y = $(cs.y))")

# Specialized printing for IndexSection
Base.isnothing(index::IndexSection) =
    all(isnothing(val) for val in (index.i, index.x, index.y))
Base.show(io::IO, index::IndexSection) =
    isnothing(index.i) ? print(io, "(x = $(index.x), y = $(index.y))") : print(io, index.i)

const log_level_map::Dict{Union{Int, String}, LogLevel} =
    Dict("debug" => Debug, "info" => Info, "warn" => Warn, "error" => Error)

Base.dirname(config::Config) = dirname(config.path)

do_water_demand(config::Config) = config.model.water_demand._was_specified
do_api(config::Config) = config.API._was_specified
do_csv(config::Config) = config.output.csv._was_specified
do_netcdf_scalar(config::Config) = config.output.netcdf_scalar._was_specified
do_netcdf_grid(config::Config) = config.output.netcdf_grid._was_specified
do_cyclic(config::Config) = !haskey(config.input.cyclic, "_was_not_specified")

do_index(index_section::IndexSection) = index_section._was_specified

subsurface_routing(config::Config) =
    config.model.type == ModelType.sbm ? RoutingType.kinematic_wave : nothing

# Pretty printing the config
function Base.show(io::IO, c::Config)
    println(io, "Wflow Config")
    for field in fieldnames(typeof(c))
        f = getfield(c, field)
        f === nothing || println(io, "\t$field\t= $f")
    end
end

"Construct a path relative to both the TOML directory and the optional `dir_input`"
input_path(config::Config, path::AbstractString) =
    normpath(dirname(config), config.dir_input, path)

"Construct a path relative to both the TOML directory and the optional `dir_output`"
output_path(config::Config, path::AbstractString) =
    normpath(dirname(config), config.dir_output, path)

function variable_info(var::InputEntry)
    (; layer, scale, offset) = var

    for (i, (scale, offset)) in enumerate(zip(scale, offset))
        if !isone(scale) || !iszero(offset)
            msg = "NetCDF parameter `$(variable_name(var))` is modified with scale `$scale` and offset `$offset`"
            !isnothing(layer) && (msg *= "at index $(layer[i])")
            @info "$msg."
        end
    end
end
