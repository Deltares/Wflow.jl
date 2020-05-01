"Parsed TOML configuration"
struct Config
    dict::Dict{String,Any}
end

# allows using getproperty, e.g. config.input.time instead of config["input"]["time"]
function Base.getproperty(config::Config, f::Symbol)
    dict = getfield(config, :dict)
    a = dict[String(f)]
    # if it is a Dict, wrap the result in Config to keep the getproperty behavior
    return a isa AbstractDict ? Config(a) : a
end

# also used in autocomplete
Base.propertynames(config::Config) = collect(keys(getfield(config, :dict)))
