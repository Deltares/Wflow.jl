using Wflow
using JSON

function to_dict(standard_name_map)
    out = Dict{String, Dict{String, Any}}()

    for (standard_name, metadata) in standard_name_map
        out[standard_name] = Dict(
            "unit" => Wflow.to_string(metadata.unit; BMI_standard = true),
            "default" => metadata.default,
            "fill" => metadata.fill,
            "type" => metadata.type,
            "description" => metadata.description,
            "allow_missing" => metadata.allow_missing,
            "allow_dynamic_input" => metadata.allow_dynamic_input,
        )
    end

    return out
end

for (name, standard_name_map) in Wflow.STANDARD_NAME_MAPS
    dict = to_dict(standard_name_map)
    path = normpath(@__DIR__, "$(name)_metadata.json")
    open(path, "w") do io
        print(io, JSON.json(dict, 2))
    end
    @info "Wrote $name parameter/variable metadata to $path."
end
