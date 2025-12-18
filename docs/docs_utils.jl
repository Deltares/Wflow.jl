using Wflow
using PrettyTables

function generate_table(model_type::Type, filter_flag::Symbol)
    standard_name_map = Wflow.standard_name_map(model_type)
    data = Vector{String}[]

    for (standard_name, metadata) in standard_name_map
        (; description, unit, default, flags) = metadata
        if filter_flag ∈ flags
            push!(
                data,
                [
                    "<code>$standard_name</code>",
                    description,
                    string(unit),
                    isnothing(default) ? "-" : string(default),
                ],
            )
        end
    end
    isempty(data) && error("No data for table.")
    data = permutedims(hcat(data...))
    pretty_table(
        data;
        column_labels = ["Standard name", "Description", "Unit", "Default"],
        alignment = :l,
        backend = :html,
        allow_html_in_cells = true,
    )
end
