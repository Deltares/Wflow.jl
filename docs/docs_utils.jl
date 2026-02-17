using Wflow
using PrettyTables

# Forcing precompilation of Wflow can help Quarto pick up changes in the
# kernel that do not automatically yield re-precompilation
using Pkg
Pkg.precompile("Wflow")

# Generate table from model type and filter tag
generate_table(model_type::Type, filter_tag::Symbol; kwargs...) =
    generate_table([model_type], filter_tag; kwargs...)

# Generate table from model types and filter tags
function generate_table(model_types::Vector{<:Type}, filter_tag::Symbol; kwargs...)
    data = Tuple{String, Wflow.ParameterMetadata}[]

    for model_type in model_types
        standard_name_map = Wflow.standard_name_map(model_type)
        for (standard_name, metadata) in standard_name_map
            if filter_tag ∈ metadata.tags
                push!(data, (standard_name, metadata))
            end
        end
    end
    return generate_table(data; kwargs...)
end

# Generate table from parameter metadata
function generate_table(
    data::Vector{Tuple{String, Wflow.ParameterMetadata}};
    relative_widths::Union{Vector{Int}, Nothing} = nothing,
)
    isempty(data) && error("No data for table.")
    data_pretty = Vector{String}[]
    for (standard_name, metadata) in data
        (; description, unit, default) = metadata
        push!(
            data_pretty,
            [
                "<code>$standard_name</code>",
                replace(description, r"`([^`]+)`" => s"<code>\1</code>"),
                string(unit),
                isnothing(default) ? "-" : string(default),
            ],
        )
    end

    data_pretty = permutedims(hcat(data_pretty...))

    # Check whether there are any defaults
    if any(!=("-"), view(data_pretty, :, 4))
        column_labels = ["Standard name", "Description", "Unit", "Default"]
    else
        column_labels = ["Standard name", "Description", "Unit"]
        data_pretty = data_pretty[:, 1:3]
        relative_widths = relative_widths[1:3]
    end

    io = IOBuffer()
    pretty_table(
        io,
        data_pretty;
        column_labels,
        alignment = :l,
        backend = :html,
        allow_html_in_cells = true,
        table_class = "table-striped table-hover caption-top table",
    )
    table_html = String(take!(io))

    # Insert column width data
    if !isnothing(relative_widths)
        total = sum(relative_widths)
        percentages = @. (100 * relative_widths) ÷ total

        colgroup =
            "<colgroup>\n" *
            join(["    <col style=\"width: $n%\">" for n in percentages], "\n") *
            "\n  </colgroup>\n  "

        colgroup_loc = findfirst("<thead>", table_html)[1]
        table_html =
            table_html[1:(colgroup_loc - 1)] * colgroup * table_html[colgroup_loc:end]
    end

    display(MIME("text/html"), HTML(table_html))
end
