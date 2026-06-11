using Statistics: mean
using SpecialFunctions: expint
using StaticArrays: SVector

"Return the first row of a Wflow output CSV file as a NamedTuple"
function csv_first_row(path)
    # silly function to avoid needing CSV.jl as a test dependency
    header, dataline = open(path) do io
        header = readline(io)
        dataline = readline(io)
        (header, dataline)
    end

    names = Tuple(Symbol.(split(header, ',')))
    ncol = length(names)
    # this assumes the first column is a time, the rest a float
    types = Tuple{DateTime, fill(Float64, ncol - 1)...}

    parts = split(dataline, ',')
    values = parse.(Float64, parts[2:end])
    row = NamedTuple{names, types}((DateTime(parts[1]), values...))
    return row
end

function run_piave(model, steps)
    q = zeros(steps)
    ssf_storage = zeros(steps)
    riv_storage = zeros(steps)
    for i in 1:steps
        Wflow.run_timestep!(model)
        ssf_storage[i] = mean(model.routing.subsurface_flow.variables.storage)
        riv_storage[i] = mean(model.routing.river_flow.variables.storage)
        q[i] = model.routing.river_flow.variables.q_average[1]
    end
    return q, riv_storage, ssf_storage
end

function initial_head(x)
    return 2 * √x
end

"""
    transient_aquifer_1d(x, time, conductivity, specific_yield, aquifer_length, beta)

Non-steady flow in an unconfined rectangular aquifer, with Dirichlet h(0, t) = 0
on the left edge, and a Neumann Boundary Condition (dh/dx = 0) on the right.
"""
function transient_aquifer_1d(x, time, conductivity, specific_yield, aquifer_length, beta)
    return initial_head(x) / 1.0 +
           (beta * conductivity * initial_head(aquifer_length) * time) /
           (specific_yield * aquifer_length * aquifer_length)
end

function homogenous_aquifer(nrow, ncol)
    shape = (nrow, ncol)
    # Domain, geometry
    domain = ones(Bool, shape)
    dx = fill(10.0, ncol)
    dy = fill(10.0, nrow)
    indices, reverse_indices = Wflow.active_indices(domain, false)
    connectivity = Wflow.Connectivity(indices, reverse_indices, dx, dy)
    ncell = connectivity.ncell

    variables = Wflow.ConstantHeadVariables(; head = Float64[])
    constanthead = Wflow.ConstantHead(; variables, index = Int64[])

    timestepping = Wflow.TimeStepping()

    parameters = Wflow.GroundwaterFlowParameters(;
        hydraulic_conductivity = fill(10.0 / 86400.0, ncell),
        top = fill(10.0, ncell),
        bottom = fill(0.0, ncell),
        area = fill(100.0, ncell),
        specific_yield = fill(0.15, ncell),
        hydraulic_conductivity_scale_parameter = fill(3.0, ncell),
    )
    variables = Wflow.GroundwaterFlowVariables(;
        n_cells = ncell,
        head = [0.0, 7.5, 20.0],
        conductance = fill(0.0, connectivity.nconnection),
        storage = fill(0.0, ncell),
        q_net = fill(0.0, ncell),
        exfiltwater_cumulative = fill(0.0, ncell),
    )

    gwf_model = Wflow.GroundwaterFlowModel(;
        timestepping,
        parameters,
        variables,
        connectivity,
        constanthead,
    )
    return gwf_model
end

function init_sbm_soil_model(n_cells, N; kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)
    kwargs[:n_cells] = n_cells

    if !haskey(kwargs, :kv_profile)
        kwargs[:kv_profile] = nothing
    end

    if !haskey(kwargs, :vegetation_parameter_set)
        kwargs[:vegetation_parameter_set] = Wflow.VegetationParameters(;
            rooting_depth = get(kwargs, :rooting_depth, []),
            leaf_area_index = nothing,
            storage_wood = nothing,
            light_extinction_coefficient = nothing,
            storage_specific_leaf = nothing,
            canopy_gap_fraction = [],
            maximum_canopy_storage = [],
            crop_coefficient = [],
        )
    end

    if !haskey(kwargs, :maximum_number_of_layers)
        kwargs[:maximum_number_of_layers] = 0
    end

    # Vectors of SVectors
    for field_name in [
        :unsaturated_layer_depth,
        :unsaturated_layer_thickness,
        :volumetric_water_content,
        :relative_volumetric_water_content,
        :actual_layer_thickness,
        :rootfraction,
        :vertical_hydraulic_conductivity_factor,
        :brooks_corey_exponent,
        :cumulative_layer_depth,
    ]
        if !haskey(kwargs, field_name)
            kwargs[field_name] = SVector{N, Float64}[]
        end
    end

    # Vectors of other types
    for field_name in [
        # Variables
        :unsaturated_store_capacity,
        :saturated_water_depth,
        :drainable_water_depth,
        :water_table_depth,
        :n_unsatlayers,
        :total_soil_water_storage,
        # Parameters
        :number_of_layers,
        :theta_s,
        :theta_r,
        :theta_fc,
        :soil_water_capacity,
        :air_entry_pressure,
        :soil_thickness,
        :infiltration_capacity_compacted_soil,
        :infiltration_capacity_soil,
        :maximum_leakage,
        :cap_hmax,
        :cap_n,
        :w_soil,
        :cf_soil,
        :compacted_soil_area_fraction,
        :wet_root_distribution_parameter,
        :h1,
        :h2,
        :h3_high,
        :h3_low,
        :h4,
        :alpha_h1,
        :soil_fraction,
    ]
        if !haskey(kwargs, field_name)
            kwargs[field_name] = []
        end
    end

    kwargs_variables =
        filter(pair -> pair.first ∈ fieldnames(Wflow.SbmSoilVariables), kwargs)
    variables = Wflow.SbmSoilVariables(; kwargs_variables...)

    kwargs_parameters =
        filter(pair -> pair.first ∈ fieldnames(Wflow.SbmSoilParameters), kwargs)
    parameters = Wflow.SbmSoilParameters(; kwargs_parameters...)

    kwargs_bc = filter(pair -> pair.first ∈ fieldnames(Wflow.SbmSoilBC), kwargs)
    boundary_conditions = Wflow.SbmSoilBC(; kwargs_bc...)

    return Wflow.SbmSoilModel(; n_cells, variables, parameters, boundary_conditions)
end

"""
River Flow Model without any restrictions on the fields, so that only the
data required in certain functions has to be supplied (e.g. in the form of NamedTuple).
"""
@kwdef struct DummyRiver{A, B, V} <: Wflow.AbstractRiverFlowModel
    allocation::A = nothing
    boundary_conditions::B = nothing
    variables::V = nothing
end

no_nan(x::Float64) = isnan(x) ? 0.0 : x
get_mean(f::Vector{Float64}) = mean(filter(!isnan, f)) |> no_nan
get_mean(f::Vector{SVector{N, Float64}}) where {N} = no_nan.(
    SVector{N}([mean(filter(!isnan, [v[i] for v in f])) for i in 1:length(first(f))]),
)

function get_means(obj)
    d = Dict{Symbol, Union{Float64, SVector{N, Float64} where N}}()
    for s in propertynames(obj)
        f = getfield(obj, s)
        if f isa Union{Vector{Float64}, Vector{SVector{N, Float64}} where {N}}
            d[s] = get_mean(f)
        end
    end
    return d
end

function test_means(obj::Any, means::Dict{Symbol})
    failed = Symbol[]
    for (s, v) in means
        v_obj = get_mean(getfield(obj, s))
        if !(v_obj ≈ v)
            push!(failed, s)
            err = v - v_obj
            fac = v ./ v_obj
            println("-"^50)
            println("$s: err = (v_expected - v_actual) = ($v - $v_obj) = $err")
            @show fac
        end
    end
    return isempty(failed)
end
