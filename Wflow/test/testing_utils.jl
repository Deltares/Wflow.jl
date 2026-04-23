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
    for timestep_idx in 1:steps
        Wflow.run_timestep!(model)
        ssf_storage[timestep_idx] = mean(model.routing.subsurface_flow.variables.storage)
        riv_storage[timestep_idx] = mean(model.routing.river_flow.variables.storage)
        q[timestep_idx] = model.routing.river_flow.variables.q_av[1]
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
    n_land_cells = connectivity.ncell

    variables = Wflow.ConstantHeadVariables(; head = Float64[])
    constanthead = Wflow.ConstantHead(; variables, index = Int64[])

    timestepping = Wflow.TimeStepping()

    parameters = Wflow.GroundwaterFlowParameters(;
        k = fill(10.0, n_land_cells),
        top = fill(10.0, n_land_cells),
        bottom = fill(0.0, n_land_cells),
        area = fill(100.0, n_land_cells),
        specific_yield = fill(0.15, n_land_cells),
        f = fill(3.0, n_land_cells),
    )
    variables = Wflow.GroundwaterFlowVariables(;
        n_land_cells,
        head = [0.0, 7.5, 20.0],
        conductance = fill(0.0, connectivity.nconnection),
        storage = fill(0.0, n_land_cells),
        q_net = fill(0.0, n_land_cells),
        q_in_av = fill(0.0, n_land_cells),
        q_av = fill(0.0, n_land_cells),
        exfiltwater = fill(0.0, n_land_cells),
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

function init_sbm_soil_model(n_land_cells, N; kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)
    kwargs[:n_land_cells] = n_land_cells

    if !haskey(kwargs, :kv_profile)
        kwargs[:kv_profile] = nothing
    end

    if !haskey(kwargs, :vegetation_parameter_set)
        kwargs[:vegetation_parameter_set] = Wflow.VegetationParameters(;
            rootingdepth = get(kwargs, :rootingdepth, []),
            leaf_area_index = nothing,
            storage_wood = nothing,
            kext = nothing,
            storage_specific_leaf = nothing,
            canopygapfraction = [],
            cmax = [],
            kc = [],
        )
    end

    if !haskey(kwargs, :maxlayers)
        kwargs[:maxlayers] = 0
    end

    # Vectors of SVectors
    for field_name in [
        :ustorelayerdepth,
        :ustorelayerthickness,
        :vwc,
        :vwc_perc,
        :act_thickl,
        :rootfraction,
        :kvfrac,
        :c,
        :sumlayers,
    ]
        if !haskey(kwargs, field_name)
            kwargs[field_name] = SVector{N, Float64}[]
        end
    end

    # Vectors of other types
    for field_name in [
        # Variables
        :ustorecapacity,
        :satwaterdepth,
        :drainable_waterdepth,
        :zi,
        :n_unsatlayers,
        :total_soilwater_storage,
        # Parameters
        :nlayers,
        :theta_s,
        :theta_r,
        :theta_fc,
        :soilwatercapacity,
        :hb,
        :soilthickness,
        :infiltcappath,
        :infiltcapsoil,
        :maxleakage,
        :cap_hmax,
        :cap_n,
        :w_soil,
        :cf_soil,
        :pathfrac,
        :rootdistpar,
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

    return Wflow.SbmSoilModel(; n_land_cells, variables, parameters, boundary_conditions)
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
