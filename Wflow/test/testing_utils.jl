using Statistics: mean
using SpecialFunctions: expint

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
        q[i] = model.routing.river_flow.variables.q_av[1]
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

"""
    drawdown_theis(distance, time, discharge, transmissivity, storativity)

Non-steady flow in a confined aquifer, using the well function of Theis.
"""
function drawdown_theis(distance, time, discharge, transmissivity, storativity)
    u = (storativity * distance^2) / (4 * transmissivity * time)
    return discharge / (4 * pi * transmissivity) * expint(u)
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

    parameters = Wflow.ConfinedAquiferParameters(;
        k = fill(10.0, ncell),
        top = fill(10.0, ncell),
        bottom = fill(0.0, ncell),
        area = fill(100.0, ncell),
        specific_storage = fill(0.1, ncell),
        storativity = fill(1.0, ncell),
    )
    variables = Wflow.AquiferVariables(;
        n = ncell,
        head = [0.0, 7.5, 20.0],
        conductance = fill(0.0, connectivity.nconnection),
        storage = fill(0.0, ncell),
        q_net = fill(0.0, ncell),
        q_in_av = fill(0.0, ncell),
        q_out_av = fill(0.0, ncell),
        exfiltwater = fill(0.0, ncell),
    )
    conf_aqf = Wflow.ConfinedAquifer(; parameters, variables)

    parameters = Wflow.UnconfinedAquiferParameters(;
        k = fill(10.0, ncell),
        top = fill(10.0, ncell),
        bottom = fill(0.0, ncell),
        area = fill(100.0, ncell),
        specific_yield = fill(0.15, ncell),
        f = fill(3.0, ncell),
    )
    variables = Wflow.AquiferVariables(;
        n = ncell,
        head = [0.0, 7.5, 20.0],
        conductance = fill(0.0, connectivity.nconnection),
        storage = fill(0.0, ncell),
        q_net = fill(0.0, ncell),
        q_in_av = fill(0.0, ncell),
        q_out_av = fill(0.0, ncell),
        exfiltwater = fill(0.0, ncell),
    )
    unconf_aqf = Wflow.UnconfinedAquifer(; parameters, variables)
    return (connectivity, conf_aqf, unconf_aqf)
end

function simple_soil(n, N; kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)
    kwargs[:n] = n

    if !haskey(kwargs, :kv_profile)
        kwargs[:kv_profile] = nothing
    end

    # Vectors of SVectors
    for field_name in [:vwc, :vwc_perc, :act_thickl, :rootfraction, :kvfrac]
        kwargs[field_name] = SVector{N, Float64}[]
    end

    # Vectors of other types
    for field_name in [
        # Variables
        :ustorecapacity,
        :satwaterdepth,
        :zi,
        :n_unsatlayers,
        :total_soilwater_storage,
        # Parameters
        :nlayers,
        :theta_s,
        :theta_r,
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

    return Wflow.SbmSoilModel(; n, variables, parameters)
end

@kwdef struct DummyRiver{A, B, V} <: Wflow.AbstractRiverFlowModel
    allocation::A = nothing
    boundary_conditions::B = nothing
    variables::V = nothing
end
