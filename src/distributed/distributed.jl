function partition_graph(g::SimpleGraph, nparts::Integer)

    G = Metis.graph(g)
    G.xadj .= G.xadj .- Cint(1)
    G.adjncy .= G.adjncy .- Cint(1)

    part = Vector{Cint}(undef, G.nvtxs)
    vwgt = isdefined(G, :vwgt) ? G.vwgt : C_NULL
    edgecut = fill(Cint(0), 1)

    options = fill(Cint(-1), Metis.METIS_NOPTIONS)
    options[Int(Metis.METIS_OPTION_CONTIG)+1] = Cint(1)
    options[Int(Metis.METIS_OPTION_OBJTYPE)+1] = Cint(Metis.METIS_OBJTYPE_CUT)

    Metis.METIS_PartGraphKway(
        Ref{Metis.idx_t}(G.nvtxs),
        Ref{Metis.idx_t}(1),
        G.xadj,
        G.adjncy,
        vwgt,
        C_NULL,
        C_NULL,
        Ref{Metis.idx_t}(nparts),
        C_NULL,
        C_NULL,
        options,
        edgecut,
        part,
    )

    return part .+ 1
end

function set_pits!(graph, ldd, subdomains, nprocs)
    for i = 1:nprocs
        inds = findall(x -> x == i, subdomains)
        g, vmap = induced_subgraph(graph, inds)
        toposort = topological_sort_by_dfs(g)
        pit_index = vmap[toposort[end]]
        ldd[pit_index] = 5
    end
end

function write_static_netcdf(ds, data, var_in, var_out, size, inds)

    buffer = zeros(Union{eltype(data),Missing}, size)
    fill!(buffer, missing)
    buffer[inds] .= data

    _, data_dim_order = read_dims(ds[var_in], (x = :, y = :))
    if first(data_dim_order) !== :x
        buffer = permutedims(buffer)
    end
    dims_increasing = dim_directions(ds, (:x, :y))
    reverse_data!(buffer, dims_increasing)
    if var_out in keys(ds)
        ds[var_out] .= buffer
    else
        defVar(ds, var_out, buffer, dimnames(ds[var_in]))
    end

end

function initialize_sbm_models(config::Config)

    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    if rank === 0
        static_path = input_path(config, config.input.path_static)

        nc = NCDataset(static_path, "a")
        subcatch_2d =
            ncread(nc, config, "subcatchment"; optional = false, allow_missing = true)
        inds, _ = active_indices(subcatch_2d, missing)
        ldd_2d = ncread(nc, config, "ldd"; optional = false, allow_missing = true)
        ldd = ldd_2d[inds]
        modelsize_2d = size(subcatch_2d)

        graph = flowgraph(ldd, inds, pcr_dir)
        subdomains = partition_graph(SimpleGraph(graph), nprocs - 1)

        var_in = config.input.subcatchment
        var_out = "subdomains"
        write_static_netcdf(nc, subdomains, var_in, var_out, modelsize_2d, inds)
        set_pits!(graph, ldd, subdomains, nprocs - 1)
        var_in = config.input.ldd
        var_out = "ldd_subdomains"
        write_static_netcdf(nc, ldd, var_in, var_out, modelsize_2d, inds)
        close(nc)
    end
    MPI.Barrier(comm)

    if rank !== 0
        config.input.subdomain = rank
        config.input.subdomains = "subdomains"
        config.input.ldd_subdomains = "ldd_subdomains"
        if haskey(config, "output") && haskey(config.output, "path")
            path = config.output.path
            config.output.path =
                string(splitext(path)[1], string("_rank_$(rank)"), splitext(path)[2])
        end
        if haskey(config, "state") && haskey(config.state, "path_output")
            path = config.state.path_output
            config.state.path_output =
                string(splitext(path)[1], string("_rank_$(rank)"), splitext(path)[2])
        end

        model = initialize_sbm_model(config)
    end
end
