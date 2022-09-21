function split_count(N::Integer, n::Integer)
    size_section, extra = divrem(N, n)
    sizes = vcat(repeat([size_section+1], extra), repeat([size_section], n-extra))
    return sizes
end

function get_start_end_indices(sizes)
    indices = Vector{NTuple{2, Int}}(undef, 0)
    index_s = similar(sizes)
    index_s[1] = 1
    push!(indices, (index_s[1], sizes[1]))
    for i=2:length(sizes)
        index_s[i] = index_s[i-1] + sizes[i-1]
        push!(indices,(index_s[i], index_s[i] + sizes[i]-1 ))
    end

    return indices
end

function broadcast_to_ranks(obj, comm::MPI.Comm)
    return MPI.bcast(obj, root, comm)
end

broadcast_to_ranks(obj, comm::Nothing) = obj

function scatter_to_ranks(data, sizes, indices, nprocs, comm::MPI.Comm, rank)

    
    remote_size = sizes[rank + 1]
    A_local = zeros(eltype(data), remote_size)
    remote_buf = MPI.Buffer(A_local)
    
    reqs = Vector{MPI.Request}(undef, 0)

    if rank == root
        for sendrank in 0:nprocs-1

            # Get the indices on the root buffer to send to the remote buffer
            ilo, ihi = indices[sendrank + 1]
            data_on_root = @view data[ilo:ihi]
            root_buf = MPI.Buffer(data_on_root)
            sendtag = sendrank + 1000
            sreq =  MPI.Isend(root_buf, sendrank, sendtag, comm)
            push!(reqs, sreq)
        end
    end

    recievetag = rank + 1000
    rreq = MPI.Irecv!(remote_buf, root, recievetag, comm)
    push!(reqs, rreq)
    
    MPI.Waitall!(reqs)
    return A_local
end

scatter_to_ranks(A::Vector{T}, sizes, indices, nprocs, comm::Nothing, rank) where {T} = A