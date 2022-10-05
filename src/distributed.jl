function broadcast_to_ranks(obj, comm::MPI.Comm)
    return MPI.bcast(obj, root, comm)
end

broadcast_to_ranks(obj, comm::Nothing) = obj

function scatter_to_ranks(
    data,
    indices,
    active_ranks::Vector{Int},
    comm::MPI.Comm,
    rank::Int,
)

    reqs = Vector{MPI.Request}(undef, 0)

    if rank == root
        for (i, sendrank) in enumerate(active_ranks)
            # Get the indices on the root buffer to send to the remote buffer
            inds = indices[i]
            data_on_root = data[inds]
            root_buf = MPI.Buffer(data_on_root)
            sendtag = sendrank + 1000
            sreq = MPI.Isend(root_buf, sendrank, sendtag, comm)
            push!(reqs, sreq)
        end
    end
    if rank in active_ranks
        i = findall(x -> x == rank, active_ranks)
        remote_size = length(indices[i])
        A_local = zeros(eltype(data), remote_size)
        remote_buf = MPI.Buffer(A_local)
        recievetag = rank + 1000
        rreq = MPI.Irecv!(remote_buf, root, recievetag, comm)
        push!(reqs, rreq)
        MPI.Waitall!(reqs)
        return A_local
    else
        return []
    end
end

function scatter_to_ranks(
    data,
    indices,
    nprocs::Int,
    comm::MPI.Comm,
    rank::Int;
    default_size = 1,
)

    reqs = Vector{MPI.Request}(undef, 0)

    size = length(indices[rank+1])
    remote_size = size == 0 ? default_size : length(indices[rank+1])
    A_local = zeros(eltype(data), remote_size)
    remote_buf = MPI.Buffer(A_local)

    if rank == root
        for sendrank = 0:nprocs-1
            # Get the indices on the root buffer to send to the remote buffer
            inds = indices[sendrank+1]
            data_on_root = data[inds]
            root_buf = MPI.Buffer(data_on_root)
            sendtag = sendrank + 1000
            sreq = MPI.Isend(root_buf, sendrank, sendtag, comm)
            push!(reqs, sreq)
        end
    end

    recievetag = rank + 1000
    rreq = MPI.Irecv!(remote_buf, root, recievetag, comm)
    push!(reqs, rreq)

    MPI.Waitall!(reqs)
    return A_local
end

scatter_to_ranks(A, indices, nprocs::Int, comm::Nothing, rank::Int) = A
