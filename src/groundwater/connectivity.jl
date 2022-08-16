"""
    Connectivity{T}

Stores connection data between cells. Connections are stored in a compressed
sparse column (CSC) adjacency matrix: only non-zero values are stored.
Primarily, this consist of two vectors:
* the row value vector holds the cell indices of neighbors
* the column pointers marks the start and end index of the row value vector

This matrix is square: n = m = ncell. nconnection is equal to nnz (the number of
non-zero values).

* ncell: the number of (active) cells in the simulation
* nconnection: the number of connections between cells
* length1: for every connection, the length in the first cell, size nconnection
* length2: for every connection, the length in the second cell, size nconnection
* width: width for every connection, (approx.) perpendicular to length1 and
  length2, size nconnection
* colptr: CSC column pointer (size ncell + 1)
* rowval: CSC row value (size nconnection)
"""
struct Connectivity{T}
    ncell::Int
    nconnection::Int
    length1::Vector{T}
    length2::Vector{T}
    width::Vector{T}
    colptr::Vector{Int}
    rowval::Vector{Int}
end


"""
    connections(C::Connectivity, id::Int)

Returns connections for a single cell, identified by ``id``.
"""
connections(C::Connectivity, id::Int) = C.colptr[id]:(C.colptr[id+1]-1)


"""
    connection_geometry(I, J, Δx, Δy)

Compute geometrical properties of connections for structured input.
"""
function connection_geometry(I, J, Δx, Δy)
    if I[1] != J[1]  # connection in y
        length1 = 0.5 * Δy[I[1]]
        length2 = 0.5 * Δy[J[1]]
        width = Δx[I[2]]
    elseif I[2] != J[2]  # connection in x
        length1 = 0.5 * Δx[I[2]]
        length2 = 0.5 * Δx[J[2]]
        width = Δy[I[1]]
    else
        # TODO: more specific exception? --> Martijn
        error("Inconsistent index")
    end
    return (length1, length2, width)
end


# Define cartesian indices for neighbors
const neighbors = (
    CartesianIndex(0, -1),
    CartesianIndex(-1, 0),
    CartesianIndex(1, 0),
    CartesianIndex(0, 1),
)

# Constructor for the Connectivity structure for structured input
function Connectivity(indices, reverse_indices, Δx::Vector{T}, Δy::Vector{T}) where {T}
    # indices: These map from the 1D internal domain to the 2D external domain.
    # reverse_indices: from the 2D external domain to the 1D internal domain,
    # providing an Int which can be used as a linear index
    nrow, ncol = size(reverse_indices)
    # Pre-allocate output, allocate for full potential number of neighbors (4)
    ncell = length(indices)
    colptr = Vector{Int}(undef, ncell + 1)
    rowval = Vector{Int}(undef, ncell * 4)
    length1 = similar(rowval, T)
    length2 = similar(rowval, T)
    width = similar(rowval, T)

    i = 1  # column index of sparse matrix
    j = 1  # row index of sparse matrix
    for I in indices # loop over active indices
        colptr[j] = i
        # Strictly increasing numbering for any row
        # (Required by a CSCSparseMatrix, if you want to convert)
        for neighbor in neighbors
            J = I + neighbor
            if (1 <= J[1] <= nrow) && (1 <= J[2] <= ncol && reverse_indices[J] != 0) # Check if it's inbounds and neighbor is active
                rowval[i] = reverse_indices[J]
                length1[i], length2[i], width[i] = connection_geometry(I, J, Δx, Δy)
                i += 1
            end
        end
        j += 1
    end
    colptr[j] = i

    nconnection = i - 1
    return Connectivity(
        ncell,
        nconnection,
        length1[1:nconnection],
        length2[1:nconnection],
        width[1:nconnection],
        colptr,
        rowval[1:nconnection],
    )
end
