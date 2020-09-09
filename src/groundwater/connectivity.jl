"""
Stores connection data between cells
"""
struct Connectivity{T,R,L}
    ncell::Int
    nconnection::Int
    # Map internal nodes to external id's
    external::Vector{CartesianIndex}
    # Map external nodes to internal id's
    internal::Array{Int,2}
    length1::Vector{T}
    length2::Vector{T}
    width::Vector{T}
    colptr::Vector{Int}
    rowval::Vector{Int}
end


"""
Returns connections for a single cell, identified by ``id``.
"""
connections(C::Connectivity, id::Int) = C._colptr[id]:(C._colptr[id + 1] - 1)


"""
Assign cell id's to the two dimensional array and 0 to inactive cells.
"""
function create_mapping(domain::Array{Bool,2})
    internal = similar(domain, Int)
    external = similar(domain, CartesianIndex)
    id = 0
    for I in CartesianIndices(domain)
        if domain[I]
            id += 1
            internal[I] = id
            external[id] = I
        else
            internal[I] = 0
        end
    end
    return internal, external[1:id]
end


"""
Compute geometrical properties of connections for structured input.
"""
function connection_geometry(I, J, Δx, Δy)
    if I[1] != J[1]  # connection in x
        length1 = 0.5 * Δx[J[1]]
        length2 = 0.5 * Δx[I[1]]
        width = Δy[I[2]]
    elseif I[2] != J[2]  # connection in y
        length1 = 0.5 * Δy[J[2]]
        length2 = 0.5 * Δy[I[2]]
        width = Δx[I[1]]
    else
        error("Inconsistent index")
    end
    return (length1, length2, width)
end


function Connectivity(domain::Array{Bool,2}, Δx::Vector{T}, Δy::Vector{T})
    nrow, ncol = size(domain)
    # Create internal <-> external id mappings
    internal, external = create_mapping(domain)

    # Pre-allocate output, allocate for full number of neighbors
    # We'll store only the part we need
    n_active = length(external)
    colptr = Vector{Int}(undef, n_active + 1)
    rowval = Vector{Int}(undef, n_active * 4)
    length1 = similar(rowval, T)
    length2 = similar(rowval, T)
    width = similar(rowval, T)

    # Define cartesian indices for neighbors
    left = CartesianIndex(-1, 0, 0)
    right = CartesianIndex(1, 0, 0)
    back = CartesianIndex(0, -1, 0)
    front = CartesianIndex(0, 1, 0)

    colptr[1] = 1
    i = 1  # column index of sparse matrix
    j = 1  # row index of sparse matrix
    for I in CartesianIndices(internal)
        row_i = i

        # Strictly increasing numbering for a row
        # (Required by a CSCSparseMatrix)
        for neighbor in (front, left, right, back)
            J = I + neighbor
            
            if (1 <= J[1] <= ncol) && (1 <= J[2] <= nrow) # Check if it's inbounds
                rowval[i] = internal[J]
                length1[i], length2[i], width[i] = connection_geometry(
                    I, J, Δx, Δy
                )
                i += 1
            end
        end

        # Store 
        j += 1
        colptr[j] = row_i
    end

    n_connection = i - 1
    
    return Connectivity(
            n_active,
            n_connection,
            external,
            internal,
            length1[1:n_connection],
            length2[1:n_connection],
            width[1:n_connection],
            colptr,
            rowval,
    )
end

        


            