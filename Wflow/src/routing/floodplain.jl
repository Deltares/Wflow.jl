abstract type AbstractFloodPlain end

"""
    FloodPlainProfile

Floodplain `storage` is a function of `depth` (flood depth intervals). Based on the
cumulative floodplain `storage` a floodplain profile as a function of `flood_depth` is
derived with floodplain area `a` (cumulative) and wetted perimeter radius `p` (cumulative).
"""
@with_kw struct FloodPlainProfile{N}
    depth::Vector{Float64}        # Flood depth [m]
    storage::Array{Float64, 2}    # Flood storage (cumulative) [m³]
    width::Array{Float64, 2}      # Flood width [m]
    a::Array{Float64, 2}          # Flow area (cumulative) [m²]
    p::Array{Float64, 2}          # Wetted perimeter (cumulative) [m]
end

"Initialize floodplain profile `FloodPlainProfile`"
function FloodPlainProfile(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver;
    index_pit::Vector{Int} = [],
)
    (; indices) = domain.network
    (; flow_width, flow_length) = domain.parameters
    storage = ncread(
        dataset,
        config,
        "floodplain_water__sum_of_volume_per_depth";
        optional = false,
        sel = indices,
        type = Float64,
        dimname = :flood_depth,
    )
    n = length(indices)

    # for convenience (interpolation) flood depth 0.0 m is added, with associated area (a),
    # storage, width (river width) and wetted perimeter (p).
    storage = vcat(fill(Float64(0), n)', storage)
    start_storage = storage
    flood_depths = Float64.(dataset["flood_depth"][:])
    pushfirst!(flood_depths, 0.0)
    n_depths = length(flood_depths)

    p = zeros(n_depths, n)
    a = zeros(n_depths, n)
    segment_storage = zeros(n_depths, n)
    width = zeros(n_depths, n)
    width[1, :] = flow_width[1:n]

    # determine flow area (a), width and wetted perimeter (p) FloodPlain
    h = diff(flood_depths)
    incorrect_vol = 0
    riv_cells = 0
    error_vol = 0
    for i in 1:n
        riv_cell = 0
        diff_storage = diff(storage[:, i])

        for j in 1:(n_depths - 1)
            # assume rectangular shape of flood depth segment
            width[j + 1, i] = diff_storage[j] / (h[j] * flow_length[i])
            # check provided flood storage (floodplain width should be constant or increasing
            # as a function of flood depth)
            if width[j + 1, i] < width[j, i]
                # raise warning only if difference is larger than rounding error of 0.01 m³
                if ((width[j, i] - width[j + 1, i]) * h[j] * flow_length[i]) > 0.01
                    incorrect_vol += 1
                    riv_cell = 1
                    error_vol =
                        error_vol +
                        ((width[j, i] - width[j + 1, i]) * h[j] * flow_length[i])
                end
                width[j + 1, i] = width[j, i]
            end
            a[j + 1, i] = width[j + 1, i] * h[j]
            p[j + 1, i] = (width[j + 1, i] - width[j, i]) + 2.0 * h[j]
            segment_storage[j + 1, i] = a[j + 1, i] * flow_length[i]
            if j == 1
                # for interpolation wetted perimeter at flood depth 0.0 is required
                p[j, i] = p[j + 1, i] - 2.0 * h[j]
            end
        end

        p[2:end, i] = cumsum(p[2:end, i])
        a[:, i] = cumsum(a[:, i])
        storage[:, i] = cumsum(segment_storage[:, i])

        riv_cells += riv_cell
    end

    if incorrect_vol > 0
        perc_riv_cells = round(100.0 * (riv_cells / n); digits = 2)
        perc_error_vol = round(100.0 * (error_vol / sum(start_storage[end, :])); digits = 2)
        @warn string(
            "The provided storage of $incorrect_vol rectangular floodplain schematization",
            " segments for $riv_cells river cells ($perc_riv_cells % of total river cells)",
            " is not correct and has been increased with $perc_error_vol % of provided storage.",
        )
    end

    # set floodplain parameters for ghost points
    storage = hcat(storage, storage[:, index_pit])
    width = hcat(width, width[:, index_pit])
    a = hcat(a, a[:, index_pit])
    p = hcat(p, p[:, index_pit])

    # initialize floodplain profile parameters
    profile = FloodPlainProfile{n_depths}(; storage, width, depth = flood_depths, a, p)
    return profile
end

"Struct to store local inertial floodplain flow model parameters"
@with_kw struct LocalInertialFloodPlainParameters{P}
    profile::P                          # floodplain profile
    mannings_n::Vector{Float64}         # manning's roughness [s m-1/3]
    mannings_n_sq::Vector{Float64}      # manning's roughness squared at edge [(s m-1/3)2]
    zb_max::Vector{Float64}             # maximum bankfull elevation at edge [m]
end

"Struct to store kinematic wave floodplain flow model parameters"
@with_kw struct KinWaveFloodPlainParameters{P}
    profile::P                          # floodplain profile
    mannings_n::Vector{Float64}         # manning's roughness [s m-1/3]
end

function get_floodplain_mannings_n(dataset, config, indices)
    mannings_n = ncread(
        dataset,
        config,
        "floodplain_water_flow__manning_n_parameter";
        sel = indices,
        defaults = 0.072,
        type = Float64,
    )
    return mannings_n
end

"Initialize local inertial floodplain flow model parameters"
function LocalInertialFloodPlainParameters(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    zb_floodplain::Vector{Float64},
    index_pit::Vector{Int},
)
    (; indices, nodes_at_edge, graph) = domain.network
    (; flow_length) = domain.parameters
    n_edges = ne(graph)
    profile = FloodPlainProfile(dataset, config, domain; index_pit)
    mannings_n = get_floodplain_mannings_n(dataset, config, indices)
    # manning roughness at edges
    append!(mannings_n, mannings_n[index_pit]) # copy to ghost nodes
    mannings_n_sq = fill(Float64(0), n_edges)
    zb_max = fill(Float64(0), n_edges)
    for i in 1:n_edges
        src_node = nodes_at_edge.src[i]
        dst_node = nodes_at_edge.dst[i]
        mannings_n_i =
            (
                mannings_n[dst_node] * flow_length[dst_node] +
                mannings_n[src_node] * flow_length[src_node]
            ) / (flow_length[dst_node] + flow_length[src_node])
        mannings_n_sq[i] = mannings_n_i * mannings_n_i
        zb_max[i] = max(zb_floodplain[src_node], zb_floodplain[dst_node])
    end
    parameters =
        LocalInertialFloodPlainParameters(profile, mannings_n, mannings_n_sq, zb_max)
    return parameters
end

"Struct to store local inertial floodplain flow model variables"
@with_kw struct LocalInertialFloodPlainVariables
    n::Int
    n_edges::Int
    storage::Vector{Float64} = zeros(n)    # storage [m³]
    h::Vector{Float64}                     # water depth [m]
    error::Vector{Float64} = zeros(n)      # error storage [m³]
    a::Vector{Float64} = zeros(n_edges)    # flow area at edge [m²]
    r::Vector{Float64} = zeros(n_edges)    # hydraulic radius at edge [m]
    hf::Vector{Float64} = zeros(n_edges)   # water depth at edge [m]
    q0::Vector{Float64} = zeros(n_edges)   # discharge at edge at previous time step
    q::Vector{Float64} = zeros(n_edges)    # discharge at edge  [m³ s⁻¹]
    q_av::Vector{Float64} = zeros(n_edges) # average river discharge at edge  [m³ s⁻¹] for model timestep Δt
    hf_index::Vector{Int} = zeros(Int, n_edges) # edge index with `hf` [-] above depth threshold
end

"Struct to store kinematic wave floodplain flow model variables"
@with_kw struct KinWaveFloodPlainVariables
    n::Int
    storage::Vector{Float64} = zeros(n) # storage [m³]
    h::Vector{Float64} = zeros(n)       # water depth [m]
    q::Vector{Float64} = zeros(n)       # discharge [m³ s⁻¹]
    q_av::Vector{Float64} = zeros(n)    # average river discharge [m³ s⁻¹] for model timestep Δt
end

"Local inertial floodplain flow model"
@with_kw struct LocalInertialFloodPlain{P} <: AbstractFloodPlain
    parameters::LocalInertialFloodPlainParameters{P}
    variables::LocalInertialFloodPlainVariables
end

"Kinematic wave floodplain flow model"
@with_kw struct KinWaveFloodPlain{P} <: AbstractFloodPlain
    parameters::KinWaveFloodPlainParameters{P}
    variables::KinWaveFloodPlainVariables
end

"Determine the initial floodplain storage"
function initialize_storage!(river, domain::Domain, nriv::Int)
    (; flow_width, flow_length) = domain.river.parameters
    (; floodplain) = river
    (; profile) = floodplain.parameters
    for i in 1:nriv
        i1, i2 = interpolation_indices(floodplain.variables.h[i], profile.depth)
        a = flow_area(
            profile.width[i2, i],
            profile.a[i1, i],
            profile.depth[i1],
            floodplain.variables.h[i],
        )
        a = max(a - (flow_width[i] * floodplain.variables.h[i]), 0.0)
        floodplain.variables.storage[i] = flow_length[i] * a
    end
    return nothing
end

"helper function to get interpolation indices"
function interpolation_indices(x, v::AbstractVector)
    i1 = 1
    for i in eachindex(v)
        if v[i] <= x
            i1 = i
        end
    end
    if i1 == length(v)
        i2 = i1
    else
        i2 = i1 + 1
    end
    return i1, i2
end

"""
    flow_area(width, area, depth, h)

Compute floodplain flow area based on flow depth `h` and floodplain `depth`, `area` and
`width` of a floodplain profile.
"""
function flow_area(width, area, depth, h)
    dh = h - depth  # depth at i1
    area = area + (width * dh) # area at i1, width at i2
    return area
end

"""
    function wetted_perimeter(p, depth, h)

Compute floodplain wetted perimeter based on flow depth `h` and floodplain `depth` and
wetted perimeter `p` of a floodplain profile.
"""
function wetted_perimeter(p, depth, h)
    dh = h - depth # depth at i1
    p += 2.0 * dh # p at i1
    return p
end

"Compute flood depth by interpolating flood storage `flood_storage` using flood depth intervals."
function flood_depth(
    profile::FloodPlainProfile,
    flood_storage::Float64,
    flow_length::Float64,
    i::Int,
)
    i1, i2 = interpolation_indices(flood_storage, @view profile.storage[:, i])
    ΔA = (flood_storage - profile.storage[i1, i]) / flow_length
    dh = ΔA / profile.width[i2, i]
    flood_depth = profile.depth[i1] + dh
    return flood_depth
end

"Initialize floodplain geometry and `LocalInerialFloodPlain` variables and parameters"
function LocalInertialFloodPlain(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    zb_floodplain::Vector{Float64},
)
    (; indices, local_drain_direction, graph) = domain.network
    n = length(indices)
    index_pit = findall(x -> x == 5, local_drain_direction)
    n_edges = ne(graph)
    parameters =
        LocalInertialFloodPlainParameters(dataset, config, domain, zb_floodplain, index_pit)
    variables =
        LocalInertialFloodPlainVariables(; n, n_edges, h = zeros(n + length(index_pit)))

    floodplain = LocalInertialFloodPlain(; parameters, variables)
    return floodplain
end

function KinWaveFloodPlain(dataset::NCDataset, config::Config, domain::DomainRiver)
    (; indices) = domain.network
    profile = FloodPlainProfile(dataset, config, domain)
    mannings_n = get_floodplain_mannings_n(dataset, config, indices)
    profile = FloodPlainProfile(dataset, config, domain)
    parameters = KinWaveFloodPlainParameters(profile, mannings_n)
    variables = KinWaveFloodPlainVariables()
    floodplain = KinWaveFloodPlain(; parameters, variables)
    return floodplain
end
