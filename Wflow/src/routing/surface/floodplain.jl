abstract type AbstractFloodPlainModel end
abstract type AbstractFloodPlainParameters end
abstract type AbstractFloodPlainVariables end

"""
    FloodPlainProfile

Floodplain `storage` is a function of `depth` (flood depth intervals). Based on the
cumulative floodplain `storage` a floodplain profile as a function of `flood_depth` is
derived with floodplain area `a` (cumulative) and wetted perimeter radius `p` (cumulative).
"""
@with_kw struct FloodPlainProfile
    depth::Vector{Float64}     # Flood depth [m]
    storage::Matrix{Float64}   # Flood storage (cumulative) [m³]
    width::Matrix{Float64}     # Flood width [m]
    a::Matrix{Float64}         # Flow area (cumulative) [m²]
    p::Matrix{Float64}         # Wetted perimeter (cumulative) [m]
end

"Initialize floodplain profile `FloodPlainProfile`"
function FloodPlainProfile(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver;
    index_pit::Vector{Int} = Int[],
)
    (; indices) = domain.network
    (; flow_width, flow_length) = domain.parameters
    storage = ncread(
        dataset,
        config,
        "floodplain_water__sum_of_volume_per_depth",
        Routing;
        sel = indices,
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
    profile = FloodPlainProfile(; storage, width, depth = flood_depths, a, p)
    return profile
end

"Struct to store floodplain flow model parameters on a staggered grid"
@with_kw struct FloodPlainStaggeredParameters <: AbstractFloodPlainParameters
    profile::FloodPlainProfile                  # floodplain profile
    # edge parameters
    mannings_n::Vector{Float64} = Float64[]     # manning's roughness at edge [s m-1/3]
    mannings_n_sq::Vector{Float64} = Float64[]  # manning's roughness squared at edge [(s m-1/3)2]
    zb_max::Vector{Float64} = Float64[]         # maximum bankfull elevation at edge [m]
    slope::Vector{Float64} = Float64[]          # slope at edge [-]
end

"Initialize floodplain flow model parameters on a staggered grid"
function FloodPlainStaggeredParameters(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    zb_floodplain::Vector{Float64},
    index_pit::Vector{Int},
)
    (; indices, nodes_at_edge, graph) = domain.network
    (; flow_length) = domain.parameters
    (; river_routing) = config.model
    n_edges = ne(graph)
    profile = FloodPlainProfile(dataset, config, domain; index_pit)
    mannings_n = ncread(
        dataset,
        config,
        "floodplain_water_flow__manning_n_parameter",
        Routing;
        sel = indices,
    )
    # manning roughness at edges
    append!(mannings_n, mannings_n[index_pit]) # copy to ghost nodes
    mannings_n_at_edge =
        compute_mannings_n_at_edge(mannings_n, flow_length, nodes_at_edge, n_edges)
    zb_max_edge = compute_value_at_edge(zb_floodplain, nodes_at_edge, n_edges, maximum)

    if river_routing == RoutingType.local_inertial
        mannings_n_sq_at_edge = mannings_n_at_edge .* mannings_n_at_edge
        slope_at_edge = []
    elseif river_routing == RoutingType.manning_staggered
        mannings_n_sq_at_edge = []
        flow_length_at_edge =
            compute_value_at_edge(flow_length, nodes_at_edge, n_edges, mean)
        slope_at_edge = compute_slope_at_edge(
            zb_floodplain,
            flow_length_at_edge,
            nodes_at_edge,
            n_edges,
        )
    end

    parameters = FloodPlainStaggeredParameters(;
        profile,
        mannings_n = mannings_n_at_edge,
        mannings_n_sq = mannings_n_sq_at_edge,
        zb_max = zb_max_edge,
        slope = slope_at_edge,
    )
    return parameters
end

"Struct to store floodplain flow model variables on a staggered grid"
@with_kw struct FloodPlainStaggeredVariables <: AbstractFloodPlainVariables
    n::Int
    n_edges::Int
    # edge variables
    a::Vector{Float64} = zeros(n_edges)    # flow area at edge [m²]
    r::Vector{Float64} = zeros(n_edges)    # hydraulic radius at edge [m]
    hf::Vector{Float64} = zeros(n_edges)   # water depth at edge [m]
    q0::Vector{Float64} = zeros(n_edges)   # discharge at edge at previous time step
    q::Vector{Float64} = zeros(n_edges)    # discharge at edge  [m³ s⁻¹]
    q_av::Vector{Float64} = zeros(n_edges) # average river discharge at edge  [m³ s⁻¹] for model timestep Δt
    hf_index::Vector{Int} = zeros(Int, n_edges) # edge index with `hf` [-] above depth threshold
    # node variables
    storage::Vector{Float64} = zeros(n)    # storage [m³]
    h::Vector{Float64}                     # water depth [m]
    error::Vector{Float64} = zeros(n)      # error storage [m³]
end

"Floodplain flow model"
@with_kw struct FloodPlainModel{
    T <: AbstractRoutingMethod,
    P <: AbstractFloodPlainParameters,
    V <: AbstractFloodPlainVariables,
} <: AbstractFloodPlainModel
    routing_method::T
    parameters::P
    variables::V
end

"Struct to store floodplain parameters"
@with_kw struct FloodPlainParameters <: AbstractFloodPlainParameters
    profile::FloodPlainProfile              # floodplain profile
    mannings_n::Vector{Float64} = Float64[] # manning's roughness[s m-1/3]
end

"Struct to store floodplain variables"
@with_kw struct FloodPlainVariables <: AbstractFloodPlainVariables
    n::Int
    q::Vector{Float64} = zeros(n)               # discharge [m³ s⁻¹]
    q_av::Vector{Float64} = zeros(n)            # average floodplain discharge [m³ s⁻¹] for model timestep Δt
    qin::Vector{Float64} = zeros(n)             # floodplain inflow from upstream cells [m³ s⁻¹]
    qin_av::Vector{Float64} = zeros(n)          # Average floodplain inflow from upstream cells [m³ s⁻¹] for model timestep Δt
    storage::Vector{Float64} = zeros(n)         # storage [m³]
    h::Vector{Float64} = zeros(n)               # water depth [m]
    flow_capacity::Vector{Float64} = zeros(n)   # flow capacity [m³ dt⁻¹]
end

"Determine the initial floodplain storage"
function initialize_storage!(river, domain::Domain, nriv::Int)
    (; flow_width, flow_length) = domain.river.parameters
    (; floodplain) = river
    (; profile) = floodplain.parameters
    for i in 1:nriv
        i1, i2 = interpolation_indices(floodplain.variables.h[i], profile.depth)
        a = compute_floodplain_flow_area(profile, floodplain.variables.h[i], i, i1, i2)
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
Compute flood flow area (including area above channel) based on flow depth `h` and
floodplain `depth`, `area` and `width` of a floodplain profile.
"""
function compute_flood_flow_area(
    profile::FloodPlainProfile,
    h::Float64,
    idx::Int,
    i1::Int,
    i2::Int,
)
    (; a, width, depth) = profile
    dh = h - depth[i1]  # depth at i1
    area = a[i1, idx] + (width[i2, idx] * dh) # area at i1, width at i2
    return area
end

"""
Compute floodplain wetted perimeter based on flow depth `h` and floodplain `depth` and
wetted perimeter `p` of a floodplain profile.
"""
function compute_wetted_perimeter(profile::FloodPlainProfile, h::Float64, idx::Int, i1::Int)
    (; p, depth) = profile
    dh = h - depth[i1] # depth at i1
    p = p[i1, idx] + 2.0 * dh # p at i1
    return p
end

"Compute flood depth by interpolating flood storage `flood_storage` using flood depth intervals."
function compute_flood_depth(
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

function compute_floodplain_flow_area(
    profile::FloodPlainProfile,
    h::Float64,
    idx::Int,
    i1::Int,
    i2::Int,
)
    channel_area = profile.width[1, idx] * h
    floodplain_flow_area = compute_flood_flow_area(profile, h, idx, i1, i2)
    floodplain_flow_area = max(floodplain_flow_area - channel_area, 0.0)
    return floodplain_flow_area
end

function active_floodplain_cells(river_flow_model)
    (; floodplain) = river_flow_model
    (; hf) = river_flow_model.variables
    (; active_e, h_thresh) = river_flow_model.parameters

    n = 0
    @inbounds for i in active_e
        @inbounds if hf[i] > h_thresh
            n += 1
            floodplain.variables.hf_index[n] = i
        else
            floodplain.variables.q[i] = 0.0
        end
    end
    return n
end

"Initialize floodplain geometry, model variables and parameters on staggered grid"
function FloodPlainModel(
    dataset::NCDataset,
    config::Config,
    domain::DomainRiver,
    zb_floodplain::Vector{Float64},
)
    (; indices, local_drain_direction, graph) = domain.network
    (; river_routing) = config.model
    n = length(indices)
    index_pit = findall(x -> x == 5, local_drain_direction)
    n_edges = ne(graph)
    parameters =
        FloodPlainStaggeredParameters(dataset, config, domain, zb_floodplain, index_pit)
    if river_routing == RoutingType.local_inertial
        routing_method = LocalInertial()
    elseif river_routing == RoutingType.manning_staggered
        routing_method = ManningStaggered()
    end
    variables = FloodPlainStaggeredVariables(; n, n_edges, h = zeros(n + length(index_pit)))

    floodplain_model = FloodPlainModel(; routing_method, parameters, variables)
    return floodplain_model
end

function FloodPlainModel(dataset::NCDataset, config::Config, domain::DomainRiver)
    (; indices) = domain.network
    n = length(indices)
    profile = FloodPlainProfile(dataset, config, domain)
    mannings_n = ncread(
        dataset,
        config,
        "floodplain_water_flow__manning_n_parameter",
        Routing;
        sel = indices,
    )
    parameters = FloodPlainParameters(; profile, mannings_n)
    variables = FloodPlainVariables(; n)
    floodplain_model =
        FloodPlainModel(; routing_method = KinematicWave(), parameters, variables)
    return floodplain_model
end
