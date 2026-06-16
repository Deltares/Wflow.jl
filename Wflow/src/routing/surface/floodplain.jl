abstract type AbstractFloodPlainModel end
abstract type AbstractFloodPlainParameters end
abstract type AbstractFloodPlainVariables end

""""
    FloodPlainModel

Floodplain flow model for each river cell (as part of a `RiverFlowModel`). The
`FloodPlainModel` can be part of three different river routing schemes:
* Local inertial equation on a staggered grid.
* Kinematic wave using Manning's equation on a staggered grid.
* Kinematic wave approach by solving the kinematic wave equation using Newton's method.

On a staggered grid the river and floodplain flow routing scheme are equal. For the
kinematic wave river flow routing solved by using Newton's method, the floodplain flow
routing is solved by using Manning's equation, as solving kinematic wave routing for near
zero flows is computationally expensive.
"""
@with_kw struct FloodPlainModel{
    T <: AbstractRoutingMethod,
    P <: AbstractFloodPlainParameters,
    V <: AbstractFloodPlainVariables,
} <: AbstractFloodPlainModel
    routing_method::T
    parameters::P
    variables::V
end

"""
    FloodPlainProfile

Floodplain `storage` is a function of `depth` (flood depth intervals). Based on the
cumulative floodplain `storage` a floodplain profile as a function of `flood_depth` is
derived with floodplain area `a` (cumulative) and wetted perimeter radius `p` (cumulative).
"""
@with_kw struct FloodPlainProfile
    # Flood depth [m]
    depth::Vector{Float64}
    # Flood storage (cumulative) [m³]
    storage::Matrix{Float64}
    # Flood width [m]
    width::Matrix{Float64}
    # Flow area (cumulative) [m²]
    flow_area::Matrix{Float64}
    # Wetted perimeter (cumulative) [m]
    wetted_perimeter::Matrix{Float64}
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

    wetted_perimeter = zeros(n_depths, n)
    flow_area = zeros(n_depths, n)
    segment_storage = zeros(n_depths, n)
    width = zeros(n_depths, n)
    width[1, :] = flow_width[1:n]

    # determine flow area (a), width and wetted perimeter (p) FloodPlainModel
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
            flow_area[j + 1, i] = width[j + 1, i] * h[j]
            wetted_perimeter[j + 1, i] = (width[j + 1, i] - width[j, i]) + 2.0 * h[j]
            segment_storage[j + 1, i] = flow_area[j + 1, i] * flow_length[i]
            if j == 1
                # for interpolation wetted perimeter at flood depth 0.0 is required
                wetted_perimeter[j, i] = wetted_perimeter[j + 1, i] - 2.0 * h[j]
            end
        end

        wetted_perimeter[2:end, i] = cumsum(wetted_perimeter[2:end, i])
        flow_area[:, i] = cumsum(flow_area[:, i])
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
    flow_area = hcat(flow_area, flow_area[:, index_pit])
    wetted_perimeter = hcat(wetted_perimeter, wetted_perimeter[:, index_pit])

    # initialize floodplain profile parameters
    profile = FloodPlainProfile(;
        storage,
        width,
        depth = flood_depths,
        flow_area,
        wetted_perimeter,
    )
    return profile
end

"Struct to store floodplain flow model parameters on a staggered grid"
@with_kw struct FloodPlainStaggeredParameters <: AbstractFloodPlainParameters
    # floodplain profile
    profile::FloodPlainProfile
    # manning's roughness [s m-1/3]
    mannings_n::Vector{Float64} = Float64[]
    # manning's roughness at edge [s m-1/3]
    mannings_n_at_edge::Vector{Float64} = Float64[]
    # manning's roughness squared at edge [(s m-1/3)2]
    mannings_n_sq_at_edge::Vector{Float64} = Float64[]
    # maximum bankfull elevation at edge [m]
    zb_max_at_edge::Vector{Float64} = Float64[]
    # slope at edge [-]
    slope_at_edge::Vector{Float64} = Float64[]
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
    zb_max_at_edge = compute_value_at_edge(zb_floodplain, nodes_at_edge, n_edges, maximum)

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
        mannings_n,
        mannings_n_at_edge,
        mannings_n_sq_at_edge,
        zb_max_at_edge,
        slope_at_edge,
    )
    return parameters
end

"Struct to store floodplain flow model variables on a staggered grid"
@with_kw struct FloodPlainStaggeredVariables <: AbstractFloodPlainVariables
    n::Int
    n_edges::Int
    # water depth at edge [m]
    water_depth_at_edge::Vector{Float64} = zeros(n_edges)
    # discharge at edge at previous time step
    q_previous::Vector{Float64} = zeros(n_edges)
    # discharge at edge  [m³ s⁻¹]
    q::Vector{Float64} = zeros(n_edges)
    # cumulative discharge at edge [m³] for model timestep dt
    q_cumulative::Vector{Float64} = zeros(n_edges)
    # average discharge at edge [m³ s⁻¹] for model timestep dt
    q_average::Vector{Float64} = zeros(n_edges)
    # edge index with `water_depth_at_edge` [-] above depth threshold
    hf_index::Vector{Int} = zeros(Int, n_edges)
    # storage [m³]
    storage::Vector{Float64} = zeros(n)
    # water depth [m]
    h::Vector{Float64}
    # error storage [m³]
    error::Vector{Float64} = zeros(n)
end

"Struct to store floodplain parameters"
@with_kw struct FloodPlainParameters <: AbstractFloodPlainParameters
    # floodplain profile
    profile::FloodPlainProfile
    # manning's roughness[s m-1/3]
    mannings_n::Vector{Float64} = Float64[]
    # slope  [m m⁻¹]
    slope::Vector{Float64} = Float64[]
end

"Struct to store floodplain variables"
@with_kw struct FloodPlainVariables <: AbstractFloodPlainVariables
    n::Int
    # discharge [m³ s⁻¹]
    q::Vector{Float64} = zeros(n)
    # cumulative discharge for model timestep dt
    q_cumulative::Vector{Float64} = zeros(n)
    # average discharge [m³ s⁻¹] for model timestep Δt
    q_average::Vector{Float64} = zeros(n)
    # inflow from upstream cells [m³ s⁻¹]
    qin::Vector{Float64} = zeros(n)
    # cumulative inflow from upstream cells for model timestep dt
    qin_cumulative::Vector{Float64} = zeros(n)
    # average inflow from upstream cells [m³ s⁻¹] for model timestep Δt
    qin_average::Vector{Float64} = zeros(n)
    # storage [m³]
    storage::Vector{Float64} = zeros(n)
    # water depth [m]
    h::Vector{Float64} = zeros(n)
    # flow capacity [m³ s⁻¹]
    flow_capacity::Vector{Float64} = zeros(n)
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
floodplain `depth`, `flow_area` and `width` of a floodplain profile.
"""
function compute_flood_flow_area(
    profile::FloodPlainProfile,
    h::Float64,
    idx::Int,
    i1::Int,
    i2::Int,
)
    dh = h - profile.depth[i1]  # depth at i1
    flow_area = profile.flow_area[i1, idx] + (profile.width[i2, idx] * dh) # area at i1, width at i2
    return flow_area
end

"""
Compute floodplain wetted perimeter based on flow depth `h` and floodplain `depth` and
`wetted_perimeter` of a floodplain profile.
"""
function compute_wetted_perimeter(profile::FloodPlainProfile, h::Float64, idx::Int, i1::Int)
    dh = h - profile.depth[i1] # depth at i1
    wetted_perimeter = profile.wetted_perimeter[i1, idx] + 2.0 * dh # p at i1
    return wetted_perimeter
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

"Compute floodplain flow area (excluding river channel area)"
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
    (; water_depth_at_edge) = river_flow_model.variables
    (; active_e, h_thresh) = river_flow_model.parameters

    n = 0
    @inbounds for i in active_e
        @inbounds if water_depth_at_edge[i] > h_thresh
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
    index_pit = findall(x -> x == LDD_PIT, local_drain_direction)
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

"""
Initialize floodplain geometry, model variables and parameters for floodplain flow routing
as part of kinematic wave river flow routing (solved using Newton's method).
"""
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
    slope = ncread(dataset, config, "floodplain__slope", Routing; sel = indices)
    clamp!(slope, 0.00001, Inf)
    parameters = FloodPlainParameters(; profile, mannings_n, slope)
    variables = FloodPlainVariables(; n)
    floodplain_model = FloodPlainModel(; routing_method = Manning(), parameters, variables)
    return floodplain_model
end
