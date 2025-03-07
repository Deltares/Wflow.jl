"Struct to store (shared) vegetation parameters"
@with_kw struct VegetationParameters
    # Leaf area index [m² m⁻²]
    leaf_area_index::Union{Vector{Float64}, Nothing}
    # Storage woody part of vegetation [mm]
    storage_wood::Union{Vector{Float64}, Nothing}
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    kext::Union{Vector{Float64}, Nothing}
    # Specific leaf storage [mm]
    storage_specific_leaf::Union{Vector{Float64}, Nothing}
    # Canopy gap fraction [-]
    canopygapfraction::Vector{Float64}
    # Maximum canopy storage [mm] 
    cmax::Vector{Float64}
    # Rooting depth [mm]
    rootingdepth::Vector{Float64}
    # Crop coefficient Kc [-]
    kc::Vector{Float64}
end

"Initialize (shared) vegetation parameters"
function VegetationParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    lens = lens_input_parameter(config, "vegetation_root__depth")
    rootingdepth =
        ncread(dataset, config, lens; sel = indices, defaults = 750.0, type = Float64)
    lens = lens_input_parameter(config, "vegetation__crop_factor")
    kc = ncread(dataset, config, lens; sel = indices, defaults = 1.0, type = Float64)
    do_cyclic = haskey(config.input, "cyclic")
    if do_cyclic && haskey(config.input.cyclic, "vegetation__leaf-area_index")
        lens = lens_input_parameter(
            config,
            "vegetation__specific-leaf_storage";
            optional = false,
        )
        storage_specific_leaf = ncread(dataset, config, lens; sel = indices, type = Float64)
        lens = lens_input_parameter(
            config,
            "vegetation_wood_water__storage_capacity";
            optional = false,
        )
        storage_wood = ncread(dataset, config, lens; sel = indices, type = Float64)
        lens = lens_input_parameter(
            config,
            "vegetation_canopy__light-extinction_coefficient";
            optional = false,
        )
        kext = ncread(dataset, config, lens; sel = indices, type = Float64)
        vegetation_parameter_set = VegetationParameters(;
            leaf_area_index = fill(MISSING_VALUE, n),
            storage_wood,
            kext,
            storage_specific_leaf,
            canopygapfraction = fill(MISSING_VALUE, n),
            cmax = fill(MISSING_VALUE, n),
            rootingdepth,
            kc,
        )
    else
        lens = lens_input_parameter(config, "vegetation_canopy__gap_fraction")
        canopygapfraction =
            ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float64)
        lens = lens_input_parameter(config, "vegetation_water__storage_capacity")
        cmax = ncread(dataset, config, lens; sel = indices, defaults = 1.0, type = Float64)
        vegetation_parameter_set = VegetationParameters(;
            leaf_area_index = nothing,
            storage_wood = nothing,
            kext = nothing,
            storage_specific_leaf = nothing,
            canopygapfraction,
            cmax,
            rootingdepth,
            kc,
        )
    end
    return vegetation_parameter_set
end

"Struct to store (shared) land parameters"
@with_kw struct LandParameters
    # cell length x direction [m]
    x_length::Vector{Float64} = Float64[]
    # cell length y direction [m]
    y_length::Vector{Float64} = Float64[]
    # cell area [m²]
    area::Vector{Float64} = Float64[]
    # flow width [m]
    flow_width::Vector{Float64} = Float64[]
    # suface flow width [m]
    surface_flow_width::Vector{Float64} = Float64[]
    # flow length [m]
    flow_length::Vector{Float64} = Float64[]
    # flow fraction [-] to river
    fraction_to_river::Vector{Float64} = Float64[]
    # slope [-]
    slope::Vector{Float64} = Float64[]
    # water body (reservoir and lake) location
    waterbody::Vector{Bool} = Bool[]
    # river location
    river::Vector{Bool} = Bool[]
    # fraction of river [-]
    river_fraction = Float64[]
    # fraction of open water (excluding rivers) [-]
    water_fraction = Float64[]
end

function get_water_fraction(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand,
    river_fraction::Vector{Float64},
)
    lens = lens_input_parameter(config, "land~water-covered__area_fraction")
    water_fraction =
        ncread(dataset, config, lens; sel = network.indices, defaults = 0.0, type = Float64)
    water_fraction = max.(water_fraction .- river_fraction, 0.0)
    return water_fraction
end

function get_river_fraction(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand,
    river::Vector{Bool},
    area::Vector{Float64},
)
    lens = lens_input_parameter(config, "river__width"; optional = false)
    river_width_2d = ncread(dataset, config, lens; type = Float64, fill = 0)
    river_width = river_width_2d[network.indices]

    lens = lens_input_parameter(config, "river__length"; optional = false)
    river_length_2d = ncread(dataset, config, lens; type = Float64, fill = 0)
    river_length = river_length_2d[network.indices]

    n = length(river)
    river_fraction = fill(MISSING_VALUE, n)
    for i in 1:n
        river_fraction[i] = if river[i]
            min((river_length[i] * river_width[i]) / (area[i]), 1.0)
        else
            0.0
        end
    end
    return river_fraction
end

function get_cell_lengths(dataset::NCDataset, config::Config, network::NetworkLand)
    y_coords = read_y_axis(dataset)
    x_coords = read_x_axis(dataset)
    y = permutedims(repeat(y_coords; outer = (1, length(x_coords))))[network.indices]
    cellength = abs(mean(diff(x_coords)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    x_length, y_length = cell_lengths(y, cellength, sizeinmetres)
    return x_length, y_length
end

function get_landsurface_slope(dataset::NCDataset, config::Config, network::NetworkLand)
    lens = lens_input_parameter(config, "land_surface__slope"; optional = false)
    slope = ncread(dataset, config, lens; sel = network.indices, type = Float64)
    clamp!(slope, 0.00001, Inf)
    return slope
end

function river_mask(dataset::NCDataset, config::Config, network::NetworkLand)
    lens = lens_input(config, "river_location__mask"; optional = false)
    river_2d = ncread(dataset, config, lens; type = Bool, fill = false)
    river = river_2d[network.indices]
    return river
end

function waterbody_mask(dataset::NCDataset, config::Config, network::NetworkLand)
    # Get waterbodies mask
    do_reservoirs = get(config.model, "reservoir", false)::Bool
    do_lakes = get(config.model, "lake", false)::Bool
    waterbodies = fill(0, length(network.indices))
    if do_reservoirs
        lens = lens_input(config, "reservoir_area__count"; optional = false)
        reservoirs =
            ncread(dataset, config, lens; sel = network.indices, type = Float64, fill = 0)
        waterbodies = waterbodies .+ reservoirs
    end
    if do_lakes
        lens = lens_input(config, "lake_area__count"; optional = false)
        lakes =
            ncread(dataset, config, lens; sel = network.indices, type = Float64, fill = 0)
        waterbodies = waterbodies .+ lakes
    end
    waterbodies = Vector{Bool}(waterbodies .> 0)
    return waterbodies
end

"Initialize (shared) land parameters for model type `SedimentModel`"
function LandParameters(dataset::NCDataset, config::Config, network::NetworkLand)
    x_length, y_length = get_cell_lengths(dataset, config, network)
    area = x_length .* y_length
    flow_width = map(get_flow_width, network.local_drain_direction, x_length, y_length)
    slope = get_landsurface_slope(dataset, config, network)
    waterbody = waterbody_mask(dataset, config, network)
    river = river_mask(dataset, config, network)

    land_parameters = LandParameters(; area, flow_width, slope, waterbody, river)
    return land_parameters
end

"Initialize (shared) land parameters for model type `SbmModel` or `SbmGwfModel`"
function LandParameters(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand,
    river_indices::Vector{Int},
)
    x_length, y_length = get_cell_lengths(dataset, config, network)
    area = x_length .* y_length
    flow_width = map(get_flow_width, network.local_drain_direction, x_length, y_length)
    flow_length = map(get_flow_length, network.local_drain_direction, x_length, y_length)
    slope = get_landsurface_slope(dataset, config, network)
    river = river_mask(dataset, config, network)
    river_fraction = get_river_fraction(dataset, config, network, river, area)

    water_fraction = get_water_fraction(dataset, config, network, river_fraction)

    land_area = @. (1.0 - river_fraction) * area
    surface_flow_width = map(get_surface_width, flow_width, flow_length, land_area, river)

    fraction_to_river = flow_fraction_to_river(
        network.graph,
        network.local_drain_direction,
        river_indices,
        slope,
    )

    waterbody = waterbody_mask(dataset, config, network)

    land_parameters = LandParameters(;
        x_length,
        y_length,
        area,
        flow_width,
        surface_flow_width,
        flow_length,
        slope,
        river,
        fraction_to_river,
        waterbody,
        river_fraction,
        water_fraction,
    )
    return land_parameters
end

"Struct to store (shared) river parameters"
@with_kw struct RiverParameters
    # river flow width [m]
    flow_width::Vector{Float64} = Float64[]
    # river flow length [m]
    flow_length::Vector{Float64} = Float64[]
    # slope [-]
    slope::Vector{Float64} = Float64[]
    # water body (reservoir and lake) location
    waterbody::Vector{Bool} = Bool[]
    # grid cell area [m²]
    cell_area::Vector{Float64} = Float64[]
end

"Initialize (shared) river parameters"
function RiverParameters(dataset::NCDataset, config::Config, network::NetworkRiver)
    (; indices) = network
    lens = lens_input_parameter(config, "river__length"; optional = false)
    flow_length = ncread(dataset, config, lens; sel = indices, type = Float64)
    minimum(flow_length) > 0 || error("river length must be positive on river cells")

    lens = lens_input_parameter(config, "river__width"; optional = false)
    flow_width = ncread(dataset, config, lens; sel = indices, type = Float64)
    minimum(flow_width) > 0 || error("river width must be positive on river cells")

    lens = lens_input_parameter(config, "river__slope"; optional = false)
    slope = ncread(dataset, config, lens; sel = indices, type = Float64)
    clamp!(slope, 0.00001, Inf)

    river_parameters = RiverParameters(; flow_width, flow_length, slope)
    return river_parameters
end

@with_kw struct SharedParameters
    river::RiverParameters
    land::LandParameters
end

function SharedParameters(dataset::NCDataset, config::Config, network::Network)
    (; land_indices) = network.river
    land_params = LandParameters(dataset, config, network.land, land_indices)
    river_params = RiverParameters(dataset, config, network.river)

    @reset river_params.cell_area = land_params.area[network.river.land_indices]
    @reset river_params.waterbody = land_params.waterbody[network.river.land_indices]

    shared_params = SharedParameters(; river = river_params, land = land_params)
    return shared_params
end