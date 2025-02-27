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

"Struct to store land geometry parameters"
@with_kw struct LandGeometry
    # cell area [m²]
    area::Vector{Float64}
    # drain width [m]
    width::Vector{Float64}
    # drain slope [-]
    slope::Vector{Float64}
end

"Initialize land geometry parameters"
function LandGeometry(nc::NCDataset, config::Config, inds::Vector{CartesianIndex{2}})
    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc; outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)
    area = xl .* yl
    lens = lens_input(config, "local_drain_direction"; optional = false)
    ldd = ncread(nc, config, lens; sel = inds, allow_missing = true)
    drain_width = map(get_flow_width, ldd, xl, yl)
    lens = lens_input_parameter(config, "land_surface__slope"; optional = false)
    landslope = ncread(nc, config, lens; sel = inds, type = Float64)
    clamp!(landslope, 0.00001, Inf)

    land_parameter_set = LandGeometry(; area = area, width = drain_width, slope = landslope)
    return land_parameter_set
end

"Struct to store river geometry parameters"
@with_kw struct RiverGeometry
    # river width [m]
    width::Vector{Float64}
    # river length [m]
    length::Vector{Float64}
    # slope [-]
    slope::Vector{Float64}
end

"Initialize river geometry parameters"
function RiverGeometry(nc::NCDataset, config::Config, inds::Vector{CartesianIndex{2}})
    lens = lens_input_parameter(config, "river__width"; optional = false)
    riverwidth = ncread(nc, config, lens; sel = inds, type = Float64)
    lens = lens_input_parameter(config, "river__length"; optional = false)
    riverlength = ncread(nc, config, lens; sel = inds, type = Float64)
    lens = lens_input_parameter(config, "river__slope"; optional = false)
    riverslope = ncread(nc, config, lens; sel = inds, type = Float64)

    minimum(riverlength) > 0 || error("river length must be positive on river cells")
    minimum(riverwidth) > 0 || error("river width must be positive on river cells")
    clamp!(riverslope, 0.00001, Inf)

    river_parameter_set =
        RiverGeometry(; width = riverwidth, length = riverlength, slope = riverslope)
    return river_parameter_set
end