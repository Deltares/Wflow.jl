"Struct to store (shared) vegetation parameters"
@get_units @grid_loc @with_kw struct VegetationParameters{T}
    # Leaf area index [m² m⁻²]
    leaf_area_index::Union{Vector{T}, Nothing} | "m2 m-2"
    # Storage woody part of vegetation [mm]
    storage_wood::Union{Vector{T}, Nothing} | "mm"
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    kext::Union{Vector{T}, Nothing} | "-"
    # Specific leaf storage [mm]
    storage_specific_leaf::Union{Vector{T}, Nothing} | "mm"
    # Canopy gap fraction [-]
    canopygapfraction::Vector{T} | "-"
    # Maximum canopy storage [mm] 
    cmax::Vector{T} | "mm"
    # Rooting depth [mm]
    rootingdepth::Vector{T} | "mm"
    # Crop coefficient Kc [-]
    kc::Vector{T} | "-"
end

"Initialize (shared) vegetation parameters"
function VegetationParameters(dataset, config, indices)
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
@get_units @grid_loc @with_kw struct LandGeometry{T}
    # cell area [m^2]
    area::Vector{T} | "m^2"
    # drain width [m]
    width::Vector{T} | "m"
    # drain slope
    slope::Vector{T} | "-"
end

"Initialize land geometry parameters"
function LandGeometry(nc, config, inds)
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

    land_parameter_set =
        LandGeometry{Float64}(; area = area, width = drain_width, slope = landslope)
    return land_parameter_set
end

"Struct to store river geometry parameters"
@get_units @grid_loc @with_kw struct RiverGeometry{T}
    # river width [m]
    width::Vector{T} | "m"
    # river length
    length::Vector{T} | "m"
    # slope
    slope::Vector{T} | "-"
end

"Initialize river geometry parameters"
function RiverGeometry(nc, config, inds)
    lens = lens_input_parameter(config, "river__width"; optional = false)
    riverwidth = ncread(nc, config, lens; sel = inds, type = Float64)
    lens = lens_input_parameter(config, "river__length"; optional = false)
    riverlength = ncread(nc, config, lens; sel = inds, type = Float64)
    lens = lens_input_parameter(config, "river__slope"; optional = false)
    riverslope = ncread(nc, config, lens; sel = inds, type = Float64)

    minimum(riverlength) > 0 || error("river length must be positive on river cells")
    minimum(riverwidth) > 0 || error("river width must be positive on river cells")
    clamp!(riverslope, 0.00001, Inf)

    river_parameter_set = RiverGeometry{Float64}(;
        width = riverwidth,
        length = riverlength,
        slope = riverslope,
    )
    return river_parameter_set
end