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
    rootingdepth = ncread(
        dataset,
        config,
        "vertical.vegetation_parameter_set.rootingdepth";
        sel = indices,
        defaults = 750.0,
        type = Float,
    )
    kc = ncread(
        dataset,
        config,
        "vertical.vegetation_parameter_set.kc";
        sel = indices,
        defaults = 1.0,
        type = Float,
    )
    if haskey(config.input.vertical.vegetation_parameter_set, "leaf_area_index")
        storage_specific_leaf = ncread(
            dataset,
            config,
            "vertical.vegetation_parameter_set.storage_specific_leaf";
            optional = false,
            sel = indices,
            type = Float,
        )
        storage_wood = ncread(
            dataset,
            config,
            "vertical.vegetation_parameter_set.storage_wood";
            optional = false,
            sel = indices,
            type = Float,
        )
        kext = ncread(
            dataset,
            config,
            "vertical.vegetation_parameter_set.kext";
            optional = false,
            sel = indices,
            type = Float,
        )
        vegetation_parameter_set = VegetationParameters(;
            leaf_area_index = fill(mv, n),
            storage_wood,
            kext,
            storage_specific_leaf,
            canopygapfraction = fill(mv, n),
            cmax = fill(mv, n),
            rootingdepth,
            kc,
        )
    else
        canopygapfraction = ncread(
            dataset,
            config,
            "vertical.vegetation_parameter_set.canopygapfraction";
            sel = indices,
            defaults = 0.1,
            type = Float,
        )
        cmax = ncread(
            dataset,
            config,
            "vertical.vegetation_parameter_set.cmax";
            sel = indices,
            defaults = 1.0,
            type = Float,
        )
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

@get_units @grid_loc @with_kw struct LandParameters{T}
    # cell area [m^2]
    area::Vector{T} | "m^2"
    # drain width [m]
    width::Vector{T} | "m"
    # drain slope
    slope::Vector{T} | "-"
end

function LandParameters(nc, config, inds)
    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc; outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)
    area = xl .* yl
    ldd = ncread(nc, config, "ldd"; optional = false, sel = inds, allow_missing = true)
    drain_width = map(get_flow_width, ldd, xl, yl)
    landslope = ncread(
        nc,
        config,
        "vertical.land_parameter_set.slope";
        optional = false,
        sel = inds,
        type = Float,
    )
    clamp!(landslope, 0.00001, Inf)

    land_parameter_set =
        LandParameters{Float}(; area = area, width = drain_width, slope = landslope)
    return land_parameter_set
end

@get_units @grid_loc @with_kw struct RiverParameters{T}
    # river width [m]
    width::Vector{T} | "m"
    # river length
    length::Vector{T} | "m"
    # slope
    slope::Vector{T} | "-"
end

function RiverParameters(nc, config, inds)
    riverwidth = ncread(
        nc,
        config,
        "lateral.river_parameter_set.width";
        optional = false,
        sel = inds,
        type = Float,
    )
    riverlength = ncread(
        nc,
        config,
        "lateral.river_parameter_set.length";
        optional = false,
        sel = inds,
        type = Float,
    )
    riverslope = ncread(
        nc,
        config,
        "lateral.river_parameter_set.slope";
        optional = false,
        sel = inds,
        type = Float,
    )
    minimum(riverlength) > 0 || error("river length must be positive on river cells")
    minimum(riverwidth) > 0 || error("river width must be positive on river cells")
    clamp!(riverslope, 0.00001, Inf)

    river_parameter_set = RiverParameters{Float}(;
        width = riverwidth,
        length = riverlength,
        slope = riverslope,
    )
    return river_parameter_set
end