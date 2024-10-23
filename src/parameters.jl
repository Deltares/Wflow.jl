"Struct to store (shared) vegetation parameters"
@get_units @grid_loc @with_kw struct VegetationParameters{T}
    # Leaf area index [m² m⁻²]
    leaf_area_index::Union{Vector{T}, Nothing} | "m2 m-2"
    # Storage woody part of vegetation [mm]
    swood::Union{Vector{T}, Nothing} | "mm"
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    kext::Union{Vector{T}, Nothing} | "-"
    # Specific leaf storage [mm]
    sl::Union{Vector{T}, Nothing} | "mm"
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
function VegetationParameters(nc, config, inds)
    n = length(inds)
    rootingdepth = ncread(
        nc,
        config,
        "vertical.vegetation_parameter_set.rootingdepth";
        sel = inds,
        defaults = 750.0,
        type = Float,
    )
    kc = ncread(
        nc,
        config,
        "vertical.vegetation_parameter_set.kc";
        sel = inds,
        defaults = 1.0,
        type = Float,
    )
    if haskey(config.input.vertical.vegetation_parameter_set, "leaf_area_index")
        sl = ncread(
            nc,
            config,
            "vertical.vegetation_parameter_set.specific_leaf";
            optional = false,
            sel = inds,
            type = Float,
        )
        swood = ncread(
            nc,
            config,
            "vertical.vegetation_parameter_set.storage_wood";
            optional = false,
            sel = inds,
            type = Float,
        )
        kext = ncread(
            nc,
            config,
            "vertical.vegetation_parameter_set.kext";
            optional = false,
            sel = inds,
            type = Float,
        )
        vegetation_parameter_set = VegetationParameters(;
            leaf_area_index = fill(mv, n),
            swood = swood,
            kext = kext,
            sl = sl,
            canopygapfraction = fill(mv, n),
            cmax = fill(mv, n),
            rootingdepth = rootingdepth,
            kc = kc,
        )
    else
        canopygapfraction = ncread(
            nc,
            config,
            "vertical.vegetation_parameter_set.canopygapfraction";
            sel = inds,
            defaults = 0.1,
            type = Float,
        )
        cmax = ncread(
            nc,
            config,
            "vertical.vegetation_parameter_set.cmax";
            sel = inds,
            defaults = 1.0,
            type = Float,
        )
        vegetation_parameter_set = VegetationParameters(;
            leaf_area_index = nothing,
            swood = nothing,
            kext = nothing,
            sl = nothing,
            canopygapfraction = canopygapfraction,
            cmax = cmax,
            rootingdepth = rootingdepth,
            kc = kc,
        )
    end
    return vegetation_parameter_set
end