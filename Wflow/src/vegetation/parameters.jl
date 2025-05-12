"Struct to store (shared) vegetation parameters"
@with_kw struct VegetationParameters
    # Leaf area index [m² m⁻²]
    leaf_area_index::Union{Vector{Float}, Nothing}
    # Storage woody part of vegetation [mm]
    storage_wood::Union{Vector{Float}, Nothing}
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    kext::Union{Vector{Float}, Nothing}
    # Specific leaf storage [mm]
    storage_specific_leaf::Union{Vector{Float}, Nothing}
    # Canopy gap fraction [-]
    canopygapfraction::Vector{Float}
    # Maximum canopy storage [mm] 
    cmax::Vector{Float}
    # Rooting depth [mm]
    rootingdepth::Vector{Float}
    # Crop coefficient Kc [-]
    kc::Vector{Float}
end

"Initialize (shared) vegetation parameters"
function VegetationParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    lens = lens_input_parameter(config, "vegetation_root__depth")
    rootingdepth = ncread(dataset, config, lens; sel = indices, type = Float)
    lens = lens_input_parameter(config, "vegetation__crop_factor")
    kc = ncread(dataset, config, lens; sel = indices, defaults = 1.0, type = Float)
    do_cyclic = haskey(config.input, "cyclic")
    if do_cyclic && haskey(config.input.cyclic, "vegetation__leaf-area_index")
        lens = lens_input_parameter(
            config,
            "vegetation__specific-leaf_storage";
            optional = false,
        )
        storage_specific_leaf = ncread(dataset, config, lens; sel = indices, type = Float)
        lens = lens_input_parameter(
            config,
            "vegetation_wood_water__storage_capacity";
            optional = false,
        )
        storage_wood = ncread(dataset, config, lens; sel = indices, type = Float)
        lens = lens_input_parameter(
            config,
            "vegetation_canopy__light-extinction_coefficient";
            optional = false,
        )
        kext = ncread(dataset, config, lens; sel = indices, type = Float)
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
        canopygapfraction = ncread(dataset, config, lens; sel = indices, type = Float)
        lens = lens_input_parameter(config, "vegetation_water__storage_capacity")
        cmax = ncread(dataset, config, lens; sel = indices, defaults = 1.0, type = Float)
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