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
    rootingdepth =
        ncread(dataset, config, "vegetation_root__depth", LandHydrologySBM; sel = indices)
    kc = ncread(dataset, config, "vegetation__crop_factor", LandHydrologySBM; sel = indices)
    if do_cyclic(config) && haskey(config.input.cyclic, "vegetation__leaf_area_index")
        storage_specific_leaf = ncread(
            dataset,
            config,
            "vegetation__specific_leaf_storage",
            LandHydrologySBM;
            sel = indices,
        )
        storage_wood = ncread(
            dataset,
            config,
            "vegetation_wood_water__storage_capacity",
            LandHydrologySBM;
            sel = indices,
        )
        kext = ncread(
            dataset,
            config,
            "vegetation_canopy__light_extinction_coefficient",
            LandHydrologySBM;
            sel = indices,
        )
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
        canopygapfraction = ncread(
            dataset,
            config,
            "vegetation_canopy__gap_fraction",
            LandHydrologySBM;
            sel = indices,
        )
        cmax = ncread(
            dataset,
            config,
            "vegetation_water__storage_capacity",
            LandHydrologySBM;
            sel = indices,
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
