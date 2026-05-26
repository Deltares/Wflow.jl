"Struct to store (shared) vegetation parameters"
@with_data_lookup struct VegetationParameters
    # Leaf area index [m² m⁻²]
    "vegetation__leaf_area_index"
    leaf_area_index::Union{Vector{Float64}, Nothing}
    # Storage woody part of vegetation [mm]
    "vegetation_wood_water__storage_capacity"
    storage_wood::Union{Vector{Float64}, Nothing}
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    "vegetation_canopy__light_extinction_coefficient"
    kext::Union{Vector{Float64}, Nothing}
    # Specific leaf storage [mm]
    "vegetation__specific_leaf_storage"
    storage_specific_leaf::Union{Vector{Float64}, Nothing}
    # Canopy gap fraction [-]
    "vegetation_canopy__gap_fraction"
    canopygapfraction::Vector{Float64}
    # Maximum canopy storage [mm]
    "vegetation_water__storage_capacity"
    cmax::Vector{Float64}
    # Rooting depth [mm]
    "vegetation_root__depth"
    rootingdepth::Vector{Float64}
    # Crop coefficient Kc [-]
    "vegetation__crop_factor"
    kc::Vector{Float64}
end

"Initialize (shared) vegetation parameters"
function VegetationParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
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
        vegetation_parameter_set = VegetationParameters(
            data_lookup;
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
        vegetation_parameter_set = VegetationParameters(
            data_lookup;
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
