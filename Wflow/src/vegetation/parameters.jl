"Struct to store (shared) vegetation parameters"
@with_kw struct VegetationParameters
    # Leaf area index [m² m⁻²]
    leaf_area_index::Union{Vector{Float64}, Nothing}
    # Storage woody part of vegetation [m]
    storage_wood::Union{Vector{Float64}, Nothing} = nothing
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    light_extinction_coefficient::Union{Vector{Float64}, Nothing} = nothing
    # Specific leaf storage [m]
    storage_specific_leaf::Union{Vector{Float64}, Nothing} = nothing
    # Canopy gap fraction [-]
    canopy_gap_fraction::Vector{Float64}
    # Maximum canopy storage [m]
    maximum_canopy_storage::Vector{Float64}
    # Rooting depth [m]
    rooting_depth::Vector{Float64} = []
    # Crop coefficient Kc [-]
    crop_coefficient::Vector{Float64} = []
    # Canopy height [m]
    canopy_height::Vector{Float64} = []
end

"Initialize (shared) vegetation parameters"
function VegetationParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    n = length(indices)
    rooting_depth =
        ncread(dataset, config, "vegetation_root__depth", LandHydrologySBM; sel = indices)
    crop_coefficient =
        ncread(dataset, config, "vegetation__crop_factor", LandHydrologySBM; sel = indices)
    canopy_height = ncread(
        dataset,
        config,
        "vegetation_canopy__height",
        LandHydrologySBM;
        sel = indices,
    )
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
        light_extinction_coefficient = ncread(
            dataset,
            config,
            "vegetation_canopy__light_extinction_coefficient",
            LandHydrologySBM;
            sel = indices,
        )
        vegetation_parameter_set = VegetationParameters(;
            leaf_area_index = fill(MISSING_VALUE, n),
            storage_wood,
            light_extinction_coefficient,
            storage_specific_leaf,
            canopy_gap_fraction = fill(MISSING_VALUE, n),
            maximum_canopy_storage = fill(MISSING_VALUE, n),
            rooting_depth,
            crop_coefficient,
            canopy_height,
        )
    else
        canopy_gap_fraction = ncread(
            dataset,
            config,
            "vegetation_canopy__gap_fraction",
            LandHydrologySBM;
            sel = indices,
        )
        maximum_canopy_storage = ncread(
            dataset,
            config,
            "vegetation_water__storage_capacity",
            LandHydrologySBM;
            sel = indices,
        )
        vegetation_parameter_set = VegetationParameters(;
            leaf_area_index = nothing,
            storage_wood = nothing,
            light_extinction_coefficient = nothing,
            storage_specific_leaf = nothing,
            canopy_gap_fraction,
            maximum_canopy_storage,
            rooting_depth,
            crop_coefficient,
            canopy_height,
        )
    end
    return vegetation_parameter_set
end
