"Struct to store (shared) vegetation parameters"
@kwdef struct VegetationParameters
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
end