abstract type AbstractSoilModel end
abstract type AbstractVerticalConductivityProfile end

"Struct for storing SBM soil model variables"
@kwdef struct SbmSoilVariables{N}
    n::Int
    maximum_number_of_layers::Int
    # Calculated soil water pressure head h3 of the root water uptake reduction function (Feddes) [m]
    h3::Vector{Float64} = fill(MISSING_VALUE, n)
    # Unsaturated store capacity [m]
    unsaturated_store_capacity::Vector{Float64} = fill(MISSING_VALUE, n)
    # Amount of water in the unsaturated store, per layer [m]
    unsaturated_layer_depth::Vector{SVector{N, Float64}} =
        fill(zeros(SVector{maximum_number_of_layers, Float64}), n)
    # Thickness of unsaturated zone, per layer [m]
    unsaturated_layer_thickness::Vector{SVector{N, Float64}} =
        fill(zeros(SVector{maximum_number_of_layers, Float64}), n)
    # Saturated store [m]
    saturated_water_depth::Vector{Float64} = fill(MISSING_VALUE, n)
    # Drainable water store [m]
    drainable_water_depth::Vector{Float64} = fill(MISSING_VALUE, n)
    # Pseudo-water table depth [m] (top of the saturated zone)
    water_table_depth::Vector{Float64} = fill(MISSING_VALUE, n)
    # Number of unsaturated soil layers
    n_unsatlayers::Vector{Int} = zeros(Int, n)
    # Transpiration [m s⁻¹]
    transpiration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual evaporation from unsaturated store [m s⁻¹]
    actual_evaporation_unsaturated_store::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil evaporation from unsaturated and saturated store [m s⁻¹]
    soil_evaporation::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil evaporation from saturated store [m s⁻¹]
    soil_evaporation_saturated_zone::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual capillary rise [m s⁻¹]
    actual_capillary_flux::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual transpiration from saturated store [m s⁻¹]
    actual_evaporation_saturated_zone::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total actual evapotranspiration [m s⁻¹]
    actual_evapotranspiration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration into the unsaturated zone [m s⁻¹]
    actual_infiltration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration non-compacted fraction [m s⁻¹]
    actual_infiltration_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration compacted fraction [m s⁻¹]
    actual_infiltration_compacted_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual infiltration (compacted and the non-compacted areas) [m s⁻¹]
    infiltration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Infiltration excess water [m s⁻¹]
    infiltration_excess::Vector{Float64} = fill(MISSING_VALUE, n)
    # Water that cannot infiltrate due to saturated soil (saturation excess) [m s⁻¹]
    saturation_excess_water::Vector{Float64} = fill(MISSING_VALUE, n)
    # Water exfiltrating during saturation excess conditions [m s⁻¹]
    exfiltration_saturated_water::Vector{Float64} = fill(MISSING_VALUE, n)
    # Excess water for non-compacted fraction [m s⁻¹]
    excess_water_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Excess water for compacted fraction [m s⁻¹]
    excess_water_compacted_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total surface runoff from infiltration and saturation excess (excluding actual open water evaporation) [m s⁻¹]
    runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Net surface runoff (surface runoff - actual open water evaporation) [m s⁻¹]
    net_runoff::Vector{Float64} = fill(MISSING_VALUE, n)
    # Volumetric water content [-] per soil layer (including theta_r and saturated zone)
    volumetric_water_content::Vector{SVector{N, Float64}} =
        fill(zeros(SVector{maximum_number_of_layers, Float64}), n)
    # Volumetric water content [%] per soil layer (including theta_r and saturated zone)
    relative_volumetric_water_content::Vector{SVector{N, Float64}} =
        fill(zeros(SVector{maximum_number_of_layers, Float64}), n)
    # Root water storage [m] in unsaturated and saturated zone (excluding theta_r)
    root_zone_storage::Vector{Float64} = fill(MISSING_VALUE, n)
    # Volumetric water content [-] in root zone (including theta_r and saturated zone)
    volumetric_water_content_root_zone::Vector{Float64} = fill(MISSING_VALUE, n)
    # Volumetric water content [%] in root zone (including theta_r and saturated zone)
    relative_volumetric_water_content_root_zone::Vector{Float64} = fill(MISSING_VALUE, n)
    # Amount of available water in the unsaturated zone [m]
    unsaturated_store_depth::Vector{Float64} = zeros(n)
    # Downward flux from unsaturated to saturated zone [m s⁻¹]
    transfer::Vector{Float64} = fill(MISSING_VALUE, n)
    # Net recharge to saturated store [m s⁻¹]
    recharge::Vector{Float64} = fill(MISSING_VALUE, n)
    # Actual leakage from saturated store [m s⁻¹]
    actual_leakage::Vector{Float64} = fill(MISSING_VALUE, n)
    # Total water storage (excluding floodplain volume and reservoirs) [m]
    total_storage::Vector{Float64} = zeros(Float64, n)
    # Total soil water storage [m]
    total_soil_water_storage::Vector{Float64} = fill(MISSING_VALUE, n)
    # Top soil temperature [K]
    soil_surface_temperature::Vector{Float64} = fill(283.15, n)
    # Soil infiltration reduction factor (when soil is frozen) [-]
    f_infiltration_reduction::Vector{Float64} = ones(n)
end

"Exponential depth profile of vertical hydraulic conductivity at the soil surface"
struct KvExponential <: AbstractVerticalConductivityProfile
    # Vertical hydraulic conductivity [m s⁻¹] at soil surface
    kv_0::Vector{Float64}
    # A scaling parameter [m⁻¹] (controls exponential decline of kv_0)
    hydraulic_conductivity_scale_parameter::Vector{Float64}
end

"Exponential constant depth profile of vertical hydraulic conductivity"
struct KvExponentialConstant <: AbstractVerticalConductivityProfile
    exponential::KvExponential
    # Depth [m] from soil surface for which exponential decline of kv_0 is valid
    z_exp::Vector{Float64}
end

"Layered depth profile of vertical hydraulic conductivity"
struct KvLayered{N} <: AbstractVerticalConductivityProfile
    # Vertical hydraulic conductivity [m s⁻¹] per soil layer
    kv::Vector{SVector{N, Float64}}
end

"Layered exponential depth profile of vertical hydraulic conductivity"
struct KvLayeredExponential{N} <: AbstractVerticalConductivityProfile
    # A scaling parameter [m⁻¹] (controls exponential decline of kv_0)
    hydraulic_conductivity_scale_parameter::Vector{Float64}
    # Vertical hydraulic conductivity [m s⁻¹] per soil layer
    kv::Vector{SVector{N, Float64}}
    # Number of soil layers [-] with vertical hydraulic conductivity value `kv`
    nlayers_kv::Vector{Int}
    # Depth [m] from soil surface for which layered profile is valid
    z_layered::Vector{Float64}
end

"Struct for storing SBM soil model parameters"
@kwdef struct SbmSoilParameters{
    N,
    M,
    Kv <: Union{AbstractVerticalConductivityProfile, Nothing},
}
    n::Int
    # Maximum number of soil layers [-]
    maximum_number_of_layers::Int
    # Number of soil layers [-]
    number_of_layers::Vector{Int} = zeros(Int, n)
    # Saturated water content (porosity) [-]
    theta_s::Vector{Float64} = fill(MISSING_VALUE, n)
    # Residual water content [-]
    theta_r::Vector{Float64} = fill(MISSING_VALUE, n)
    # Field capacity water content [-]
    theta_fc::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soilwater capacity [m]
    soil_water_capacity::Vector{Float64} = fill(MISSING_VALUE, n)
    # Multiplication factor [-] applied to kv_z (vertical flow)
    vertical_hydraulic_conductivity_factor::Vector{SVector{N, Float64}} =
        fill(ones(SVector{maximum_number_of_layers, Float64}), n)
    # Air entry pressure [m] of soil (Brooks-Corey)
    air_entry_pressure::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil thickness [m]
    soil_thickness::Vector{Float64} = fill(MISSING_VALUE, n)
    # Thickness of soil layers [m]
    actual_layer_thickness::Vector{SVector{N, Float64}} =
        fill(zeros(SVector{maximum_number_of_layers, Float64}), n)
    # Cumulative sum of soil layers [m], starting at soil surface (0)
    cumulative_layer_depth::Vector{SVector{M, Float64}} =
        fill(zeros(SVector{maximum_number_of_layers + 1, Float64}), n)
    # Infiltration capacity of the compacted areas [m s⁻¹]
    infiltration_capacity_compacted_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil infiltration capacity [m s⁻¹]
    infiltration_capacity_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Maximum leakage [m s⁻¹] from saturated zone
    maximum_leakage::Vector{Float64} = fill(MISSING_VALUE, n)
    # Parameter [m] controlling capillary rise
    cap_hmax::Vector{Float64} = fill(MISSING_VALUE, n)
    # Coefficient [-] controlling capillary rise
    cap_n::Vector{Float64} = fill(MISSING_VALUE, n)
    # Brooks-Corey power coefficient [-] for each soil layer
    brooks_corey_exponent::Vector{SVector{N, Float64}} =
        fill(zeros(SVector{maximum_number_of_layers, Float64}), n)
    # Soil temperature smooth factor [-]
    w_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Controls soil infiltration reduction factor when soil is frozen [-]
    cf_soil::Vector{Float64} = fill(MISSING_VALUE, n)
    # Fraction of compacted area  [-]
    compacted_soil_area_fraction::Vector{Float64} = fill(MISSING_VALUE, n)
    # Controls how roots are linked to water table [m⁻¹]
    wet_root_distribution_parameter::Vector{Float64} = fill(MISSING_VALUE, n)
    # Fraction of the root length density in each soil layer [-]
    rootfraction::Vector{SVector{N, Float64}} =
        fill(zeros(SVector{maximum_number_of_layers, Float64}), n)
    # Soil water pressure head h1 of the root water uptake reduction function (Feddes) [m]
    h1::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil water pressure head h2 of the root water uptake reduction function (Feddes) [m]
    h2::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil water pressure head h3_high of the root water uptake reduction function (Feddes) [m]
    h3_high::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil water pressure head h3_low of the root water uptake reduction function (Feddes) [m]
    h3_low::Vector{Float64} = fill(MISSING_VALUE, n)
    # Soil water pressure head h4 of the root water uptake reduction function (Feddes) [m]
    h4::Vector{Float64} = fill(MISSING_VALUE, n)
    # Root water uptake reduction at soil water pressure head h1 (0.0 or 1.0) [-]
    alpha_h1::Vector{Float64} = fill(MISSING_VALUE, n)
    # Bare soil fraction [-]
    bare_soil_fraction::Vector{Float64} = fill(MISSING_VALUE, n)
    # Vertical hydraulic conductivity profile type
    kv_profile::Kv = nothing
    # Vegetation parameter set
    vegetation_parameter_set::VegetationParameters = VegetationParameters(; n)
end

"Struct for storing SBM soil model boundary conditions"
@kwdef struct SbmSoilBC
    n::Int
    # Water flux at the soil surface [m s⁻¹]
    water_flux_surface::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential transpiration rate [m s⁻¹]
    potential_transpiration::Vector{Float64} = fill(MISSING_VALUE, n)
    # Potential soil evaporation rate [m s⁻¹]
    potential_soilevaporation::Vector{Float64} = fill(MISSING_VALUE, n)
end

"SBM soil model"
@kwdef struct SbmSoilModel{N, M, Kv} <: AbstractSoilModel
    n::Int
    maximum_number_of_layers::Int
    boundary_conditions::SbmSoilBC = SbmSoilBC(; n)
    parameters::SbmSoilParameters{N, M, Kv} =
        SbmSoilParameters(; n, maximum_number_of_layers)
    variables::SbmSoilVariables{N} = SbmSoilVariables(; n, maximum_number_of_layers)
end
