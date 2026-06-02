abstract type AbstractSedimentToRiverModel end

"Struct to store total sediment reaching the river model variables"
@with_data_lookup struct SedimentToRiverVariables
    n::Int
    # Total sediment rate to the river [kg s⁻¹]
    "land_surface_water_sediment__to_river_mass_flow_rate"
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store total sediment reaching the river model boundary conditions"
@kwdef struct SedimentToRiverBC
    n::Int
    # Deposition material rate [kg s⁻¹]
    deposition::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store total sediment reaching the river model"
@kwdef struct SedimentToRiverModel <: AbstractSedimentToRiverModel
    n::Int
    boundary_conditions::SedimentToRiverBC = SedimentToRiverBC(; n)
    variables::SedimentToRiverVariables = SedimentToRiverVariables(; n)
end

"Initialize total sediment reaching the river model"
function SedimentToRiverModel(
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
)
    n = length(indices)
    variables = SedimentToRiverVariables(data_lookup; n)
    sediment_to_river_model = SedimentToRiverModel(; n, variables)
    return sediment_to_river_model
end

"Update total sediment reaching the river model boundary conditions"
function update_bc_sediment_to_river_model!(
    sediment_to_river_model::SedimentToRiverModel,
    sediment_transport_model::SedimentLandTransportModel,
)
    (; deposition) = sediment_to_river_model.boundary_conditions
    @. deposition = sediment_transport_model.variables.deposition
end

"Update total sediment reaching the river model for a single timestep"
function update_sediment_to_river_model!(
    sediment_to_river_model::SedimentToRiverModel,
    rivers::Vector{Bool},
    dt::Float64,
)
    (; deposition) = sediment_to_river_model.boundary_conditions
    (; sediment_rate) = sediment_to_river_model.variables

    for (i, river) in enumerate(rivers)
        sediment_rate[i] = river ? deposition[i] : 0.0
    end
end

"Struct to store differentiated sediment reaching the river model variables"
@with_data_lookup struct SedimentToRiverDifferentiationVariables
    n::Int
    # Total sediment rate [kg s⁻¹]
    "land_surface_water_sediment__to_river_mass_flow_rate"
    sediment_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Clay rate [kg s⁻¹]
    "land_surface_water_clay__to_river_mass_flow_rate"
    clay_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt rate [kg s⁻¹]
    "land_surface_water_silt__to_river_mass_flow_rate"
    silt_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand rate [kg s⁻¹]
    "land_surface_water_sand__to_river_mass_flow_rate"
    sand_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates rate [kg s⁻¹]
    "land_surface_water_small_aggregates__to_river_mass_flow_rate"
    sagg_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates rate [kg s⁻¹]
    "land_surface_water_large_aggregates__to_river_mass_flow_rate"
    lagg_rate::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store differentiated sediment reaching the river model boundary conditions"
@kwdef struct SedimentToRiverDifferentiationBC
    n::Int
    # Clay deposition rate [kg s⁻¹]
    deposition_clay::Vector{Float64} = fill(MISSING_VALUE, n)
    # Silt deposition rate [kg s⁻¹]
    deposition_silt::Vector{Float64} = fill(MISSING_VALUE, n)
    # Sand deposition rate [kg s⁻¹]
    deposition_sand::Vector{Float64} = fill(MISSING_VALUE, n)
    # Small aggregates deposition rate [kg s⁻¹]
    deposition_sagg::Vector{Float64} = fill(MISSING_VALUE, n)
    # Large aggregates deposition rate [kg s⁻¹]
    deposition_lagg::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct to store differentiated sediment reaching the river model"
@kwdef struct SedimentToRiverDifferentiationModel <: AbstractSedimentToRiverModel
    n::Int
    boundary_conditions::SedimentToRiverDifferentiationBC =
        SedimentToRiverDifferentiationBC(; n)
    variables::SedimentToRiverDifferentiationVariables =
        SedimentToRiverDifferentiationVariables(; n)
end

"Initialize differentiated sediment reaching the river model"
function SedimentToRiverDifferentiationModel(
    indices::Vector{CartesianIndex{2}};
    data_lookup::DataLookup = DataLookup(),
)
    n = length(indices)
    variables = SedimentToRiverDifferentiationVariables(data_lookup; n)
    sediment_to_river_model = SedimentToRiverDifferentiationModel(; n, variables)
    return sediment_to_river_model
end

"Update differentiated sediment reaching the river model boundary conditions"
function update_bc_sediment_to_river_model!(
    sediment_to_river_model::SedimentToRiverDifferentiationModel,
    sediment_transport_model::SedimentLandTransportDifferentiationModel,
)
    (;
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_sagg,
        deposition_lagg,
    ) = sediment_to_river_model.boundary_conditions
    @. deposition_clay = sediment_transport_model.variables.deposition_clay
    @. deposition_silt = sediment_transport_model.variables.deposition_silt
    @. deposition_sand = sediment_transport_model.variables.deposition_sand
    @. deposition_sagg = sediment_transport_model.variables.deposition_sagg
    @. deposition_lagg = sediment_transport_model.variables.deposition_lagg
end

"Update differentiated sediment reaching the river model for a single timestep"
function update_sediment_to_river_model!(
    sediment_to_river_model::SedimentToRiverDifferentiationModel,
    rivers::Vector{Bool},
    dt::Float64,
)
    (;
        deposition_clay,
        deposition_silt,
        deposition_sand,
        deposition_sagg,
        deposition_lagg,
    ) = sediment_to_river_model.boundary_conditions
    (; sediment_rate, clay_rate, silt_rate, sand_rate, sagg_rate, lagg_rate) =
        sediment_to_river_model.variables

    for (i, river) in enumerate(rivers)
        if river
            clay_rate[i] = deposition_clay[i]
            silt_rate[i] = deposition_silt[i]
            sand_rate[i] = deposition_sand[i]
            sagg_rate[i] = deposition_sagg[i]
            lagg_rate[i] = deposition_lagg[i]
        else
            clay_rate[i] = 0.0
            silt_rate[i] = 0.0
            sand_rate[i] = 0.0
            sagg_rate[i] = 0.0
            lagg_rate[i] = 0.0
        end
    end
    @. sediment_rate = clay_rate + silt_rate + sand_rate + sagg_rate + lagg_rate
end
