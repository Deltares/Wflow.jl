abstract type AbstractInterceptionModel end

"Struct for storing interception model variables"
@with_kw struct InterceptionVariables
    n::Int
    # Canopy potential evaporation [m s⁻¹]
    canopy_potevap::Vector{Float64} = fill(MISSING_VALUE, n)
    # Interception loss by evaporation [m s⁻¹]
    interception_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Canopy storage [m]
    canopy_storage::Vector{Float64} = zeros(n)
    # Stemflow [m s⁻¹]
    stemflow::Vector{Float64} = fill(MISSING_VALUE, n)
    # Throughfall [m s⁻¹]
    throughfall::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Struct for storing Gash interception model parameters"
@with_kw struct GashParameters
    # ratio [-] of wet canopy [m s⁻¹] and the average precipitation intensity [m s⁻¹] on a saturated canopy
    evaporation_to_precipitation_ratio::Vector{Float64}
    vegetation_parameter_set::VegetationParameters
end

"Gash interception model"
@with_kw struct GashInterceptionModel <: AbstractInterceptionModel
    n::Int
    parameters::GashParameters
    variables::InterceptionVariables = InterceptionVariables(; n)
end

"Initialize Gash interception model"
function GashInterceptionModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    vegetation_parameter_set::VegetationParameters,
)
    evaporation_to_precipitation_ratio = ncread(
        dataset,
        config,
        "vegetation_canopy_water__mean_evaporation_to_mean_precipitation_ratio",
        LandHydrologySBM;
        sel = indices,
    )
    n = length(indices)
    parameters =
        GashParameters(; evaporation_to_precipitation_ratio, vegetation_parameter_set)
    model = GashInterceptionModel(; n, parameters)
    return model
end

"Update Gash interception model for a single timestep"
function update_interception_model!(
    interception_model::GashInterceptionModel,
    atmospheric_forcing::AtmosphericForcing,
    dt::Float64,
)
    (; evaporation_to_precipitation_ratio, vegetation_parameter_set) =
        interception_model.parameters
    (; leaf_area_index, canopy_gap_fraction, maximum_canopy_storage, crop_coefficient) =
        vegetation_parameter_set
    (; canopy_potevap, throughfall, interception_rate, stemflow, canopy_storage) =
        interception_model.variables
    (; precipitation, potential_evaporation) = atmospheric_forcing
    n = length(precipitation)
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(interception_model.parameters.vegetation_parameter_set)
        threaded_foreach(1:n; basesize = 1000) do i
            canopyfraction = 1.0 - canopy_gap_fraction[i]
            ewet = canopyfraction * potential_evaporation[i] * crop_coefficient[i]
            evaporation_to_precipitation_ratio[i] =
                precipitation[i] > 0.0 ?
                min(
                    0.25,
                    ewet / max(
                        to_SI(1e-4, MM_PER_DT; dt_val = dt),
                        canopyfraction * precipitation[i],
                    ),
                ) : 0.0
        end
    end
    threaded_foreach(1:n; basesize = 1000) do i
        canopy_potevap[i] =
            crop_coefficient[i] * potential_evaporation[i] * (1.0 - canopy_gap_fraction[i])
        throughfall[i], interception_rate[i], stemflow[i], canopy_storage[i] =
            rainfall_interception_gash(
                maximum_canopy_storage[i],
                evaporation_to_precipitation_ratio[i],
                canopy_gap_fraction[i],
                precipitation[i],
                canopy_storage[i],
                canopy_potevap[i],
                dt,
            )
    end
    return nothing
end

"Rutter interception model"
@with_kw struct RutterInterceptionModel <: AbstractInterceptionModel
    parameters::VegetationParameters
    variables::InterceptionVariables
end

"Initialize Rutter interception model"
function RutterInterceptionModel(vegetation_parameter_set::VegetationParameters, n::Int)
    variables = InterceptionVariables(; n)
    interception_model =
        RutterInterceptionModel(; parameters = vegetation_parameter_set, variables)
    return interception_model
end

"Update Rutter interception model for a single timestep"
function update_interception_model!(
    interception_model::RutterInterceptionModel,
    atmospheric_forcing::AtmosphericForcing,
    dt::Float64,
)
    (; leaf_area_index, canopy_gap_fraction, maximum_canopy_storage, crop_coefficient) =
        interception_model.parameters
    (; canopy_potevap, throughfall, interception_rate, stemflow, canopy_storage) =
        interception_model.variables
    (; precipitation, potential_evaporation) = atmospheric_forcing
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(interception_model.parameters)
    end
    n = length(precipitation)
    threaded_foreach(1:n; basesize = 1000) do i
        canopy_potevap[i] =
            crop_coefficient[i] * potential_evaporation[i] * (1.0 - canopy_gap_fraction[i])
        throughfall[i], interception_rate[i], stemflow[i], canopy_storage[i] =
            rainfall_interception_modrut(
                precipitation[i],
                canopy_potevap[i],
                canopy_storage[i],
                canopy_gap_fraction[i],
                maximum_canopy_storage[i],
                dt,
            )
    end
    return nothing
end

"Update canopy parameters `maximum_canopy_storage` and `canopy_gap_fraction` based on `leaf_area_index` for a single timestep"
function update_canopy_parameters!(parameters::VegetationParameters)
    (;
        leaf_area_index,
        storage_wood,
        light_extinction_coefficient,
        storage_specific_leaf,
        canopy_gap_fraction,
        maximum_canopy_storage,
    ) = parameters

    n = length(leaf_area_index)
    threaded_foreach(1:n; basesize = 1000) do i
        maximum_canopy_storage[i] =
            storage_specific_leaf[i] * leaf_area_index[i] + storage_wood[i]
        canopy_gap_fraction[i] = exp(-light_extinction_coefficient[i] * leaf_area_index[i])
    end
    return nothing
end

"Return potential transpiration rate based on the interception rate"
get_potential_transpiration(interception_model::AbstractInterceptionModel) = @. max(
    0.0,
    interception_model.variables.canopy_potevap -
    interception_model.variables.interception_rate,
)
