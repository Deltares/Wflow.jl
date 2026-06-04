abstract type AbstractInterceptionModel end

"Struct for storing interception model variables"
@with_kw struct InterceptionVariables
    n_cells::Int
    # Canopy potential evaporation [m s⁻¹]
    canopy_potevap::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Interception loss by evaporation [m s⁻¹]
    interception_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Canopy storage [m]
    canopy_storage::Vector{Float64} = zeros(n_cells)
    # Stemflow [m s⁻¹]
    stemflow::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Throughfall [m s⁻¹]
    throughfall::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct for storing Gash interception model parameters"
@with_kw struct GashParameters
    # ratio [-] of wet canopy [m s⁻¹] and the average precipitation intensity [m s⁻¹] on a saturated canopy
    evaporation_to_precipitation_ratio::Vector{Float64}
    vegetation_parameter_set::VegetationParameters
end

"Gash interception model"
@with_kw struct GashInterceptionModel <: AbstractInterceptionModel
    n_cells::Int
    parameters::GashParameters
    variables::InterceptionVariables = InterceptionVariables(; n_cells)
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
    n_cells = length(indices)
    parameters =
        GashParameters(; evaporation_to_precipitation_ratio, vegetation_parameter_set)
    model = GashInterceptionModel(; n_cells, parameters)
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
    n_cells = length(precipitation)
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(interception_model.parameters.vegetation_parameter_set)
        threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
            canopyfraction = 1.0 - canopy_gap_fraction[cell_idx]
            ewet =
                canopyfraction *
                potential_evaporation[cell_idx] *
                crop_coefficient[cell_idx]
            evaporation_to_precipitation_ratio[cell_idx] =
                precipitation[cell_idx] > 0.0 ?
                min(
                    0.25,
                    ewet / max(
                        to_SI(1e-4, MM_PER_DT; dt_val = dt),
                        canopyfraction * precipitation[cell_idx],
                    ),
                ) : 0.0
        end
    end
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        canopy_potevap[cell_idx] =
            crop_coefficient[cell_idx] *
            potential_evaporation[cell_idx] *
            (1.0 - canopy_gap_fraction[cell_idx])
        throughfall[cell_idx],
        interception_rate[cell_idx],
        stemflow[cell_idx],
        canopy_storage[cell_idx] = rainfall_interception_gash(
            maximum_canopy_storage[cell_idx],
            evaporation_to_precipitation_ratio[cell_idx],
            canopy_gap_fraction[cell_idx],
            precipitation[cell_idx],
            canopy_storage[cell_idx],
            canopy_potevap[cell_idx],
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
function RutterInterceptionModel(
    vegetation_parameter_set::VegetationParameters,
    n_cells::Int,
)
    variables = InterceptionVariables(; n_cells)
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
    n_cells = length(precipitation)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        canopy_potevap[cell_idx] =
            crop_coefficient[cell_idx] *
            potential_evaporation[cell_idx] *
            (1.0 - canopy_gap_fraction[cell_idx])
        throughfall[cell_idx],
        interception_rate[cell_idx],
        stemflow[cell_idx],
        canopy_storage[cell_idx] = rainfall_interception_modrut(
            precipitation[cell_idx],
            canopy_potevap[cell_idx],
            canopy_storage[cell_idx],
            canopy_gap_fraction[cell_idx],
            maximum_canopy_storage[cell_idx],
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

    n_cells = length(leaf_area_index)
    threaded_foreach(1:n_cells; basesize = 1000) do cell_idx
        maximum_canopy_storage[cell_idx] =
            storage_specific_leaf[cell_idx] * leaf_area_index[cell_idx] +
            storage_wood[cell_idx]
        canopy_gap_fraction[cell_idx] =
            exp(-light_extinction_coefficient[cell_idx] * leaf_area_index[cell_idx])
    end
    return nothing
end

"Return potential transpiration rate based on the interception rate"
get_potential_transpiration(interception_model::AbstractInterceptionModel) = @. max(
    0.0,
    interception_model.variables.canopy_potevap -
    interception_model.variables.interception_rate,
)
