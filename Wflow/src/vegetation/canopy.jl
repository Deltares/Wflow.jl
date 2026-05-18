abstract type AbstractInterceptionModel end

"Struct for storing interception model variables"
@with_kw struct InterceptionVariables
    n_cells::Int
    # Canopy potential evaporation [mm Δt⁻¹]
    canopy_potevap::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Interception loss by evaporation [mm Δt⁻¹]
    interception_rate::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Canopy storage [mm]
    canopy_storage::Vector{Float64} = zeros(n_cells)
    # Stemflow [mm Δt⁻¹]
    stemflow::Vector{Float64} = fill(MISSING_VALUE, n_cells)
    # Throughfall [mm Δt⁻¹]
    throughfall::Vector{Float64} = fill(MISSING_VALUE, n_cells)
end

"Struct for storing Gash interception model parameters"
@with_kw struct GashParameters
    # ratio [-] of wet canopy [mm Δt⁻¹] and the average precipitation intensity [mm Δt⁻¹] on a saturated canopy
    e_r::Vector{Float64}
    vegetation_parameter_set::VegetationParameters
end

"Gash interception model"
@with_kw struct GashInterceptionModel <: AbstractInterceptionModel
    parameters::GashParameters
    variables::InterceptionVariables
end

"Initialize Gash interception model"
function GashInterceptionModel(
    dataset::NCDataset,
    config::Config,
    land_indices_2d::Vector{CartesianIndex{2}},
    vegetation_parameter_set::VegetationParameters,
)
    e_r = ncread(
        dataset,
        config,
        "vegetation_canopy_water__mean_evaporation_to_mean_precipitation_ratio",
        LandHydrologySBM;
        sel = land_indices_2d,
    )
    n_cells = length(land_indices_2d)
    parameters = GashParameters(; e_r, vegetation_parameter_set)
    variables = InterceptionVariables(; n_cells)
    interception_model = GashInterceptionModel(; parameters, variables)
    return interception_model
end

"Update Gash interception model for a single timestep"
function update_interception_model!(
    interception_model::GashInterceptionModel,
    atmospheric_forcing::AtmosphericForcing,
)
    (; parameters, variables) = interception_model
    (; leaf_area_index, canopygapfraction, cmax, kc) = parameters.vegetation_parameter_set
    (; canopy_potevap, throughfall, interception_rate, stemflow, canopy_storage, n_cells) =
        variables
    e_r = parameters.e_r
    (; precipitation, potential_evaporation) = atmospheric_forcing
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(interception_model.parameters.vegetation_parameter_set)
        threaded_foreach(1:n_cells; basesize = 1000) do land_cell_idx
            canopyfraction = 1.0 - canopygapfraction[land_cell_idx]
            ewet = canopyfraction * potential_evaporation[land_cell_idx] * kc[land_cell_idx]
            e_r[land_cell_idx] =
                precipitation[land_cell_idx] > 0.0 ?
                min(
                    0.25,
                    ewet / max(0.0001, canopyfraction * precipitation[land_cell_idx]),
                ) : 0.0
        end
    end
    threaded_foreach(1:n_cells; basesize = 1000) do land_cell_idx
        canopy_potevap[land_cell_idx] =
            kc[land_cell_idx] *
            potential_evaporation[land_cell_idx] *
            (1.0 - canopygapfraction[land_cell_idx])
        throughfall[land_cell_idx],
        interception_rate[land_cell_idx],
        stemflow[land_cell_idx],
        canopy_storage[land_cell_idx] = rainfall_interception_gash(
            cmax[land_cell_idx],
            e_r[land_cell_idx],
            canopygapfraction[land_cell_idx],
            precipitation[land_cell_idx],
            canopy_storage[land_cell_idx],
            canopy_potevap[land_cell_idx],
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
)
    (; leaf_area_index, canopygapfraction, cmax, kc) = interception_model.parameters
    (; canopy_potevap, throughfall, interception_rate, stemflow, canopy_storage) =
        interception_model.variables
    (; precipitation, potential_evaporation, n_cells) = atmospheric_forcing
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(interception_model.parameters)
    end
    threaded_foreach(1:n_cells; basesize = 1000) do land_cell_idx
        canopy_potevap[land_cell_idx] =
            kc[land_cell_idx] *
            potential_evaporation[land_cell_idx] *
            (1.0 - canopygapfraction[land_cell_idx])
        throughfall[land_cell_idx],
        interception_rate[land_cell_idx],
        stemflow[land_cell_idx],
        canopy_storage[land_cell_idx] = rainfall_interception_modrut(
            precipitation[land_cell_idx],
            canopy_potevap[land_cell_idx],
            canopy_storage[land_cell_idx],
            canopygapfraction[land_cell_idx],
            cmax[land_cell_idx],
        )
    end
    return nothing
end

"Update canopy parameters `cmax` and `canopygapfraction` based on `leaf_area_index` for a single timestep"
function update_canopy_parameters!(parameters::VegetationParameters)
    (;
        leaf_area_index,
        storage_wood,
        kext,
        storage_specific_leaf,
        canopygapfraction,
        cmax,
    ) = parameters

    n_cells = length(leaf_area_index)
    threaded_foreach(1:n_cells; basesize = 1000) do land_cell_idx
        cmax[land_cell_idx] =
            storage_specific_leaf[land_cell_idx] * leaf_area_index[land_cell_idx] +
            storage_wood[land_cell_idx]
        canopygapfraction[land_cell_idx] =
            exp(-kext[land_cell_idx] * leaf_area_index[land_cell_idx])
    end
    return nothing
end

"Return potential transpiration rate based on the interception rate"
get_potential_transpiration(interception_model::AbstractInterceptionModel) = @. max(
    0.0,
    interception_model.variables.canopy_potevap -
    interception_model.variables.interception_rate,
)
