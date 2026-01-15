abstract type AbstractInterceptionModel end

"Struct for storing interception model variables"
@with_kw struct InterceptionVariables
    n::Int
    # Canopy potential evaporation [mm Δt⁻¹]
    canopy_potevap::Vector{Float64} = fill(MISSING_VALUE, n)
    # Interception loss by evaporation [mm Δt⁻¹]
    interception_rate::Vector{Float64} = fill(MISSING_VALUE, n)
    # Canopy storage [mm]
    canopy_storage::Vector{Float64} = zeros(n)
    # Stemflow [mm Δt⁻¹]
    stemflow::Vector{Float64} = fill(MISSING_VALUE, n)
    # Throughfall [mm Δt⁻¹]
    throughfall::Vector{Float64} = fill(MISSING_VALUE, n)
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
    indices::Vector{CartesianIndex{2}},
    vegetation_parameter_set::VegetationParameters,
)
    e_r = ncread(
        dataset,
        config,
        "vegetation_canopy_water__mean_evaporation_to_mean_precipitation_ratio",
        LandHydrologySBM;
        sel = indices,
    )
    n = length(indices)
    parameters = GashParameters(; e_r, vegetation_parameter_set)
    variables = InterceptionVariables(; n)
    model = GashInterceptionModel(; parameters, variables)
    return model
end

"Update Gash interception model for a single timestep"
function update!(model::GashInterceptionModel, atmospheric_forcing::AtmosphericForcing)
    (; leaf_area_index, canopygapfraction, cmax, kc) =
        model.parameters.vegetation_parameter_set
    (; canopy_potevap, throughfall, interception_rate, stemflow, canopy_storage) =
        model.variables
    (; precipitation, potential_evaporation) = atmospheric_forcing
    e_r = model.parameters.e_r
    n = length(precipitation)
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(model.parameters.vegetation_parameter_set)
        threaded_foreach(1:n; basesize = 1000) do i
            canopyfraction = 1.0 - canopygapfraction[i]
            ewet = canopyfraction * potential_evaporation[i] * kc[i]
            e_r[i] =
                precipitation[i] > 0.0 ?
                min(0.25, ewet / max(0.0001, canopyfraction * precipitation[i])) : 0.0
        end
    end
    threaded_foreach(1:n; basesize = 1000) do i
        canopy_potevap[i] = kc[i] * potential_evaporation[i] * (1.0 - canopygapfraction[i])
        throughfall[i], interception_rate[i], stemflow[i], canopy_storage[i] =
            rainfall_interception_gash(
                cmax[i],
                e_r[i],
                canopygapfraction[i],
                precipitation[i],
                canopy_storage[i],
                canopy_potevap[i],
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
    model = RutterInterceptionModel(; parameters = vegetation_parameter_set, variables)
    return model
end

"Update Rutter interception model for a single timestep"
function update!(model::RutterInterceptionModel, atmospheric_forcing::AtmosphericForcing)
    (; leaf_area_index, canopygapfraction, cmax, kc) = model.parameters
    (; canopy_potevap, throughfall, interception_rate, stemflow, canopy_storage) =
        model.variables
    (; precipitation, potential_evaporation) = atmospheric_forcing
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(model.parameters)
    end
    n = length(precipitation)
    threaded_foreach(1:n; basesize = 1000) do i
        canopy_potevap[i] = kc[i] * potential_evaporation[i] * (1.0 - canopygapfraction[i])
        throughfall[i], interception_rate[i], stemflow[i], canopy_storage[i] =
            rainfall_interception_modrut(
                precipitation[i],
                canopy_potevap[i],
                canopy_storage[i],
                canopygapfraction[i],
                cmax[i],
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

    n = length(leaf_area_index)
    threaded_foreach(1:n; basesize = 1000) do i
        cmax[i] = storage_specific_leaf[i] * leaf_area_index[i] + storage_wood[i]
        canopygapfraction[i] = exp(-kext[i] * leaf_area_index[i])
    end
    return nothing
end

"Return potential transpiration rate based on the interception rate"
get_potential_transpiration(model::AbstractInterceptionModel) =
    @. max(0.0, model.variables.canopy_potevap - model.variables.interception_rate)
