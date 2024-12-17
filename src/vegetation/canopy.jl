abstract type AbstractInterceptionModel{T} end

"Struct for storing interception model variables"
@get_units @grid_loc @with_kw struct InterceptionVariables{T}
    # Canopy potential evaporation [mm Δt⁻¹]
    canopy_potevap::Vector{T}
    # Interception loss by evaporation [mm Δt⁻¹]
    interception_rate::Vector{T}
    # Canopy storage [mm]
    canopy_storage::Vector{T} | "mm"
    # Stemflow [mm Δt⁻¹]
    stemflow::Vector{T}
    # Throughfall [mm Δt⁻¹]
    throughfall::Vector{T}
end

"Initialize interception model variables"
function InterceptionVariables(T::Type{<:AbstractFloat}, n::Int)
    return InterceptionVariables(;
        canopy_potevap = fill(mv, n),
        interception_rate = fill(mv, n),
        canopy_storage = zeros(T, n),
        stemflow = fill(mv, n),
        throughfall = fill(mv, n),
    )
end

"Struct for storing Gash interception model parameters"
@get_units @grid_loc @with_kw struct GashParameters{T}
    # wet canopy [mm Δt⁻¹] and the average precipitation intensity [mm Δt⁻¹] on a saturated canopy
    e_r::Vector{T} | "-"
    vegetation_parameter_set::VegetationParameters{T}
end

"Gash interception model"
@with_kw struct GashInterceptionModel{T} <: AbstractInterceptionModel{T}
    parameters::GashParameters{T}
    variables::InterceptionVariables{T}
end

"Initialize Gash interception model"
function GashInterceptionModel(dataset, config, indices, vegetation_parameter_set)
    lens = lens_input_parameter(
        "vegetation_canopy_water__mean_evaporation-to-mean_precipitation_ratio",
    )
    e_r = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)
    n = length(indices)
    params =
        GashParameters(; e_r = e_r, vegetation_parameter_set = vegetation_parameter_set)
    vars = InterceptionVariables(Float, n)
    model = GashInterceptionModel(; parameters = params, variables = vars)
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
        update_canopy_parameters!(model)
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
@with_kw struct RutterInterceptionModel{T} <: AbstractInterceptionModel{T}
    parameters::VegetationParameters{T}
    variables::InterceptionVariables{T}
end

"Initialize Rutter interception model"
function RutterInterceptionModel(vegetation_parameter_set, n)
    vars = InterceptionVariables(n)
    model =
        RutterInterceptionModel(; parameters = vegetation_parameter_set, variables = vars)
    return model
end

"Update Rutter interception model for a single timestep"
function update!(model::RutterInterceptionModel, atmospheric_forcing::AtmosphericForcing)
    (; leaf_area_index, canopygapfraction, cmax, kc) =
        model.parameters.vegetation_parameter_set
    (; canopy_potevap, throughfall, interception_rate, stemflow, canopy_storage) =
        model.variables
    (; precipitation, potential_evaporation) = atmospheric_forcing
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(model)
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
function update_canopy_parameters!(model::AbstractInterceptionModel)
    (;
        leaf_area_index,
        storage_wood,
        kext,
        storage_specific_leaf,
        canopygapfraction,
        cmax,
    ) = model.parameters.vegetation_parameter_set

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