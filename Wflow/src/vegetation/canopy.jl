abstract type AbstractInterceptionModel end

"Struct for storing interception model variables"
@with_kw struct InterceptionVariables
    # Canopy potential evaporation [mm Δt⁻¹]
    canopy_potevap::Vector{Float}
    # Interception loss by evaporation [mm Δt⁻¹]
    interception_rate::Vector{Float}
    # Canopy storage [mm]
    canopy_storage::Vector{Float}
    # Stemflow [mm Δt⁻¹]
    stemflow::Vector{Float}
    # Throughfall [mm Δt⁻¹]
    throughfall::Vector{Float}
end

"Initialize interception model variables"
function InterceptionVariables(n::Int)
    return InterceptionVariables(;
        canopy_potevap = fill(MISSING_VALUE, n),
        interception_rate = fill(MISSING_VALUE, n),
        canopy_storage = zeros(Float, n),
        stemflow = fill(MISSING_VALUE, n),
        throughfall = fill(MISSING_VALUE, n),
    )
end

"Struct for storing Gash interception model parameters"
@with_kw struct GashParameters
    # ratio [-] of wet canopy [mm Δt⁻¹] and the average precipitation intensity [mm Δt⁻¹] on a saturated canopy
    e_r::Vector{Float}
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
    lens = lens_input_parameter(
        config,
        "vegetation_canopy_water__mean_evaporation-to-mean_precipitation_ratio",
    )
    e_r = ncread(dataset, config, lens; sel = indices, defaults = 0.1, type = Float)
    n = Int(length(indices))
    params =
        GashParameters(; e_r = e_r, vegetation_parameter_set = vegetation_parameter_set)
    vars = InterceptionVariables(n)
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
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(model)
        AK.foreachindex(precipitation; scheduler = :polyester, min_elems = 1000) do i
            canopyfraction = 1.0 - canopygapfraction[i]
            ewet = canopyfraction * potential_evaporation[i] * kc[i]
            e_r[i] =
                precipitation[i] > 0.0 ?
                min(0.25, ewet / max(0.0001, canopyfraction * precipitation[i])) : 0.0
        end
    end
    AK.foreachindex(throughfall; scheduler = :polyester, min_elems = 1000) do i
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
    vars = InterceptionVariables(n)
    model =
        RutterInterceptionModel(; parameters = vegetation_parameter_set, variables = vars)
    return model
end

"Update Rutter interception model for a single timestep"
function update!(model::RutterInterceptionModel, atmospheric_forcing::AtmosphericForcing)
    (; leaf_area_index, canopygapfraction, cmax, kc) = model.parameters
    (; canopy_potevap, throughfall, interception_rate, stemflow, canopy_storage) =
        model.variables
    (; precipitation, potential_evaporation) = atmospheric_forcing
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(model)
    end
    AK.foreachindex(precipitation; scheduler = :polyester, min_elems = 1000) do i
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
    ) = model.parameters

    AK.foreachindex(cmax; scheduler = :polyester, min_elems = 1000) do i
        cmax[i] = storage_specific_leaf[i] * leaf_area_index[i] + storage_wood[i]
        canopygapfraction[i] = exp(-kext[i] * leaf_area_index[i])
    end
    return nothing
end

"Return potential transpiration rate based on the interception rate"
get_potential_transpiration(model::AbstractInterceptionModel) =
    @. max(0.0, model.variables.canopy_potevap - model.variables.interception_rate)