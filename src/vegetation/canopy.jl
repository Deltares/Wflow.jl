@get_units @grid_loc @with_kw struct InterceptionModelVars{T}
    # Canopy potential evaporation [mm Δt⁻¹]
    canopy_potevap::Vector{T}
    # Interception loss by evaporation [mm Δt⁻¹]
    interception_flux::Vector{T}
    # Canopy storage [mm]
    canopy_storage::Vector{T} | "mm"
    # Stemflow [mm Δt⁻¹]
    stemflow::Vector{T}
    # Throughfall [mm Δt⁻¹]
    throughfall::Vector{T}
end

function interception_model_vars(n)
    vars = InterceptionModelVars(;
        canopy_potevap = fill(mv, n),
        interception_flux = fill(mv, n),
        canopy_storage = fill(0.0, n),
        stemflow = fill(mv, n),
        throughfall = fill(mv, n),
    )
    return vars
end

abstract type AbstractInterceptionModel end

@get_units @grid_loc @with_kw struct GashParameters{T}
    # wet canopy [mm Δt⁻¹] and the average precipitation intensity [mm Δt⁻¹] on a saturated canopy
    e_r::Vector{T} | "-"
    vegetation_parameter_set::VegetationParameters{T}
end

@with_kw struct GashInterceptionModel{T} <: AbstractInterceptionModel
    parameters::GashParameters{T}
    variables::InterceptionModelVars{T}
end

@with_kw struct RutterInterceptionModel{T} <: AbstractInterceptionModel
    parameters::VegetationParameters{T}
    variables::InterceptionModelVars{T}
end

function initialize_gash_interception_model(nc, config, inds, vegetation_parameter_set)
    e_r = ncread(
        nc,
        config,
        "vertical.bucket.parameters.eoverr";
        sel = inds,
        defaults = 0.1,
        type = Float,
    )

    params =
        GashParameters(; e_r = e_r, vegetation_parameter_set = vegetation_parameter_set)
    vars = interception_model_vars(length(inds))
    model = GashInterceptionModel(; parameters = params, variables = vars)
    return model
end

function initialize_rutter_interception_model(vegetation_parameter_set, n)
    vars = interception_model_vars(n)
    model =
        RutterInterceptionModel(; parameters = vegetation_parameter_set, variables = vars)
    return model
end

function update_canopy_parameters!(model::I) where {I <: AbstractInterceptionModel}
    (; leaf_area_index, swood, kext, sl, canopygapfraction, cmax) =
        model.parameters.vegetation_parameter_set

    n = length(leaf_area_index)
    threaded_foreach(1:n; basesize = 1000) do i
        cmax[i] = sl[i] * leaf_area_index[i] + swood[i]
        canopygapfraction[i] = exp(-kext[i] * leaf_area_index[i])
    end
end

get_canopygapfraction(model::GashInterceptionModel) =
    model.parameters.vegetation_parameter_set.canopygapfraction

get_canopygapfraction(model::RutterInterceptionModel) = model.parameters.canopygapfraction

get_potential_transpiration(model::AbstractInterceptionModel) =
    @. max(0.0, model.variables.canopy_potevap - model.variables.interception_flux)

function update!(model::GashInterceptionModel, atmospheric_forcing::AtmosphericForcing)
    (; leaf_area_index, canopygapfraction, cmax, kc) =
        model.parameters.vegetation_parameter_set
    (; canopy_potevap, throughfall, interception_flux, stemflow, canopy_storage) =
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
        throughfall[i], interception_flux[i], stemflow[i], canopy_storage[i] =
            rainfall_interception_gash(
                cmax[i],
                e_r[i],
                canopygapfraction[i],
                precipitation[i],
                canopy_storage[i],
                canopy_potevap[i],
            )
    end
end

function update!(model::RutterInterceptionModel, atmospheric_forcing::AtmosphericForcing)
    (; leaf_area_index, canopygapfraction, cmax, kc) =
        model.parameters.vegetation_parameter_set
    (; canopy_potevap, throughfall, interception_flux, stemflow, canopy_storage) =
        model.variables
    (; precipitation, potential_evaporation) = atmospheric_forcing
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(model)
    end
    n = length(precipitation)
    threaded_foreach(1:n; basesize = 1000) do i
        canopy_potevap[i] = kc[i] * potential_evaporation[i] * (1.0 - canopygapfraction[i])
        throughfall[i], interception_flux[i], stemflow[i], canopy_storage[i] =
            rainfall_interception_modrut(
                precipitation[i],
                canopy_potevap[i],
                canopy_storage[i],
                canopygapfraction[i],
                cmax[i],
            )
    end
end
