@get_units @with_kw struct InterceptionModelVars{T}
    # Canopy potential evaporation [mm Δt⁻¹]
    canopy_potevap::Vector{T}
    # Interception loss by evaporation [mm Δt⁻¹]
    interception::Vector{T}
    # Canopy storage [mm]
    canopystorage::Vector{T} | "mm"
    # Stemflow [mm Δt⁻¹]
    stemflow::Vector{T}
    # Throughfall [mm Δt⁻¹]
    throughfall::Vector{T}
end

function interception_model_vars(n)
    vars = InterceptionModelVars(
        canopy_potevap = fill(mv, n),
        interception = fill(mv, n),
        canopystorage = fill(0.0, n),
        stemflow = fill(mv, n),
        throughfall = fill(mv, n),
    )
    return vars
end

@get_units @with_kw struct VegetationParameters{T}
    # Leaf area index [m² m⁻²]
    leaf_area_index::Union{Vector{T},Nothing} | "m2 m-2"
    # Storage woody part of vegetation [mm]
    swood::Union{Vector{T},Nothing} | "mm"
    # Extinction coefficient [-] (to calculate canopy gap fraction)
    kext::Union{Vector{T},Nothing} | "-"
    # Specific leaf storage [mm]
    sl::Union{Vector{T},Nothing} | "mm"
    # Canopy gap fraction [-]
    canopygapfraction::Vector{T} | "-"
    # Maximum canopy storage [mm] 
    cmax::Vector{T} | "mm"
    # Rooting depth [mm]
    rootingdepth::Vector{T} | "mm"
    # Crop coefficient Kc [-]
    kc::Vector{T} | "-"
end

abstract type AbstractInterceptionModel{T} end

@get_units @with_kw struct GashParameters{T}
    # wet canopy [mm Δt⁻¹] and the average precipitation intensity [mm Δt⁻¹] on a saturated canopy
    e_r::Vector{T} | "-"
    veg_param_set::VegetationParameters{T} | "-"
end

@get_units @with_kw struct GashInterceptionModel{T} <: AbstractInterceptionModel{T}
    parameters::GashParameters{T} | "-"
    variables::InterceptionModelVars{T} | "-"
end

@get_units @with_kw struct RutterInterceptionModel{T} <: AbstractInterceptionModel{T}
    parameters::VegetationParameters{T} | "-"
    variables::InterceptionModelVars{T} | "-"
end

function initialize_vegetation_params(nc, config, inds)
    n = length(inds)
    rootingdepth = ncread(
        nc,
        config,
        "vertical.rootingdepth";
        sel = inds,
        defaults = 750.0,
        type = Float,
    )
    kc = ncread(nc, config, "vertical.kc"; sel = inds, defaults = 1.0, type = Float)
    if haskey(
        config.input.vertical.interception_model.parameters.veg_param_set,
        "leaf_area_index",
    )
        sl = ncread(
            nc,
            config,
            "vertical.specific_leaf";
            optional = false,
            sel = inds,
            type = Float,
        )
        swood = ncread(
            nc,
            config,
            "vertical.storage_wood";
            optional = false,
            sel = inds,
            type = Float,
        )
        kext =
            ncread(nc, config, "vertical.kext"; optional = false, sel = inds, type = Float)
        vegetation_params = VegetationParameters(
            leaf_area_index = fill(mv, n),
            swood = swood,
            kext = kext,
            sl = sl,
            canopygapfraction = fill(mv, n),
            cmax = fill(mv, n),
            rootingdepth = rootingdepth,
            kc = kc,
        )
    else
        canopygapfraction = ncread(
            nc,
            config,
            "vertical.canopygapfraction";
            sel = inds,
            defaults = 0.1,
            type = Float,
        )
        cmax = ncread(nc, config, "vertical.cmax"; sel = inds, defaults = 1.0, type = Float)
        vegetation_params = VegetationParameters(
            leaf_area_index = nothing,
            swood = nothing,
            kext = nothing,
            sl = nothing,
            canopygapfraction = canopygapfraction,
            cmax = cmax,
            rootingdepth = rootingdepth,
            kc = kc,
        )
    end
    return vegetation_params
end

function initialize_gash_interception_model(nc, config, inds, vegetation_params)
    e_r = ncread(nc, config, "vertical.eoverr"; sel = inds, defaults = 0.1, type = Float)

    params = GashParameters(e_r = e_r, veg_param_set = vegetation_params)
    vars = interception_model_vars(length(inds))
    model = GashInterceptionModel(parameters = params, variables = vars)
    return model
end

function initialize_rutter_interception_model(vegetation_params, n)
    vars = interception_model_vars(n)
    model = RutterInterceptionModel(parameters = vegetation_params, variables = vars)
    return model
end

function update_canopy_parameters!(
    interception_model::I,
) where {I<:AbstractInterceptionModel}
    (; leaf_area_index, swood, kext, sl, canopygapfraction, cmax) =
        interception_model.parameters.veg_param_set

    n = length(leaf_area_index)
    threaded_foreach(1:n, basesize = 1000) do i
        cmax[i] = sl[i] * leaf_area_index[i] + swood[i]
        canopygapfraction[i] = exp(-kext[i] * leaf_area_index[i])
    end
end

function update(
    interception_model::GashInterceptionModel,
    precipitation,
    potential_evaporation,
)
    (; leaf_area_index, canopygapfraction, cmax, kc) =
        interception_model.parameters.veg_param_set
    (; canopy_potevap, throughfall, interception, stemflow, canopystorage) =
        interception_model.variables
    e_r = interception_model.parameters.e_r
    n = length(precipitation)
    if !isnothing(leaf_area_index)
        update_canopy_parameters!(interception_model)
        threaded_foreach(1:n, basesize = 1000) do i
            canopyfraction = 1.0 - canopygapfraction[i]
            ewet = canopyfraction * potential_evaporation[i] * kc[i]
            e_r[i] =
                precipitation[i] > 0.0 ?
                min(0.25, ewet / max(0.0001, canopyfraction * precipitation[i])) : 0.0
        end
    end
    threaded_foreach(1:n, basesize = 1000) do i
        canopy_potevap[i] = kc[i] * potential_evaporation[i] * (1.0 - canopygapfraction[i])
        throughfall[i], interception[i], stemflow[i], canopystorage[i] =
            rainfall_interception_gash(
                cmax[i],
                e_r[i],
                canopygapfraction[i],
                precipitation[i],
                canopystorage[i],
                canopy_potevap[i],
            )
    end
end

function update(interception::RutterInterceptionModel) end
