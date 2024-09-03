
@get_units @with_kw struct GlacierModelVars{T}
    # Water within the glacier [mm]
    glacier_store::Vector{T} | "mm"
    # Glacier melt [mm Δt⁻¹]  
    glacier_melt::Vector{T}
end

function glacier_model_vars(glacier_store, n)
    vars = GlacierModelVars(; glacier_store = glacier_store, glacier_melt = fill(mv, n))
    return vars
end

@get_units @with_kw struct SnowStateBC{T}
    # Snow storage [mm]
    snow_storage::Vector{T} | "mm"
end

function glacier_model_bc(snow_storage)
    bc = SnowStateBC(; snow_storage = snow_storage)
    return bc
end

@get_units @with_kw struct GlacierHbvParameters{T}
    # Threshold temperature for snowfall above glacier [ᵒC]
    g_tt::Vector{T} | "ᵒC"
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹] for glacier
    g_cfmax::Vector{T} | "mm ᵒC-1 dt-1"
    # Fraction of the snowpack on top of the glacier converted into ice [Δt⁻¹]
    g_sifrac::Vector{T} | "dt-1"
    # Fraction covered by a glacier [-]
    glacierfrac::Vector{T} | "-"
    # Maximum snow to glacier conversion rate [mm Δt⁻¹]
    max_snow_to_glacier::T
end

abstract type AbstractGlacierModel end

@get_units @with_kw struct GlacierHbvModel{T} <: AbstractGlacierModel
    boundary_conditions::SnowStateBC{T} | "-"
    parameters::GlacierHbvParameters{T} | "-"
    variables::GlacierModelVars{T} | "-"
end

struct NoGlacierModel <: AbstractGlacierModel end

function initialize_glacier_hbv_params(nc, config, inds, dt)
    g_tt = ncread(
        nc,
        config,
        "vertical.g_tt";
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    g_cfmax =
        ncread(
            nc,
            config,
            "vertical.g_cfmax";
            sel = inds,
            defaults = 3.0,
            type = Float,
            fill = 0.0,
        ) .* (dt / basetimestep)
    g_sifrac =
        ncread(
            nc,
            config,
            "vertical.g_sifrac";
            sel = inds,
            defaults = 0.001,
            type = Float,
            fill = 0.0,
        ) .* (dt / basetimestep)
    glacierfrac = ncread(
        nc,
        config,
        "vertical.glacierfrac";
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    max_snow_to_glacier = 8.0 * (dt / basetimestep)
    glacier_hbv_params = GlacierHbvParameters(;
        g_tt = g_tt,
        g_cfmax = g_cfmax,
        g_sifrac = g_sifrac,
        glacierfrac = glacierfrac,
        max_snow_to_glacier = max_snow_to_glacier,
    )
    return glacier_hbv_params
end

function initialize_glacier_hbv_model(nc, config, inds, dt, bc)
    n = length(inds)
    params = initialize_glacier_hbv_params(nc, config, inds, dt)
    glacier_store = ncread(
        nc,
        config,
        "vertical.glacierstore";
        sel = inds,
        defaults = 5500.0,
        type = Float,
        fill = 0.0,
    )
    vars = glacier_model_vars(glacier_store, n)
    model =
        GlacierHbvModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

function update!(model::GlacierHbvModel, atmospheric_forcing::AtmosphericForcing)
    (; temperature) = atmospheric_forcing
    (; glacier_store, glacier_melt) = model.variables
    (; snow_storage) = model.boundary_conditions
    (; g_tt, g_cfmax, g_sifrac, glacierfrac, max_snow_to_glacier) = model.parameters

    n = length(temperature)

    threaded_foreach(1:n; basesize = 1000) do i
        snow_storage[i], _, glacier_store[i], glacier_melt[i] = glacier_hbv(
            glacierfrac[i],
            glacier_store[i],
            snow_storage[i],
            temperature[i],
            g_tt[i],
            g_cfmax[i],
            g_sifrac[i],
            max_snow_to_glacier,
        )
    end
end

function update!(model::NoGlacierModel, atmospheric_forcing::AtmosphericForcing)
    return nothing
end

get_glacier_melt(model::NoGlacierModel) = 0.0
get_glacier_melt(model::AbstractGlacierModel) =
    @. model.variables.glacier_melt * model.variables.glacier_frac
get_glacier_frac(model::NoGlacierModel) = 0.0
get_glacier_frac(model::AbstractGlacierModel) = model.variables.glacier_frac
get_glacier_store(model::NoGlacierModel) = 0.0
get_glacier_store(model::AbstractGlacierModel) =
    @. model.variables.glacier_store * model.variables.glacier_frac