
@get_units @with_kw struct GlacierModelVars{T}
    # Water within the glacier [mm]
    glacierstore::Vector{T} | "mm"
    # Glacier melt [mm Δt⁻¹]  
    glaciermelt::Vector{T}
end

function glacier_model_vars(glacierstore, n)
    vars = GlacierModelVars(; glacierstore = glacierstore, glaciermelt = fill(mv, n))
    return vars
end

@get_units @with_kw struct SnowStateBC{T}
    # Snow storage [mm]
    snow::Vector{T} | "mm"
end

function glacier_model_bc(snow)
    bc = SnowStateBC(; snow = snow)
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
    glacierstore = ncread(
        nc,
        config,
        "vertical.glacierstore";
        sel = inds,
        defaults = 5500.0,
        type = Float,
        fill = 0.0,
    )
    vars = glacier_model_vars(glacierstore, n)
    model =
        GlacierHbvModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

function update(model::GlacierHbvModel, atmospheric_forcing::AtmosphericForcing)
    (; temperature) = atmospheric_forcing
    (; glacierstore, glaciermelt) = model.variables
    (; snow) = model.boundary_conditions
    (; g_tt, g_cfmax, g_sifrac, glacierfrac, max_snow_to_glacier) = model.parameters

    n = length(temperature)

    threaded_foreach(1:n; basesize = 1000) do i
        snow[i], _, glacierstore[i], glaciermelt[i] = glacier_hbv(
            glacierfrac[i],
            glacierstore[i],
            snow[i],
            temperature[i],
            g_tt[i],
            g_cfmax[i],
            g_sifrac[i],
            max_snow_to_glacier,
        )
    end
end

function update(model::NoGlacierModel, atmospheric_forcing::AtmosphericForcing)
    return nothing
end

glaciermelt(model::NoGlacierModel) = 0.0
glaciermelt(model::AbstractGlacierModel) =
    @. model.variables.glaciermelt * model.variables.glacierfrac