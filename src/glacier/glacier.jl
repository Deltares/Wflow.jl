abstract type AbstractGlacierModel{T} end

"Struct for storing glacier model variables"
@get_units @grid_loc @with_kw struct GlacierVariables{T}
    # Water within the glacier [mm]
    glacier_store::Vector{T} | "mm"
    # Glacier melt [mm Δt⁻¹]  
    glacier_melt::Vector{T}
end

"Initialize glacier model variables"
function GlacierVariables(nc, config, inds)
    glacier_store = ncread(
        nc,
        config,
        "vertical.glacier.variables.glacier_store";
        sel = inds,
        defaults = 5500.0,
        type = Float,
        fill = 0.0,
    )
    n = length(glacier_store)
    vars = GlacierVariables(; glacier_store = glacier_store, glacier_melt = fill(mv, n))
    return vars
end

"Struct for storing boundary condition (snow storage from a snow model) of a glacier model"
@get_units @grid_loc @with_kw struct SnowStateBC{T}
    # Snow storage [mm]
    snow_storage::Vector{T} | "mm"
end

"Struct for storing glacier HBV model parameters"
@get_units @grid_loc @with_kw struct GlacierHbvParameters{T}
    # Threshold temperature for glacier melt [ᵒC]
    g_ttm::Vector{T} | "ᵒC"
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹] for glacier
    g_cfmax::Vector{T} | "mm ᵒC-1 dt-1"
    # Fraction of the snowpack on top of the glacier converted into ice [Δt⁻¹]
    g_sifrac::Vector{T} | "dt-1"
    # Fraction covered by a glacier [-]
    glacier_frac::Vector{T} | "-"
    # Maximum snow to glacier conversion rate [mm Δt⁻¹]
    max_snow_to_glacier::T
end

"Glacier HBV model"
@with_kw struct GlacierHbvModel{T} <: AbstractGlacierModel{T}
    boundary_conditions::SnowStateBC{T}
    parameters::GlacierHbvParameters{T}
    variables::GlacierVariables{T}
end

struct NoGlacierModel{T} <: AbstractGlacierModel{T} end

"Initialize glacier HBV model parameters"
function GlacierHbvParameters(nc, config, inds, dt)
    g_ttm = ncread(
        nc,
        config,
        "vertical.glacier.parameters.g_ttm";
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    g_cfmax =
        ncread(
            nc,
            config,
            "vertical.glacier.parameters.g_cfmax";
            sel = inds,
            defaults = 3.0,
            type = Float,
            fill = 0.0,
        ) .* (dt / basetimestep)
    g_sifrac =
        ncread(
            nc,
            config,
            "vertical.glacier.parameters.g_sifrac";
            sel = inds,
            defaults = 0.001,
            type = Float,
            fill = 0.0,
        ) .* (dt / basetimestep)
    glacier_frac = ncread(
        nc,
        config,
        "vertical.glacier.parameters.glacier_frac";
        sel = inds,
        defaults = 0.0,
        type = Float,
        fill = 0.0,
    )
    max_snow_to_glacier = 8.0 * (dt / basetimestep)
    glacier_hbv_params =
        GlacierHbvParameters(; g_ttm, g_cfmax, g_sifrac, glacier_frac, max_snow_to_glacier)
    return glacier_hbv_params
end

"Initialize glacier HBV model"
function GlacierHbvModel(nc, config, inds, dt, bc)
    params = GlacierHbvParameters(nc, config, inds, dt)
    vars = GlacierVariables(nc, config, inds)
    model =
        GlacierHbvModel(; boundary_conditions = bc, parameters = params, variables = vars)
    return model
end

"Update glacier HBV model for a single timestep"
function update!(model::GlacierHbvModel, atmospheric_forcing::AtmosphericForcing)
    (; temperature) = atmospheric_forcing
    (; glacier_store, glacier_melt) = model.variables
    (; snow_storage) = model.boundary_conditions
    (; g_ttm, g_cfmax, g_sifrac, glacier_frac, max_snow_to_glacier) = model.parameters

    n = length(temperature)

    threaded_foreach(1:n; basesize = 1000) do i
        snow_storage[i], _, glacier_store[i], glacier_melt[i] = glacier_hbv(
            glacier_frac[i],
            glacier_store[i],
            snow_storage[i],
            temperature[i],
            g_ttm[i],
            g_cfmax[i],
            g_sifrac[i],
            max_snow_to_glacier,
        )
    end
    return nothing
end

function update!(model::NoGlacierModel, atmospheric_forcing::AtmosphericForcing)
    return nothing
end

# wrapper methods
get_glacier_melt(model::NoGlacierModel) = 0.0
get_glacier_melt(model::AbstractGlacierModel) = model.variables.glacier_melt
get_glacier_fraction(model::NoGlacierModel) = 0.0
get_glacier_fraction(model::AbstractGlacierModel) = model.parameters.glacier_frac
get_glacier_store(model::NoGlacierModel) = 0.0
get_glacier_store(model::AbstractGlacierModel) = model.variables.glacier_store