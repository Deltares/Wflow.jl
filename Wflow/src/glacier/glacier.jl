abstract type AbstractGlacierModel end

"Struct for storing glacier model variables"
@with_kw struct GlacierVariables
    n::Int
    # Water within the glacier [mm]
    glacier_store::Vector{Float64}
    # Glacier melt [mm Δt⁻¹]
    glacier_melt::Vector{Float64} = fill(MISSING_VALUE, n)
end

"Initialize glacier model variables"
function GlacierVariables(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
)
    glacier_store = ncread(
        dataset,
        config,
        "glacier_ice__initial_leq_depth",
        LandHydrologySBM;
        sel = indices,
    )
    n = length(glacier_store)
    vars = GlacierVariables(; n, glacier_store)
    return vars
end

"Struct for storing boundary condition (snow storage from a snow model) of a glacier model"
@with_kw struct SnowStateBC
    # Snow storage [mm]
    snow_storage::Vector{Float64}
end

"Struct for storing glacier HBV model parameters"
@with_kw struct GlacierHbvParameters
    # Threshold temperature for glacier melt [ᵒC]
    g_ttm::Vector{Float64}
    # Degree-day factor [mm ᵒC⁻¹ Δt⁻¹] for glacier
    g_cfmax::Vector{Float64}
    # Fraction of the snowpack on top of the glacier converted into ice [Δt⁻¹]
    g_sifrac::Vector{Float64}
    # Fraction covered by a glacier [-]
    glacier_frac::Vector{Float64}
    # Maximum snow to glacier conversion rate [mm Δt⁻¹]
    max_snow_to_glacier::Float64
end

"Glacier HBV model"
@with_kw struct GlacierHbvModel <: AbstractGlacierModel
    boundary_conditions::SnowStateBC
    parameters::GlacierHbvParameters
    variables::GlacierVariables
end

struct NoGlacierModel <: AbstractGlacierModel
    n::Int
end

"Initialize glacier HBV model parameters"
function GlacierHbvParameters(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
)
    g_ttm = ncread(
        dataset,
        config,
        "glacier_ice__melting_temperature_threshold",
        LandHydrologySBM;
        sel = indices,
    )
    g_cfmax =
        ncread(
            dataset,
            config,
            "glacier_ice__degree_day_coefficient",
            LandHydrologySBM;
            sel = indices,
        ) .* (dt / BASETIMESTEP)
    g_sifrac =
        ncread(
            dataset,
            config,
            "glacier_firn_accumulation__snowpack_dry_snow_leq_depth_fraction",
            LandHydrologySBM;
            sel = indices,
        ) .* (dt / BASETIMESTEP)
    glacier_frac = ncread(
        dataset,
        config,
        "glacier_surface__area_fraction",
        LandHydrologySBM;
        sel = indices,
    )
    max_snow_to_glacier = 8.0 * (dt / BASETIMESTEP)
    glacier_hbv_params =
        GlacierHbvParameters(; g_ttm, g_cfmax, g_sifrac, glacier_frac, max_snow_to_glacier)
    return glacier_hbv_params
end

"Initialize glacier HBV model"
function GlacierHbvModel(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second,
    boundary_conditions::SnowStateBC,
)
    parameters = GlacierHbvParameters(dataset, config, indices, dt)
    variables = GlacierVariables(dataset, config, indices)
    model = GlacierHbvModel(; boundary_conditions, parameters, variables)
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
get_glacier_melt(model::NoGlacierModel) = Zeros(model.n)
get_glacier_melt(model::AbstractGlacierModel) = model.variables.glacier_melt
get_glacier_fraction(model::NoGlacierModel) = Zeros(model.n)
get_glacier_fraction(model::AbstractGlacierModel) = model.parameters.glacier_frac
get_glacier_store(model::NoGlacierModel) = Zeros(model.n)
get_glacier_store(model::AbstractGlacierModel) = model.variables.glacier_store
