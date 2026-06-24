abstract type AbstractGlacierModel end

"Struct for storing glacier model variables"
@with_kw struct GlacierVariables
    n::Int
    # Water within the glacier [m]
    glacier_store::Vector{Float64}
    # Glacier melt [m s⁻¹]
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
    # Snow storage [m]
    snow_storage::Vector{Float64}
end

"Struct for storing glacier HBV model parameters"
@with_kw struct GlacierHbvParameters
    # Threshold temperature for glacier melt [K]
    temperature_threshold_melt::Vector{Float64}
    # Degree-day factor [m K⁻¹ s⁻¹] for glacier
    degree_day_factor::Vector{Float64}
    # Fraction of the snowpack on top of the glacier converted into ice [s⁻¹]
    snow_to_ice_fraction::Vector{Float64}
    # Fraction covered by a glacier [-]
    glacier_fraction::Vector{Float64}
    # Maximum snow to glacier conversion rate [m s⁻¹]
    maximum_snow_to_ice_rate::Float64
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
    temperature_threshold_melt = ncread(
        dataset,
        config,
        "glacier_ice__melting_temperature_threshold",
        LandHydrologySBM;
        sel = indices,
    )
    degree_day_factor = ncread(
        dataset,
        config,
        "glacier_ice__degree_day_coefficient",
        LandHydrologySBM;
        sel = indices,
    )
    snow_to_ice_fraction = ncread(
        dataset,
        config,
        "glacier_firn_accumulation__snowpack_dry_snow_leq_depth_fraction",
        LandHydrologySBM;
        sel = indices,
    )
    glacier_fraction = ncread(
        dataset,
        config,
        "glacier_surface__area_fraction",
        LandHydrologySBM;
        sel = indices,
    )
    # 8 mm day⁻¹
    maximum_snow_to_ice_rate = to_SI(8.0, MM_PER_DAY)
    return GlacierHbvParameters(;
        temperature_threshold_melt,
        degree_day_factor,
        snow_to_ice_fraction,
        glacier_fraction,
        maximum_snow_to_ice_rate,
    )
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
    glacier_model = GlacierHbvModel(; boundary_conditions, parameters, variables)
    return glacier_model
end

"Update glacier HBV model for a single timestep"
function update_glacier_model!(
    glacier_model::GlacierHbvModel,
    atmospheric_forcing::AtmosphericForcing,
    dt::Float64,
)
    (; temperature) = atmospheric_forcing
    (; glacier_store, glacier_melt) = glacier_model.variables
    (; snow_storage) = glacier_model.boundary_conditions
    (;
        temperature_threshold_melt,
        degree_day_factor,
        snow_to_ice_fraction,
        glacier_fraction,
        maximum_snow_to_ice_rate,
    ) = glacier_model.parameters

    n = length(temperature)

    threaded_foreach(1:n; basesize = 1000) do idx
        snow_storage[idx], _, glacier_store[idx], glacier_melt[idx] = glacier_hbv(
            glacier_fraction[idx],
            glacier_store[idx],
            snow_storage[idx],
            temperature[idx],
            temperature_threshold_melt[idx],
            degree_day_factor[idx],
            snow_to_ice_fraction[idx],
            maximum_snow_to_ice_rate,
            dt,
        )
    end
    return nothing
end

function update_glacier_model!(
    glacier_model::NoGlacierModel,
    atmospheric_forcing::AtmosphericForcing,
    dt::Float64,
)
    return nothing
end

# wrapper methods
get_glacier_melt(glacier_model::NoGlacierModel) = Zeros(glacier_model.n)
get_glacier_melt(glacier_model::AbstractGlacierModel) = glacier_model.variables.glacier_melt
get_glacier_fraction(glacier_model::NoGlacierModel) = Zeros(glacier_model.n)
get_glacier_fraction(glacier_model::AbstractGlacierModel) =
    glacier_model.parameters.glacier_fraction
get_glacier_store(glacier_model::NoGlacierModel) = Zeros(glacier_model.n)
get_glacier_store(glacier_model::AbstractGlacierModel) =
    glacier_model.variables.glacier_store
