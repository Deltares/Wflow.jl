
struct NoMassBalance <: AbstractMassBalance end

@with_kw struct MassBalance
    n::Int
    storage_prev::Vector{Float64} = fill(MISSING_VALUE, n)
    error::Vector{Float64} = fill(MISSING_VALUE, n)
    relative_error::Vector{Float64} = fill(MISSING_VALUE, n)
end

struct HydrologicalMassBalance <: AbstractMassBalance
    land::MassBalance
end

function HydrologicalMassBalance(n::Int, modelsettings::NamedTuple)
    do_water_mass_balance = modelsettings.water_mass_balance
    if do_water_mass_balance
        HydrologicalMassBalance(MassBalance(; n))
    else
        NoMassBalance()
    end
end

function total_storage(model::LandHydrologySBM, i::Int)
    (; total_soilwater_storage) = model.soil.variables
    (; canopy_storage) = model.interception.variables
    (; snow, glacier, demand) = model

    snow_storage = get_snow_storage(snow)[i] + get_snow_water(snow)[i]
    glacier_storage = get_glacier_store(glacier)[i] * get_glacier_fraction(glacier)[i]
    paddy_storage = get_water_depth(demand.paddy)[i]

    total_storage =
        total_soilwater_storage[i] +
        canopy_storage[i] +
        snow_storage +
        glacier_storage +
        paddy_storage

    return total_storage
end

function mass_storage_prev!(model, ::HydrologicalMassBalance)
    (; storage_prev) = model.mass_balance.land
    for i in eachindex(storage_prev)
        storage_prev[i] = total_storage(model.land, i)
    end
end

function mass_storage_prev!(model, ::NoMassBalance)
    return nothing
end

function compute_land_hydrology_balance!(
    model::AbstractModel{<:Union{SbmModel, SbmGwfModel}},
)
    (; storage_prev, error, relative_error) = model.mass_balance.land
    (; allocation, snow) = model.land
    (; net_runoff, actevap, actleakage) = model.land.soil.variables
    (; net_runoff_river) = model.land.runoff.variables
    (; precipitation) = model.land.atmospheric_forcing
    (; area) = model.domain.land.parameters
    (; subsurface_flow) = model.routing

    for i in eachindex(storage_prev)
        subsurface_flow_in = get_inflow(subsurface_flow)[i]
        subsurface_flux_in =
            1000.0 * subsurface_flow_in / area[i] * (model.clock.dt / BASETIMESTEP)
        total_input =
            subsurface_flux_in +
            precipitation[i] +
            get_irrigation_allocated(allocation)[i] +
            get_snow_in(snow)[i]

        subsurface_flow_out = get_outflow(subsurface_flow)[i]
        subsurface_flux_out =
            1000.0 * subsurface_flow_out / area[i] * (model.clock.dt / BASETIMESTEP)
        vertical_flux_out =
            net_runoff[i] +
            actevap[i] +
            net_runoff_river[i] +
            actleakage[i] +
            get_groundwater_abstraction_flux(allocation)[i]
        total_output = subsurface_flux_out + vertical_flux_out + get_snow_in(snow)[i]

        storage = total_storage(model.land, i)

        error[i] = (total_input - total_output - (storage - storage_prev[i]))
        relative_error[i] = error[i] / ((total_input + total_output) / 2)
    end
    return nothing
end

function compute_mass_balance!(model, ::HydrologicalMassBalance)
    compute_land_hydrology_balance!(model)
    return nothing
end

function compute_mass_balance!(model, ::NoMassBalance)
    return nothing
end
