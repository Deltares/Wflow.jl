
struct NoMassBalance <: AbstractMassBalance end

@with_kw struct MassBalance
    n::Int
    storage_prev::Vector{Float64} = fill(MISSING_VALUE, n)
    error::Vector{Float64} = fill(MISSING_VALUE, n)
    relative_error::Vector{Float64} = fill(MISSING_VALUE, n)
end

struct HydrologicalMassBalance <: AbstractMassBalance
    soil_water::MassBalance
end

function HydrologicalMassBalance(n::Int, modelsettings::NamedTuple)
    do_water_mass_balance = modelsettings.water_mass_balance
    if do_water_mass_balance
        HydrologicalMassBalance(MassBalance(; n))
    else
        NoMassBalance()
    end
end

function mass_storage_prev!(model, ::HydrologicalMassBalance)
    (; soil_water) = model.mass_balance
    soil_water.storage_prev .= model.land.soil.variables.total_soilwater_storage
end

function mass_storage_prev!(model, ::NoMassBalance)
    return nothing
end

function compute_soil_water_balance!(model::AbstractModel{<:SbmModel})
    (; storage_prev, error, relative_error) = model.mass_balance.soil_water
    (;
        actinfilt,
        exfiltsatwater,
        exfiltustore,
        transpiration,
        soilevap,
        actleakage,
        total_soilwater_storage,
    ) = model.land.soil.variables
    (; area) = model.domain.land.parameters

    (; ssf, ssfin) = model.routing.subsurface_flow.variables

    for i in eachindex(storage_prev)
        ssf_flux_in = 1000.0 * ssfin[i] / area[i] * (model.clock.dt / BASETIMESTEP)
        total_input = ssf_flux_in + actinfilt[i]

        ssf_flux_out = 1000.0 * ssf[i] / area[i] * (model.clock.dt / BASETIMESTEP)
        total_output =
            ssf_flux_out +
            exfiltsatwater[i] +
            exfiltustore[i] +
            soilevap[i] +
            actleakage[i] +
            transpiration[i]
        storage = total_soilwater_storage[i]
        error[i] = (total_input - total_output - (storage - storage_prev[i]))
        relative_error[i] = error[i] / ((total_input + total_output) / 2)
    end
    return nothing
end

function compute_mass_balance!(model, ::HydrologicalMassBalance)
    compute_soil_water_balance!(model)
    return nothing
end

function compute_mass_balance!(model, ::NoMassBalance)
    return nothing
end
