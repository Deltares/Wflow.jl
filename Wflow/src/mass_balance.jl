
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

function compute_soil_water_balance!(model::AbstractModel{<:Union{SbmModel, SbmGwfModel}})
    (; storage_prev, error, relative_error) = model.mass_balance.soil_water
    (; soil) = model.land
    (; area) = model.domain.land.parameters
    (; subsurface_flow) = model.routing

    for i in eachindex(storage_prev)
        subsurface_flow_in = get_inflow(subsurface_flow, i)
        subsurface_flux_in =
            1000.0 * subsurface_flow_in / area[i] * (model.clock.dt / BASETIMESTEP)
        total_input = subsurface_flux_in + get_vertical_flux_in(soil, i)

        subsurface_flow_out = get_outflow(subsurface_flow, i)
        subsurface_flux_out =
            1000.0 * subsurface_flow_out / area[i] * (model.clock.dt / BASETIMESTEP)
        total_output = subsurface_flux_out + get_vertical_flux_out(soil, i)

        storage = soil.variables.total_soilwater_storage[i]

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
