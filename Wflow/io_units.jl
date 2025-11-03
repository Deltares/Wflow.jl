"""
Store a unit as a product of powers, for instance:
m/s -> IOUnitPowers(; m = (0,1), s = (1, 0))

Positive and negative powers are split to support e.g.
m²m⁻² -> IOUnitPowers(; m = (2, 2))
"""
@kwdef struct IOUnitPowers
    # Time
    dt::Tuple{Int, Int} = (0, 0) # time step
    s::Tuple{Int, Int} = (0, 0) # second
    # Length
    m::Tuple{Int, Int} = (0, 0) # meter
    mm::Tuple{Int, Int} = (0, 0) # millimeter
    # Temperature
    degC::Tuple{Int, Int} = (0, 0) # Degree Celcius
    # Mass
    g::Tuple{Int, Int} = (0, 0) # gram
    t::Tuple{Int, Int} = (0, 0) # tonne
end

"""
Represent the IOUnitPowers as a string,
following the BMI standard:
https://bmi.csdms.io/en/stable/bmi.var_funcs.html#get-var-units
"""
function Base.String(io_unit_powers::IOUnitPowers)
    (; dt, s, m, mm, degC, g, t) = io_unit_powers
    symbols = ("dt", "s", "m", "mm", "°C", "g", "t")
    powers = (dt, s, m, mm, degC, g, t)
    out = String[]

    # Positive powers
    for (symbol, powers_) in zip(symbols, powers)
        power = powers_[2]
        if !iszero(power)
            term = isone(power) ? symbol : symbol * String(power)
            push!(out, term)
        end
    end

    # Negative powers
    for (symbol, powers_) in zip(symbols, powers)
        power = powers_[1]
        if !iszero(power)
            push!(out, "$symbol-$power")
        end
    end

    return join(out, " ")
end

"""
Obtain the conversion factor from an IOUnitPowers into
SI units, e.g.:
mm/dt -> 1e-3 / dt_val 
"""
function to_SI_factor(
    io_unit_powers::IOUnitPowers;
    dt_val::Union{Nothing, Float64} = nothing,
)
    (; dt, mm, g, t) = io_unit_powers

    isnothing(dt_val) &&
        !all(iszero(dt)) &&
        error("Unit depends on dt but `dt_val` was not passed.")

    factor = 1.0

    # s, m, degC are already SI units
    for (powers, factor_) in ((dt, dt_val), (mm, 1e-3), (g, 1e-3), (t, 1e3))
        factor *= factor_^(powers[2] - powers[1])
    end

    return factor
end

const input_unit_map = Dict{String, IOUnitPowers}(
    ### LandHydrologySBM
    ## AtmosphericForcing
    "atmosphere_water__precipitation_volume_flux" =>
        IOUnitPowers(; mm = (0, 1), dt = (1, 0)),
    "land_surface_water__potential_evaporation_volume_flux" =>
        IOUnitPowers(; mm = (0, 1), dt = (1, 0)),
    "atmosphere_air__temperature" => IOUnitPowers(; degC = (0, 1)),
    ## VegetationParameters
)