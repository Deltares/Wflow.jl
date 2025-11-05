"""
Store a unit as a product of powers, for instance:
m/s -> Unit(; m = 1, s = -1)

Positive and negative powers can be split to support e.g.
m²m⁻² -> Unit(; m = (2, 2))

The `absolute_temperature` flag indicates whether °C => K
requires a +273.15 shift.
"""
struct Unit
    # Temperature
    absolute_temperature::Bool # Expected to be first field!
    K::SVector{2, Rational{Int}} # Kelvin, SI standard
    degC::SVector{2, Rational{Int}} # Degree Celcius
    # Time
    s::SVector{2, Rational{Int}} # second, SI standard
    d::SVector{2, Rational{Int}} # day
    dt::SVector{2, Rational{Int}} # time step
    # Length
    m::SVector{2, Rational{Int}} # meter, SI standard
    cm::SVector{2, Rational{Int}} # centimeter
    mm::SVector{2, Rational{Int}} # millimeter
    # Mass
    kg::SVector{2, Rational{Int}} # kilogram, SI standard
    g::SVector{2, Rational{Int}} # gram
    t::SVector{2, Rational{Int}} # tonne
    # Fraction
    percentage::SVector{2, Rational{Int}} # percentage, converted to unitless fraction
end

const Units = fieldnames(Unit)[2:end]

function Unit(; absolute_temperature = false, kwargs...)
    out = Dict(unit => SVector(0 // 1, 0 // 1) for unit in Units)
    for (unit, val) in pairs(kwargs)
        @assert unit ∈ fieldnames(Unit) "Unrecognized unit $unit."
        powers = if val isa Number
            val > 0 ? (0, val) : (-val, 0)
        else
            val
        end
        out[unit] = powers
    end
    return Unit(
        absolute_temperature,
        ntuple(i -> SVector{2, Rational{Int}}(out[Units[i]]), Val(length(Units)))...,
    )
end

"""
Represent the Unit as a string,
following the BMI standard:
https://bmi.csdms.io/en/stable/bmi.var_funcs.html#get-var-units
"""
function Base.String(unit::Unit)
    n_units = length(Units)
    symbols = ntuple(i -> begin
        symb = Units[i]
        if symb == :degC
            "°C"
        elseif symb == :percentage
            "%"
        else
            String(symb)
        end
    end, Val(n_units))
    powers = ntuple(i -> getfield(unit, Units[i]), Val(n_units))
    out = String[]

    # Positive powers
    for (symbol, powers_) in zip(symbols, powers)
        power = powers_[2]
        (; num, den) = power
        power_str = isinteger(power) ? string(Int(power)) : "$num/$den"
        if !iszero(power)
            term = isone(power) ? symbol : "$symbol$power_str"
            push!(out, term)
        end
    end

    # Negative powers
    for (symbol, powers_) in zip(symbols, powers)
        power = powers_[1]
        (; num, den) = power
        power_str = isinteger(power) ? string(Int(power)) : "$num/$den"
        if !iszero(power)
            push!(out, "$symbol-$power_str")
        end
    end

    out = join(out, " ")
    # Dash if unitless
    return isempty(out) ? "-" : out
end

Base.show(io::IO, unit::Unit) = print(io, String(unit))

"""
Obtain the conversion factor from an Unit into
SI units, e.g.:
mm/dt -> 1e-3 / dt_val
"""
function to_SI_factor(unit::Unit; dt_val::Union{Nothing, Float64} = nothing)
    (; d, dt, mm, g, t, percentage) = unit

    if isnothing(dt_val)
        if !all(iszero(dt))
            throw(ArgumentError("Unit depends on `dt` but `dt_val` was not passed."))
        else
            dt_val = 1.0
        end
    end

    factor = 1.0

    # Only units that contribute to the factor are relevant
    for (powers, factor_) in
        ((d, 86400.0), (dt, dt_val), (mm, 1e-3), (g, 1e-3), (t, 1e3), (percentage, 1e-2))
        factor *= factor_^(powers[2] - powers[1])
    end

    return factor
end

"""
Convert a general unit to an SI standard unit
"""
function to_SI(unit::Unit)
    unit_SI = Unit(;
        s = unit.s + unit.d + unit.dt,
        m = unit.m + +unit.cm + unit.mm,
        K = unit.K + unit.degC,
        kg = unit.kg + unit.g + unit.t,
    )
    @assert to_SI_factor(unit_SI) |> isone
    return unit_SI
end

"""
Convert the given value of the given unit to the value of the corresponding standard SI unit
"""
function to_SI(x::AbstractFloat, unit::Unit; dt_val::Union{Nothing, Float64} = nothing)
    return if unit == Unit(; degC = 1, absolute_temperature = true)
        # Special case for absolute temperatures in °C
        x + 273.15
    else
        x * to_SI_factor(unit; dt_val)
    end
end

# Fallback method for non-Floats, e.g. nothing, missing, Bool, Int
to_SI(x, unit::Unit; kwargs...) = x

"""
Convert an array of values to the values in the corresponding SI unit in-place.
"""
function to_SI!(x::AbstractArray, unit::Unit; dt_val::Union{Nothing, Float64} = nothing)
    unit_ref = Ref(unit)
    @. x = to_SI(x, unit_ref; dt_val)
    return x
end
