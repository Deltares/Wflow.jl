"""
Store a unit as a product of powers, for instance:
m/s -> Unit(; m = 1, s = -1)

Positive and negative powers can be split to support e.g.
m²m⁻² -> Unit(; m = (2, 2))
"""
struct Unit
    # Time
    s::SVector{2, Int} # second, SI standard
    d::SVector{2, Int} # day
    dt::SVector{2, Int} # time step
    # Length
    m::SVector{2, Int} # meter, SI standard
    mm::SVector{2, Int} # millimeter
    # Temperature
    K::SVector{2, Int} # Kelvin, SI standard
    degC::SVector{2, Int} # Degree Celcius
    # Mass
    kg::SVector{2, Int} # kilogram, SI standard
    g::SVector{2, Int} # gram
    t::SVector{2, Int} # tonne
    # Fraction
    percentage::SVector{2, Int} # percentage, converted to unitless fraction
end

function Unit(; kwargs...)
    out = Dict(unit => SVector(0, 0) for unit in fieldnames(Unit))
    for (unit, val) in pairs(kwargs)
        @assert unit ∈ fieldnames(Unit) "Unrecognized unit $unit."
        powers = if val isa Int
            val > 0 ? (0, val) : (-val, 0)
        elseif (val isa Union{Tuple{Int, Int}, SVector{2, Int}} && val[1] ≥ 0 && val[2] ≥ 0)
            SVector{2, Int}(val)
        else
            throw(ArgumentError("Invalid input for $unit: $val."))
        end
        out[unit] = powers
    end
    return Unit(
        ntuple(
            i -> SVector{2, Int}(out[fieldnames(Unit)[i]]),
            Val(length(fieldnames(Unit))),
        )...,
    )
end

"""
Represent the Unit as a string,
following the BMI standard:
https://bmi.csdms.io/en/stable/bmi.var_funcs.html#get-var-units
"""
function Base.String(unit::Unit)
    units = fieldnames(Unit)
    n_units = length(units)
    symbols = ntuple(i -> begin
        symb = units[i]
        if symb == :degC
            "°C"
        elseif symb == :percentage
            "%"
        else
            String(symb)
        end
    end, Val(n_units))
    powers = ntuple(i -> getfield(unit, units[i]), Val(n_units))
    out = String[]

    # Positive powers
    for (symbol, powers_) in zip(symbols, powers)
        power = powers_[2]
        if !iszero(power)
            term = isone(power) ? symbol : "$symbol$power"
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
            error("Unit depends on dt but `dt_val` was not passed.")
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
        m = unit.m + unit.mm,
        K = unit.K + unit.degC,
        kg = unit.kg + unit.g + unit.t,
    )
    @assert to_SI_factor(unit_SI) |> isone
    return unit_SI
end

"""
Convert the given value of the given unit to the value of the corresponding standard SI unit
"""
function to_SI(x::Number, unit::Unit; dt_val::Union{Nothing, Float64} = nothing)
    return if unit == Unit(; degC = 1)
        # Special case: if the unit is "°C" the temperature is assumed to be absolute
        x + 274.15
    else
        x * to_SI_factor(unit; dt_val)
    end
end
