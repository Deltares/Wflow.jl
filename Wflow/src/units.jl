const PowersType = SVector{2, Rational{Int}}

"""
Store a unit as a product of powers, for instance:
m/s -> Unit(; m = 1, s = -1)

Positive and negative powers can be split to support e.g.
m²m⁻² -> Unit(; m = (2, 2))

The `absolute_temperature` flag indicates whether °C => K
requires a +273.15 shift.

Adding support for a new unit is easy:
- Add a field to the `Unit` struct
- Specify how it translates to SI standard units in `to_SI_data`
- If the symbol for the unit is different from its field name in `Unit`,
    add it to `UnitStrings`
"""
struct Unit
    # Temperature
    absolute_temperature::Bool # Expected to be first field!
    K::PowersType # Kelvin, SI standard
    degC::PowersType # degree Celcius
    # Time
    s::PowersType # second, SI standard
    ms::PowersType # millisecond
    min::PowersType # minute
    h::PowersType # hour
    d::PowersType # day
    dt::PowersType # time step
    # Length
    m::PowersType # meter, SI standard
    cm::PowersType # centimeter
    mm::PowersType # millimeter
    μm::PowersType # micrometer
    # Volume
    L::PowersType # liter
    # Mass
    kg::PowersType # kilogram, SI standard
    g::PowersType # gram
    t::PowersType # tonne
    # Energy
    J::PowersType # joule
    # Fraction
    percentage::PowersType # percentage, converted to unitless fraction in the SI standard
    ppm::PowersType # parts per million, converted to unitless fraction in the SI standard
    # Factor for converting a value in this unit to the above mentioned standard SI units (apart from dt)
    to_SI_factor_without_dt::Float64 # Expected to be last field!
end

const Units = fieldnames(Unit)[2:(end - 1)]
const N_UNITS = length(Units)
const STANDARD_UNITS = [:K, :s, :m, :kg]

function Unit(absolute_temperature, powers_all...)
    to_SI_factor_without_dt = 1.0
    for (i, powers) in enumerate(powers_all)
        unit = Units[i]
        (unit == :dt) && continue
        @assert all(≥(0), powers) "Expected non-negative input, got $powers for $unit."
        net_power = powers[2] - powers[1]
        if (unit ∉ STANDARD_UNITS) && !iszero(net_power)
            factor = to_SI_data[i].factor
            to_SI_factor_without_dt *= factor^net_power
        end
    end
    return Unit(absolute_temperature, powers_all..., to_SI_factor_without_dt)
end

# Get a tuple of all power values
@generated function get_powers_tuple(unit::Unit)
    exprs = [:(getfield(unit, $(i + 1))) for i in 1:N_UNITS]  # +1 to skip absolute_temperature
    return :(tuple($(exprs...)))
end

function Unit(; absolute_temperature = false, kwargs...)
    powers = ntuple(
        i -> begin
            unit = Units[i]
            powers = return if unit in keys(kwargs)
                val = kwargs[unit]

                if val isa Number
                    val > 0 ? (0, val) : (-val, 0)
                else
                    val
                end
            else
                (0 // 1, 0 // 1)
            end
            PowersType(powers)
        end,
        N_UNITS,
    )
    for unit in keys(kwargs)
        (unit ∉ Units) && argument_error("Unrecognized unit $unit.")
    end
    return Unit(absolute_temperature, powers...)
end

const to_SI_data = @NamedTuple{factor::Float64, unit_SI::Unit}[
    (factor = 1.0, unit_SI = Unit(; K = 1)), # K
    (factor = 1.0, unit_SI = Unit(; K = 1)), # degC
    (factor = 1.0, unit_SI = Unit(; s = 1)), # s
    (factor = 1e-3, unit_SI = Unit(; s = 1)), # ms
    (factor = 60.0, unit_SI = Unit(; s = 1)), # min
    (factor = 3600, unit_SI = Unit(; s = 1)), # h
    (factor = 86400.0, unit_SI = Unit(; s = 1)), # d
    (factor = NaN, unit_SI = Unit(; s = 1)), # dt
    (factor = 1.0, unit_SI = Unit(; m = 1)), # m
    (factor = 1e-2, unit_SI = Unit(; m = 1)), # cm
    (factor = 1e-3, unit_SI = Unit(; m = 1)), # mm
    (factor = 1e-6, unit_SI = Unit(; m = 1)), # μm
    (factor = 1e-3, unit_SI = Unit(; m = 3)), # L
    (factor = 1.0, unit_SI = Unit(; kg = 1)), # kg
    (factor = 1e-3, unit_SI = Unit(; kg = 1)), # g
    (factor = 1e3, unit_SI = Unit(; kg = 1)), # t
    (factor = 1.0, unit_SI = Unit(; kg = 1, m = 2, s = -2)), # J
    (factor = 1e-2, unit_SI = Unit()), # percentage
    (factor = 1e-6, unit_SI = Unit()), # ppm
]

# Predefined units used within the code
const EMPTY_UNIT = Unit()
const ABSOLUTE_DEGREES = Unit(; degC = 1, absolute_temperature = true)
const MM_PER_MIN = Unit(; mm = 1, min = -1)
const MM_PER_DAY = Unit(; mm = 1, d = -1)
const MM_PER_HOUR = Unit(; mm = 1, h = -1)
const HOUR = Unit(; h = 1)
const MM = Unit(; mm = 1)
const KG_PER_MIN = Unit(; kg = 1, min = -1)
const M3_PER_MIN = Unit(; m = 3, min = -1)
const TON_PER_DT = Unit(; t = 1, dt = -1)
const GRAM_PER_L = Unit(; g = 1, L = 1)
const CM_PER_S = Unit(; cm = 1, s = -1)
const TON_PER_M3 = Unit(; t = 1, m = -3)
const PPM = Unit(; ppm = 1)
const M3_PER_DAY = Unit(; m = 3, d = -1)

function Base.:*(u1::Unit, u2::Unit)
    absolute_temperature_new = u1.absolute_temperature || u2.absolute_temperature
    powers1 = get_powers_tuple(u1)
    powers2 = get_powers_tuple(u2)
    powers_new = ntuple(i -> powers1[i] + powers2[i], N_UNITS)
    return Unit(absolute_temperature_new, powers_new...)
end

function Base.:^(u::Unit, n::Rational{Int})
    powers = get_powers_tuple(u)
    powers_new = if (n > 0)
        ntuple(i -> powers[i] * n, N_UNITS)
    else
        ntuple(i -> reverse(powers[i] * -n), N_UNITS)
    end
    return Unit(u.absolute_temperature, powers_new...)
end

# Unit strings that are not the same as their field name in
# the Unit struct
const UnitStrings = Dict{Symbol, String}(
    :degC => "°C",
    :percentage => "%",
    :d => "day",
    :dt => "Δt",
    :t => "ton",
)

const SUPERSCRIPT_NUMBERS = (
    '0' => '⁰',
    '1' => '¹',
    '2' => '²',
    '3' => '³',
    '4' => '⁴',
    '5' => '⁵',
    '6' => '⁶',
    '7' => '⁷',
    '8' => '⁸',
    '9' => '⁹',
    '-' => '⁻',
)

function power_string(power::Rational{Int}, BMI_standard::Bool)
    (; num, den) = power
    return if isinteger(power)
        n = string(Int(power))
        BMI_standard ? n : replace(n, SUPERSCRIPT_NUMBERS...)
    else
        if BMI_standard
            "$num/$den"
        else
            num_str = replace(string(num), SUPERSCRIPT_NUMBERS...)
            den_str = replace(string(den), SUPERSCRIPT_NUMBERS...)
            "$(num_str)ᐟ$den_str"
        end
    end
end

"""
Represent the Unit as a string,
following the BMI standard if `BMI_standard = true`:
https://bmi.csdms.io/en/stable/bmi.var_funcs.html#get-var-units
"""
function to_string(unit::Unit; BMI_standard = false)
    n_units = length(Units)
    symbols = ntuple(i -> get(UnitStrings, Units[i], String(Units[i])), Val(n_units))
    powers = ntuple(i -> getfield(unit, Units[i]), Val(n_units))
    out = String[]

    # Positive powers
    for (symbol, powers_) in zip(symbols, powers)
        if (symbol == :dt) && !BMI_standard
            symbol = :Δt
        end
        power = powers_[2]
        if !iszero(power)
            term = isone(power) ? symbol : "$symbol$(power_string(power, BMI_standard))"
            push!(out, term)
        end
    end

    # Negative powers
    for (symbol, powers_) in zip(symbols, powers)
        power = powers_[1]
        if !iszero(power)
            push!(out, "$symbol$(power_string(-power, BMI_standard))")
        end
    end

    out = join(out, " ")
    # Dash if unitless
    return isempty(out) ? "-" : out
end

Base.string(unit::Unit) = to_string(unit)
Base.show(io::IO, unit::Unit) = print(io, to_string(unit))

"""
Convert a general unit to an SI standard unit
"""
function to_SI(unit::Unit)
    powers_all = get_powers_tuple(unit)
    unit_new = Unit(
        unit.absolute_temperature,
        ntuple(Returns(SVector(0 // 1, 0 // 1)), N_UNITS)...,
    )

    for (i, powers) in enumerate(powers_all)
        net_power = powers[2] - powers[1]
        if !iszero(net_power)
            unit_new *= to_SI_data[i].unit_SI^net_power
        end
    end

    return unit_new
end

function to_SI_factor(unit::Unit; dt_val::Union{Nothing, Number} = nothing)
    powers_dt = unit.dt
    net_power_dt = powers_dt[2] - powers_dt[1]
    factor = unit.to_SI_factor_without_dt
    if !iszero(net_power_dt)
        isnothing(dt_val) &&
            throw(error("Unit $unit depends on dt but dt_val was not supplied."))
        factor *= dt_val^net_power_dt
    end
    return factor
end

"""
Convert the given value of the given unit to the value of the corresponding standard SI unit
"""
function to_SI(x::AbstractFloat, unit::Unit; dt_val::Union{Nothing, Number} = nothing)
    return if unit == ABSOLUTE_DEGREES
        # Special case for absolute temperatures in °C
        x + 273.15
    else
        x * to_SI_factor(unit; dt_val)
    end
end

# Fallback method for non-Floats, e.g. nothing, missing, Bool, Int
to_SI(x, unit::Unit; kwargs...) = x

"""
Convert an array of values to the values in the corresponding SI unit in-place
"""
function to_SI!(x::AbstractArray, unit::Unit; dt_val::Union{Nothing, Number} = nothing)
    unit_ref = Ref(unit)
    @. x = to_SI(x, unit_ref; dt_val)
    return x
end

"""
Convert the given value of the SI equivalent of the given unit to the value in the given unit
"""
function from_SI(x::AbstractFloat, unit::Unit; dt_val::Union{Nothing, Number} = nothing)
    return if unit == ABSOLUTE_DEGREES
        # Special case for absolute temperatures in °C
        x - 273.15
    else
        x / to_SI_factor(unit; dt_val)
    end
end

# Fallback method for non-Floats, e.g. nothing, missing, Bool, Int
from_SI(x, unit::Unit; kwargs...) = x

"""
Convert the given values of the SI equivalent of the given unit to the values in the given unit in-place
"""
function from_SI!(
    x::Union{AbstractArray},
    unit::Unit;
    dt_val::Union{Nothing, Number} = nothing,
)
    unit_ref = Ref(unit)
    @. x = from_SI(x, unit_ref; dt_val)
end

"""
Get the unit by checking all `*_standard_name_map` dictionaries
and make sure that there is no ambiguity.
"""
function get_unit(variable_name::AbstractString)::Unit
    unit = nothing

    for standard_name_map in (
        sbm_standard_name_map,
        sediment_standard_name_map,
        domain_standard_name_map,
        routing_standard_name_map,
    )
        nt = get(standard_name_map, variable_name, nothing)
        if !isnothing(nt)
            if isnothing(unit)
                unit = nt.unit
            else
                @assert nt.unit == unit "Unit ambiguity found for variable name `$variable_name`: $unit, $(nt.unit)."
            end
        end
    end
    return isnothing(unit) ? Unit() : unit
end
