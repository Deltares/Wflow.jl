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
    μm::SVector{2, Rational{Int}} # micrometer
    # Mass
    kg::SVector{2, Rational{Int}} # kilogram, SI standard
    g::SVector{2, Rational{Int}} # gram
    t::SVector{2, Rational{Int}} # tonne
    # Energy
    J::SVector{2, Rational{Int}} # Joules
    # Fraction
    percentage::SVector{2, Rational{Int}} # percentage, converted to unitless fraction in the SI standard
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

function Base.:*(u1::Unit, u2::Unit)
    # Multiply units by adding their powers
    Unit(
        u1.absolute_temperature || u2.absolute_temperature,
        ntuple(
            i -> getfield(u1, Units[i]) .+ getfield(u2, Units[i]),
            Val(length(Units)),
        )...,
    )
end

function Base.:^(u::Unit, n::Rational{Int})
    # Raise unit to a power by multiplying all powers by n
    Unit(
        u.absolute_temperature,
        ntuple(i -> getfield(u, Units[i]) .* n, Val(length(Units)))...,
    )
end

# Unit strings that are not the same as their field name in
# the Unit struct
const UnitStrings = Dict{Symbol, String}(:degC => "°C", :percentage => "%")

function power_string(power::Rational{Int}, BMI_standard::Bool)
    (; num, den) = power
    return if isinteger(power)
        n = string(Int(power))
        BMI_standard ? n : Subscripts.super(n)
    else
        if BMI_standard
            "$num/$den"
        else
            num_str = Subscripts.super(string(num))
            den_str = Subscripts.super(string(den))
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

const to_SI_data = Dict{Symbol, @NamedTuple{factor::Float64, unit_SI::Unit}}(
    :K => (factor = 1.0, unit_SI = Unit(; K = 1)),
    :degC => (factor = 1.0, unit_SI = Unit(; K = 1)),
    :s => (factor = 1.0, unit_SI = Unit(; s = 1)),
    :d => (factor = 86400.0, unit_SI = Unit(; s = 1)),
    :dt => (factor = NaN, unit_SI = Unit(; s = 1)),
    :m => (factor = 1.0, unit_SI = Unit(; m = 1)),
    :cm => (factor = 1e-2, unit_SI = Unit(; m = 1)),
    :mm => (factor = 1e-3, unit_SI = Unit(; m = 1)),
    :μm => (factor = 1e-6, unit_SI = Unit(; m = 1)),
    :kg => (factor = 1.0, unit_SI = Unit(; kg = 1)),
    :g => (factor = 1e-3, unit_SI = Unit(; kg = 1)),
    :t => (factor = 1e3, unit_SI = Unit(; kg = 1)),
    :J => (factor = 1.0, unit_SI = Unit(; kg = 1, m = 2, s = -2)),
    :percentage => (factor = 1e-2, unit_SI = Unit()),
)

"""
Obtain the conversion factor from an Unit into
SI units, e.g.:
mm/dt -> 1e-3 / dt_val
"""
function to_SI_factor(unit::Unit; dt_val::Union{Nothing, Float64} = nothing)
    # Model dependent units
    (; dt) = unit

    if isnothing(dt_val)
        if !all(iszero(dt))
            throw(ArgumentError("Unit depends on `dt` but `dt_val` was not passed."))
        else
            dt_val = 1.0
        end
    end

    factor = 1.0

    for base_unit in Units
        factor_unit = to_SI_data[base_unit].factor
        if isnan(factor_unit)
            # Model dependent conversion factors
            factor_unit = if base_unit == :dt
                dt_val
            end
        end
        powers = getfield(unit, base_unit)::SVector{2, Rational{Int}}
        factor *= factor_unit^(powers[2] - powers[1])
    end
    return factor
end

"""
Convert a general unit to an SI standard unit
"""
function to_SI(unit::Unit)
    # Start with unitless
    unit_SI = Unit()

    # Convert each unit field using to_SI_data
    for unit_field in Units
        powers = getfield(unit, unit_field)
        # Get the SI equivalent for this unit
        si_info = to_SI_data[unit_field]
        # Calculate net power (positive - negative)
        net_power = powers[2] - powers[1]
        # Raise the SI unit to the net power and multiply
        unit_contribution = si_info.unit_SI^net_power
        unit_SI = unit_SI * unit_contribution
    end

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
