# Not all special functions have been implemented yet:
# https://github.com/JuliaMath/SpecialFunctions.jl/issues/19
# (MIT licensed)
#
# This expint implementation is copied from:
# https://github.com/mschauer/Bridge.jl/blob/a719e8a57aaa77fa4ca289c25978a7e080908826/src/expint.jl#L1
########################################################################
# Inlined, optimized code for the exponential integral E₁ in double precison
# by Steven G. Johnson (@stevengj), Code under MIT license.
# For more explanations, see course notes by Steven G. Johnson at
#     https://github.com/stevengj/18S096-iap17/blob/master/pset3/pset3-solutions.ipynb

# n coefficients of the Taylor series of E₁(z) + log(z), in type T:
function E1_taylor_coefficients(::Type{T}, n::Integer) where {T <: Number}
    n < 0 && throw(ArgumentError("$n ≥ 0 is required"))
    n == 0 && return T[]
    n == 1 && return T[-eulergamma]
    # iteratively compute the terms in the series, starting with k=1
    term::T = 1
    terms = T[-eulergamma, term]
    for k in 2:n
        term = -term * (k - 1) / (k * k)
        push!(terms, term)
    end
    return terms
end

# inline the Taylor expansion for a given order n, in double precision
macro E1_taylor64(z, n::Integer)
    c = E1_taylor_coefficients(Float64, n)
    taylor = :(@evalpoly zz)
    append!(taylor.args, c)
    quote
        let zz = $(esc(z))
            $taylor - log(zz)
        end
    end
end

# for numeric-literal coefficients: simplify to a ratio of two polynomials:
# return (p,q): the polynomials p(x) / q(x) corresponding to E1_cf(x, a...),
# but without the exp(-x) term
function E1_cfpoly(n::Integer, ::Type{T} = BigInt) where {T <: Real}
    q = Polynomials.Polynomial(T[1])
    p = x = Polynomials.Polynomial(T[0, 1])
    for i in n:-1:1
        p, q = x * p + (1 + i) * q, p # from cf = x + (1+i)/cf = x + (1+i)*q/p
        p, q = p + i * q, p     # from cf = 1 + i/cf = 1 + i*q/p
    end
    # do final 1/(x + inv(cf)) = 1/(x + q/p) = p/(x*p + q)
    return p, x * p + q
end

macro E1_cf64(z, n::Integer)
    p, q = E1_cfpoly(n, BigInt)
    num_expr = :(@evalpoly zz)
    append!(num_expr.args, Float64.(Polynomials.coeffs(p)))
    den_expr = :(@evalpoly zz)
    append!(den_expr.args, Float64.(Polynomials.coeffs(q)))
    quote
        let zz = $(esc(z))
            exp(-zz) * $num_expr / $den_expr
        end
    end
end

# exponential integral function E₁(z)
function expint(z::Union{Float64, Complex{Float64}})
    xSq = real(z)^2
    ySq = imag(z)^2
    if real(z) > 0 && xSq + 0.233 * ySq ≥ 7.84 # use cf expansion, ≤ 30 terms
        if (xSq ≥ 546121) & (real(z) > 0) # underflow
            return zero(z)
        elseif xSq + 0.401 * ySq ≥ 58.0 # ≤ 15 terms
            if xSq + 0.649 * ySq ≥ 540.0 # ≤ 8 terms
                xSq + ySq ≥ 4e4 && return @E1_cf64 z 4
                return @E1_cf64 z 8
            end
            return @E1_cf64 z 15
        end
        return @E1_cf64 z 30
    else # use Taylor expansion, ≤ 37 terms
        rSq = xSq + ySq
        return if rSq ≤ 0.36
            (
                if rSq ≤ 2.8e-3
                    (rSq ≤ 2e-7 ? @E1_taylor64(z, 4) : @E1_taylor64(z, 8))
                else
                    @E1_taylor64(z, 15)
                end
            )
        else
            @E1_taylor64(z, 37)
        end
    end
end
function expint(
    z::Union{T, Complex{T}, Rational{T}, Complex{Rational{T}}},
) where {T <: Integer}
    return expint(float(z))
end

######################################################################
# exponential integral Eₙ(z)

function expint(n::Integer, z)
    if n == 1
        return expint(z)
    elseif n < 1
        # backwards recurrence from E₀ = e⁻ᶻ/z
        zinv = inv(z)
        exp_minus_z = exp(-z)
        Ei = zinv * exp_minus_z
        for i in 1:(-n)
            Ei = zinv * (exp_minus_z + i * Ei)
        end
        return Ei
    else
        # forwards recurrence from E₁
        exp_minus_z = exp(-z)
        Ei = expint(z)
        Ei *= !isinf(Ei)
        for i in 2:n
            Ei = (exp_minus_z - z * Ei) / (i - 1)
        end
        return Ei
    end
end

"Return the first row of a Wflow output CSV file as a NamedTuple"
function csv_first_row(path)
    # silly function to avoid needing CSV.jl as a test dependency
    header, dataline = open(path) do io
        header = readline(io)
        dataline = readline(io)
        (header, dataline)
    end

    names = Tuple(Symbol.(split(header, ',')))
    ncol = length(names)
    # this assumes the first column is a time, the rest a float
    types = Tuple{DateTime, fill(Float64, ncol - 1)...}

    parts = split(dataline, ',')
    values = parse.(Float64, parts[2:end])
    row = NamedTuple{names, types}((DateTime(parts[1]), values...))
    return row
end
