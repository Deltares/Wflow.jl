"""
    scurve(x, a, b, c)

Sigmoid "S"-shaped curve.

# Arguments
- `x::Real`: input
- `a::Real`: determines the center level
- `b::Real`: determines the amplitude of the curve (range: (0, b⁻¹))
- `c::Real`: determines the steepness or "stepwiseness" of the curve.
             The higher c the sharper the function. A negative c reverses the function.
"""
function scurve(x::Real, a::Real, b::Real, c::Real)::Real
    s = inv(b + exp(-c * (x - a)))
    return s
end

# 2.5x faster power method
"Faster method for exponentiation"
pow(x::Real, y::Real)::Real = exp(y * log(x))

"""
    bounded_power(base, power)

Computes min(base^power, 1) without computing the power if the result is known to be larger
than 1. Assumes base, power > 0.
"""
function bounded_power(base::T, power) where {T}
    return if base > 1
        one(T)
    else
        pow(base, power)
    end
end

"Partition indices with at least size `basesize`"
function _partition(xs::Integer, basesize::Integer)
    n = Int(max(1, xs ÷ basesize))
    return (Int(1 + ((i - 1) * xs) ÷ n):Int((i * xs) ÷ n) for i in 1:n)
end

"""
    threaded_foreach(f, x::AbstractArray; basesize::Integer)

Run function `f` in parallel by spawning tasks, each task iterates over a chunk of size
`basesize`. Falls back to a serial loop when running single-threaded.
"""
function threaded_foreach(f::Function, x::AbstractArray; basesize::Integer)::Nothing
    len = length(x)
    partitions = _partition(len, basesize)
    if length(partitions) > 1 && Threads.nthreads() > 1
        @sync for p in partitions
            Threads.@spawn begin
                for i in eachindex(p)
                    f(@inbounds p[i])
                end
            end
        end
    else
        for i in eachindex(x)
            f(@inbounds x[i])
        end
    end
    return nothing
end

"""
    set_layerthickness(reference_depth, cum_depth, thickness)

Calculate actual soil thickness of layers based on a reference depth (e.g. soil depth or
water table depth) `reference_depth`, a SVector `cum_depth` with cumulative soil depth
starting at soil surface (0), and a SVector `thickness` with thickness per soil layer.
"""
function set_layerthickness(
    reference_depth::Real,
    cum_depth::SVector,
    thickness::SVector{N, Float64},
)::SVector{N, Float64} where {N}
    thicknesslayers = thickness .* MISSING_VALUE
    for i in 1:length(thicknesslayers)
        if reference_depth > cum_depth[i + 1]
            thicknesslayers = setindex(thicknesslayers, thickness[i], i)
        elseif reference_depth - cum_depth[i] > 0.0
            thicknesslayers = setindex(thicknesslayers, reference_depth - cum_depth[i], i)
        end
    end
    return thicknesslayers
end

function number_of_active_layers(thickness::SVector)::Int
    number_of_layers = length(thickness) - sum(isnan.(thickness))
    return number_of_layers
end

"Set lower bound for drainable porosity"
function lower_bound_drainable_porosity(theta_s, theta_fc; lower_bound = 0.02)
    return max(theta_s - theta_fc, lower_bound)
end

"""
    hydraulic_conductivity_at_depth(p, vertical_hydraulic_conductivity_factor, z, i, n)

Return vertical hydraulic conductivity `kv_z` at depth `z` for index `i` using multiplication
factor at soil layer `n` and vertical hydraulic conductivity profile `p`.
"""
function hydraulic_conductivity_at_depth(
    p::KvExponential,
    vertical_hydraulic_conductivity_factor,
    z,
    i,
    n,
)
    kv_z =
        vertical_hydraulic_conductivity_factor[i][n] *
        p.kv_0[i] *
        exp(-p.hydraulic_conductivity_scale_parameter[i] * z)
    return kv_z
end

function hydraulic_conductivity_at_depth(
    p::KvExponentialConstant,
    vertical_hydraulic_conductivity_factor,
    z,
    i,
    n,
)
    (; kv_0, hydraulic_conductivity_scale_parameter) = p.exponential
    if z < p.z_exp[i]
        kv_z =
            vertical_hydraulic_conductivity_factor[i][n] *
            kv_0[i] *
            exp(-hydraulic_conductivity_scale_parameter[i] * z)
    else
        kv_z =
            vertical_hydraulic_conductivity_factor[i][n] *
            kv_0[i] *
            exp(-hydraulic_conductivity_scale_parameter[i] * p.z_exp[i])
    end
    return kv_z
end

function hydraulic_conductivity_at_depth(
    p::KvLayered,
    vertical_hydraulic_conductivity_factor,
    z,
    i,
    n,
)
    kv_z = vertical_hydraulic_conductivity_factor[i][n] * p.kv[i][n]
    return kv_z
end

function hydraulic_conductivity_at_depth(
    p::KvLayeredExponential,
    vertical_hydraulic_conductivity_factor,
    z,
    i,
    n,
)
    return if z < p.z_layered[i]
        vertical_hydraulic_conductivity_factor[i][n] * p.kv[i][n]
    else
        n = p.nlayers_kv[i]
        vertical_hydraulic_conductivity_factor[i][n] *
        p.kv[i][n] *
        exp(-p.hydraulic_conductivity_scale_parameter[i] * (z - p.z_layered[i]))
    end
end
