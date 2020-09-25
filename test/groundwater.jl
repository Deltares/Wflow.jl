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

using Base.MathConstants: eulergamma

# n coefficients of the Taylor series of E₁(z) + log(z), in type T:
function E₁_taylor_coefficients(::Type{T}, n::Integer) where T <: Number
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
macro E₁_taylor64(z, n::Integer)
    c = E₁_taylor_coefficients(Float64, n)
    taylor = :(@evalpoly zz)
    append!(taylor.args, c)
    quote
        let zz = $(esc(z))
            $taylor - log(zz)
        end
    end
end

# for numeric-literal coefficients: simplify to a ratio of two polynomials:
import Polynomials
# return (p,q): the polynomials p(x) / q(x) corresponding to E₁_cf(x, a...),
# but without the exp(-x) term
function E₁_cfpoly(n::Integer, ::Type{T}=BigInt) where T <: Real
    q = Polynomials.Polynomial(T[1])
    p = x = Polynomials.Polynomial(T[0,1])
    for i in n:-1:1
        p, q = x * p + (1 + i) * q, p # from cf = x + (1+i)/cf = x + (1+i)*q/p
        p, q = p + i * q, p     # from cf = 1 + i/cf = 1 + i*q/p
    end
    # do final 1/(x + inv(cf)) = 1/(x + q/p) = p/(x*p + q)
    return p, x * p + q
end

macro E₁_cf64(z, n::Integer)
    p, q = E₁_cfpoly(n, BigInt)
    num_expr =  :(@evalpoly zz)
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
function expint(z::Union{Float64,Complex{Float64}})
    x² = real(z)^2
    y² = imag(z)^2
    if real(z) > 0 && x² + 0.233 * y² ≥ 7.84 # use cf expansion, ≤ 30 terms
        if (x² ≥ 546121) & (real(z) > 0) # underflow
            return zero(z)
        elseif x² + 0.401 * y² ≥ 58.0 # ≤ 15 terms
            if x² + 0.649 * y² ≥ 540.0 # ≤ 8 terms
                x² + y² ≥ 4e4 && return @E₁_cf64 z 4
                return @E₁_cf64 z 8
            end
            return @E₁_cf64 z 15
        end
        return @E₁_cf64 z 30
    else # use Taylor expansion, ≤ 37 terms
        r² = x² + y²
        return r² ≤ 0.36 ? (r² ≤ 2.8e-3 ? (r² ≤ 2e-7 ? @E₁_taylor64(z,4) :
                                                        @E₁_taylor64(z,8)) :
                                                        @E₁_taylor64(z,15)) :
                                                        @E₁_taylor64(z,37)
    end
end
expint(z::Union{T,Complex{T},Rational{T},Complex{Rational{T}}}) where {T <: Integer} = expint(float(z))

######################################################################
# exponential integral Eₙ(z)

function expint(n::Integer, z)
    if n == 1
        return expint(z)
    elseif n < 1
        # backwards recurrence from E₀ = e⁻ᶻ/z
        zinv = inv(z)
        e⁻ᶻ = exp(-z)
        Eᵢ = zinv * e⁻ᶻ
        for i in 1:-n
            Eᵢ = zinv * (e⁻ᶻ + i * Eᵢ)
        end
        return Eᵢ
    elseif n > 1
        # forwards recurrence from E₁
        e⁻ᶻ = exp(-z)
        Eᵢ = expint(z)
        Eᵢ *= !isinf(Eᵢ)
        for i in 2:n
            Eᵢ = (e⁻ᶻ - z * Eᵢ) / (i - 1)
        end
        return Eᵢ
    end
end

"""
    drawdown_theis(distance, time, discharge, transmissivity, storativity)

Non-steady flow in a confined aquifer, using the well function of Theis.
"""
function drawdown_theis(distance, time, discharge, transmissivity, storativity)
    u = (storativity * distance^2) / (4 * transmissivity * time)
    return discharge / (4 * π  * transmissivity) * expint(1, u)
end


# Avoid 0.0, since it's a singularity
# X = -9995.0:10.0:9995.0
# Y = -9995.0:10.0:9995.0
# time = 0.0:100.0:10_000.0
# ϕ = [-drawdown_theis(√(x^2 + y^2), t, discharge, kD, S) for t in time, x in X, y in Y]


@testset "groundwater" begin
    ncol = 2
    nrow = 3
    # TODO: this seems weird to me, I guess x is firstdim due to column major?
    # which happens because of spatial raster layout?
    # relates to how active_indices works
    shape = (ncol, nrow)
    Δx = [10.0, 20.0]
    Δy = [5.0, 15.0, 25.0]
    collect_connections(con, cell_id) = [con.rowval[nzi] for nzi in Wflow.connections(con, cell_id)]

    @testset "unit: connectivity" begin
        @testset "connection_geometry: y" begin
            I = CartesianIndex(1, 1)
            J = CartesianIndex(2, 1)
            @test Wflow.connection_geometry(I, J, Δx, Δy) == (2.5, 7.5, 10.0)
            @test_throws Exception Wflow.connection_geometry(I, I, Δx, Δy)
        end

        @testset "connection_geometry: x" begin
            I = CartesianIndex(1, 1)
            J = CartesianIndex(1, 2)
            @test Wflow.connection_geometry(I, J, Δx, Δy) == (5.0, 10.0, 5.0)   
            @test_throws Exception Wflow.connection_geometry(I, I, Δx, Δy)
        end

        @testset "Connectivity 1D(x)" begin
            # +---+---+
            # | 1 | 2 |
            # +---+---+
            domain = ones(Bool, (1, 2))
            indices, reverse_indices = Wflow.active_indices(domain, false)
            conn = Wflow.Connectivity(indices, reverse_indices, Δx, [5.0])
            @test conn.ncell == 2 
            @test conn.nconnection == 2
            @test conn.length1 == [5.0, 10.0]
            @test conn.length2 == [10.0, 5.0]
            @test conn.width == [5.0, 5.0]
            @test conn.colptr == [1, 2, 3]
            @test conn.rowval == [2, 1]
            @test collect_connections(conn, 1) == [2]
            @test collect_connections(conn, 2) == [1]
        end

        @testset "Connectivity 1D(y)" begin
            # +---+
            # | 1 |
            # +---+
            # | 2 |
            # +---+
            # | 3 |
            # +---+
            domain = ones(Bool, (3, 1))
            indices, reverse_indices = Wflow.active_indices(domain, false)
            conn = Wflow.Connectivity(indices, reverse_indices, [10.0], Δy)
            @test conn.ncell == 3 
            @test conn.nconnection == 4
            @test conn.length1 == [2.5, 7.5, 7.5, 12.5]
            @test conn.length2 == [7.5, 2.5, 12.5, 7.5]
            @test conn.width == [10.0, 10.0, 10.0, 10.0]
            @test conn.colptr == [1, 2, 4, 5]
            @test conn.rowval == [2, 1, 3, 2]
            @test collect_connections(conn, 1) == [2]
            @test collect_connections(conn, 2) == [1, 3]
            @test collect_connections(conn, 3) == [2]
        end

        @testset "Connectivity 2D" begin
            # +---+---+
            # | 1 | 4 |
            # +---+---+
            # | 2 | 5 |
            # +---+---+
            # | 3 | 6 |
            # +---+---+
            domain = ones(Bool, (nrow, ncol))
            indices, reverse_indices = Wflow.active_indices(domain, false)
            conn = Wflow.Connectivity(indices, reverse_indices, Δx, Δy)
            @test conn.ncell == 6
            @test conn.nconnection == 14
            @test conn.colptr == [1, 3, 6, 8, 10, 13, 15]
            @test conn.rowval == [2, 4, 1, 3, 5, 2, 6, 1, 5, 2, 4, 6, 3, 5] 
            @test collect_connections(conn, 1) == [2, 4]
            @test collect_connections(conn, 2) == [1, 3, 5]
            @test collect_connections(conn, 3) == [2, 6]
            @test collect_connections(conn, 4) == [1, 5]
            @test collect_connections(conn, 5) == [2, 4, 6]
            @test collect_connections(conn, 6) == [3, 5]
        end
    end

    @testset "unit: aquifer, boundary conditions" begin
        @testset "harmonicmean_conductance" begin
            # harmonicmean_conductance(k1, k2, H1, H2, l1, l2, width)
            @test Wflow.harmonicmean_conductance(10.0, 10.0, 5.0, 5.0, 0.5, 0.5, 1.0) == 50.0
            @test Wflow.harmonicmean_conductance(10.0, 10.0, 0.0, 5.0, 0.5, 0.5, 1.0) == 0.0
            @test Wflow.harmonicmean_conductance(10.0, 10.0, 5.0, 0.0, 0.5, 0.5, 1.0) == 0.0
            # kD of 10 and 20 -> harmonicmean = 1/(1/10 + 1/20)
            @test Wflow.harmonicmean_conductance(10.0, 10.0, 1.0, 2.0, 1.0, 1.0, 1.0) ≈ (6.0 + 2.0 / 3.0)
        end

        nrow = 1
        ncol = 3
        shape = (nrow, ncol)
        # Domain, geometry
        domain = ones(Bool, shape)
        Δx = fill(10.0, ncol)
        Δy = fill(10.0, nrow)
        indices, reverse_indices = Wflow.active_indices(domain, false)
        connectivity = Wflow.Connectivity(indices, reverse_indices, Δx, Δy)

        ncell = connectivity.ncell

        conf_aqf = Wflow.ConfinedAquifer(
            [0.0, 7.5, 20.0],  # head
            fill(10.0, ncell),  # k
            fill(10.0, ncell),  # top
            fill(0.0, ncell),  # bottom
            fill(100.0, ncell),  # area
            fill(1.0e-5, ncell), # specific storage
            fill(1.0e-4, ncell),  # storativity 
            fill(0.0, connectivity.nconnection),  # conductance
        )
        unconf_aqf = Wflow.UnconfinedAquifer(
            [0.0, 7.5, 20.0],
            fill(10.0, ncell),
            fill(10.0, ncell),
            fill(0.0, ncell),
            fill(100.0, ncell),
            fill(0.15, ncell),
        )

        @testset "saturated_thickness-confined" begin
            @test (
                Wflow.saturated_thickness(conf_aqf, 1) == 
                Wflow.saturated_thickness(conf_aqf, 2) == 
                Wflow.saturated_thickness(conf_aqf, 3) == 
                10.0
                )
        end

        @testset "saturated_thickness-unconfined" begin
            @test Wflow.saturated_thickness(unconf_aqf, 1) == 0.0
            @test Wflow.saturated_thickness(unconf_aqf, 2) == 7.5
            @test Wflow.saturated_thickness(unconf_aqf, 3) == 10.0
        end

        @testset "horizontal_conductance" begin
            @test (
                Wflow.horizontal_conductance(1, 2, 1, conf_aqf, connectivity) == 
                Wflow.harmonicmean_conductance(10.0, 10.0, 10.0, 10.0, 5.0, 5.0, 10.0)
            )
        end

        @testset "conductance" begin
            @test Wflow.conductance(conf_aqf, connectivity, 2, 3, 3) == 0.0  # just returns a set value
            @test Wflow.conductance(unconf_aqf, connectivity, 2, 3, 3) > 0.0   # computes a value based on head
            @test Wflow.conductance(unconf_aqf, connectivity, 1, 2, 1) == 0.0
        end

        @testset "flux-confined" begin
            Q = zeros(3)
            Wflow.flux!(Q, conf_aqf, connectivity)
            @test all(Q .== 0.0)
        end

        @testset "flux-unconfined" begin
            Q = zeros(3)
            Wflow.flux!(Q, unconf_aqf, connectivity)
            @test Q[1] == 0.0 # no flow 1 <-> 2
            @test Q[2] > 0.0 # flow 3 -> 2
            @test Q[3] < 0.0
            @test -Q[2] == Q[3]
        end

        @testset "river" begin
            river = Wflow.River(
                [2.0, 2.0], [100.0, 100.0], [200.0, 200.0], [1.0, 1.0], [1, 3]
            )
            Q = zeros(3)
            Wflow.flux!(Q, river, conf_aqf)
            # infiltration, below bottom, flux is (stage - bottom) * inf_cond
            @test Q[1] == (2.0 - 1.0) * 100.0
            # drainage, flux is () * exf_cond
            @test Q[3] == (2.0 - 20.0) * 200.0
        end

        @testset "drainage" begin
            drainage = Wflow.Drainage([2.0, 2.0], [100.0, 100.0], [1, 2])
            Q = zeros(3)
            Wflow.flux!(Q, drainage, conf_aqf)
            @test Q[1] == 0.0
            @test Q[2] == 100.0 * (2.0 - 7.5)
        end

        @testset "headboundary" begin
            headboundary = Wflow.HeadBoundary([2.0, 2.0], [100.0, 100.0], [1, 2])
            Q = zeros(3)
            Wflow.flux!(Q, headboundary, conf_aqf)
            @test Q[1] == 100.0 * (2.0 - 0.0)
            @test Q[2] == 100.0 * (2.0 - 7.5)
        end

        @testset "recharge" begin
            recharge = Wflow.Recharge([1.0e-3, 1.0e-3, 1.0e-3], [1, 2, 3])
            Q = zeros(3)
            Wflow.flux!(Q, recharge, conf_aqf)
            @test all(Q .== 1.0e-3 * 100.0)
        end

        @testset "well" begin
            well = Wflow.Well([-1000.0], [1])
            Q = zeros(3)
            Wflow.flux!(Q, well, conf_aqf)
            @test Q[1] == -1000.0
        end
    end
    

#    @testset "integration: flow 1D" begin
#        aquifer, connectivity = homogenous_aquifer(3, 1)
#        constanthead = ConstantHead([5.0 10.0], [1 3])
#        gwf = GroundwaterFlow(
#            aquifer,
#            connectivity,
#            constanthead,
#            [],
#        )
#    
#        # A hundred timesteps
#        Δt = 0.5 # d
#        for _ in 1:100
#            Wflow.update(gwf, Δt)
#        end
#    
#        @test gwf.head ≈ [5.0 7.5 10.0]
#    end

end