
function initial_head(x)
    2 * √x
end

"""
    transient_aquifer_1d(x, time, conductivity, specific_yield, aquifer_length, Β)

Non-steady flow in an unconfined rectangular aquifer, with Dirichlet h(0, t) = 0
on the left edge, and a Neumann Boundary Condition (dh/dx = 0) on the right.
"""
function transient_aquifer_1d(x, time, conductivity, specific_yield, aquifer_length, Β)
    initial_head(x) / 1.0 +
    (Β * conductivity * initial_head(aquifer_length) * time) /
    (specific_yield * aquifer_length * aquifer_length)
end


"""
    drawdown_theis(distance, time, discharge, transmissivity, storativity)

Non-steady flow in a confined aquifer, using the well function of Theis.
"""
function drawdown_theis(distance, time, discharge, transmissivity, storativity)
    u = (storativity * distance^2) / (4 * transmissivity * time)
    return discharge / (4 * π * transmissivity) * expint(u)
end


function homogenous_aquifer(nrow, ncol)
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
        fill(0.1, ncell), # specific storage
        fill(1.0, ncell),  # storativity
        fill(0.0, connectivity.nconnection),  # conductance
        false,
    )
    unconf_aqf = Wflow.UnconfinedAquifer(
        [0.0, 7.5, 20.0],  # head
        fill(10.0, ncell),   # k
        fill(10.0, ncell),  # top
        fill(0.0, ncell),  # bottom
        fill(100.0, ncell),  # area
        fill(0.15, ncell),  # specific yield
        fill(0.0, connectivity.nconnection),  # conductance
        false, # toggle for reduction factor
        fill(3.0, ncell) # conductance reduction factor
    )
    return (connectivity, conf_aqf, unconf_aqf)
end


@testset "groundwater" begin
    ncol = 2
    nrow = 3
    shape = (ncol, nrow)
    Δx = [10.0, 20.0]
    Δy = [5.0, 15.0, 25.0]
    collect_connections(con, cell_id) =
        [con.rowval[nzi] for nzi in Wflow.connections(con, cell_id)]

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

        @testset "Connectivity 2D - partially inactive" begin
            # +---+---+
            # | 1 | 3 |
            # +---+---+
            # | X | 4 |
            # +---+---+
            # | 2 | 5 |
            # +---+---+
            domain = ones(Bool, (nrow, ncol))
            domain[2, 1] = false
            indices, reverse_indices = Wflow.active_indices(domain, false)
            conn = Wflow.Connectivity(indices, reverse_indices, Δx, Δy)
            @test conn.ncell == 5
            @test conn.nconnection == 8
            @test conn.colptr == [1, 2, 3, 5, 7, 9]
            @test conn.rowval == [3, 5, 1, 4, 3, 5, 2, 4]
            @test collect_connections(conn, 1) == [3]
            @test collect_connections(conn, 2) == [5]
            @test collect_connections(conn, 3) == [1, 4]
            @test collect_connections(conn, 4) == [3, 5]
            @test collect_connections(conn, 5) == [2, 4]
        end
    end

    @testset "unit: aquifer, boundary conditions" begin
        @testset "harmonicmean_conductance" begin
            # harmonicmean_conductance(k1, k2, H1, H2, l1, l2, width)
            @test Wflow.harmonicmean_conductance(10.0, 10.0, 5.0, 5.0, 0.5, 0.5, 1.0) ==
                  50.0
            @test Wflow.harmonicmean_conductance(10.0, 10.0, 0.0, 5.0, 0.5, 0.5, 1.0) == 0.0
            @test Wflow.harmonicmean_conductance(10.0, 10.0, 5.0, 0.0, 0.5, 0.5, 1.0) == 0.0
            # kD of 10 and 20 -> harmonicmean = 1/(1/10 + 1/20)
            @test Wflow.harmonicmean_conductance(10.0, 10.0, 1.0, 2.0, 1.0, 1.0, 1.0) ≈
                  (6.0 + 2.0 / 3.0)
        end

        nrow = 1
        ncol = 3
        connectivity, conf_aqf, unconf_aqf = homogenous_aquifer(nrow, ncol)
        Wflow.initialize_conductance!(conf_aqf, connectivity)
        Wflow.initialize_conductance!(unconf_aqf, connectivity)
        ncell = connectivity.ncell

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
            @test Wflow.conductance(conf_aqf, 2, 3, 3, false, connectivity) == 100.0
            @test Wflow.conductance(unconf_aqf, 2, 3, 3, false, connectivity) == 100.0  # upstream sat. thickness
            @test Wflow.conductance(unconf_aqf, 1, 2, 1, false, connectivity) == 75.0  # upstream sat. thickness
        end

        @testset "minimum_head-confined" begin
            original_head = copy(conf_aqf.head)
            conf_aqf.head[1] = -10.0
            @test Wflow.check_flux(-1.0, conf_aqf, 1) == -1.0
            @test Wflow.minimum_head(conf_aqf)[1] == -10.0
            conf_aqf.head .= original_head
        end

        @testset "minimum_head-unconfined" begin
            original_head = copy(unconf_aqf.head)
            unconf_aqf.head[1] = -10.0
            @test Wflow.check_flux(-1.0, unconf_aqf, 1) == 0.0
            @test Wflow.minimum_head(conf_aqf)[1] == 0.0
            unconf_aqf.head .= original_head
        end

        @testset "stable_timestep" begin
            @test Wflow.stable_timestep(conf_aqf, false) == 0.25
        end

        # Parametrization in setup is as follows:
        # [0.0, 7.5, 20.0],  # head
        # fill(10.0, ncell),  # k
        # fill(10.0, ncell),  # top
        # fill(0.0, ncell),  # bottom

        @testset "flux-confined" begin
            Q = zeros(3)
            Wflow.flux!(Q, conf_aqf, connectivity, false)
            # kD = 10 * 10 = 100
            # dH = 7.5, 12.5
            @test Q == [750.0, 500.0, -1250.0]  # TODO
        end

        @testset "flux-unconfined" begin
            Q = zeros(3)
            Wflow.flux!(Q, unconf_aqf, connectivity, false)
            # KD is based on upstream saturated thickness, i.e. 7.5 m and 20.0 m (which is capped to 10.0)
            @test Q == [562.5, 687.5, -1250.0]
        end

        @testset "river" begin
            river = Wflow.River(
                [2.0, 2.0],
                [100.0, 100.0],
                [200.0, 200.0],
                [1.0, 1.0],
                [0.0, 0.0],
                [1, 3],
            )
            Q = zeros(3)
            Wflow.flux!(Q, river, conf_aqf)
            # infiltration, below bottom, flux is (stage - bottom) * inf_cond
            @test Q[1] == (2.0 - 1.0) * 100.0
            # drainage, flux is () * exf_cond
            @test Q[3] == (2.0 - 20.0) * 200.0
        end

        @testset "drainage" begin
            drainage = Wflow.Drainage([2.0, 2.0], [100.0, 100.0], [0.0, 0.0], [1, 2])
            Q = zeros(3)
            Wflow.flux!(Q, drainage, conf_aqf)
            @test Q[1] == 0.0
            @test Q[2] == 100.0 * (2.0 - 7.5)
        end

        @testset "headboundary" begin
            headboundary =
                Wflow.HeadBoundary([2.0, 2.0], [100.0, 100.0], [0.0, 0.0], [1, 2])
            Q = zeros(3)
            Wflow.flux!(Q, headboundary, conf_aqf)
            @test Q[1] == 100.0 * (2.0 - 0.0)
            @test Q[2] == 100.0 * (2.0 - 7.5)
        end

        @testset "recharge" begin
            recharge = Wflow.Recharge([1.0e-3, 1.0e-3, 1.0e-3], [0.0, 0.0, 0.0], [1, 2, 3])
            Q = zeros(3)
            Wflow.flux!(Q, recharge, conf_aqf)
            @test all(Q .== 1.0e-3 * 100.0)
        end

        @testset "well" begin
            well = Wflow.Well([-1000.0], [0.0], [1])
            Q = zeros(3)
            Wflow.flux!(Q, well, conf_aqf)
            @test Q[1] == -1000.0
        end
    end

    @testset "integration: steady 1D" begin
        connectivity, aquifer, _ = homogenous_aquifer(3, 1)
        constanthead = Wflow.ConstantHead([2.0, 4.0], [1, 3])
        gwf = Wflow.GroundwaterFlow(
            aquifer,
            connectivity,
            constanthead,
            Wflow.AquiferBoundaryCondition[],
        )
        # Set constant head (dirichlet) boundaries
        gwf.aquifer.head[gwf.constanthead.index] .= gwf.constanthead.head

        Q = zeros(3)
        Δt = 0.25 # days
        for _ = 1:50
            Wflow.update(gwf, Q, Δt)
        end

        @test gwf.aquifer.head ≈ [2.0, 3.0, 4.0]
    end

    @testset "integration: unconfined transient 1D" begin
        nrow = 1
        ncol = 9
        shape = (nrow, ncol)
        conductivity = 200.0
        top = 150.0
        bottom = 0.0
        specific_yield = 0.15
        cellsize = 500.0
        Β = 1.12
        aquifer_length = cellsize * ncol
        gwf_f = 3.0
        exp_conductivity = false

        # Domain, geometry
        domain = ones(Bool, shape)
        Δx = fill(cellsize, ncol)
        Δy = fill(cellsize, nrow)
        indices, reverse_indices = Wflow.active_indices(domain, false)
        connectivity = Wflow.Connectivity(indices, reverse_indices, Δx, Δy)
        ncell = connectivity.ncell
        xc = collect(range(0.0, stop = aquifer_length - cellsize, step = cellsize))
        aquifer = Wflow.UnconfinedAquifer(
            initial_head.(xc),
            fill(conductivity, ncell),
            fill(top, ncell),
            fill(bottom, ncell),
            fill(cellsize * cellsize, ncell),
            fill(specific_yield, ncell),
            fill(0.0, connectivity.nconnection),
            exp_conductivity,
            fill(gwf_f, ncell),
        )
        # constant head on left boundary, 0 at 0
        constanthead = Wflow.ConstantHead([0.0], [1])
        gwf = Wflow.GroundwaterFlow(
            aquifer,
            connectivity,
            constanthead,
            Wflow.AquiferBoundaryCondition[],
        )

        Δt = Wflow.stable_timestep(gwf.aquifer, exp_conductivity)
        Q = zeros(ncell)
        time = 20.0
        nstep = Int(ceil(time / Δt))
        time = nstep * Δt

        for i = 1:nstep
            Wflow.update(gwf, Q, Δt)
            # Gradient dh/dx is positive, all flow to the left
            @test all(diff(gwf.aquifer.head) .> 0.0)
        end

        ϕ_analytical = [
            transient_aquifer_1d(x, time, conductivity, specific_yield, aquifer_length, Β) for x in xc
        ]
        difference = gwf.aquifer.head .- ϕ_analytical
        # @test all(difference .< ?)  #TODO
    end

    @testset "integration: confined transient radial 2D" begin
        halfnrow = 20
        wellrow = halfnrow + 1
        nrow = halfnrow * 2 + 1
        ncol = nrow
        shape = (nrow, ncol)
        conductivity = 5.0
        top = 10.0
        bottom = 0.0
        transmissivity = (top - bottom) * conductivity
        cellsize = 10.0
        startinghead = top
        specific_storage = 0.015
        storativity = 0.15
        aquifer_length = cellsize * ncol
        discharge = -50.0
        exp_conductivity = false

        # Domain, geometry
        domain = ones(Bool, shape)
        Δx = fill(cellsize, ncol)
        Δy = fill(cellsize, nrow)
        indices, reverse_indices = Wflow.active_indices(domain, false)
        connectivity = Wflow.Connectivity(indices, reverse_indices, Δx, Δy)
        ncell = connectivity.ncell
        aquifer = Wflow.ConfinedAquifer(
            fill(startinghead, ncell),
            fill(conductivity, ncell),
            fill(top, ncell),
            fill(bottom, ncell),
            fill(cellsize * cellsize, ncell),
            fill(specific_storage, ncell),
            fill(storativity, ncell),
            fill(0.0, connectivity.nconnection), # conductance, to be set
            exp_conductivity,
        )

        cell_index = reshape(collect(range(1, ncell, step = 1)), shape)
        indices = vcat(cell_index[1, :], cell_index[end, :])# , cell_index[:, 1], cell_index[:, end],)
        constanthead = Wflow.ConstantHead(fill(10.0, size(indices)), indices)
        # Place a well in the middle of the domain
        well = Wflow.Well([discharge], [0.0], [reverse_indices[wellrow, wellrow]])
        gwf = Wflow.GroundwaterFlow(aquifer, connectivity, constanthead, [well])

        Δt = Wflow.stable_timestep(gwf.aquifer, exp_conductivity)
        Q = zeros(ncell)
        time = 20.0
        nstep = Int(ceil(time / Δt))
        time = nstep * Δt

        for i = 1:nstep
            Wflow.update(gwf, Q, Δt)
        end

        # test for symmetry on x and y axes
        ϕ = reshape(gwf.aquifer.head, shape)
        @test ϕ[1:halfnrow, :] ≈ ϕ[end:-1:halfnrow+2, :]
        @test ϕ[:, 1:halfnrow] ≈ ϕ[:, end:-1:halfnrow+2]

        # compare with analytical solution
        start = -0.5 * aquifer_length + 0.5 * cellsize
        stop = 0.5 * aquifer_length - 0.5 * cellsize
        X = collect(range(start, stop = stop, step = cellsize))
        ϕ_analytical =
            [drawdown_theis(x, time, discharge, transmissivity, storativity) for x in X] .+ 10.0
        # compare left-side, since it's symmetric anyway. Skip the well cell, and its first neighbor
        difference = ϕ[1:halfnrow-1, halfnrow] - ϕ_analytical[1:halfnrow-1]
        @test all(difference .< 0.02)
    end
end
