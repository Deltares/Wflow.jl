
function initial_head(x)
    return 2 * √x
end

"""
    transient_aquifer_1d(x, time, conductivity, specific_yield, aquifer_length, beta)

Non-steady flow in an unconfined rectangular aquifer, with Dirichlet h(0, t) = 0
on the left edge, and a Neumann Boundary Condition (dh/dx = 0) on the right.
"""
function transient_aquifer_1d(x, time, conductivity, specific_yield, aquifer_length, beta)
    return initial_head(x) / 1.0 +
           (beta * conductivity * initial_head(aquifer_length) * time) /
           (specific_yield * aquifer_length * aquifer_length)
end

"""
    drawdown_theis(distance, time, discharge, transmissivity, storativity)

Non-steady flow in a confined aquifer, using the well function of Theis.
"""
function drawdown_theis(distance, time, discharge, transmissivity, storativity)
    u = (storativity * distance^2) / (4 * transmissivity * time)
    return discharge / (4 * pi * transmissivity) * expint(u)
end

function homogenous_aquifer(nrow, ncol)
    shape = (nrow, ncol)
    # Domain, geometry
    domain = ones(Bool, shape)
    dx = fill(10.0, ncol)
    dy = fill(10.0, nrow)
    indices, reverse_indices = Wflow.active_indices(domain, false)
    connectivity = Wflow.Connectivity(indices, reverse_indices, dx, dy)
    ncell = connectivity.ncell

    parameters = Wflow.ConfinedAquiferParameters(;
        k = fill(10.0, ncell),
        top = fill(10.0, ncell),
        bottom = fill(0.0, ncell),
        area = fill(100.0, ncell),
        specific_storage = fill(0.1, ncell),
        storativity = fill(1.0, ncell),
    )
    variables = Wflow.ConfinedAquiferVariables(;
        head = [0.0, 7.5, 20.0],
        conductance = fill(0.0, connectivity.nconnection),
        storage = fill(0.0, ncell),
    )
    conf_aqf = Wflow.ConfinedAquifer(; parameters, variables)

    parameters = Wflow.UnconfinedAquiferParameters(;
        k = fill(10.0, ncell),
        top = fill(10.0, ncell),
        bottom = fill(0.0, ncell),
        area = fill(100.0, ncell),
        specific_yield = fill(0.15, ncell),
        f = fill(3.0, ncell),
    )
    variables = Wflow.UnconfinedAquiferVariables(;
        head = [0.0, 7.5, 20.0],
        conductance = fill(0.0, connectivity.nconnection),
        storage = fill(0.0, ncell),
    )
    unconf_aqf = Wflow.UnconfinedAquifer(; parameters, variables)
    return (connectivity, conf_aqf, unconf_aqf)
end

@testset "groundwater" begin
    ncol = 2
    nrow = 3
    shape = (ncol, nrow)
    dx = [10.0, 20.0]
    dy = [5.0, 15.0, 25.0]
    collect_connections(con, cell_id) =
        [con.rowval[nzi] for nzi in Wflow.connections(con, cell_id)]

    @testset "unit: connectivity" begin
        @testset "connection_geometry: y" begin
            I = CartesianIndex(1, 1)
            J = CartesianIndex(2, 1)
            @test Wflow.connection_geometry(I, J, dx, dy) == (2.5, 7.5, 10.0)
            @test_throws Exception Wflow.connection_geometry(I, I, dx, dy)
        end

        @testset "connection_geometry: x" begin
            I = CartesianIndex(1, 1)
            J = CartesianIndex(1, 2)
            @test Wflow.connection_geometry(I, J, dx, dy) == (5.0, 10.0, 5.0)
            @test_throws Exception Wflow.connection_geometry(I, I, dx, dy)
        end

        @testset "Connectivity 1D(x)" begin
            # +---+---+
            # | 1 | 2 |
            # +---+---+
            domain = ones(Bool, (1, 2))
            indices, reverse_indices = Wflow.active_indices(domain, false)
            conn = Wflow.Connectivity(indices, reverse_indices, dx, [5.0])
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
            conn = Wflow.Connectivity(indices, reverse_indices, [10.0], dy)
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
            conn = Wflow.Connectivity(indices, reverse_indices, dx, dy)
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
            conn = Wflow.Connectivity(indices, reverse_indices, dx, dy)
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
            # harmonicmean_conductance(kH1, kH2, l1, l2, width)
            @test Wflow.harmonicmean_conductance(10.0 * 5.0, 10.0 * 5.0, 0.5, 0.5, 1.0) ==
                  50.0
            @test Wflow.harmonicmean_conductance(10.0 * 0.0, 10.0 * 5.0, 0.5, 0.5, 1.0) ==
                  0.0
            @test Wflow.harmonicmean_conductance(10.0 * 5.0, 10.0 * 0.0, 0.5, 0.5, 1.0) ==
                  0.0
            # kD of 10 and 20 -> harmonicmean = 1/(1/10 + 1/20)
            @test Wflow.harmonicmean_conductance(10.0 * 1.0, 10.0 * 2.0, 1.0, 1.0, 1.0) ≈
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
                Wflow.harmonicmean_conductance(10.0 * 10.0, 10.0 * 10.0, 5.0, 5.0, 10.0)
            )
        end

        @testset "conductance" begin
            conductivity_profile = "uniform"
            @test Wflow.conductance(
                conf_aqf,
                2,
                3,
                3,
                conductivity_profile,
                connectivity,
            ) == 100.0
            @test Wflow.conductance(
                unconf_aqf,
                2,
                3,
                3,
                conductivity_profile,
                connectivity,
            ) == 100.0  # upstream sat. thickness
            @test Wflow.conductance(
                unconf_aqf,
                1,
                2,
                1,
                conductivity_profile,
                connectivity,
            ) == 75.0  # upstream sat. thickness
        end

        @testset "minimum_head-confined" begin
            original_head = copy(conf_aqf.variables.head)
            conf_aqf.variables.head[1] = -10.0
            @test Wflow.check_flux(-1.0, conf_aqf, 1) == -1.0
            @test Wflow.minimum_head(conf_aqf)[1] == -10.0
            conf_aqf.variables.head .= original_head
        end

        @testset "minimum_head-unconfined" begin
            original_head = copy(unconf_aqf.variables.head)
            unconf_aqf.variables.head[1] = -10.0
            @test Wflow.check_flux(-1.0, unconf_aqf, 1) == 0.0
            @test Wflow.minimum_head(conf_aqf)[1] == 0.0
            unconf_aqf.variables.head .= original_head
        end

        @testset "stable_timestep" begin
            conductivity_profile = "uniform"
            @test Wflow.stable_timestep(conf_aqf, conductivity_profile) == 0.25
        end

        # Parametrization in setup is as follows:
        # [0.0, 7.5, 20.0],  # head
        # fill(10.0, ncell),  # k
        # fill(10.0, ncell),  # top
        # fill(0.0, ncell),  # bottom

        @testset "flux-confined" begin
            Q = zeros(3)
            conductivity_profile = "uniform"
            Wflow.flux!(Q, conf_aqf, connectivity, conductivity_profile)
            # kD = 10 * 10 = 100
            # dH = 7.5, 12.5
            @test Q == [750.0, 500.0, -1250.0]  # TODO
        end

        @testset "flux-unconfined" begin
            Q = zeros(3)
            conductivity_profile = "uniform"
            Wflow.flux!(Q, unconf_aqf, connectivity, conductivity_profile)
            # KD is based on upstream saturated thickness, i.e. 7.5 m and 20.0 m (which is capped to 10.0)
            @test Q == [562.5, 687.5, -1250.0]
        end

        @testset "river" begin
            parameters = Wflow.GwfRiverParameters(;
                infiltration_conductance = [100.0, 100.0],
                exfiltration_conductance = [200.0, 200.0],
                bottom = [1.0, 1.0],
            )
            variables = Wflow.GwfRiverVariables(;
                stage = [2.0, 2.0],
                storage = [20.0, 20.0],
                flux = [0.0, 0.0],
            )
            river = Wflow.GwfRiver(; parameters, variables, index = [1, 3])
            Q = zeros(3)
            Wflow.flux!(Q, river, conf_aqf)
            # infiltration, below bottom, flux is (stage - bottom) * inf_cond, limited by
            # river storage (20.0)
            @test Q[1] == 20.0
            # drainage, flux is (stage - head) * exf_cond
            @test Q[3] == (2.0 - 20.0) * 200.0
        end

        @testset "drainage" begin
            parameters = Wflow.DrainageParameters(;
                elevation = [2.0, 2.0],
                conductance = [100.0, 100.0],
            )
            variables = Wflow.DrainageVariables(; flux = [0.0, 0.0])
            drainage = Wflow.Drainage(; parameters, variables, index = [1, 2])
            Q = zeros(3)
            Wflow.flux!(Q, drainage, conf_aqf)
            @test Q[1] == 0.0
            @test Q[2] == 100.0 * (2.0 - 7.5)
        end

        @testset "headboundary" begin
            parameters = Wflow.HeadBoundaryParameters(; conductance = [100.0, 100.0])
            variables = Wflow.HeadBoundaryVariables(; head = [2.0, 2.0], flux = [0.0, 0.0])

            headboundary = Wflow.HeadBoundary(; parameters, variables, index = [1, 2])
            Q = zeros(3)
            Wflow.flux!(Q, headboundary, conf_aqf)
            @test Q[1] == 100.0 * (2.0 - 0.0)
            @test Q[2] == 100.0 * (2.0 - 7.5)
        end

        @testset "recharge" begin
            variables = Wflow.RechargeVariables(;
                rate = [1.0e-3, 1.0e-3, 1.0e-3],
                flux = [0.0, 0.0, 0.0],
            )
            recharge = Wflow.Recharge(; variables, index = [1, 2, 3])
            Q = zeros(3)
            Wflow.flux!(Q, recharge, conf_aqf)
            @test all(Q .== 1.0e-3 * 100.0)
        end

        @testset "well" begin
            variables = Wflow.WellVariables(; volumetric_rate = [-1000.0], flux = [0.0])
            well = Wflow.Well(; variables, index = [1])
            Q = zeros(3)
            Wflow.flux!(Q, well, conf_aqf)
            @test Q[1] == -1000.0
        end
    end

    @testset "integration: steady 1D" begin
        connectivity, aquifer, _ = homogenous_aquifer(3, 1)
        variables = Wflow.ConstantHeadVariables(; head = [2.0, 4.0])
        constanthead = Wflow.ConstantHead(; variables, index = [1, 3])
        conductivity_profile = "uniform"
        gwf = Wflow.GroundwaterFlow(;
            aquifer = aquifer,
            connectivity = connectivity,
            constanthead = constanthead,
            boundaries = NamedTuple(),
        )
        # Set constant head (dirichlet) boundaries
        gwf.aquifer.variables.head[gwf.constanthead.index] .=
            gwf.constanthead.variables.head

        Q = zeros(3)
        dt = 0.25 # days
        for _ in 1:50
            Wflow.update!(gwf, Q, dt, conductivity_profile)
        end

        @test gwf.aquifer.variables.head ≈ [2.0, 3.0, 4.0]
    end

    @testset "integration: steady 1D, exponential conductivity" begin
        connectivity, aquifer, _ = homogenous_aquifer(3, 1)
        variables = Wflow.ConstantHeadVariables(; head = [2.0, 4.0])
        constanthead = Wflow.ConstantHead(; variables, index = [1, 3])
        conductivity_profile = "exponential"
        gwf = Wflow.GroundwaterFlow(;
            aquifer = aquifer,
            connectivity = connectivity,
            constanthead = constanthead,
            boundaries = NamedTuple(),
        )
        # Set constant head (dirichlet) boundaries
        gwf.aquifer.variables.head[gwf.constanthead.index] .=
            gwf.constanthead.variables.head

        Q = zeros(3)
        dt = 0.25 # days
        for _ in 1:50
            Wflow.update!(gwf, Q, dt, conductivity_profile)
        end

        @test gwf.aquifer.variables.head ≈ [2.0, 3.0, 4.0]
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
        beta = 1.12
        aquifer_length = cellsize * ncol
        gwf_f = 3.0
        conductivity_profile = "uniform"

        # Domain, geometry
        domain = ones(Bool, shape)
        dx = fill(cellsize, ncol)
        dy = fill(cellsize, nrow)
        indices, reverse_indices = Wflow.active_indices(domain, false)
        connectivity = Wflow.Connectivity(indices, reverse_indices, dx, dy)
        ncell = connectivity.ncell
        xc = collect(range(0.0; stop = aquifer_length - cellsize, step = cellsize))

        variables = Wflow.UnconfinedAquiferVariables(;
            head = initial_head.(xc),
            conductance = fill(0.0, connectivity.nconnection),
            storage = fill(0.0, ncell),
        )
        parameters = Wflow.UnconfinedAquiferParameters(;
            k = fill(conductivity, ncell),
            top = fill(top, ncell),
            bottom = fill(bottom, ncell),
            area = fill(cellsize * cellsize, ncell),
            specific_yield = fill(specific_yield, ncell),
            f = fill(gwf_f, ncell),
        )

        aquifer = Wflow.UnconfinedAquifer(; parameters, variables)
        # constant head on left boundary, 0 at 0
        variables = Wflow.ConstantHeadVariables(; head = [0.0])
        constanthead = Wflow.ConstantHead(; variables, index = [1])
        gwf = Wflow.GroundwaterFlow(;
            aquifer = aquifer,
            connectivity = connectivity,
            constanthead = constanthead,
            boundaries = NamedTuple(),
        )

        dt = Wflow.stable_timestep(gwf.aquifer, conductivity_profile)
        Q = zeros(ncell)
        time = 20.0
        nstep = Int(ceil(time / dt))
        time = nstep * dt

        for i in 1:nstep
            Wflow.update!(gwf, Q, dt, conductivity_profile)
            # Gradient dh/dx is positive, all flow to the left
            @test all(diff(gwf.aquifer.variables.head) .> 0.0)
        end

        head_analytical = [
            transient_aquifer_1d(
                x,
                time,
                conductivity,
                specific_yield,
                aquifer_length,
                beta,
            ) for x in xc
        ]
        difference = gwf.aquifer.variables.head .- head_analytical
        # @test all(difference .< ?)  #TODO
    end

    @testset "integration: unconfined transient 1D, exponential conductivity" begin
        nrow = 1
        ncol = 9
        shape = (nrow, ncol)
        conductivity = 200.0
        top = 150.0
        bottom = 0.0
        specific_yield = 0.15
        cellsize = 500.0
        beta = 1.12
        aquifer_length = cellsize * ncol
        gwf_f = 3.0
        conductivity_profile = "exponential"

        # Domain, geometry
        domain = ones(Bool, shape)
        dx = fill(cellsize, ncol)
        dy = fill(cellsize, nrow)
        indices, reverse_indices = Wflow.active_indices(domain, false)
        connectivity = Wflow.Connectivity(indices, reverse_indices, dx, dy)
        ncell = connectivity.ncell
        xc = collect(range(0.0; stop = aquifer_length - cellsize, step = cellsize))

        variables = Wflow.UnconfinedAquiferVariables(;
            head = initial_head.(xc),
            conductance = fill(0.0, connectivity.nconnection),
            storage = fill(0.0, ncell),
        )
        parameters = Wflow.UnconfinedAquiferParameters(;
            k = fill(conductivity, ncell),
            top = fill(top, ncell),
            bottom = fill(bottom, ncell),
            area = fill(cellsize * cellsize, ncell),
            specific_yield = fill(specific_yield, ncell),
            f = fill(gwf_f, ncell),
        )

        aquifer = Wflow.UnconfinedAquifer(; parameters, variables)
        # constant head on left boundary, 0 at 0
        variables = Wflow.ConstantHeadVariables(; head = [0.0])
        constanthead = Wflow.ConstantHead(; variables, index = [1])
        gwf = Wflow.GroundwaterFlow(;
            aquifer = aquifer,
            connectivity = connectivity,
            constanthead = constanthead,
            boundaries = NamedTuple(),
        )

        dt = Wflow.stable_timestep(gwf.aquifer, conductivity_profile)
        Q = zeros(ncell)
        time = 20.0
        nstep = Int(ceil(time / dt))
        time = nstep * dt

        for i in 1:nstep
            Wflow.update!(gwf, Q, dt, conductivity_profile)
            # Gradient dh/dx is positive, all flow to the left
            @test all(diff(gwf.aquifer.variables.head) .> 0.0)
        end

        head_analytical = [
            transient_aquifer_1d(
                x,
                time,
                conductivity,
                specific_yield,
                aquifer_length,
                beta,
            ) for x in xc
        ]
        difference = gwf.aquifer.variables.head .- head_analytical
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
        conductivity_profile = "uniform"

        # Domain, geometry
        domain = ones(Bool, shape)
        dx = fill(cellsize, ncol)
        dy = fill(cellsize, nrow)
        indices, reverse_indices = Wflow.active_indices(domain, false)
        connectivity = Wflow.Connectivity(indices, reverse_indices, dx, dy)
        ncell = connectivity.ncell

        parameters = Wflow.ConfinedAquiferParameters(;
            k = fill(conductivity, ncell),
            top = fill(top, ncell),
            bottom = fill(bottom, ncell),
            area = fill(cellsize * cellsize, ncell),
            specific_storage = fill(specific_storage, ncell),
            storativity = fill(storativity, ncell),
        )
        variables = Wflow.ConfinedAquiferVariables(;
            head = fill(startinghead, ncell),
            conductance = fill(0.0, connectivity.nconnection),
            storage = fill(0.0, ncell),
        )
        aquifer = Wflow.ConfinedAquifer(; parameters, variables)

        cell_index = reshape(collect(range(1, ncell; step = 1)), shape)
        indices = vcat(cell_index[1, :], cell_index[end, :])# , cell_index[:, 1], cell_index[:, end],)
        variables = Wflow.ConstantHeadVariables(; head = fill(10.0, size(indices)))
        constanthead = Wflow.ConstantHead(; variables, index = indices)
        # Place a well in the middle of the domain
        variables = Wflow.WellVariables(; volumetric_rate = [discharge], flux = [0.0])
        well = Wflow.Well(; variables, index = [reverse_indices[wellrow, wellrow]])
        gwf = Wflow.GroundwaterFlow(;
            aquifer = aquifer,
            connectivity = connectivity,
            constanthead = constanthead,
            boundaries = (; well,),
        )

        dt = Wflow.stable_timestep(gwf.aquifer, conductivity_profile)
        Q = zeros(ncell)
        time = 20.0
        nstep = Int(ceil(time / dt))
        time = nstep * dt

        for i in 1:nstep
            Wflow.update!(gwf, Q, dt, conductivity_profile)
        end

        # test for symmetry on x and y axes
        head = reshape(gwf.aquifer.variables.head, shape)
        @test head[1:halfnrow, :] ≈ head[end:-1:(halfnrow + 2), :]
        @test head[:, 1:halfnrow] ≈ head[:, end:-1:(halfnrow + 2)]

        # compare with analytical solution
        start = -0.5 * aquifer_length + 0.5 * cellsize
        stop = 0.5 * aquifer_length - 0.5 * cellsize
        X = collect(range(start; stop = stop, step = cellsize))
        head_analytical =
            [drawdown_theis(x, time, discharge, transmissivity, storativity) for x in X] .+ 10.0
        # compare left-side, since it's symmetric anyway. Skip the well cell, and its first neighbor
        difference = head[1:(halfnrow - 1), halfnrow] - head_analytical[1:(halfnrow - 1)]
        @test all(difference .< 0.02)
    end
end
