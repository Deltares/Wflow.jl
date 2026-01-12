@testitem "unit: connectivity" begin
    ncol = 2
    nrow = 3
    shape = (ncol, nrow)
    dx = [10.0, 20.0]
    dy = [5.0, 15.0, 25.0]
    collect_connections(con, cell_id) =
        [con.rowval[nzi] for nzi in Wflow.connections(con, cell_id)]

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

@testitem "unit: aquifer, boundary conditions" begin
    include("testing_utils.jl")
    @testset "harmonicmean_conductance" begin
        # harmonicmean_conductance(kH1, kH2, l1, l2, width)
        @test Wflow.harmonicmean_conductance(10.0 * 5.0, 10.0 * 5.0, 0.5, 0.5, 1.0) == 50.0
        @test Wflow.harmonicmean_conductance(10.0 * 0.0, 10.0 * 5.0, 0.5, 0.5, 1.0) == 0.0
        @test Wflow.harmonicmean_conductance(10.0 * 5.0, 10.0 * 0.0, 0.5, 0.5, 1.0) == 0.0
        # kD of 10 and 20 -> harmonicmean = 1/(1/10 + 1/20)
        @test Wflow.harmonicmean_conductance(10.0 * 1.0, 10.0 * 2.0, 1.0, 1.0, 1.0) ≈
              (6.0 + 2.0 / 3.0)
    end

    nrow = 1
    ncol = 3
    connectivity, unconf_aqf = homogenous_aquifer(nrow, ncol)
    Wflow.initialize_conductance!(unconf_aqf, connectivity)
    ncell = connectivity.ncell

    @testset "saturated_thickness-unconfined" begin
        @test Wflow.saturated_thickness(unconf_aqf, 1) == 0.0
        @test Wflow.saturated_thickness(unconf_aqf, 2) == 7.5
        @test Wflow.saturated_thickness(unconf_aqf, 3) == 10.0
    end

    @testset "conductance" begin
        conductivity_profile = Wflow.GwfConductivityProfileType.uniform
        @test Wflow.conductance(unconf_aqf, 2, 3, 3, conductivity_profile, connectivity) ==
              100.0  # upstream sat. thickness
        @test Wflow.conductance(unconf_aqf, 1, 2, 1, conductivity_profile, connectivity) ==
              75.0  # upstream sat. thickness
    end

    @testset "minimum_head-unconfined" begin
        original_head = copy(unconf_aqf.variables.head)
        unconf_aqf.variables.head[1] = -10.0
        @test Wflow.check_flux(-1.0, unconf_aqf, 1) == 0.0
        @test Wflow.minimum_head(unconf_aqf)[1] == 0.0
        unconf_aqf.variables.head .= original_head
    end

    @testset "stable_timestep" begin
        conductivity_profile = Wflow.GwfConductivityProfileType.uniform
        cfl = 0.25
        @test Wflow.stable_timestep(unconf_aqf, conductivity_profile, cfl) == 0.0375
    end

    # Parametrization in setup is as follows:
    # [0.0, 7.5, 20.0],  # head
    # fill(10.0, ncell),  # k
    # fill(10.0, ncell),  # top
    # fill(0.0, ncell),  # bottom

    @testset "flux-unconfined" begin
        dt = 1.0
        unconf_aqf.variables.q_net .= 0.0
        conductivity_profile = Wflow.GwfConductivityProfileType.uniform
        Wflow.flux!(unconf_aqf, connectivity, conductivity_profile, dt)
        # KD is based on upstream saturated thickness, i.e. 7.5 m and 20.0 m (which is capped to 10.0)
        @test unconf_aqf.variables.q_net == [562.5, 687.5, -1250.0]
    end

    @testset "river" begin
        dt = 1.0
        n = 2
        parameters = Wflow.GwfRiverParameters(;
            infiltration_conductance = [100.0, 100.0],
            exfiltration_conductance = [200.0, 200.0],
            bottom = [1.0, 1.0],
        )
        variables = Wflow.GwfRiverVariables(;
            n,
            stage = [2.0, 2.0],
            storage = [20.0, 20.0],
            flux = [0.0, 0.0],
            flux_av = [0.0, 0.0],
        )
        river = Wflow.GwfRiver(; parameters, variables, index = [1, 3])
        unconf_aqf.variables.q_net .= 0.0
        Wflow.flux!(river, unconf_aqf, dt)
        # infiltration, below bottom, flux is (stage - bottom) * inf_cond, limited by
        # river storage (20.0)
        @test unconf_aqf.variables.q_net[1] == 20.0
        # drainage, flux is (stage - head) * exf_cond
        @test unconf_aqf.variables.q_net[3] == (2.0 - 20.0) * 200.0
    end

    @testset "drainage" begin
        dt = 1.0
        n = 2
        parameters =
            Wflow.DrainageParameters(; elevation = [2.0, 2.0], conductance = [100.0, 100.0])
        variables = Wflow.DrainageVariables(; n, flux = [0.0, 0.0], flux_av = [0.0, 0.0])
        drainage = Wflow.Drainage(; parameters, variables, index = [1, 2])
        unconf_aqf.variables.q_net .= 0.0
        Wflow.flux!(drainage, unconf_aqf, dt)
        @test unconf_aqf.variables.q_net[1] == 0.0
        @test unconf_aqf.variables.q_net[2] == 100.0 * (2.0 - 7.5)
    end

    @testset "headboundary" begin
        dt = 1.0
        parameters = Wflow.HeadBoundaryParameters(; conductance = [100.0, 100.0])
        variables = Wflow.HeadBoundaryVariables(;
            head = [2.0, 2.0],
            flux = [0.0, 0.0],
            flux_av = [0.0, 0.0],
        )

        headboundary = Wflow.HeadBoundary(; parameters, variables, index = [1, 2])
        unconf_aqf.variables.q_net .= 0.0
        Wflow.flux!(headboundary, unconf_aqf, dt)
        @test unconf_aqf.variables.q_net[1] == 100.0 * (2.0 - 0.0)
        @test unconf_aqf.variables.q_net[2] == 100.0 * (2.0 - 7.5)
    end

    @testset "recharge" begin
        dt = 1.0
        n = 3
        variables = Wflow.RechargeVariables(;
            n,
            rate = [1.0e-3, 1.0e-3, 1.0e-3],
            flux = [0.0, 0.0, 0.0],
            flux_av = [0.0, 0.0, 0.0],
        )
        recharge = Wflow.Recharge(; n, variables, index = [1, 2, 3])
        unconf_aqf.variables.q_net .= 0.0
        Wflow.flux!(recharge, unconf_aqf, dt)
        @test all(unconf_aqf.variables.q_net .== 1.0e-3 * 100.0)
    end

    @testset "well" begin
        dt = 1.0
        variables = Wflow.WellVariables(;
            volumetric_rate = [-1000.0],
            flux = [0.0],
            flux_av = [0.0],
        )
        well = Wflow.Well(; variables, index = [2])
        unconf_aqf.variables.q_net .= 0.0
        Wflow.flux!(well, unconf_aqf, dt)
        @test unconf_aqf.variables.q_net[2] == -1000.0
    end
end

@testitem "integration: unconfined transient 1D" begin
    include("testing_utils.jl")

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
    conductivity_profile = Wflow.GwfConductivityProfileType.uniform

    # Domain, geometry
    domain = ones(Bool, shape)
    dx = fill(cellsize, ncol)
    dy = fill(cellsize, nrow)
    indices, reverse_indices = Wflow.active_indices(domain, false)
    connectivity = Wflow.Connectivity(indices, reverse_indices, dx, dy)
    ncell = connectivity.ncell
    xc = collect(range(0.0; stop = aquifer_length - cellsize, step = cellsize))

    variables = Wflow.AquiferVariables(;
        n = ncell,
        head = initial_head.(xc),
        conductance = fill(0.0, connectivity.nconnection),
        storage = fill(0.0, ncell),
        q_net = fill(0.0, ncell),
        q_in_av = fill(0.0, ncell),
        q_out_av = fill(0.0, ncell),
        exfiltwater = fill(0.0, ncell),
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
    timestepping = Wflow.TimeStepping(; cfl = 0.25)
    gwf = Wflow.GroundwaterFlow(; timestepping, aquifer, connectivity, constanthead)

    time = 20.0
    t = 0.0
    (; cfl) = gwf.timestepping
    while t < time
        global t
        gwf.aquifer.variables.q_net .= 0.0
        dt_s = Wflow.stable_timestep(gwf.aquifer, conductivity_profile, cfl)
        dt_s = Wflow.check_timestepsize(dt_s, t, time)
        Wflow.update_fluxes!(gwf, conductivity_profile, dt_s)
        Wflow.update_head!(gwf, dt_s)
        t = t + dt_s
        # Gradient dh/dx is positive, all flow to the left
        @test all(diff(gwf.aquifer.variables.head) .> 0.0)
    end

    head_analytical = [
        transient_aquifer_1d(x, time, conductivity, specific_yield, aquifer_length, beta) for x in xc
    ]
    difference = gwf.aquifer.variables.head .- head_analytical
    # @test all(difference .< ?)  #TODO
end

@testitem "integration: unconfined transient 1D, exponential conductivity" begin
    include("testing_utils.jl")
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
    conductivity_profile = Wflow.GwfConductivityProfileType.exponential

    # Domain, geometry
    domain = ones(Bool, shape)
    dx = fill(cellsize, ncol)
    dy = fill(cellsize, nrow)
    indices, reverse_indices = Wflow.active_indices(domain, false)
    connectivity = Wflow.Connectivity(indices, reverse_indices, dx, dy)
    ncell = connectivity.ncell
    xc = collect(range(0.0; stop = aquifer_length - cellsize, step = cellsize))

    variables = Wflow.AquiferVariables(;
        n = ncell,
        head = initial_head.(xc),
        conductance = fill(0.0, connectivity.nconnection),
        storage = fill(0.0, ncell),
        q_net = fill(0.0, ncell),
        q_in_av = fill(0.0, ncell),
        q_out_av = fill(0.0, ncell),
        exfiltwater = fill(0.0, ncell),
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
    timestepping = Wflow.TimeStepping(; cfl = 0.25)
    gwf = Wflow.GroundwaterFlow(; timestepping, aquifer, connectivity, constanthead)

    time = 20.0
    t = 0.0
    (; cfl) = gwf.timestepping
    while t < time
        global t
        gwf.aquifer.variables.q_net .= 0.0
        dt_s = Wflow.stable_timestep(gwf.aquifer, conductivity_profile, cfl)
        dt_s = Wflow.check_timestepsize(dt_s, t, time)
        Wflow.update_fluxes!(gwf, conductivity_profile, dt_s)
        Wflow.update_head!(gwf, dt_s)
        t = t + dt_s
        # Gradient dh/dx is positive, all flow to the left
        @test all(diff(gwf.aquifer.variables.head) .> 0.0)
    end

    head_analytical = [
        transient_aquifer_1d(x, time, conductivity, specific_yield, aquifer_length, beta) for x in xc
    ]
    difference = gwf.aquifer.variables.head .- head_analytical
    # @test all(difference .< ?)  #TODO
end
