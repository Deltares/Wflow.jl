@testitem "Run SBM" begin
    using Dates
    using Statistics: mean
    include("testing_utils.jl")

    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    model = Wflow.Model(config)
    (; domain) = model

    Wflow.run_timestep!(model)

    # test if the first timestep was written to the CSV file
    flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
    @testset "CSV output" begin
        row = csv_first_row(model.writer.csv_path)

        @test row.time == DateTime("2000-01-02T00:00:00")
        @test row.Q ≈ 6.850474762647178
        @test row.storage ≈ 2.753500258633802e7
        @test row.temp_bycoord ≈ 2.390000104904175
        @test row.vwc_layer2_bycoord ≈ 0.25938707367146907
        @test row.temp_byindex ≈ 2.390000104904175
        @test row.Q_6336050 ≈ 0.003175570694335172
        @test row.Q_6336510 ≈ 0.012785802015566708
        @test row.Q_6836100 ≈ 0.008740434542207666
        @test row.Q_6336500 ≈ 0.00240582054420714
        @test row.Q_6836190 ≈ 0.002534862969007946
        @test row.Q_6336800 ≈ 0.004367717681790628
        @test row.Q_6336900 ≈ 0.003934027688070696
        @test row.Q_6336930 ≈ 0.013810820648092513
        @test row.Q_6336910 ≈ 0.003073485680460027
        @test row.Q_6136500 ≈ 0.00043538067672607716
        @test row.Q_6136520 ≈ 0.0004119095811616616
        @test row.Q_6136150 ≈ 0.0029657574491729643
        @test row.Q_6136151 ≈ 0.0027844732219454186
        @test row.Q_6136160 ≈ 3.472843959583444
        @test row.Q_6136202 ≈ 1.306258882043624
        @test row.recharge_1 ≈ -0.002257181032501202
    end

    @testset "NetCDF scalar output" begin
        ds = model.writer.dataset_scalar
        @test ds["time"][1] == DateTime("2000-01-02T00:00:00")
        @test ds["Q"][:][1:20] ≈ [
            0.7411495,
            1.7804804,
            1.7792363,
            1.5534992,
            13.25166,
            4.7680607,
            0.09143825,
            4.354929,
            0.010488547,
            6.0434036,
            6.0332,
            0.0059943455,
            0.00660558,
            0.006399016,
            0.0041626506,
            1.2043761,
            0.0022155696,
            1.7879109,
            8.354377,
            2.735113,
        ]
        @test ds["river_gauge__count"].attrib["cf_role"] == "timeseries_id"
        @test ds["temp_index"][:] ≈ [2.39]
        @test ds["temp_coord"][:] ≈ [2.39]
        @test keys(ds.dim) ==
              ["time", "layer", "river_gauge__count", "temp_bycoord", "temp_byindex"]
    end

    @testset "first timestep" begin
        sbm = model.land.soil
        snow = model.land.snow
        @test snow.parameters.tt[50063] ≈ 0.0

        @test model.clock.iteration == 1

        @test sbm.parameters.theta_s[50063] ≈ 0.48755401372909546
        @test sbm.parameters.theta_r[50063] ≈ 0.15943120419979095
        @test mean(sbm.variables.runoff) ≈ 0.04033701711287029
        @test mean(sbm.variables.soilevap) ≈ 0.02122698830889417
        @test mean(sbm.variables.actevap) ≈ 0.33545623834952143
        @test mean(sbm.variables.actinfilt) ≈ 1.6444774688444848
        @test snow.variables.snow_storage[5] ≈ 3.768513390588815
        @test mean(snow.variables.snow_storage) ≈ 0.038019723676094325
        @test sbm.variables.total_storage[50063] ≈ 559.886196102425
        @test sbm.variables.total_storage[429] ≈ 599.6375723706255 # river cell
    end

    # run the second timestep
    Wflow.run_timestep!(model)

    @testset "second timestep" begin
        sbm = model.land.soil
        snow = model.land.snow
        @test sbm.parameters.theta_s[50063] ≈ 0.48755401372909546
        @test sbm.parameters.theta_r[50063] ≈ 0.15943120419979095
        @test mean(sbm.variables.net_runoff) ≈ 0.22092252199097845
        @test mean(sbm.variables.runoff) ≈ 0.22129098393293603
        @test mean(sbm.variables.soilevap) ≈ 0.018738226205510525
        @test mean(sbm.variables.actevap) ≈ 0.2744312454440506
        @test mean(sbm.variables.actinfilt) ≈ 0.08861928614420159
        @test snow.variables.snow_storage[5] ≈ 3.843412524052313
        @test mean(snow.variables.snow_storage) ≈ 0.03461317061870949
        @test sbm.variables.total_storage[50063] ≈ 559.9633189802773
        @test sbm.variables.total_storage[429] ≈ 620.6174448020174  # river cell
    end

    @testset "subsurface flow" begin
        ssf = model.routing.subsurface_flow.variables.ssf
        @test sum(ssf) ≈ 6.250079949202134e7
        @test ssf[domain.land.network.order[1]] ≈ 699.3636285243076
        @test ssf[domain.land.network.order[end - 100]] ≈ 2395.6159482448143
        @test ssf[domain.land.network.order[end]] ≈ 287.61501877867994
    end

    @testset "overland flow" begin
        q = model.routing.overland_flow.variables.q_av
        @test sum(q) ≈ 264.6944359620204
        @test q[26625] ≈ 0.0
        @test q[39308] ≈ 0.0
        @test q[domain.land.network.order[end]] ≈ 1.0e-30
    end

    @testset "river flow" begin
        q = model.routing.river_flow.variables.q_av
        @test sum(q) ≈ 3696.5789619579173
        @test q[1622] ≈ 0.0007502850311515928
        @test q[43] ≈ 11.458528971675033
        @test q[domain.river.network.order[end]] ≈ 0.04397648668235761
    end

    @testset "reservoir simple" begin
        res = model.routing.river_flow.boundary_conditions.reservoir
        @test res.variables.outflow[1] ≈ 0.2174998614438593
        @test res.variables.outflow_av[1] ≈ 0.21749986282401396
        @test res.boundary_conditions.inflow[1] ≈ 0.0005130607130643568
        @test res.variables.storage[1] ≈ 2.751299001489657f7
        @test res.variables.actevap[1] ≈ 0.5400000810623169
        @test res.boundary_conditions.precipitation[1] ≈ 0.17999997735023499
        @test res.boundary_conditions.evaporation[1] ≈ 0.5400000810623169
    end

    Wflow.close_files(model; delete_output = false)

    # test without lateral snow transport
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.snow_gravitational_transport__flag = false

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)

    @testset "lateral snow transport off" begin
        snow = model.land.snow
        @test snow.variables.snow_storage[5] ≈ 3.7686103651001375
        @test mean(snow.variables.snow_storage) ≈ 0.03801972367609432
        @test mean(snow.variables.snow_water) ≈ 0.0025756728488273866
        @test mean(snow.variables.swe) ≈ 0.0405953965249217
    end

    # test without snow model
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.snow__flag = false
    pop!(config.output.netcdf_grid.variables, "snowpack_dry_snow__leq_depth")
    pop!(config.output.netcdf_grid.variables, "snowpack_liquid_water__depth")
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)

    @testset "snow model not included" begin
        snow = model.land.snow
        @test typeof(model.land.snow) == Wflow.NoSnowModel
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "Set pit and multithreading multiple basins" begin
    using Dates
    # test for setting a pit and multithreading multiple basins (by setting 2 extra pits
    # resulting in 3 basins)
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.pit__flag = true
    config.input.basin_pit_location__mask = "wflow_pits"
    config.time.endtime = DateTime(2000, 1, 9)
    config.logging.loglevel = "info"

    model = Wflow.run(config)

    @testset "timing" begin
        # clock has been reset
        @test model.clock.time == Wflow.cftime(config.time.starttime, config.time.calendar)
        @test model.clock.iteration == 0
    end

    @testset "river flow at basin outlets and downstream of one pit" begin
        q = model.routing.river_flow.variables.q_av
        @test q[4009] ≈ 8.426694842173548 # pit/ outlet, CartesianIndex(141, 228)
        @test q[4020] ≈ 0.006370691658310787 # downstream of pit 4009, CartesianIndex(141, 229)
        @test q[2508] ≈ 131.40631419288573 # pit/ outlet
        @test q[5808] ≈ 0.11941498848157915  # pit/ outlet
    end
end

@testitem "Changing forcing and cyclic LAI parameter" begin
    # Run unchanged
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)

    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    precip = copy(model.land.atmospheric_forcing.precipitation)
    evap = copy(model.land.atmospheric_forcing.potential_evaporation)
    lai = copy(model.land.vegetation_parameters.leaf_area_index)
    res_evap = copy(
        model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.evaporation,
    )

    # Run changed
    config.input.forcing["atmosphere_water__precipitation_volume_flux"] =
        Wflow.init_config_section(
            Wflow.InputEntry,
            Dict("scale" => 2.0, "netcdf_variable_name" => "precip"),
        )
    config.input.forcing["land_surface_water__potential_evaporation_volume_flux"] =
        Wflow.init_config_section(
            Wflow.InputEntry,
            Dict("scale" => 3.0, "offset" => 1.50, "netcdf_variable_name" => "pet"),
        )
    config.input.cyclic["vegetation__leaf_area_index"] = Wflow.init_config_section(
        Wflow.InputEntry,
        Dict("scale" => 1.6, "netcdf_variable_name" => "LAI"),
    )
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "changed dynamic parameters" begin
        res = model.routing.river_flow.boundary_conditions.reservoir
        land = model.land
        @test land.atmospheric_forcing.precipitation[2] / precip[2] ≈ 2.0f0
        @test (land.atmospheric_forcing.potential_evaporation[100] - 1.50) / evap[100] ≈
              3.0f0
        @test land.vegetation_parameters.leaf_area_index[100] / lai[100] ≈ 1.6f0
        @test (res.boundary_conditions.evaporation[2] - 1.50) / res_evap[2] ≈ 3.0f0
    end
end

@testitem "Cyclic river and reservoir external inflow (kinematic wave routing)" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    config.input.cyclic["river_water__external_inflow_volume_flow_rate"] = "inflow"
    config.input.cyclic["reservoir_water__external_inflow_volume_flow_rate"] = "reservoir_inflow"

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "kinematic wave routing: river and reservoir external inflow (cyclic)" begin
        (; reservoir) = model.routing.river_flow.boundary_conditions
        @test model.routing.river_flow.boundary_conditions.external_inflow[44] ≈ 0.75
        @test model.routing.river_flow.boundary_conditions.actual_external_abstraction_av[44] ==
              0.0
        @test model.routing.river_flow.variables.q_av[44] ≈ 10.156053077341845
        @test reservoir.boundary_conditions.external_inflow[2] == -1.0
        @test reservoir.boundary_conditions.actual_external_abstraction_av[2] == 1.0
        @test reservoir.boundary_conditions.inflow[2] ≈ -0.9054049318713108
        @test reservoir.variables.outflow_av[2] ≈ 3.000999922024245
    end
end

@testitem "Cyclic river and reservoir external inflow (local inertial routing)" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    config.input.cyclic["river_water__external_inflow_volume_flow_rate"] = "inflow"
    config.input.cyclic["reservoir_water__external_inflow_volume_flow_rate"] = "reservoir_inflow"

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "local inertial routing: river and reservoir external inflow (cyclic)" begin
        (; reservoir) = model.routing.river_flow.boundary_conditions
        @test model.routing.river_flow.boundary_conditions.external_inflow[44] ≈ 0.75
        @test model.routing.river_flow.boundary_conditions.actual_external_abstraction_av[44] ==
              0.0
        @test model.routing.river_flow.variables.q_av[44] ≈ 10.119602100591411
        @test reservoir.boundary_conditions.external_inflow[2] == -1.0
        @test reservoir.boundary_conditions.actual_external_abstraction_av[2] == 1.0
        @test reservoir.boundary_conditions.inflow[2] ≈ -0.9090882985601229
        @test reservoir.variables.outflow_av[2] ≈ 3.000999922022744
    end
end

@testitem "External negative inflow" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    config.dir_output = mktempdir()
    model.routing.river_flow.boundary_conditions.external_inflow[44] = -10.0
    (; actual_external_abstraction_av, external_inflow) =
        model.routing.river_flow.boundary_conditions
    (; q_av) = model.routing.river_flow.variables
    @testset "river external negative inflow" begin
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_av[44] ≈ 1.5965227273142157
        @test q_av[44] ≈ 1.4328200140906533
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_av[44] ≈ 5.416747895349456
        @test q_av[44] ≈ 4.000700150630984
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_av[44] ≈ 9.756847289612736
        @test external_inflow[44] == -10.0
        @test q_av[44] ≈ 7.275452569433593
    end
end

@testitem "Fixed forcing (precipitation = 2.5)" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.input.forcing["atmosphere_water__precipitation_volume_flux"] = 2.5
    model = Wflow.Model(config)
    Wflow.load_fixed_forcing!(model)

    @testset "fixed precipitation forcing (initialize)" begin
        @test maximum(model.land.atmospheric_forcing.precipitation) ≈ 2.5
        @test minimum(model.land.atmospheric_forcing.precipitation) ≈ 0.0
        @test all(
            isapprox.(
                model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation,
                2.5,
            ),
        )
    end

    Wflow.run_timestep!(model)

    @testset "fixed precipitation forcing (first timestep)" begin
        @test maximum(model.land.atmospheric_forcing.precipitation) ≈ 2.5
        @test minimum(model.land.atmospheric_forcing.precipitation) ≈ 0.0
        @test all(
            isapprox.(
                model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.precipitation,
                2.5,
            ),
        )
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "Local-inertial option for river flow river_routing" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "river flow and depth (local inertial)" begin
        q = model.routing.river_flow.variables.q_av
        @test sum(q) ≈ 3692.1984174051863
        @test q[1622] ≈ 7.260572774917902e-5
        @test q[43] ≈ 11.248731903973153
        @test q[501] ≈ 3.163313530437266
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.00191277920611667
        @test h[43] ≈ 0.447403557220214
        @test h[501] ≈ 0.37810388857945015
        q_channel = model.routing.river_flow.variables.q_channel_av
        @test q ≈ q_channel
    end
end

@testitem "External negative inflow local-inertial river flow" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)
    model.routing.river_flow.boundary_conditions.external_inflow[44] = -10.0
    (; actual_external_abstraction_av, external_inflow, reservoir) =
        model.routing.river_flow.boundary_conditions
    (; q_av) = model.routing.river_flow.variables
    @testset "river external negative inflow (local inertial)" begin
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_av[44] ≈ 3.0038267868094692
        @test q_av[44] ≈ 0.0007694313441619584
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_av[44] ≈ 9.420226369186521
        @test q_av[44] ≈ 0.009169625606576883
        Wflow.run_timestep!(model)
        @test actual_external_abstraction_av[44] ≈ 9.999999999999991
        @test external_inflow[44] == -10.0
        @test q_av[44] ≈ 6.917023458482667
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "Local-inertial option for river and overland flow" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-land-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()

    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "river and overland flow and depth (local inertial)" begin
        q = model.routing.river_flow.variables.q_av
        @test sum(q) ≈ 2415.8565080484427
        @test q[1622] ≈ 7.289953869913158e-5
        @test q[43] ≈ 5.235679976913814
        @test q[501] ≈ 1.395242045016401
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.001915431360203216
        @test h[43] ≈ 0.30691170109894245
        @test h[501] ≈ 0.3039452016829932
        qx = model.routing.overland_flow.variables.qx
        qy = model.routing.overland_flow.variables.qy
        @test qx[[26, 35, 631]] ≈ [0.17466039055941732, 0.002077783779618663, 0.0]
        @test qy[[26, 35, 631]] ≈ [0.12287975269080074, 0.019011432839373375, 0.0]
        h = model.routing.overland_flow.variables.h
        @test h[[26, 35, 631]] ≈ [0.07163764409112827, 0.009015112665828464, 0.0]
    end

    Wflow.close_files(model; delete_output = false)
end

@testitem "Local-inertial option for river flow including 1D floodplain schematization" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-floodplain-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    model = Wflow.Model(config)

    (; flow_length, flow_length) = model.domain.river.parameters
    fp = model.routing.river_flow.floodplain.parameters.profile
    river = model.routing.river_flow
    dh = diff(fp.depth)
    Δv = diff(fp.storage[:, 3])
    Δa = diff(fp.a[:, 3])

    @testset "river flow (local inertial) floodplain schematization" begin
        # floodplain geometry checks (index 3)
        @test fp.storage[:, 3] ≈ [0.0, 8641.0, 19011.0, 31685.0, 51848.0, 80653.0]
        @test fp.width[:, 3] ≈ [
            30.0,
            99.28617594254938,
            119.15260323159785,
            145.6258527827648,
            231.6754039497307,
            330.9730700179533,
        ]
        @test fp.p[:, 3] ≈ [
            69.28617594254938,
            70.28617594254938,
            91.15260323159785,
            118.62585278276481,
            205.6754039497307,
            305.9730700179533,
        ]
        @test fp.a[:, 3] ≈ [
            0.0,
            49.64308797127469,
            109.21938958707361,
            182.032315978456,
            297.8700179533214,
            463.35655296229805,
        ]
        @test dh .* fp.width[2:end, 3] * flow_length[3] ≈ Δv
        @test fp.a[:, 3] * flow_length[3] ≈ fp.storage[:, 3]
        # flood depth from flood storage (8000.0)
        flood_vol = 8000.0
        river.variables.storage[3] = flood_vol + river.parameters.bankfull_storage[3]
        i1, i2 = Wflow.interpolation_indices(flood_vol, fp.storage[:, 3])
        @test (i1, i2) == (1, 2)
        flood_depth = Wflow.flood_depth(fp, flood_vol, flow_length[3], 3)
        @test flood_depth ≈ 0.46290938548779076
        @test (flood_depth - fp.depth[i1]) * fp.width[i2, 3] * flow_length[3] +
              fp.storage[i1, 3] ≈ flood_vol
        # flood depth from flood storage (12000.0)
        flood_vol = 12000.0
        river.variables.storage[3] = flood_vol + river.parameters.bankfull_storage[3]
        i1, i2 = Wflow.interpolation_indices(flood_vol, fp.storage[:, 3])
        @test (i1, i2) == (2, 3)
        flood_depth = Wflow.flood_depth(fp, flood_vol, flow_length[3], 3)
        @test flood_depth ≈ 0.6619575699132112
        @test (flood_depth - fp.depth[i1]) * fp.width[i2, 3] * flow_length[3] +
              fp.storage[i1, 3] ≈ flood_vol
        # test extrapolation of segment
        flood_vol = 95000.0
        river.variables.storage[3] = flood_vol + river.parameters.bankfull_storage[3]
        i1, i2 = Wflow.interpolation_indices(flood_vol, fp.storage[:, 3])
        @test (i1, i2) == (6, 6)
        flood_depth = Wflow.flood_depth(fp, flood_vol, flow_length[3], 3)
        @test flood_depth ≈ 2.749036625585836
        @test (flood_depth - fp.depth[i1]) * fp.width[i2, 3] * flow_length[3] +
              fp.storage[i1, 3] ≈ flood_vol
        river.variables.storage[3] = 0.0 # reset storage
        # flow area and wetted perimeter based on hf
        h = 0.5
        i1, i2 = Wflow.interpolation_indices(h, fp.depth)
        @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
              49.64308797127469
        @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 70.28617594254938
        h = 1.5
        i1, i2 = Wflow.interpolation_indices(h, fp.depth)
        @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
              182.032315978456
        @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 118.62585278276481
        h = 1.7
        i1, i2 = Wflow.interpolation_indices(h, fp.depth)
        @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
              228.36739676840216
        @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 119.02585278276482
        h = 3.2
        i1, i2 = Wflow.interpolation_indices(h, fp.depth)
        @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
              695.0377019748654
        @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 307.3730700179533
        h = 4.0
        i1, i2 = Wflow.interpolation_indices(h, fp.depth)
        @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
              959.816157989228
        @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 308.9730700179533
        @test Wflow.flow_area(fp.width[i2, 4], fp.a[i1, 4], fp.depth[i1], h) ≈
              407.6395313908081
        @test Wflow.wetted_perimeter(fp.p[i1, 4], fp.depth[i1], h) ≈ 90.11775307900271
    end

    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "river flow (local inertial) with floodplain schematization simulation" begin
        q = model.routing.river_flow.variables.q_av
        @test sum(q) ≈ 3682.763689938466
        @test q[1622] ≈ 7.260572782240141e-5
        @test q[43] ≈ 11.248731903973168
        @test q[501] ≈ 3.121950941199987
        @test q[5808] ≈ 0.0022169051490145645
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.0019127792060060795
        @test h[43] ≈ 0.44740355722021324
        @test h[501] ≈ 0.37317174589196656
        @test h[5808] ≈ 0.0073332218703272825
    end

    # set boundary condition local inertial routing from netCDF file
    config.input.static["model_boundary_condition_river__length"] = "riverlength_bc"
    config.input.static["model_boundary_condition_river_bank_water__depth"] = "riverdepth_bc"
    model = Wflow.Model(config)
    Wflow.run_timestep!(model)
    Wflow.run_timestep!(model)

    @testset "change boundary condition for local inertial routing (including floodplain)" begin
        q = model.routing.river_flow.variables.q_av
        @test sum(q) ≈ 3682.950549280907
        @test q[1622] ≈ 7.260572782240141e-5
        @test q[43] ≈ 11.248731903973168
        @test q[501] ≈ 3.121950941199987
        @test q[5808] ≈ 0.05466707293198311
        h = model.routing.river_flow.variables.h
        @test h[1622] ≈ 0.0019127792060060795
        @test h[43] ≈ 0.44740355722021324
        @test h[501] ≈ 0.37317174589196656
        @test h[5808] ≈ 2.0000266345482336
    end
    Wflow.close_files(model; delete_output = false)

    # test different ksat profiles
    @testset "ksat profiles (SBM)" begin
        tomlpath = joinpath(@__DIR__, "sbm_config.toml")
        function get_config(profile)
            config = Wflow.Config(tomlpath)
            config.dir_output = mktempdir()
            config.model.saturated_hydraulic_conductivity_profile = profile
            config.input.static["soil_layer_water__vertical_saturated_hydraulic_conductivity"] = "kv"
            config.input.static["soil_exponential_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth"] =
                400.0
            config.input.static["soil_layered_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth"] =
                400.0
            config
        end

        i = 100

        @testset "exponential profile" begin
            config = get_config(Wflow.VerticalConductivityProfile.exponential)
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.zi[i]
            kvfrac = soil.parameters.kvfrac
            kv_z = Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2)
            @test kv_z ≈ kvfrac[i][2] * kv_profile.kv_0[i] * exp(-kv_profile.f[i] * z)
            @test subsurface_flow.variables.ssfmax[i] ≈ 28.32720603576582
            @test subsurface_flow.variables.ssf[i] ≈ 11683.330684556406
        end

        @testset "exponential constant profile" begin
            config = get_config(Wflow.VerticalConductivityProfile.exponential_constant)
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.zi[i]
            kvfrac = soil.parameters.kvfrac
            kv_z = Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2)
            @test kv_z ≈
                  kvfrac[i][2] *
                  kv_profile.exponential.kv_0[i] *
                  exp(-kv_profile.exponential.f[i] * z)
            kv_400 = Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, 400.0, i, 2)
            kv_1000 =
                Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, 1000.0, i, 3)
            @test kv_400 ≈ kv_1000
            @test all(kv_profile.z_exp .== 400.0)
            @test subsurface_flow.variables.ssfmax[i] ≈ 49.38558575188426
            @test subsurface_flow.variables.ssf[i] ≈ 24810.460986497365
        end

        @testset "layered profile" begin
            config = get_config(Wflow.VerticalConductivityProfile.layered)
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.zi[i]
            kvfrac = soil.parameters.kvfrac
            @test Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2) ≈
                  kv_profile.kv[i][2]
            Wflow.kh_layered_profile!(soil, subsurface_flow, kv_profile, 86400.0)
            @test subsurface_flow.parameters.kh_profile.kh[i] ≈ 47.508932674632355
            @test subsurface_flow.variables.ssfmax[i] ≈ 30.237094380100316
            @test subsurface_flow.variables.ssf[i] ≈ 14546.518932613191
        end

        config = get_config(Wflow.VerticalConductivityProfile.layered_exponential)

        @testset "layered exponential profile" begin
            model = Wflow.Model(config)
            (; soil) = model.land
            (; kv_profile) = soil.parameters
            (; subsurface_flow) = model.routing
            z = soil.variables.zi[i]
            kvfrac = soil.parameters.kvfrac
            @test Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2) ≈
                  kv_profile.kv[i][2]
            @test kv_profile.nlayers_kv[i] == 2
            Wflow.kh_layered_profile!(soil, subsurface_flow, kv_profile, 86400.0)
            @test subsurface_flow.parameters.kh_profile.kh[i] ≈ 33.76026208801769
            @test all(kv_profile.z_layered[1:10] .== 400.0)
            @test subsurface_flow.variables.ssfmax[i] ≈ 23.4840490395906
            @test subsurface_flow.variables.ssf[i] ≈ 10336.88327617503
        end

        @testset "river flow layered exponential profile" begin
            config.dir_output = mktempdir()
            model = Wflow.Model(config)
            Wflow.run_timestep!(model)
            Wflow.run_timestep!(model)
            q = model.routing.river_flow.variables.q_av
            @test sum(q) ≈ 3032.617045063436
            @test q[1622] ≈ 0.0006987386860043929
            @test q[43] ≈ 8.710529495056752
        end

        Wflow.close_files(model; delete_output = false)
    end
end

@testitem "run wflow sbm" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.time.endtime = "2000-01-05"
    Wflow.run(config)
end

@testitem "water balance sbm (kinematic wave routing)" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; land_water_balance, routing) = model.mass_balance
    (; overland_water_balance, river_water_balance, subsurface_water_balance) = routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 1e-9, land_water_balance.error)
        @test all(re -> abs(re) < 1e-9, land_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, overland_water_balance.error)
        @test all(re -> abs(re) < 6.6e11, overland_water_balance.relative_error)
        inds = findall(x -> x > 1e-3, model.routing.overland_flow.variables.q_av)
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error[inds])
        @test all(e -> abs(e) < 1.e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
        @test all(e -> abs(e) < 1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, land_water_balance.error)
        @test all(re -> abs(re) < 1e-9, land_water_balance.relative_error)
        @test all(e -> abs(e) < 1.e-9, routing.overland_water_balance.error)
        @test all(re -> abs(re) < 5.4e11, routing.overland_water_balance.relative_error)
        inds = findall(x -> x > 1e-3, model.routing.overland_flow.variables.q_av)
        @test all(re -> abs(re) < 1e-9, routing.overland_water_balance.relative_error[inds])
        @test all(e -> abs(e) < 3e-5, river_water_balance.error)
        @test all(re -> abs(re) < 12.2, river_water_balance.relative_error)
        inds = findall(x -> x > 1e-3, model.routing.river_flow.variables.q_av)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error[inds])
        @test all(e -> abs(e) < 1e-9, subsurface_water_balance.error)
        @test all(re -> abs(re) < 1e-9, subsurface_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "water balance river local inertial routing" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; river_water_balance) = model.mass_balance.routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "water balance river local inertial routing with floodplain" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-floodplain-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; river_water_balance) = model.mass_balance.routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, river_water_balance.error)
        @test all(re -> abs(re) < 1e-9, river_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "water balance river and land local inertial routing" begin
    tomlpath = joinpath(@__DIR__, "sbm_river-land-local-inertial_config.toml")
    config = Wflow.Config(tomlpath)
    config.dir_output = mktempdir()
    config.model.water_mass_balance__flag = true
    model = Wflow.Model(config)
    (; overland_water_balance) = model.mass_balance.routing
    Wflow.run_timestep!(model)
    @testset "water balance first timestep" begin
        @test all(e -> abs(e) < 1e-9, overland_water_balance.error)
        @test all(re -> abs(re) < 1e-9, overland_water_balance.relative_error)
    end
    Wflow.run_timestep!(model)
    @testset "water balance second timestep" begin
        @test all(e -> abs(e) < 1e-9, overland_water_balance.error)
        @test all(re -> abs(re) < 1e-9, overland_water_balance.relative_error)
    end
    Wflow.close_files(model; delete_output = false)
end
