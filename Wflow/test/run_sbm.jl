using Dates

tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.Model(config)
(; domain) = model

Wflow.run_timestep!(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2000-01-02T00:00:00")
    @test row.Q ≈ 6.873232832284681
    @test row.storage ≈ 2.7543324658385325e7
    @test row.temp_bycoord ≈ 2.390000104904175
    @test row.vwc_layer2_bycoord ≈ 0.25938809638672006
    @test row.temp_byindex ≈ 2.390000104904175
    @test row.Q_6336050 ≈ 0.003199581692927209
    @test row.Q_6336510 ≈ 0.012862611828173352
    @test row.Q_6836100 ≈ 0.00875995439640751
    @test row.Q_6336500 ≈ 0.0024066827521196437
    @test row.Q_6836190 ≈ 0.002527417212072159
    @test row.Q_6336800 ≈ 0.004389783779801617
    @test row.Q_6336900 ≈ 0.003937679133692149
    @test row.Q_6336930 ≈ 0.013838197736695087
    @test row.Q_6336910 ≈ 0.0030751269759244286
    @test row.Q_6136500 ≈ 0.0004352908448117654
    @test row.Q_6136520 ≈ 0.0004120247552775616
    @test row.Q_6136150 ≈ 0.002963633332122439
    @test row.Q_6136151 ≈ 0.0027878183390109136
    @test row.Q_6136160 ≈ 3.5796610032021494
    @test row.Q_6136202 ≈ 1.3584953239771165
    @test row.recharge_1 ≈ -0.002257181032501202
end

@testset "NetCDF scalar output" begin
    ds = model.writer.dataset_scalar
    @test ds["time"][1] == DateTime("2000-01-02T00:00:00")
    @test ds["Q"][:][1:20] ≈ [
        0.79544294,
        1.8472399,
        1.8459954,
        1.5499605,
        13.334657,
        4.865496,
        0.09143141,
        4.437221,
        0.010515818,
        6.1864867,
        6.1764274,
        0.0060072034,
        0.0066146287,
        0.006407934,
        0.00416659,
        1.2095966,
        0.002213421,
        1.832764,
        8.601765,
        2.7730224,
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
    @test mean(sbm.variables.runoff) ≈ 0.04177459898728149
    @test mean(sbm.variables.soilevap) ≈ 0.02122698830889417
    @test mean(sbm.variables.actevap) ≈ 0.33545623834952143
    @test mean(sbm.variables.actinfilt) ≈ 1.6444774688444848
    @test snow.variables.snow_storage[5] ≈ 3.768513390588815
    @test mean(snow.variables.snow_storage) ≈ 0.038019723676094325
    @test sbm.variables.total_storage[50063] ≈ 559.9035608052374
    @test sbm.variables.total_storage[429] ≈ 597.1794136370914 # river cell
end

# run the second timestep
Wflow.run_timestep!(model)

@testset "second timestep" begin
    sbm = model.land.soil
    snow = model.land.snow
    @test sbm.parameters.theta_s[50063] ≈ 0.48755401372909546
    @test sbm.parameters.theta_r[50063] ≈ 0.15943120419979095
    @test mean(sbm.variables.net_runoff) ≈ 0.23711326799150073
    @test mean(sbm.variables.runoff) ≈ 0.23748172993345834
    @test mean(sbm.variables.soilevap) ≈ 0.018750808322054
    @test mean(sbm.variables.actevap) ≈ 0.27366312063874443
    @test mean(sbm.variables.actinfilt) ≈ 0.08863102527394363
    @test snow.variables.snow_storage[5] ≈ 3.843412524052313
    @test mean(snow.variables.snow_storage) ≈ 0.03461317061870949
    @test sbm.variables.total_storage[50063] ≈ 560.0152135062889
    @test sbm.variables.total_storage[429] ≈ 616.8916185820142  # river cell
end

@testset "subsurface flow" begin
    ssf = model.routing.subsurface_flow.variables.ssf
    @test sum(ssf) ≈ 6.371790004604776e7
    @test ssf[domain.land.network.order[1]] ≈ 717.0520437557641
    @test ssf[domain.land.network.order[end - 100]] ≈ 2333.80105657076
    @test ssf[domain.land.network.order[end]] ≈ 288.19428729403984
end

@testset "overland flow" begin
    q = model.routing.overland_flow.variables.q_av
    @test sum(q) ≈ 285.29655590226673
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[domain.land.network.order[end]] ≈ 1.0e-30
end

@testset "river flow" begin
    q = model.routing.river_flow.variables.q_av
    @test sum(q) ≈ 3846.1068741546355
    @test q[1622] ≈ 0.0007514485505472345
    @test q[43] ≈ 11.928086691686428
    @test q[domain.river.network.order[end]] ≈ 0.044189018390543434
end

@testset "reservoir simple" begin
    res = model.routing.river_flow.boundary_conditions.reservoir
    @test res.variables.outflow[1] ≈ 0.2174998614438593
    @test res.variables.outflow_av[1] ≈ 0.21749986282401396
    @test res.boundary_conditions.inflow[1] ≈ 0.00051287944327482
    @test res.variables.storage[1] ≈ 2.751299001489657f7
    @test res.variables.storage_av[1] ≈ 2.752388968718314f7
    @test res.variables.actevap[1] ≈ 0.5400000810623169
    @test res.boundary_conditions.precipitation[1] ≈ 0.17999997735023499
    @test res.boundary_conditions.evaporation[1] ≈ 0.5400000810623169
end

# set these variables for comparison in "changed dynamic parameters"
precip = copy(model.land.atmospheric_forcing.precipitation)
evap = copy(model.land.atmospheric_forcing.potential_evaporation)
lai = copy(model.land.vegetation_parameters.leaf_area_index)
res_evap = copy(
    model.routing.river_flow.boundary_conditions.reservoir.boundary_conditions.evaporation,
)

Wflow.close_files(model; delete_output = false)

# test without lateral snow transport
tomlpath = joinpath(@__DIR__, "sbm_config.toml")

config = Wflow.Config(tomlpath)
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

# test for setting a pit and multithreading multiple basins (by setting 2 extra pits
# resulting in 3 basins)
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)
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
    @test q[4009] ≈ 8.533285784984278 # pit/ outlet, CartesianIndex(141, 228)
    @test q[4020] ≈ 0.006779712743030111 # downstream of pit 4009, CartesianIndex(141, 229)
    @test q[2508] ≈ 150.2672574147025 # pit/ outlet
    @test q[5808] ≈ 0.12319287597731672 # pit/ outlet
end

# test changing forcing and cyclic LAI parameter
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

config.input.forcing["atmosphere_water__precipitation_volume_flux"] =
    Wflow.init_config_section(
        Wflow.InputEntry,
        Dict("scale" => 2.0, "netcdf" => Dict("variable" => Dict("name" => "precip"))),
    )
config.input.forcing["land_surface_water__potential_evaporation_volume_flux"] =
    Wflow.init_config_section(
        Wflow.InputEntry,
        Dict(
            "scale" => 3.0,
            "offset" => 1.50,
            "netcdf" => Dict("variable" => Dict("name" => "pet")),
        ),
    )
config.input.cyclic["vegetation__leaf_area_index"] = Wflow.init_config_section(
    Wflow.InputEntry,
    Dict("scale" => 1.6, "netcdf" => Dict("variable" => Dict("name" => "LAI"))),
)

model = Wflow.Model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "changed dynamic parameters" begin
    res = model.routing.river_flow.boundary_conditions.reservoir
    land = model.land
    @test land.atmospheric_forcing.precipitation[2] / precip[2] ≈ 2.0f0
    @test (land.atmospheric_forcing.potential_evaporation[100] - 1.50) / evap[100] ≈ 3.0f0
    @test land.vegetation_parameters.leaf_area_index[100] / lai[100] ≈ 1.6f0
    @test (res.boundary_conditions.evaporation[2] - 1.50) / res_evap[2] ≈ 3.0f0
end

# test cyclic river and reservoir external inflow (kinematic wave routing)
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

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
    @test model.routing.river_flow.variables.q_av[44] ≈ 10.545108098407255
    @test reservoir.boundary_conditions.external_inflow[2] == -1.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[2] == 1.0
    @test reservoir.boundary_conditions.inflow[2] ≈ -0.9034940832803586
    @test reservoir.variables.storage_av[2] ≈ 4.621644103995512e7
    @test reservoir.variables.outflow_av[2] ≈ 3.000999922024245
end

# test cyclic river and reservoir external inflow (local inertial routing)
tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
config = Wflow.Config(tomlpath)

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
    @test model.routing.river_flow.variables.q_av[44] ≈ 10.50644055015246
    @test reservoir.boundary_conditions.external_inflow[2] == -1.0
    @test reservoir.boundary_conditions.actual_external_abstraction_av[2] == 1.0
    @test reservoir.boundary_conditions.inflow[2] ≈ -0.9071721542108732
    @test reservoir.variables.storage_av[2] ≈ 4.6218243987631395e7
    @test reservoir.variables.outflow_av[2] ≈ 3.000999922022744
end

# test external negative inflow
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
model.routing.river_flow.boundary_conditions.external_inflow[44] = -10.0
(; actual_external_abstraction_av, external_inflow) =
    model.routing.river_flow.boundary_conditions
(; q_av) = model.routing.river_flow.variables
@testset "river external negative inflow" begin
    Wflow.run_timestep!(model)
    @test actual_external_abstraction_av[44] ≈ 1.5977234689598292
    @test q_av[44] ≈ 1.4346381048266201
    Wflow.run_timestep!(model)
    @test actual_external_abstraction_av[44] ≈ 5.6461557885428135
    @test q_av[44] ≈ 4.173559689914583
    Wflow.run_timestep!(model)
    @test actual_external_abstraction_av[44] ≈ 9.85199447859229
    @test external_inflow[44] == -10.0
    @test q_av[44] ≈ 8.159414181574551
end

# test fixed forcing (precipitation = 2.5)
config = Wflow.Config(tomlpath)
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

# test local-inertial option for river flow river_routing"
tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.Model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "river flow and depth (local inertial)" begin
    q = model.routing.river_flow.variables.q_av
    @test sum(q) ≈ 3851.882468453605
    @test q[1622] ≈ 7.289980625194041e-5
    @test q[43] ≈ 11.708212586581881
    @test q[501] ≈ 3.477737734513234
    h = model.routing.river_flow.variables.h_av
    @test h[1622] ≈ 0.001985933166759342
    @test h[43] ≈ 0.4329730866991719
    @test h[501] ≈ 0.0563033301985948
    q_channel = model.routing.river_flow.variables.q_channel_av
    @test q ≈ q_channel
end

# test external negative inflow local-inertial river flow
tomlpath = joinpath(@__DIR__, "sbm_river-local-inertial_config.toml")
config = Wflow.Config(tomlpath)
model = Wflow.Model(config)
model.routing.river_flow.boundary_conditions.external_inflow[44] = -10.0
(; actual_external_abstraction_av, external_inflow, reservoir) =
    model.routing.river_flow.boundary_conditions
(; q_av) = model.routing.river_flow.variables
@testset "river external negative inflow (local inertial)" begin
    Wflow.run_timestep!(model)
    @test actual_external_abstraction_av[44] ≈ 3.0125200444690616
    @test q_av[44] ≈ 0.0007487510598495983
    Wflow.run_timestep!(model)
    @test actual_external_abstraction_av[44] ≈ 9.73523034718189
    @test q_av[44] ≈ 0.07687497983468243
    Wflow.run_timestep!(model)
    @test actual_external_abstraction_av[44] ≈ 10.0
    @test external_inflow[44] == -10.0
    @test q_av[44] ≈ 7.904126692597465
end
model = Wflow.Model(config)

Wflow.close_files(model; delete_output = false)

# test local-inertial option for river and overland flow
tomlpath = joinpath(@__DIR__, "sbm_river-land-local-inertial_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.Model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "river and overland flow and depth (local inertial)" begin
    q = model.routing.river_flow.variables.q_av
    @test sum(q) ≈ 2495.9569197351266
    @test q[1622] ≈ 7.300041557915818e-5
    @test q[43] ≈ 5.361194661421843
    @test q[501] ≈ 1.6008404903133953
    h = model.routing.river_flow.variables.h_av
    @test h[1622] ≈ 0.001986911257475209
    @test h[43] ≈ 0.30038092757954926
    @test h[501] ≈ 0.031910646103667804
    qx = model.routing.overland_flow.variables.qx
    qy = model.routing.overland_flow.variables.qy
    @test qx[[26, 35, 631]] ≈ [0.19370097955175714, 0.026749121396015898, 0.0]
    @test qy[[26, 35, 631]] ≈ [0.1284818298277494, 1.7209662222437763, 0.0]
    h = model.routing.overland_flow.variables.h
    @test h[[26, 35, 631]] ≈
          [0.07360961858558272, 0.009134592417481706, 0.000726454163107432]
end

Wflow.close_files(model; delete_output = false)

# test local-inertial option for river flow including 1D floodplain schematization
tomlpath = joinpath(@__DIR__, "sbm_river-floodplain-local-inertial_config.toml")
config = Wflow.Config(tomlpath)
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
    @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈ 49.64308797127469
    @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 70.28617594254938
    h = 1.5
    i1, i2 = Wflow.interpolation_indices(h, fp.depth)
    @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈ 182.032315978456
    @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 118.62585278276481
    h = 1.7
    i1, i2 = Wflow.interpolation_indices(h, fp.depth)
    @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
          228.36739676840216
    @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 119.02585278276482
    h = 3.2
    i1, i2 = Wflow.interpolation_indices(h, fp.depth)
    @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈ 695.0377019748654
    @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 307.3730700179533
    h = 4.0
    i1, i2 = Wflow.interpolation_indices(h, fp.depth)
    @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈ 959.816157989228
    @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 308.9730700179533
    @test Wflow.flow_area(fp.width[i2, 4], fp.a[i1, 4], fp.depth[i1], h) ≈ 407.6395313908081
    @test Wflow.wetted_perimeter(fp.p[i1, 4], fp.depth[i1], h) ≈ 90.11775307900271
end

Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "river flow (local inertial) with floodplain schematization simulation" begin
    q = model.routing.river_flow.variables.q_av
    @test sum(q) ≈ 3841.1358969211033
    @test q[1622] ≈ 7.289980623447371e-5
    @test q[43] ≈ 11.708212586581894
    @test q[501] ≈ 3.4204477934334836
    @test q[5808] ≈ 0.0022226684528054044
    h = model.routing.river_flow.variables.h_av
    @test h[1622] ≈ 0.0019859331038689374
    @test h[43] ≈ 0.4329730866994521
    @test h[501] ≈ 0.05578682441669217
    @test h[5808] ≈ 0.005930509404131392
end

# set boundary condition local inertial routing from netCDF file
config.input.static["model_boundary_condition_river__length"] = "riverlength_bc"
config.input.static["model_boundary_condition_river_bank_water__depth"] = "riverdepth_bc"
model = Wflow.Model(config)
Wflow.run_timestep!(model)
Wflow.run_timestep!(model)

@testset "change boundary condition for local inertial routing (including floodplain)" begin
    q = model.routing.river_flow.variables.q_av
    @test sum(q) ≈ 3841.325252009684
    @test q[1622] ≈ 7.289980623447371e-5
    @test q[43] ≈ 11.708212586581894
    @test q[501] ≈ 3.4204477934334836
    @test q[5808] ≈ 0.055273361881645354
    h = model.routing.river_flow.variables.h_av
    @test h[1622] ≈ 0.0019859331038689374
    @test h[43] ≈ 0.4329730866994521
    @test h[501] ≈ 0.05578682441669217
    @test h[5808] ≈ 1.9999993313276971
end
Wflow.close_files(model; delete_output = false)

# test different ksat profiles
@testset "ksat profiles (SBM)" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")

    function get_config(profile)
        config = Wflow.Config(tomlpath)
        config.model.saturated_hydraulic_conductivity_profile = profile
        config.input.static["soil_layer_water__vertical_saturated_hydraulic_conductivity"] = "kv"
        config.input.static["soil_exponential_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth"] =
            400.0
        config.input.static["soil_layered_vertical_saturated_hydraulic_conductivity_profile_below_surface__depth"] =
            400.0
        config
    end

    i = 100
    config = get_config("exponential")

    @testset "exponential profile" begin
        config = get_config("exponential")
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
        config = get_config("exponential_constant")
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
        kv_1000 = Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, 1000.0, i, 3)
        @test kv_400 ≈ kv_1000
        @test all(kv_profile.z_exp .== 400.0)
        @test subsurface_flow.variables.ssfmax[i] ≈ 49.38558575188426
        @test subsurface_flow.variables.ssf[i] ≈ 24810.460986497365
    end

    @testset "layered profile" begin
        config = get_config("layered")
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

    config = get_config("layered_exponential")

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
        model = Wflow.Model(config)
        Wflow.run_timestep!(model)
        Wflow.run_timestep!(model)
        q = model.routing.river_flow.variables.q_av
        @test sum(q) ≈ 3300.5679346067054
        @test q[1622] ≈ 0.000698946772256762
        @test q[43] ≈ 9.667763405724566
    end

    Wflow.close_files(model; delete_output = false)
end
