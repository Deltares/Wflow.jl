using Dates

tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
(; network) = model

model = Wflow.run_timestep(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2000-01-02T00:00:00")
    @test row.Q ≈ 8.15299947254324f0
    @test row.volume ≈ 2.7535003939625636f7
    @test row.temp_bycoord ≈ 2.390000104904175f0
    @test row.vwc_layer2_bycoord ≈ 0.25938809638672006f0
    @test row.temp_byindex ≈ 2.390000104904175f0
    @test row.Q_6336050 ≈ 0.006583064321841488f0
    @test row.Q_6336510 ≈ 0.029864230092642642f0
    @test row.Q_6836100 ≈ 0.19995488963854305f0
    @test row.Q_6336500 ≈ 0.006277726622788425f0
    @test row.Q_6836190 ≈ 0.0031262850749354237f0
    @test row.Q_6336800 ≈ 0.008278375560053742f0
    @test row.Q_6336900 ≈ 0.0066141980189014385f0
    @test row.Q_6336930 ≈ 0.09141703511009937f0
    @test row.Q_6336910 ≈ 0.007475453481320056f0
    @test row.Q_6136500 ≈ 0.001834989281902289f0
    @test row.Q_6136520 ≈ 0.0022266031120691397f0
    @test row.Q_6136150 ≈ 0.006310361139139334f0
    @test row.Q_6136151 ≈ 0.007946301730645885f0
    @test row.Q_6136160 ≈ 3.927719795530719f0
    @test row.Q_6136202 ≈ 1.4162246003743886f0
    @test row.recharge_1 ≈ -0.0020800523945940217f0
end

@testset "NetCDF scalar output" begin
    ds = model.writer.dataset_scalar
    @test ds["time"][1] == DateTime("2000-01-02T00:00:00")
    @test ds["Q"][:][1:20] ≈ [
        0.7425387f0,
        1.4162246f0,
        1.4425076f0,
        1.4044669f0,
        5.738109f0,
        2.7616737f0,
        2.1128905f0,
        4.105428f0,
        0.008651769f0,
        3.9277198f0,
        4.069447f0,
        0.006356805f0,
        0.007946302f0,
        0.008135906f0,
        0.0037393502f0,
        0.70888275f0,
        0.0024000728f0,
        1.3347782f0,
        3.8374817f0,
        1.676597f0,
    ]
    @test ds["Q_gauges"].attrib["cf_role"] == "timeseries_id"
    @test ds["temp_index"][:] ≈ [2.39f0]
    @test ds["temp_coord"][:] ≈ [2.39f0]
    @test keys(ds.dim) == ["time", "layer", "Q_gauges", "temp_bycoord", "temp_byindex"]
end

@testset "first timestep" begin
    sbm = model.vertical.soil
    snow = model.vertical.snow
    @test snow.parameters.tt[50063] ≈ 0.0f0

    @test model.clock.iteration == 1

    @test sbm.parameters.theta_s[50063] ≈ 0.48755401372909546f0
    @test sbm.parameters.theta_r[50063] ≈ 0.15943120419979095f0
    @test mean(sbm.variables.runoff) ≈ 0.04177459898728149f0
    @test mean(sbm.variables.soilevap) ≈ 0.02122698830889417f0
    @test mean(sbm.variables.actevap) ≈ 0.3353001180202587f0
    @test mean(sbm.variables.actinfilt) ≈ 1.6444774688444848f0
    @test snow.variables.snow_storage[5] ≈ 3.768513390588815f0
    @test mean(snow.variables.snow_storage) ≈ 0.038019723676094325f0
    @test sbm.variables.total_storage[50063] ≈ 559.9035608052374f0
    @test sbm.variables.total_storage[429] ≈ 597.4578475404879f0 # river cell
end

# run the second timestep
model = Wflow.run_timestep(model)

@testset "second timestep" begin
    sbm = model.vertical.soil
    snow = model.vertical.snow
    @test sbm.parameters.theta_s[50063] ≈ 0.48755401372909546f0
    @test sbm.parameters.theta_r[50063] ≈ 0.15943120419979095f0
    @test mean(sbm.variables.net_runoff) ≈ 0.23734052031823816f0
    @test mean(sbm.variables.runoff) ≈ 0.23770898226019577f0
    @test mean(sbm.variables.soilevap) ≈ 0.018750808322054897f0
    @test mean(sbm.variables.actevap) ≈ 0.14545276216428166f0
    @test mean(sbm.variables.actinfilt) ≈ 0.08863102527394363f0
    @test snow.variables.snow_storage[5] ≈ 3.843412524052313f0
    @test mean(snow.variables.snow_storage) ≈ 0.03461317061870949f0
    @test sbm.variables.total_storage[50063] ≈ 560.0152135062889f0
    @test sbm.variables.total_storage[429] ≈ 617.2238533241972f0 # river cell
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 6.3761585406186976f7
    @test ssf[network.land.order[1]] ≈ 718.2802566393531f0
    @test ssf[network.land.order[end - 100]] ≈ 2337.771227118579f0
    @test ssf[network.land.order[end]] ≈ 288.19428729403984f0
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 291.4923871784623f0
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 1.0f-30
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 3625.0013368279815f0
    @test q[1622] ≈ 0.0006503254947860838f0
    @test q[43] ≈ 12.06416878694095f0
    @test q[network.river.order[end]] ≈ 0.039200124520463835f0
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    @test res.outflow[1] ≈ 0.21750000119148086f0
    @test res.inflow[1] ≈ 43.18479982574888f0
    @test res.volume[1] ≈ 2.751299001489657f7
    @test res.precipitation[1] ≈ 0.17999997735023499f0
    @test res.evaporation[1] ≈ 0.5400000810623169f0
end

# set these variables for comparison in "changed dynamic parameters"
precip = copy(model.vertical.atmospheric_forcing.precipitation)
evap = copy(model.vertical.atmospheric_forcing.potential_evaporation)
lai = copy(model.vertical.vegetation_parameter_set.leaf_area_index)
res_evap = copy(model.lateral.river.reservoir.evaporation)

Wflow.close_files(model; delete_output = false)

# test for setting a pit and multithreading multiple basins (by setting 2 extra pits
# resulting in 3 basins)
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)
config["model"]["pits"] = true
config["input"]["pits"] = "wflow_pits"
config.endtime = DateTime(2000, 1, 9)
config.loglevel = "info"

model = Wflow.run(config)

@testset "timing" begin
    # clock has been reset
    calendar = get(config, "calendar", "standard")::String
    @test model.clock.time == Wflow.cftime(config.starttime, calendar)
    @test model.clock.iteration == 0
end

@testset "river flow at basin outlets and downstream of one pit" begin
    q = model.lateral.river.q_av
    @test q[4009] ≈ 8.543145028037452f0 # pit/ outlet, CartesianIndex(141, 228)
    @test q[4020] ≈ 0.006779014715290862f0 # downstream of pit 4009, CartesianIndex(141, 229)
    @test q[2508] ≈ 150.5640617045796f0 # pit/ outlet
    @test q[5808] ≈ 0.12438899196040518f0 # pit/ outlet
end

# test changing forcing and cyclic LAI parameter
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

config.input.vertical.atmospheric_forcing.precipitation =
    Dict("scale" => 2.0, "netcdf" => Dict("variable" => Dict("name" => "precip")))
config.input.vertical.atmospheric_forcing.potential_evaporation = Dict(
    "scale" => 3.0,
    "offset" => 1.50,
    "netcdf" => Dict("variable" => Dict("name" => "pet")),
)
config.input.vertical.vegetation_parameter_set.leaf_area_index =
    Dict("scale" => 1.6, "netcdf" => Dict("variable" => Dict("name" => "LAI")))

model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "changed dynamic parameters" begin
    res = model.lateral.river.reservoir
    vertical = model.vertical
    @test vertical.atmospheric_forcing.precipitation[2] / precip[2] ≈ 2.0f0
    @test (vertical.atmospheric_forcing.potential_evaporation[100] - 1.50) / evap[100] ≈
          3.0f0
    @test vertical.vegetation_parameter_set.leaf_area_index[100] / lai[100] ≈ 1.6f0
    @test (res.evaporation[2] - 1.50) / res_evap[2] ≈ 3.0000012203408635f0
end

# test cyclic river inflow
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

config.input.cyclic =
    ["vertical.vegetation_parameter_set.leaf_area_index", "lateral.river.inflow"]
config.input.lateral.river.inflow = "inflow"

model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "river inflow (cyclic)" begin
    @test model.lateral.river.inflow[44] ≈ 0.75
    @test model.lateral.river.q_av[44] ≈ 10.723729440690567f0
end

# test fixed forcing (precipitation = 2.5)
config = Wflow.Config(tomlpath)
config.input.vertical.atmospheric_forcing.precipitation = Dict("value" => 2.5)
model = Wflow.initialize_sbm_model(config)
Wflow.load_fixed_forcing(model)

@testset "fixed precipitation forcing (initialize)" begin
    @test maximum(model.vertical.atmospheric_forcing.precipitation) ≈ 2.5
    @test minimum(model.vertical.atmospheric_forcing.precipitation) ≈ 0.0
    @test all(isapprox.(model.lateral.river.reservoir.precipitation, 2.5))
end

model = Wflow.run_timestep(model)

@testset "fixed precipitation forcing (first timestep)" begin
    @test maximum(model.vertical.atmospheric_forcing.precipitation) ≈ 2.5
    @test minimum(model.vertical.atmospheric_forcing.precipitation) ≈ 0.0
    @test all(isapprox.(model.lateral.river.reservoir.precipitation, 2.5))
end

Wflow.close_files(model; delete_output = false)

# test local-inertial option for river flow river_routing
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)
config.model.river_routing = "local-inertial"

model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "river flow and depth (local inertial)" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 3922.0644366362544f0
    @test q[1622] ≈ 7.315676375562105f-5
    @test q[43] ≈ 11.92787156357907f0
    @test q[501] ≈ 3.57855182713785f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.001987887644883981f0
    @test h[43] ≈ 0.4366415244811759f0
    @test h[501] ≈ 0.057317706869865745f0
    q_channel = model.lateral.river.q_channel_av
    @test q ≈ q_channel
end
Wflow.close_files(model; delete_output = false)

# test local-inertial option for river and overland flow
tomlpath = joinpath(@__DIR__, "sbm_swf_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "river and overland flow and depth (local inertial)" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 2380.64389229669f0
    @test q[1622] ≈ 7.328535246760549f-5
    @test q[43] ≈ 5.3566292152594155f0
    @test q[501] ≈ 1.6042388573126602f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.0019891342000364796f0
    @test h[43] ≈ 0.30026439683630496f0
    @test h[501] ≈ 0.03195324587192846f0
    qx = model.lateral.land.qx
    qy = model.lateral.land.qy
    @test qx[[26, 35, 631]] ≈ [0.1939736998417174f0, 0.026579954465883678f0, 0.0f0]
    @test qy[[26, 35, 631]] ≈ [0.12906530420401777f0, 1.7225115950614904f0, 0.0f0]
    h = model.lateral.land.h
    @test h[[26, 35, 631]] ≈
          [0.07367301172613304f0, 0.009139882310161706f0, 0.0007482998926237368f0]
end

Wflow.close_files(model; delete_output = false)

# test local-inertial option for river flow including 1D floodplain schematization
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

config.model.floodplain_1d = true
config.model.river_routing = "local-inertial"
config.model.land_routing = "kinematic-wave"
Dict(config.input.lateral.river)["floodplain"] = Dict("volume" => "floodplain_volume")
Dict(config.state.lateral.river)["floodplain"] =
    Dict("q" => "q_floodplain", "h" => "h_floodplain")

model = Wflow.initialize_sbm_model(config)

fp = model.lateral.river.floodplain.profile
river = model.lateral.river
dh = diff(fp.depth)
Δv = diff(fp.volume[:, 3])
Δa = diff(fp.a[:, 3])

@testset "river flow (local inertial) floodplain schematization" begin
    # floodplain geometry checks (index 3)
    @test fp.volume[:, 3] ≈ [0.0f0, 8641.0f0, 19011.0f0, 31685.0f0, 51848.0f0, 80653.0f0]
    @test fp.width[:, 3] ≈ [
        30.0f0,
        99.28617594254938f0,
        119.15260323159785f0,
        145.6258527827648f0,
        231.6754039497307f0,
        330.9730700179533f0,
    ]
    @test fp.p[:, 3] ≈ [
        69.28617594254938f0,
        70.28617594254938f0,
        91.15260323159785f0,
        118.62585278276481f0,
        205.6754039497307f0,
        305.9730700179533f0,
    ]
    @test fp.a[:, 3] ≈ [
        0.0f0,
        49.64308797127469f0,
        109.21938958707361f0,
        182.032315978456f0,
        297.8700179533214f0,
        463.35655296229805f0,
    ]
    @test dh .* fp.width[2:end, 3] * river.dl[3] ≈ Δv
    @test fp.a[:, 3] * river.dl[3] ≈ fp.volume[:, 3]
    # flood depth from flood volume (8000.0)
    flood_vol = 8000.0f0
    river.volume[3] = flood_vol + river.bankfull_volume[3]
    i1, i2 = Wflow.interpolation_indices(flood_vol, fp.volume[:, 3])
    @test (i1, i2) == (1, 2)
    flood_depth = Wflow.flood_depth(fp, flood_vol, river.dl[3], 3)
    @test flood_depth ≈ 0.46290938548779076f0
    @test (flood_depth - fp.depth[i1]) * fp.width[i2, 3] * river.dl[3] + fp.volume[i1, 3] ≈
          flood_vol
    # flood depth from flood volume (12000.0)
    flood_vol = 12000.0f0
    river.volume[3] = flood_vol + river.bankfull_volume[3]
    i1, i2 = Wflow.interpolation_indices(flood_vol, fp.volume[:, 3])
    @test (i1, i2) == (2, 3)
    flood_depth = Wflow.flood_depth(fp, flood_vol, river.dl[3], 3)
    @test flood_depth ≈ 0.6619575699132112f0
    @test (flood_depth - fp.depth[i1]) * fp.width[i2, 3] * river.dl[3] + fp.volume[i1, 3] ≈
          flood_vol
    # test extrapolation of segment
    flood_vol = 95000.0f0
    river.volume[3] = flood_vol + river.bankfull_volume[3]
    i1, i2 = Wflow.interpolation_indices(flood_vol, fp.volume[:, 3])
    @test (i1, i2) == (6, 6)
    flood_depth = Wflow.flood_depth(fp, flood_vol, river.dl[3], 3)
    @test flood_depth ≈ 2.749036625585836f0
    @test (flood_depth - fp.depth[i1]) * fp.width[i2, 3] * river.dl[3] + fp.volume[i1, 3] ≈
          flood_vol
    river.volume[3] = 0.0 # reset volume
    # flow area and wetted perimeter based on hf
    h = 0.5
    i1, i2 = Wflow.interpolation_indices(h, fp.depth)
    @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
          49.64308797127469f0
    @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 70.28617594254938f0
    h = 1.5
    i1, i2 = Wflow.interpolation_indices(h, fp.depth)
    @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
          182.032315978456f0
    @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 118.62585278276481f0
    h = 1.7
    i1, i2 = Wflow.interpolation_indices(h, fp.depth)
    @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
          228.36739676840216f0
    @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 119.02585278276482f0
    h = 3.2
    i1, i2 = Wflow.interpolation_indices(h, fp.depth)
    @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
          695.0377019748654f0
    @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 307.3730700179533f0
    h = 4.0
    i1, i2 = Wflow.interpolation_indices(h, fp.depth)
    @test Wflow.flow_area(fp.width[i2, 3], fp.a[i1, 3], fp.depth[i1], h) ≈
          959.816157989228f0
    @test Wflow.wetted_perimeter(fp.p[i1, 3], fp.depth[i1], h) ≈ 308.9730700179533f0
    @test Wflow.flow_area(fp.width[i2, 4], fp.a[i1, 4], fp.depth[i1], h) ≈
          407.6395313908081f0
    @test Wflow.wetted_perimeter(fp.p[i1, 4], fp.depth[i1], h) ≈ 90.11775307900271f0
end

model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "river flow (local inertial) with floodplain schematization simulation" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 3910.4728717811836f0
    @test q[1622] ≈ 7.315676384849305f-5
    @test q[43] ≈ 11.92787156357908f0
    @test q[501] ≈ 3.510668846752431f0
    @test q[5808] ≈ 0.002223993845806248f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.001987887580593841f0
    @test h[43] ≈ 0.436641524481545f0
    @test h[501] ≈ 0.05670770509802258f0
    @test h[5808] ≈ 0.005929945681367346f0
end

# set boundary condition local inertial routing from netCDF file
config.input.lateral.river.riverlength_bc = "riverlength_bc"
config.input.lateral.river.riverdepth_bc = "riverdepth_bc"
model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "change boundary condition for local inertial routing (including floodplain)" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 3910.683449719468f0
    @test q[1622] ≈ 7.315757521099307f-5
    @test q[43] ≈ 11.927871563591228f0
    @test q[501] ≈ 3.5106678593721496f0
    @test q[5808] ≈ 0.060518234525259465f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.0019878952928530183f0
    @test h[43] ≈ 0.4366415249636809f0
    @test h[501] ≈ 0.056707564314724804f0
    @test h[5808] ≈ 2.0000006940603936f0
end
Wflow.close_files(model; delete_output = false)

# test different ksat profiles
@testset "ksat profiles (SBM)" begin
    i = 100
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.input.vertical.soil.parameters.kv = "kv"
    config.input.vertical.soil.parameters.z_exp = Dict("value" => 400.0)
    config.input.vertical.soil.parameters.z_layered = Dict("value" => 400.0)

    @testset "exponential profile" begin
        model = Wflow.initialize_sbm_model(config)
        (; soil) = model.vertical
        (; kv_profile) = soil.parameters
        (; subsurface) = model.lateral
        z = soil.variables.zi[i]
        kvfrac = soil.parameters.kvfrac
        kv_z = Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2)
        @test kv_z ≈ kvfrac[i][2] * kv_profile.kv_0[i] * exp(-kv_profile.f[i] * z)
        @test subsurface.ssfmax[i] ≈ 28.32720603576582f0
        @test subsurface.ssf[i] ≈ 11683.330684556406f0
    end

    @testset "exponential constant profile" begin
        config.input.vertical.ksat_profile = "exponential_constant"
        model = Wflow.initialize_sbm_model(config)
        (; soil) = model.vertical
        (; kv_profile) = soil.parameters
        (; subsurface) = model.lateral
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
        @test subsurface.ssfmax[i] ≈ 49.38558575188426f0
        @test subsurface.ssf[i] ≈ 24810.460986497365f0
    end

    @testset "layered profile" begin
        config.input.vertical.ksat_profile = "layered"
        model = Wflow.initialize_sbm_model(config)
        (; soil) = model.vertical
        (; kv_profile) = soil.parameters
        (; subsurface) = model.lateral
        z = soil.variables.zi[i]
        kvfrac = soil.parameters.kvfrac
        @test Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2) ≈
              kv_profile.kv[i][2]
        Wflow.kh_layered_profile!(soil, subsurface, kv_profile, 86400.0)
        @test subsurface.kh_profile.kh[i] ≈ 47.508932674632355f0
        @test subsurface.ssfmax[i] ≈ 30.237094380100316f0
        @test subsurface.ssf[i] ≈ 14546.518932613191f0
    end

    @testset "layered exponential profile" begin
        config.input.vertical.ksat_profile = "layered_exponential"
        model = Wflow.initialize_sbm_model(config)
        (; soil) = model.vertical
        (; kv_profile) = soil.parameters
        (; subsurface) = model.lateral
        z = soil.variables.zi[i]
        kvfrac = soil.parameters.kvfrac
        @test Wflow.hydraulic_conductivity_at_depth(kv_profile, kvfrac, z, i, 2) ≈
              kv_profile.kv[i][2]
        @test kv_profile.nlayers_kv[i] == 2
        Wflow.kh_layered_profile!(soil, subsurface, kv_profile, 86400.0)
        @test subsurface.kh_profile.kh[i] ≈ 33.76026208801769f0
        @test all(kv_profile.z_layered[1:10] .== 400.0)
        @test subsurface.ssfmax[i] ≈ 23.4840490395906f0
        @test subsurface.ssf[i] ≈ 10336.88327617503f0
    end

    model = Wflow.run_timestep(model)
    model = Wflow.run_timestep(model)
    @testset "river flow layered exponential profile" begin
        q = model.lateral.river.q_av
        @test sum(q) ≈ 3159.38300016008f0
        @test q[1622] ≈ 0.0005972577112819149f0
        @test q[43] ≈ 10.017642376280731f0
    end

    Wflow.close_files(model; delete_output = false)
end
