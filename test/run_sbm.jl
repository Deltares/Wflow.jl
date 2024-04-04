using Dates

tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
@unpack network = model

model = Wflow.run_timestep(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2000-01-02T00:00:00")
    @test row.Q ≈ 8.1396354649197f0
    @test row.volume ≈ 2.7535003939625636f7
    @test row.temp_bycoord ≈ 2.390000104904175f0
    @test row.vwc_layer2_bycoord ≈ 0.25901943991019094f0
    @test row.temp_byindex ≈ 2.390000104904175f0
    @test row.Q_6336050 ≈ 0.006160282939850074f0
    @test row.Q_6336510 ≈ 0.029177309556712334f0
    @test row.Q_6836100 ≈ 0.19621282450614713f0
    @test row.Q_6336500 ≈ 0.006089112638001381f0
    @test row.Q_6836190 ≈ 0.0031262850749354237f0
    @test row.Q_6336800 ≈ 0.007770868657277307f0
    @test row.Q_6336900 ≈ 0.006403194169947582f0
    @test row.Q_6336930 ≈ 0.08888787154163148f0
    @test row.Q_6336910 ≈ 0.007071851236520184f0
    @test row.Q_6136500 ≈ 0.0016367337487926633f0
    @test row.Q_6136520 ≈ 0.002084670434294102f0
    @test row.Q_6136150 ≈ 0.006095549758915344f0
    @test row.Q_6136151 ≈ 0.007643634432992056f0
    @test row.Q_6136160 ≈ 3.9199531493174726f0
    @test row.Q_6136202 ≈ 1.4125847550988493f0
    @test row.recharge_1 ≈ -0.05653226176238641f0
end

@testset "NetCDF scalar output" begin
    ds = model.writer.dataset_scalar
    @test ds["time"][1] == DateTime("2000-01-02T00:00:00")
    @test ds["Q"][:][1:20] ≈ [
        0.73993874f0,
        1.4125848f0,
        1.4389194f0,
        1.4036233f0,
        5.7276225f0,
        2.7586424f0,
        2.1080604f0,
        4.1026044f0,
        0.008365918f0,
        3.919953f0,
        4.0615087f0,
        0.006050462f0,
        0.0076436345f0,
        0.00780229f0,
        0.00346879f0,
        0.7066179f0,
        0.0022446192f0,
        1.3323202f0,
        3.8271446f0,
        1.6729931f0,
    ]
    @test ds["Q_gauges"].attrib["cf_role"] == "timeseries_id"
    @test ds["temp_index"][:] ≈ [2.39f0]
    @test ds["temp_coord"][:] ≈ [2.39f0]
    @test keys(ds.dim) == ["time", "layer", "Q_gauges", "temp_bycoord", "temp_byindex"]
end

@testset "first timestep" begin
    sbm = model.vertical

    @test sbm.tt[50063] ≈ 0.0f0

    @test model.clock.iteration == 1

    @test sbm.θₛ[50063] ≈ 0.48755401372909546f0
    @test sbm.θᵣ[50063] ≈ 0.15943120419979095f0
    @test mean(sbm.runoff) ≈ 0.0416027711985709f0
    @test mean(sbm.soilevap) ≈ 0.00616346765710303f0
    @test mean(sbm.actevap) ≈ 0.5393128846181061f0
    @test mean(sbm.actinfilt) ≈ 1.4860556930083888f0
    @test sbm.snow[5] ≈ 3.592840840467347f0
    @test mean(sbm.snow) ≈ 0.035197394854698104f0
    @test sbm.total_storage[50063] ≈ 559.70849973344f0
    @test sbm.total_storage[429] ≈ 597.4578475404879f0 # river cell
end

# run the second timestep
model = Wflow.run_timestep(model)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.θₛ[50063] ≈ 0.48755401372909546f0
    @test sbm.θᵣ[50063] ≈ 0.15943120419979095f0
    @test mean(sbm.net_runoff) ≈ 0.23657364759232993f0
    @test mean(sbm.runoff) ≈ 0.2369017651524159f0
    @test mean(sbm.soilevap) ≈ 0.015275916086112940
    @test mean(sbm.actevap) ≈ 0.29716147685525746f0
    @test mean(sbm.actinfilt) ≈ 0.08829600284076534f0
    @test sbm.snow[5] ≈ 3.667748983774868f0
    @test mean(sbm.snow) ≈ 0.03209680141455693f0
    @test sbm.total_storage[50063] ≈ 559.7935411649405f0
    @test sbm.total_storage[429] ≈ 617.0062092646873f0 # river cell
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 6.370399148012509f7
    @test ssf[network.land.order[1]] ≈ 7.169036749244327f2
    @test ssf[network.land.order[end-100]] ≈ 2335.2465707069578f0
    @test ssf[network.land.order[end]] ≈ 288.19428729403944f0
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 290.5520014030802f0
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 1.0f-30
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 3615.2084127447865f0
    @test q[1622] ≈ 0.0005986272622329333f0
    @test q[43] ≈ 12.036342425160155f0
    @test q[network.river.order[end]] ≈ 0.03835913312643948f0
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
precip = copy(model.vertical.precipitation)
evap = copy(model.vertical.potential_evaporation)
lai = copy(model.vertical.leaf_area_index)
res_evap = copy(model.lateral.river.reservoir.evaporation)

Wflow.close_files(model, delete_output = false)

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
    @test q[4009] ≈ 8.60480399680283f0 # pit/ outlet, CartesianIndex(141, 228)
    @test q[4020] ≈ 0.006779014715290862f0 # downstream of pit 4009, CartesianIndex(141, 229)
    @test q[2508] ≈ 150.15321264134985f0 # pit/ outlet
    @test q[5808] ≈ 0.12625654862968252f0 # pit/ outlet
end

# test changing forcing and cyclic LAI parameter
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

config.input.vertical.precipitation =
    Dict("scale" => 2.0, "netcdf" => Dict("variable" => Dict("name" => "precip")))
config.input.vertical.potential_evaporation = Dict(
    "scale" => 3.0,
    "offset" => 1.50,
    "netcdf" => Dict("variable" => Dict("name" => "pet")),
)
config.input.vertical.leaf_area_index =
    Dict("scale" => 1.6, "netcdf" => Dict("variable" => Dict("name" => "LAI")))

model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "changed dynamic parameters" begin
    res = model.lateral.river.reservoir
    vertical = model.vertical
    @test vertical.precipitation[2] / precip[2] ≈ 2.0f0
    @test (vertical.potential_evaporation[100] - 1.50) / evap[100] ≈ 3.0f0
    @test vertical.leaf_area_index[100] / lai[100] ≈ 1.6f0
    @test (res.evaporation[2] - 1.50) / res_evap[2] ≈ 3.0000012203408635f0
end

# test cyclic river inflow
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

config.input.cyclic = ["vertical.leaf_area_index", "lateral.river.inflow"]
config.input.lateral.river.inflow = "inflow"

model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "river inflow (cyclic)" begin
    @test model.lateral.river.inflow[44] ≈ 0.75
    @test model.lateral.river.q_av[44] ≈ 10.698591283662008
end

# test fixed forcing (precipitation = 2.5)
config = Wflow.Config(tomlpath)
config.input.vertical.precipitation = Dict("value" => 2.5)
model = Wflow.initialize_sbm_model(config)
Wflow.load_fixed_forcing(model)

@testset "fixed precipitation forcing (initialize)" begin
    @test maximum(model.vertical.precipitation) ≈ 2.5
    @test minimum(model.vertical.precipitation) ≈ 0.0
    @test all(isapprox.(model.lateral.river.reservoir.precipitation, 2.5))
end

model = Wflow.run_timestep(model)

@testset "fixed precipitation forcing (first timestep)" begin
    @test maximum(model.vertical.precipitation) ≈ 2.5
    @test minimum(model.vertical.precipitation) ≈ 0.0
    @test all(isapprox.(model.lateral.river.reservoir.precipitation, 2.5))
end

Wflow.close_files(model, delete_output = false)

# test local-inertial option for river flow river_routing
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)
config.model.river_routing = "local-inertial"

model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "river flow and depth (local inertial)" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 3910.2095742376546f0
    @test q[1622] ≈ 6.0094181857060604f-5
    @test q[43] ≈ 11.900372477232786f0
    @test q[501] ≈ 3.536628093804679f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.001809965063947279f0
    @test h[43] ≈ 0.43627044208669874f0
    @test h[501] ≈ 0.05669956233680719f0
    q_channel = model.lateral.river.q_channel_av
    @test q ≈ q_channel
end
Wflow.close_files(model, delete_output = false)

# test local-inertial option for river and overland flow
tomlpath = joinpath(@__DIR__, "sbm_swf_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "river and overland flow and depth (local inertial)" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 2375.1186861861243f0
    @test q[1622] ≈ 6.011407534125278f-5
    @test q[43] ≈ 5.358152280519331f0
    @test q[501] ≈ 1.5878151534724314f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.0018101771426228585f0
    @test h[43] ≈ 0.3003008810153667f0
    @test h[501] ≈ 0.03162351626631113f0
    qx = model.lateral.land.qx
    qy = model.lateral.land.qy
    @test qx[[26, 35, 631]] ≈ [0.18613687016733824f0, 0.0004519163131931592f0, 0.0f0]
    @test qy[[26, 35, 631]] ≈ [0.12681702046955443f0, 1.7210193779889194f0, 0.0f0]
    h = model.lateral.land.h
    @test h[[26, 35, 631]] ≈
          [0.07341443653334193f0, 0.009152294150993293f0, 0.0006875940563996746f0]
end
Wflow.close_files(model, delete_output = false)

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
Δh = diff(fp.depth)
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
    @test Δh .* fp.width[2:end, 3] * river.dl[3] ≈ Δv
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
    @test sum(q) ≈ 3898.719057830299f0
    @test q[1622] ≈ 6.0094627478450016f-5
    @test q[43] ≈ 11.900372477232796f0
    @test q[501] ≈ 3.470259878228359f0
    @test q[5808] ≈ 0.002199141435542889f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.0018099697988149294f0
    @test h[43] ≈ 0.4362704420867342f0
    @test h[501] ≈ 0.05610231297517167f0
    @test h[5808] ≈ 0.005901095907129176f0
end

# set boundary condition local inertial routing from NetCDF file
config.input.lateral.river.riverlength_bc = "riverlength_bc"
config.input.lateral.river.riverdepth_bc = "riverdepth_bc"
model = Wflow.initialize_sbm_model(config)
model = Wflow.run_timestep(model)
model = Wflow.run_timestep(model)

@testset "change boundary condition for local inertial routing (including floodplain)" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 3899.01831081148f0
    @test q[1622] ≈ 6.008529317397127f-5
    @test q[43] ≈ 11.900394994314782f0
    @test q[501] ≈ 3.4702817507735366f0
    @test q[5808] ≈ 0.060081151400640784f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.001809870969211623f0
    @test h[43] ≈ 0.4362709200147433f0
    @test h[501] ≈ 0.05610440991535934f0
    @test h[5808] ≈ 2.0000006940603936f0
end
Wflow.close_files(model, delete_output = false)

# test different ksat profiles
@testset "ksat profiles (SBM)" begin
    i = 100
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.input.vertical.kv = "kv"
    config.input.vertical.z_exp = Dict("value" => 400.0)
    config.input.vertical.z_layered = Dict("value" => 400.0)

    @testset "exponential profile" begin
        model = Wflow.initialize_sbm_model(config)
        @unpack vertical = model
        z = vertical.zi[i]
        kv_z = Wflow.hydraulic_conductivity_at_depth(vertical, z, i, 2, "exponential")
        @test kv_z ≈ vertical.kvfrac[i][2] * vertical.kv₀[i] * exp(-vertical.f[i] * z)
        @test vertical.z_exp == vertical.soilthickness
        @test_throws ErrorException Wflow.kh_layered_profile(
            vertical,
            100.0,
            i,
            "exponential",
        )
        @test all(isnan.(vertical.z_layered))
        @test all(isnan.(vertical.kv[i]))
        @test all(vertical.nlayers_kv .== 0)
    end

    @testset "exponential constant profile" begin
        config.input.vertical.ksat_profile = "exponential_constant"
        model = Wflow.initialize_sbm_model(config)
        @unpack vertical = model
        z = vertical.zi[i]
        kv_z =
            Wflow.hydraulic_conductivity_at_depth(vertical, z, i, 2, "exponential_constant")
        @test kv_z ≈ vertical.kvfrac[i][2] * vertical.kv₀[i] * exp(-vertical.f[i] * z)
        kv_400 = Wflow.hydraulic_conductivity_at_depth(
            vertical,
            400.0,
            i,
            2,
            "exponential_constant",
        )
        kv_1000 = Wflow.hydraulic_conductivity_at_depth(
            vertical,
            1000.0,
            i,
            3,
            "exponential_constant",
        )
        @test kv_400 ≈ kv_1000
        @test_throws ErrorException Wflow.kh_layered_profile(
            vertical,
            100.0,
            i,
            "exponential_constant",
        )
        @test all(isnan.(vertical.z_layered))
        @test all(isnan.(vertical.kv[i]))
        @test all(vertical.nlayers_kv .== 0)
        @test all(vertical.z_exp .== 400.0)
    end

    @testset "layered profile" begin
        config.input.vertical.ksat_profile = "layered"
        model = Wflow.initialize_sbm_model(config)
        @unpack vertical = model
        z = vertical.zi[i]
        @test Wflow.hydraulic_conductivity_at_depth(vertical, z, i, 2, "layered") ≈
              vertical.kv[100][2]
        @test Wflow.kh_layered_profile(vertical, 100.0, i, "layered") ≈ 47.508932674632355f0
        @test vertical.nlayers_kv[i] == 4
        @test vertical.z_layered == vertical.soilthickness
        @test all(isnan.(vertical.z_exp))
    end

    @testset "layered exponential profile" begin
        config.input.vertical.ksat_profile = "layered_exponential"
        model = Wflow.initialize_sbm_model(config)
        @unpack vertical = model
        z = vertical.zi[i]
        @test Wflow.hydraulic_conductivity_at_depth(
            vertical,
            z,
            i,
            2,
            "layered_exponential",
        ) ≈ vertical.kv[i][2]
        @test vertical.nlayers_kv[i] == 2
        @test Wflow.kh_layered_profile(vertical, 100.0, i, "layered_exponential") ≈
              33.76026208801769f0
        @test all(vertical.z_layered[1:10] .== 400.0)
        @test all(isnan.(vertical.z_exp))
    end

    model = Wflow.run_timestep(model)
    model = Wflow.run_timestep(model)
    @testset "river flow layered exponential profile" begin
        q = model.lateral.river.q_av
        @test sum(q) ≈ 3118.8690178033266f0
        @test q[1622] ≈ 0.000548447582354063f0
        @test q[43] ≈ 9.860543811678328f0
    end

    Wflow.close_files(model, delete_output = false)
end
