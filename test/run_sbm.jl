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
    @test row.recharge_1 ≈ -0.002257181032501202f0
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
    sbm = model.vertical

    @test sbm.tt[50063] ≈ 0.0f0

    @test model.clock.iteration == 1

    @test sbm.θₛ[50063] ≈ 0.48755401372909546f0
    @test sbm.θᵣ[50063] ≈ 0.15943120419979095f0
    @test sbm.runoff[50063] == 0.0
    @test sbm.soilevap[50063] ≈ 0.011096671644901038f0
    @test sbm.snow[5] ≈ 3.7685133905888146f0
    @test sbm.total_storage[50063] ≈ 559.9035608052374f0
    @test sbm.total_storage[429] ≈ 597.4578475404879f0 # river cell
end

# run the second timestep
model = Wflow.run_timestep(model)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.θₛ[50063] ≈ 0.48755401372909546f0
    @test sbm.θᵣ[50063] ≈ 0.15943120419979095f0
    @test sbm.runoff[50063] == 0.0
    @test sbm.soilevap[50063] ≈ 0.008718333439094138f0
    @test sbm.snow[5] ≈ 3.8434125240523125f0
    @test sbm.total_storage[50063] ≈ 560.0152135062889f0
    @test sbm.total_storage[429] ≈ 617.2238533241972f0 # river cell
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 6.370399148012509f7
    @test ssf[network.land.order[1]] ≈ 7.169036749244327f2
    @test ssf[network.land.order[end-100]] ≈ 2333.801056570759f0
    @test ssf[network.land.order[end]] ≈ 288.19428729403944f0
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 291.27639107427285f0
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 1.0f-30
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 3622.7369292570543f0
    @test q[1622] ≈ 0.0006497468064774366f0
    @test q[43] ≈ 12.05767242907667f0
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
    @test q[4009] ≈ 8.51907041734622f0 # pit/ outlet, CartesianIndex(141, 228)
    @test q[4020] ≈ 0.006779014715290862f0 # downstream of pit 4009, CartesianIndex(141, 229)
    @test q[2508] ≈ 150.28398167251638f0 # pit/ outlet
    @test q[5808] ≈ 0.12419895007970105f0 # pit/ outlet
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
    @test model.lateral.river.q_av[44] ≈ 10.71846407068599
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
    @test sum(q) ≈ 3919.6025219496014f0
    @test q[1622] ≈ 7.31010246736994f-5
    @test q[43] ≈ 11.92153120707289f0
    @test q[501] ≈ 3.5736389982451895f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.001987887644883981f0
    @test h[43] ≈ 0.4366415244811759f0
    @test h[501] ≈ 0.057265962518284294f0
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
    @test sum(q) ≈ 2380.64389229669f0
    @test q[1622] ≈ 7.322956970529551f-5
    @test q[43] ≈ 5.361283165612762f0
    @test q[501] ≈ 1.6021771576366957f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.0019891342000364796f0
    @test h[43] ≈ 0.3003008810153667f0
    @test h[501] ≈ 0.031925992442532f0
    qx = model.lateral.land.qx
    qy = model.lateral.land.qy
    @test qx[[26, 35, 631]] ≈ [0.18614776104106373f0, 0.029502872625766417f0, 0.0f0]
    @test qy[[26, 35, 631]] ≈ [0.12757214437549858f0, 1.7212079599401755f0, 0.0f0]
    h = model.lateral.land.h
    @test h[[26, 35, 631]] ≈
          [0.07361854999908582f0, 0.009155393111676267f0, 0.0007258741013439351f0]
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
    @test sum(q) ≈ 3908.039208613999f0
    @test q[1622] ≈ 7.310102468091527f-5
    @test q[43] ≈ 11.921531207072922f0
    @test q[501] ≈ 3.5061516913374717f0
    h = model.lateral.river.h_av
    @test h[1622] ≈ 0.001987887580593841f0
    @test h[43] ≈ 0.436641524481545f0
    @test h[501] ≈ 0.05665942153713204f0
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
        @test sum(q) ≈ 3126.3509477318844f0
        @test q[1622] ≈ 0.0005972577112819149f0
        @test q[43] ≈ 9.880641908157857f0
    end

    Wflow.close_files(model, delete_output = false)
end
