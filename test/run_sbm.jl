using Dates

tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
@unpack network = model

Wflow.load_dynamic_input!(model)
model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2000-01-02T00:00:00")
    @test row.Q ≈ 9.88374809067884f0
    @test row.volume ≈ 3.298264381636329f7
    @test row.temp_bycoord ≈ 2.390000104904175f0
    @test row.temp_byindex ≈ 2.390000104904175f0
    @test row.Q_6336050 ≈ 0.0066621442402248245f0
    @test row.Q_6336510 ≈ 0.04727879981652821f0
    @test row.Q_6836100 ≈ 0.9339829956644901f0
    @test row.Q_6336500 ≈ 0.006084683242965108f0
    @test row.Q_6836190 ≈ 0.003038597498366973f0
    @test row.Q_6336800 ≈ 0.009077059837780043f0
    @test row.Q_6336900 ≈ 0.004167896211543442f0
    @test row.Q_6336930 ≈ 0.2352989751274897f0
    @test row.Q_6336910 ≈ 0.0076100268263675965f0
    @test row.Q_6136500 ≈ 0.0020577743577728067f0
    @test row.Q_6136520 ≈ 0.006701871174165091f0
    @test row.Q_6136150 ≈ 0.0077104586193382715f0
    @test row.Q_6136151 ≈ 0.011294064937339442f0
    @test row.Q_6136160 ≈ 4.793470200175769f0
    @test row.Q_6136202 ≈ 1.5231027963956818f0
    @test row.recharge_1 ≈ -0.05653226176238641f0
end

@testset "NetCDF scalar output" begin
    ds = model.writer.dataset_scalar
    @test ds["time"][1] == DateTime("2000-01-02T00:00:00")
    @test ds["Q"][:][1:20] ≈ [
        0.7561489f0,
        1.5231028f0,
        1.5432047f0,
        1.4594529f0,
        7.992092f0,
        3.2193222f0,
        4.8701606f0,
        4.2510085f0,
        0.011138658f0,
        4.7934704f0,
        4.9047155f0,
        0.0074881334f0,
        0.011294065f0,
        0.011618339f0,
        0.003813113f0,
        0.9464258f0,
        0.0025089371f0,
        1.4951783f0,
        5.3202457f0,
        2.2412596f0,
    ]
    @test ds["Q_gauges"].attrib["cf_role"] == "timeseries_id"
    @test ds["temp_index"][:] ≈ [2.39f0]
    @test ds["temp_coord"][:] ≈ [2.39f0]
    @test keys(ds.dim) == ["time", "Q_gauges", "temp_bycoord", "temp_byindex"]
end

@testset "first timestep" begin
    sbm = model.vertical

    @test sbm.tt[50063] ≈ 0.0f0

    @test model.clock.iteration == 2

    @test sbm.θₛ[50063] ≈ 0.48755401372909546f0
    @test sbm.runoff[50063] == 0.0
    @test sbm.soilevap[50063] == 0.0
    @test sbm.snow[5] ≈ 3.592840840467347f0
end

# run the second timestep
Wflow.load_dynamic_input!(model)
model = Wflow.update(model)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.θₛ[50063] ≈ 0.48755401372909546f0
    @test sbm.runoff[50063] == 0.0
    @test sbm.soilevap[50063] ≈ 0.006358004660566856f0
    @test sbm.snow[5] ≈ 3.667748983774868f0
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 6.370399148012509f7
    @test ssf[network.land.order[1]] ≈ 7.169036749244327f2
    @test ssf[network.land.order[end-100]] ≈ 2334.2186871917656f0
    @test ssf[network.land.order[end]] ≈ 288.19428729403944f0
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 290.2334461480497f0
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 1.0f-30
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 5072.751072834332f0
    @test q[1622] ≈ 0.0005484213345345085f0
    @test q[43] ≈ 12.544889508137787f0
    @test q[network.river.order[end]] ≈ 0.04295401190822074f0
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    @test res.outflow[1] ≈ 0.21750000119148086f0
    @test res.inflow[1] ≈ 45.53267637607559f0
    @test res.volume[1] ≈ 3.29605565479739f7
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

model = Wflow.run(config)

@testset "timing" begin
    # clock has been reset
    calendar = get(config, "calendar", "standard")::String
    @test model.clock.time == Wflow.cftime(config.starttime, calendar)
    @test model.clock.iteration == 1
end

@testset "river flow at basin outlets and downstream of one pit" begin
    q = model.lateral.river.q_av
    @test q[4009] ≈ 8.550171257466186f0 # pit/ outlet, CartesianIndex(141, 228)
    @test q[4020] ≈ 0.006513001839317472f0 # downstream of pit 4009, CartesianIndex(141, 229)
    @test q[2508] ≈ 155.82827037176472f0 # pit/ outlet
    @test q[5808] ≈ 0.18596231352540007f0 # pit/ outlet
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
Wflow.load_dynamic_input!(model)
model = Wflow.update(model)
Wflow.load_dynamic_input!(model)
model = Wflow.update(model)

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
Wflow.load_dynamic_input!(model)
model = Wflow.update(model)
Wflow.load_dynamic_input!(model)
model = Wflow.update(model)

@testset "river inflow (cyclic)" begin
    @test model.lateral.river.inflow[43] ≈ 0.75
    @test model.lateral.river.q_av[43] ≈ 8.08521095942811
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

Wflow.load_dynamic_input!(model)
model = Wflow.update(model)

@testset "fixed precipitation forcing (first timestep)" begin
    @test maximum(model.vertical.precipitation) ≈ 2.5
    @test minimum(model.vertical.precipitation) ≈ 0.0
    @test all(isapprox.(model.lateral.river.reservoir.precipitation, 2.5))
end

Wflow.close_files(model, delete_output = false)

# test local-inertial option for river flow river_routing
#tomlpath = joinpath(@__DIR__, "sbm_config.toml")
#config = Wflow.Config(tomlpath)
#config.model.river_routing = "local-inertial"

#model = Wflow.initialize_sbm_model(config)
#model = Wflow.update(model)
#model = Wflow.update(model)
