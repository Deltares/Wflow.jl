using Dates

tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

model = Wflow.initialize_sbm_model(config)
@unpack network = model

model = Wflow.update(model)

# test if the first timestep was written to the CSV file
flush(model.writer.csv_io)  # ensure the buffer is written fully to disk
@testset "CSV output" begin
    row = csv_first_row(model.writer.csv_path)

    @test row.time == DateTime("2000-01-01T00:00:00")
    @test row.Q ≈ 7.815587720999557f0
    @test row.volume ≈ 2.7783643457819197f7
    @test row.temp_bycoord ≈ 2.3279826641082764f0
    @test row.temp_byindex ≈ 2.3279826641082764f0
    @test row.Q_6336050 ≈ 0.02388490515704142f0
    @test row.Q_6336510 ≈ 0.012411069754330572f0
    @test row.Q_6836100 ≈ 0.004848609260578462f0
    @test row.Q_6336500 ≈ 0.011474801875419624f0
    @test row.Q_6836190 ≈ 0.0005262529349807199f0
    @test row.Q_6336800 ≈ 0.013467698526109831f0
    @test row.Q_6336900 ≈ 0.003408280699908697f0
    @test row.Q_6336930 ≈ 0.09773275295006544f0
    @test row.Q_6336910 ≈ 0.002147610203630289f0
    @test row.Q_6336920 ≈ 0.002649393436607859f0
    @test row.Q_6136100 ≈ 0.0008708128761381345f0
    @test row.Q_6136500 ≈ 0.000729148906480041f0
    @test row.Q_6136520 ≈ 0.002155395279410574f0
    @test row.Q_6136150 ≈ 0.0022298329229975987f0
    @test row.Q_6136151 ≈ 0.0031045027476524524f0
    @test row.Q_6136160 ≈ 3.3423894540713786f0
    @test row.Q_6136200 ≈ 1.358270076420503f0
    @test row.Q_6136201 ≈ 5.942330218662147f0
    @test row.Q_6136202 ≈ 1.680997341533662f0
    @test row.recharge_1 ≈ -0.027398093386017383f0
end

@testset "NetCDF scalar output" begin
    ds = model.writer.dataset_scalar
    @test ds["time"][1] == DateTime("2000-01-01T00:00:00")
    @test ds["Q"][:] ≈ [
        1.6809974f0,
        5.9423304f0,
        1.35827f0,
        3.3423896f0,
        0.0031045028f0,
        0.0022298328f0,
        0.0021553952f0,
        0.0007291489f0,
        0.00087081286f0,
        0.0026493934f0,
        0.0021476103f0,
        0.09773275f0,
        0.0034082807f0,
        0.013467698f0,
        0.00052625296f0,
        0.011474802f0,
        0.004848609f0,
        0.01241107f0,
        0.023884906f0,
    ]
    @test ds["Q_gauges"].attrib["cf_role"] == "timeseries_id"
    @test ds["temp_index"][:] ≈ [2.3279827f0]
    @test ds["temp_coord"][:] ≈ [2.3279827f0]
    @test keys(ds.dim) == ["time", "Q_gauges", "temp_bycoord", "temp_byindex"]
end

@testset "first timestep" begin
    sbm = model.vertical

    @test sbm.tt[50069] ≈ 1.2999999523162842f0

    @test model.clock.iteration == 2

    @test sbm.θₛ[50069] ≈ 0.48343977332115173f0
    @test sbm.runoff[50069] == 0.0
    @test sbm.soilevap[50069] == 0.0
    @test sbm.snow[50069] ≈ 0.6029989752244306f0
end

# run the second timestep
model = Wflow.update(model)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.θₛ[50069] ≈ 0.48343977332115173f0
    @test sbm.runoff[50069] == 0.0
    @test sbm.soilevap[50069] ≈ 0.005865651540305367f0
    @test sbm.snow[50069] ≈ 0.009696763863612956f0
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 6.368140761295825f16
    @test ssf[network.land.order[1]] ≈ 5.324295565183691f11
    @test ssf[network.land.order[end-100]] ≈ 7.855716879739626f11
    @test ssf[network.land.order[end]] ≈ 2.161246841709492f11
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 319.69538082148307f0
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 1.0978547520221912f-5
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 2807.884971209105f0
    @test q[1622] ≈ 0.0016288040314320486f0
    @test q[43] ≈ 7.338169165884175f0
    @test q[network.river.order[end]] ≈ 0.00610520650626283f0
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    @test res.outflow[1] ≈ 0.2174998592483153f0
    @test res.inflow[1] ≈ 50.170880189190626f0
    @test res.volume[1] ≈ 2.776162917050312f7
    @test res.precipitation[1] ≈ 0.1765228509902954f0
    @test res.evaporation[1] ≈ 0.5372688174247742f0
end

# set these variables for comparison in "changed dynamic parameters"
precip = copy(model.vertical.precipitation)
evap = copy(model.vertical.potential_evaporation)
lai = copy(model.vertical.leaf_area_index)
res_evap = copy(model.lateral.river.reservoir.evaporation)

benchmark = @benchmark Wflow.update(model)
trialmin = BenchmarkTools.minimum(benchmark)

println("SBM Model update")
print_benchmark(trialmin)
# @profview Wflow.update(model)
Wflow.close_files(model, delete_output = false)

# test for setting a pit
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

# pit is at river index 1765, CartesianIndex(141, 86)
# downstream from pit is at river index 1739, CartesianIndex(142, 85)
@testset "river flow at and downstream of pit" begin
    q = model.lateral.river.q_av
    @test q[3908] ≈ 8.133435566601927f0
    @test q[3920] ≈ 0.008996045895155833f0
end

# test changing forcing and cyclic LAI parameter
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

config.input.vertical.precipitation =
    Dict("scale" => 2.0, "netcdf" => Dict("variable" => Dict("name" => "P")))
config.input.vertical.potential_evaporation = Dict(
    "scale" => 3.0,
    "offset" => 1.50,
    "netcdf" => Dict("variable" => Dict("name" => "PET")),
)
config.input.vertical.leaf_area_index =
    Dict("scale" => 1.6, "netcdf" => Dict("variable" => Dict("name" => "LAI")))

model = Wflow.initialize_sbm_model(config)
model = Wflow.update(model)
model = Wflow.update(model)

@testset "changed dynamic parameters" begin
    res = model.lateral.river.reservoir
    vertical = model.vertical
    @test vertical.precipitation[2] / precip[2] ≈ 2.0f0
    @test (vertical.potential_evaporation[100] - 1.50) / evap[100] ≈ 3.0f0
    @test vertical.leaf_area_index[100] / lai[100] ≈ 1.6f0
    @test (res.evaporation[2] - 1.50) / res_evap[2] ≈ 3.0000012203408635f0
end
Wflow.close_files(model, delete_output = false)
