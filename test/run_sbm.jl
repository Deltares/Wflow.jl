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
    @test row.Q ≈ 7.815587720999557
    @test row.volume ≈ 4.364251626782245e7
    @test row.temp_bycoord ≈ 2.3279826641082764
    @test row.temp_byindex ≈ 2.3279826641082764
    @test row.Q_6336050 ≈ 0.02388490515704142
    @test row.Q_6336510 ≈ 0.012411069754330572
    @test row.Q_6836100 ≈ 0.004848609260578462
    @test row.Q_6336500 ≈ 0.011474801875419624
    @test row.Q_6836190 ≈ 0.0005262529349807199
    @test row.Q_6336800 ≈ 0.013467698526109831
    @test row.Q_6336900 ≈ 0.003408280699908697
    @test row.Q_6336930 ≈ 0.09773275295006544
    @test row.Q_6336910 ≈ 0.002147610203630289
    @test row.Q_6336920 ≈ 0.002649393436607859
    @test row.Q_6136100 ≈ 0.0008708128761381345
    @test row.Q_6136500 ≈ 0.000729148906480041
    @test row.Q_6136520 ≈ 0.002155395279410574
    @test row.Q_6136150 ≈ 0.0022298329229975987
    @test row.Q_6136151 ≈ 0.0031045027476524524
    @test row.Q_6136160 ≈ 3.3423894540713786
    @test row.Q_6136200 ≈ 1.358270076420503
    @test row.Q_6136201 ≈ 5.942330218662147
    @test row.Q_6136202 ≈ 1.680997341533662
    @test row.recharge_1 ≈ -0.027398093386017383
end

@testset "NetCDF scalar output" begin
    ds = model.writer.dataset_scalar
    @test ds["time"][1] == DateTime("2000-01-01T00:00:00")
    @test ds["Q"][:] ≈ [
        0.023884907,
        0.012411069,
        0.004848609,
        0.011474802,
        0.0005262529,
        0.013467698,
        0.003408281,
        0.09773275,
        0.0021476103,
        0.0026493936,
        0.0008708129,
        0.0007291489,
        0.0021553952,
        0.002229833,
        0.003104503,
        3.3423893,
        1.35827,
        5.9423304,
        1.6809973,
    ]
    @test ds["Q_gauges"].attrib["cf_role"] == "timeseries_id"
    @test ds["temp_index"][:] ≈ [2.3279827]
    @test ds["temp_coord"][:] ≈ [2.3279827]
    @test keys(ds.dim) == ["time", "Q_gauges", "temp_bycoord", "temp_byindex"]
end

@testset "first timestep" begin
    sbm = model.vertical

    @test sbm.tt[1] ≈ 1.2999999523162842

    @test model.clock.iteration == 2

    @test sbm.θₛ[1] ≈ 0.48343977332115173
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] == 0.0
    @test sbm.snow[1] ≈ 0.6029989752244306
end

# run the second timestep
model = Wflow.update(model)

@testset "second timestep" begin
    sbm = model.vertical
    @test sbm.θₛ[1] ≈ 0.48343977332115173
    @test sbm.runoff[1] == 0.0
    @test sbm.soilevap[1] ≈ 0.005865651540305367
    @test sbm.snow[1] ≈ 0.009696763863612956
end

@testset "subsurface flow" begin
    ssf = model.lateral.subsurface.ssf
    @test sum(ssf) ≈ 6.368140761295825e16
    @test ssf[network.land.order[1]] ≈ 3.0449782003445332e13
    @test ssf[network.land.order[end-100]] ≈ 7.855716879739626e11
    @test ssf[network.land.order[end]] ≈ 2.161246841709492e11
end

@testset "overland flow" begin
    q = model.lateral.land.q_av
    @test sum(q) ≈ 319.69538082148307
    @test q[26625] ≈ 0.0
    @test q[39308] ≈ 0.0
    @test q[network.land.order[end]] ≈ 1.0978547520221912e-5
end

@testset "river flow" begin
    q = model.lateral.river.q_av
    @test sum(q) ≈ 2807.884971209105
    @test q[4061] ≈ 0.0016288040314320486
    @test q[5617] ≈ 7.338169165884175
    @test q[network.river.order[end]] ≈ 0.00610520650626283
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    @test res.outflow[2] ≈ 0.2174998592483153
    @test res.inflow[2] ≈ 50.170880189190626
    @test res.volume[2] ≈ 2.776162917050312e7
    @test res.precipitation[2] ≈ 0.1765228509902954
    @test res.evaporation[2] ≈ 0.5372688174247742
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
    @test q[1765] ≈ 8.133435566601927
    @test q[1739] ≈ 0.008996045895155833
end

# test changing forcing and cyclic LAI parameter
tomlpath = joinpath(@__DIR__, "sbm_config.toml")
config = Wflow.Config(tomlpath)

Dict(config)["input"]["vertical"]["precipitation"] =
    Dict("scale" => 2.0, "netcdf" => Dict("variable" => Dict("name" => "P")))
Dict(config)["input"]["vertical"]["potential_evaporation"] = Dict(
    "scale" => 3.0,
    "offset" => 1.50,
    "netcdf" => Dict("variable" => Dict("name" => "PET")),
)
Dict(config)["input"]["vertical"]["leaf_area_index"] =
    Dict("scale" => 1.6, "netcdf" => Dict("variable" => Dict("name" => "LAI")))

model = Wflow.initialize_sbm_model(config)
model = Wflow.update(model)
model = Wflow.update(model)

@testset "changed dynamic parameters" begin
    res = model.lateral.river.reservoir
    vertical = model.vertical
    @test vertical.precipitation[2] / precip[2] ≈ 2.0
    @test (vertical.potential_evaporation[100] - 1.50) / evap[100] ≈ 3.0
    @test vertical.leaf_area_index[100] / lai[100] ≈ 1.6
    @test (res.evaporation[2] - 1.50) / res_evap[2] ≈ 3.0000012203408635
end
Wflow.close_files(model, delete_output = false)
