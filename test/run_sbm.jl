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

    @test row.time == DateTime("2000-01-01T00:00:00")
    @test row.Q ≈ 7.680739298259266f0
    @test row.volume ≈ 2.7783643457819197f7
    @test row.temp_bycoord ≈ 2.3279826641082764f0
    @test row.temp_byindex ≈ 2.3279826641082764f0
    @test row.Q_6336050 ≈ 0.02314152846188506f0
    @test row.Q_6336510 ≈ 0.011829340015725808f0
    @test row.Q_6836100 ≈ 0.004285433507347682f0
    @test row.Q_6336500 ≈ 0.011434655691502343f0
    @test row.Q_6836190 ≈ 0.0004903061800226998f0
    @test row.Q_6336800 ≈ 0.013376156569068251f0
    @test row.Q_6336900 ≈ 0.0033244344482295423f0
    @test row.Q_6336930 ≈ 0.08966505345584352f0
    @test row.Q_6336910 ≈ 0.0020127260170180474f0
    @test row.Q_6336920 ≈ 0.002580440335740291f0
    @test row.Q_6136100 ≈ 0.0008485772628065181f0
    @test row.Q_6136500 ≈ 0.0006749692699499864f0
    @test row.Q_6136520 ≈ 0.0019351916317406491f0
    @test row.Q_6136150 ≈ 0.002126215608742851f0
    @test row.Q_6136151 ≈ 0.0028804236661004205f0
    @test row.Q_6136160 ≈ 3.279985818360295f0
    @test row.Q_6136200 ≈ 1.259358422161174f0
    @test row.Q_6136201 ≈ 5.803861903796768f0
    @test row.Q_6136202 ≈ 1.6624100470042156f0
    @test row.recharge_1 ≈ -0.027398093386017383f0
end

@testset "NetCDF scalar output" begin
    ds = model.writer.dataset_scalar
    @test ds["time"][1] == DateTime("2000-01-01T00:00:00")
    @test ds["Q"][:] ≈ [
        1.66241f0,
        5.803862f0,
        1.2593584f0,
        3.279986f0,
        0.0028804236f0,
        0.0021262157f0,
        0.0019351917f0,
        0.0006749693f0,
        0.00084857724f0,
        0.0025804404f0,
        0.002012726f0,
        0.089665055f0,
        0.0033244344f0,
        0.013376157f0,
        0.0004903062f0,
        0.011434656f0,
        0.0042854333f0,
        0.01182934f0,
        0.023141528f0,
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
Wflow.load_dynamic_input!(model)
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
    @test sum(ssf) ≈ 6.368140761295825f7
    @test ssf[network.land.order[1]] ≈ 5.324295565183691f2
    @test ssf[network.land.order[end-100]] ≈ 7.855716879739626f2
    @test ssf[network.land.order[end]] ≈ 2.161246841709492f2
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
    @test sum(q) ≈ 2805.266510229795f0
    @test q[1622] ≈ 0.0015498389046670307f0
    @test q[43] ≈ 7.335279907301223f0
    @test q[network.river.order[end]] ≈ 0.006068754125930134f0
end

@testset "reservoir simple" begin
    res = model.lateral.river.reservoir
    @test res.outflow[1] ≈ 0.2174998592483153f0
    @test res.inflow[1] ≈ 49.52026446039592f0
    @test res.volume[1] ≈ 2.776162917050312f7
    @test res.precipitation[1] ≈ 0.1765228509902954f0
    @test res.evaporation[1] ≈ 0.5372688174247742f0
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
    @test q[3908] ≈ 8.107901948181988f0 # pit/ outlet
    @test q[3920] ≈ 0.008992215117496462f0 # downstream of pit 3908
    @test q[2474] ≈ 156.71758046363857f0 # pit/ outlet
    @test q[5650] ≈ 0.10849551740603759f0 # pit/ outlet
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
