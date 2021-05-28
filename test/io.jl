using NCDatasets
using Dates
using TOML
using CFTime
using Random
using UnPack

tomlpath = joinpath(@__DIR__, "sbm_config.toml")
parsed_toml = TOML.parsefile(tomlpath)
config = Wflow.Config(tomlpath)

@testset "configuration file" begin
    @test parsed_toml isa Dict{String,Any}
    @test config isa Wflow.Config
    @test Dict(config) == parsed_toml
    @test pathof(config) == tomlpath
    @test dirname(config) == dirname(tomlpath)

    # test if the values are parsed as expected
    @test config.starttime === DateTime(2000)
    @test config.endtime === DateTime(2000, 2)
    @test config.output.path == "data/output_moselle.nc"
    @test config.output isa Wflow.Config
    @test collect(keys(config.output)) == ["lateral", "vertical", "path"]
end

@testset "Clock constructor" begin
    config = Wflow.Config(tomlpath)

    # mock a NCReader object
    ncpath = joinpath(dirname(pathof(config)), config.input.path_forcing)
    ds = NCDataset(ncpath)
    reader = (; dataset = ds)

    # if these keys are missing, they are derived from the NetCDF
    pop!(Dict(config), "starttime")
    pop!(Dict(config), "endtime")
    pop!(Dict(config), "timestepsecs")
    clock = Wflow.Clock(config, reader)

    @test clock.time == DateTimeProlepticGregorian(2000, 1, 1)
    @test clock.iteration == 1
    @test clock.Δt == Second(Day(1))
    # test that the missing keys have been added to the config
    @test config.starttime == DateTime(2000, 1, 1)
    @test config.endtime == DateTime(2000, 12, 31)
    @test config.timestepsecs == 86400

    # replace the keys with different values
    config.starttime = "2003-04-05"
    config.endtime = "2003-04-06"
    config.timestepsecs = 3600
    config.calendar = "standard"

    clock = Wflow.Clock(config, reader)
    @test clock.time == DateTimeStandard(2003, 4, 5)
    @test clock.iteration == 1
    @test clock.Δt == Second(Hour(1))

    close(ds)
    config = Wflow.Config(tomlpath)  # restore the config
end

@testset "Clock{DateTimeStandard}" begin
    # 29 days in this February due to leap year
    starttime = DateTimeStandard(2000, 2, 28)
    Δt = Day(1)
    clock = Wflow.Clock(starttime, 1, Second(Δt))

    Wflow.advance!(clock)
    Wflow.advance!(clock)
    @test clock.time == DateTimeStandard(2000, 3, 1)
    @test clock.iteration == 3
    @test clock.Δt == Δt

    Wflow.rewind!(clock)
    @test clock.time == DateTimeStandard(2000, 2, 29)
    @test clock.iteration == 2
    @test clock.Δt == Δt

    config = Wflow.Config(
        Dict("starttime" => starttime, "timestepsecs" => Dates.value(Second(Δt))),
    )
    Wflow.reset_clock!(clock, config)
    @test clock.time == starttime
    @test clock.iteration == 1
    @test clock.Δt == Δt
end

@testset "Clock{DateTime360Day}" begin
    # 30 days in each month
    starttime = DateTime360Day(2000, 2, 29)
    Δt = Day(1)
    clock = Wflow.Clock(starttime, 1, Second(Δt))

    Wflow.advance!(clock)
    Wflow.advance!(clock)
    @test clock.time == DateTime360Day(2000, 3, 1)
    @test clock.iteration == 3
    @test clock.Δt == Δt

    Wflow.rewind!(clock)
    @test clock.time == DateTime360Day(2000, 2, 30)
    @test clock.iteration == 2
    @test clock.Δt == Δt

    config = Wflow.Config(
        Dict(
            "starttime" => "2020-02-29",
            "calendar" => "360_day",
            "timestepsecs" => Dates.value(Second(Δt)),
        ),
    )
    Wflow.reset_clock!(clock, config)
    @test clock.time isa DateTime360Day
    @test string(clock.time) == "2020-02-29T00:00:00"
    @test clock.iteration == 1
    @test clock.Δt == Δt
end

@testset "CFTime" begin
    @test Wflow.cftime("2006-01-02T15:04:05", "standard") ==
          DateTimeStandard(2006, 1, 2, 15, 4, 5)
    @test Wflow.cftime("2006-01-02", "proleptic_gregorian") ==
          DateTimeProlepticGregorian(2006, 1, 2)
    @test Wflow.cftime("2006-01-02T15:04:05", "360_day") ==
          DateTime360Day(2006, 1, 2, 15, 4, 5)
    @test Wflow.cftime(DateTime("2006-01-02T15:04:05"), "360_day") ==
          DateTime360Day(2006, 1, 2, 15, 4, 5)
    @test Wflow.cftime(Date("2006-01-02"), "360_day") == DateTime360Day(2006, 1, 2)
end

@testset "timecycles" begin
    @test Wflow.timecycles([Date(2020, 4, 21), Date(2020, 10, 21)]) == [(4, 21), (10, 21)]
    @test_throws ErrorException Wflow.timecycles([Date(2020, 4, 21), Date(2021, 10, 21)])
    @test_throws ErrorException Wflow.timecycles(collect(1:400))
    @test Wflow.timecycles(collect(1:12)) == collect(zip(1:12, fill(1, 12)))
    @test Wflow.timecycles(collect(1:366)) ==
          monthday.(Date(2000, 1, 1):Day(1):Date(2000, 12, 31))

    @test Wflow.monthday_passed((1, 1), (1, 1))  # same day
    @test Wflow.monthday_passed((1, 2), (1, 1))  # day later
    @test Wflow.monthday_passed((2, 1), (1, 1))  # month later
    @test Wflow.monthday_passed((2, 2), (1, 1))  # month and day later
    @test !Wflow.monthday_passed((2, 1), (2, 2))  # day before
    @test !Wflow.monthday_passed((1, 2), (2, 2))  # month before
    @test !Wflow.monthday_passed((1, 1), (2, 2))  # day and month before
end

# test reading and setting of warm states (reinit=false)
# modify existing config and initialize model with warm states
@test config.model.reinit
config["model"]["reinit"] = false
@test !config.model.reinit

model = Wflow.initialize_sbm_model(config)

@unpack vertical, clock, reader, writer = model

@testset "output and state names" begin
    ncdims = ("lon", "lat", "layer", "time")
    @test dimnames(writer.dataset["ustorelayerdepth"]) == ncdims
    ncvars = [k for k in keys(writer.dataset) if !in(k, ncdims)]
    @test "snow" in ncvars
    @test "q_river" in ncvars
    @test "q_land" in ncvars
    @test length(writer.state_parameters) == 14
end

@testset "warm states" begin
    @test Wflow.param(model, "lateral.river.reservoir.volume")[1] ≈ 2.7393418e7
    @test Wflow.param(model, "vertical.satwaterdepth")[9115] ≈ 201.51429748535156
    @test Wflow.param(model, "vertical.snow")[9115] ≈ 4.21874475479126
    @test Wflow.param(model, "vertical.tsoil")[9115] ≈ -1.9285825490951538
    @test Wflow.param(model, "vertical.ustorelayerdepth")[50069][1] ≈ 16.73013687133789
    @test Wflow.param(model, "vertical.snowwater")[9115] ≈ 0.42188167572021484
    @test Wflow.param(model, "vertical.canopystorage")[50069] ≈ 0.0
    @test Wflow.param(model, "vertical.zi")[50069] ≈ 380.67793082060416
    @test Wflow.param(model, "lateral.subsurface.ssf")[10606] ≈ 1.1614302208e11
    @test Wflow.param(model, "lateral.river.q")[149] ≈ 111.46229553222656
    @test Wflow.param(model, "lateral.river.h")[149] ≈ 8.555977821350098
    @test Wflow.param(model, "lateral.river.volume")[149] ≈ 249163.33754433296
    @test Wflow.param(model, "lateral.land.q")[10584] ≈ 1.0575231313705444
    @test Wflow.param(model, "lateral.land.h")[10584] ≈ 0.03456519544124603
    @test Wflow.param(model, "lateral.land.volume")[10584] ≈ 19379.23274404319
end

@testset "reducer" begin
    V = [6, 5, 4, 1]
    @test Wflow.reducerfunction("maximum")(V) == 6
    @test Wflow.reducerfunction("mean")(V) == 4
    @test Wflow.reducerfunction("median")(V) == 4.5
    @test Wflow.reducerfunction("first")(V) == 6
    @test Wflow.reducerfunction("last")(V) == 1
    @test_throws ErrorException Wflow.reducerfunction("other")
end

@testset "network" begin
    @unpack network = model
    @unpack indices, reverse_indices = model.network.land
    # test if the reverse index reverses the index
    linear_index = 100
    cartesian_index = indices[linear_index]
    @test cartesian_index === CartesianIndex(178, 7)
    @test reverse_indices[cartesian_index] === linear_index
end

@testset "initial parameter values" begin
    @unpack vertical = model
    @test vertical.cfmax[1] ≈ 3.7565300464630127
    @test vertical.soilthickness[1] ≈ 2000.0
    @test vertical.precipitation[49951] ≈ 2.197660446166992
end

config.input.vertical.cfmax = Dict("value" => 2.0)
config.input.vertical.soilthickness = Dict(
    "scale" => 3.0,
    "offset" => 100.0,
    "netcdf" => Dict("variable" => Dict("name" => "SoilThickness")),
)
config.input.vertical.precipitation =
    Dict("scale" => 1.5, "netcdf" => Dict("variable" => Dict("name" => "P")))

model = Wflow.initialize_sbm_model(config)

@testset "changed parameter values" begin
    @unpack vertical = model
    @test vertical.cfmax[1] == 2.0
    @test vertical.soilthickness[1] ≈ 2000.0 * 3.0 + 100.0
    @test vertical.precipitation[49951] ≈ 1.5 * 2.197660446166992
end

Wflow.close_files(model, delete_output = false)

@testset "NetCDF creation" begin
    path = Base.Filesystem.tempname()
    _ = Wflow.create_tracked_netcdf(path)
    # safe to open the same path twice
    ds = Wflow.create_tracked_netcdf(path)
    close(ds)  # path is removed on process exit
end

@testset "NetCDF read variants" begin
    NCDataset(staticmaps_moselle_path) do ds

        @test Wflow.is_increasing(ds[:x])
        @test !Wflow.is_increasing(ds[:y])

        @test Wflow.nc_dim_name(ds, :time) == :time
        @test Wflow.nc_dim_name([:longitude], :x) == :longitude
        @test Wflow.nc_dim_name([:lat], :y) == :lat
        @test_throws ErrorException Wflow.nc_dim_name(ds, :not_present)

        x = collect(Wflow.nc_dim(ds, :x))
        @test length(x) == 291
        @test x isa Vector{Union{Float64,Missing}}

        @test Wflow.internal_dim_name(:lon) == :x
        @test Wflow.internal_dim_name(:latitude) == :y
        @test Wflow.internal_dim_name(:time) == :time

        @test_throws ArgumentError Wflow.read_dims(ds["c"], (x = :, y = :))
        @test_throws ArgumentError Wflow.read_dims(ds["LAI"], (x = :, y = :))
        data, data_dim_order = Wflow.read_dims(ds["wflow_dem"], (x = :, y = :))
        @test data isa Matrix{Union{Float32,Missing}}
        @test data[end, end] === missing
        @test data[125, 1] ≈ 643.547f0
        @test data_dim_order == (:x, :y)

        @test Wflow.dim_directions(ds, (:x, :y)) === (x = true, y = false)
        @test Wflow.dim_directions(ds, (:y, :x, :layer)) ===
              (y = false, x = true, layer = true)

        data, dims = Wflow.permute_data(zeros(1, 2, 3), (:layer, :y, :x))
        @test size(data) == (3, 2, 1)
        @test dims == (:x, :y, :layer)
        data, dims = Wflow.permute_data(zeros(1, 2), (:x, :y))
        @test size(data) == (1, 2)
        @test dims == (:x, :y)
        @test_throws AssertionError size(Wflow.permute_data(zeros(1, 2, 3), (:x, :y)))

        data = collect(reshape(1:6, (2, 3)))
        # flip y, which is the second dimension
        @test Wflow.reverse_data!(data, (y = false, x = true))[1, :] == [5, 3, 1]
        # and mutate it back, the NamedTuple order should not matter
        @test Wflow.reverse_data!(data, (x = true, y = false))[1, :] == [1, 3, 5]
        # flip both dimensions at the same time
        data = Wflow.reverse_data!(data, (x = false, y = false))
        @test data[1, :] == [6, 4, 2]
        @test data[:, 1] == [6, 5]

        data = Wflow.read_standardized(ds, "wflow_dem", (x = :, y = :))
        # since in this case only the second dimension needs reversing, we can easily do it manually
        manual_fix = reverse(ds["wflow_dem"]; dims = 2)
        @test all(data .=== manual_fix)
    end
end
