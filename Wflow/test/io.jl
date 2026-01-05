@testitem "Configuration file" begin
    using TOML
    using Dates: DateTime

    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    parsed_toml = TOML.parsefile(tomlpath)
    config = Wflow.Config(tomlpath)

    @test parsed_toml isa Dict{String, Any}
    @test config isa Wflow.Config
    @test dirname(config) == dirname(tomlpath)

    # test if the values are parsed as expected
    @test config.time.starttime === DateTime(2000, 1, 1)
    @test config.time.endtime === DateTime(2000, 2)
    @test config.output.netcdf_grid.path == "output_moselle.nc"
    @test config.output.netcdf_scalar.path == "output_scalar_moselle.nc"
    @test config.output.csv.path == "output_moselle.csv"
    @test string(config) isa String

    # modifiers can also be applied
    var =
        config.input.static["soil_surface_water__vertical_saturated_hydraulic_conductivity"]
    @test var.netcdf_variable_name === "KsatVer"
    @test var.scale == [1.0]
    @test var.offset == [0.0]
    @test var.value === nothing
    @test var.layer === nothing

    # test the optional "dir_input" and "dir_output" keys
    @test config.dir_input != "."
    @test config.dir_output != "."
    @test Wflow.input_path(config, config.state.path_input) ==
          joinpath(@__DIR__, "data", "input", "instates-moselle.nc")
    @test Wflow.output_path(config, config.state.path_output) ==
          joinpath(@__DIR__, "data", "output", "outstates-moselle.nc")

    # test error is thrown for wrong non-optional model parameter
    @test_throws ErrorException Wflow.get_var(config, "not_set_in_TOML"; optional = false)
end

@testitem "Clock constructor" begin
    using NCDatasets: NCDataset
    using CFTime: DateTimeProlepticGregorian, DateTimeStandard
    using Dates: Second, Hour, Day, DateTime

    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    config.time.starttime = nothing
    config.time.endtime = nothing
    config.time.timestepsecs = nothing

    # mock a NCReader object
    ncpath = Wflow.input_path(config, config.input.path_forcing)
    ds = NCDataset(ncpath)
    reader = (; dataset = ds)

    clock = Wflow.Clock(config, reader)

    @test clock.time == DateTimeProlepticGregorian(2000, 1, 1)
    @test clock.iteration == 0
    @test clock.dt == Second(Day(1))
    # test that the missing keys have been added to the config
    @test config.time.starttime == DateTime(2000, 1, 1)
    @test config.time.endtime == DateTime(2001, 1, 1)
    @test config.time.timestepsecs == 86400

    # replace the keys with different values
    config = Wflow.Config(tomlpath)
    config.time.starttime = "2003-04-05"
    config.time.endtime = "2003-04-06"
    config.time.timestepsecs = 3600
    config.time.calendar = "standard"

    clock = Wflow.Clock(config, reader)
    @test clock.time == DateTimeStandard(2003, 4, 5)
    @test clock.iteration == 0
    @test clock.dt == Second(Hour(1))
end

@testitem "Clock{DateTimeStandard}" begin
    using CFTime: DateTimeStandard
    using Dates: Dates, Second, Day

    tomlpath = joinpath(@__DIR__, "sbm_config.toml")

    # 29 days in this February due to leap year
    starttime = DateTimeStandard(2000, 2, 28)
    dt = Day(1)
    clock = Wflow.Clock(starttime, 0, Second(dt))

    Wflow.advance!(clock)
    Wflow.advance!(clock)
    @test clock.time == DateTimeStandard(2000, 3, 1)
    @test clock.iteration == 2
    @test clock.dt == dt

    Wflow.rewind!(clock)
    @test clock.time == DateTimeStandard(2000, 2, 29)
    @test clock.iteration == 1
    @test clock.dt == dt

    config = Wflow.Config(tomlpath)
    config.time.starttime = starttime
    config.time.timestepsecs = Dates.value(Second(dt))

    Wflow.reset_clock!(clock, config)
    @test clock.time == starttime
    @test clock.iteration == 0
    @test clock.dt == dt
end

@testitem "Clock{DateTime360Day}" begin
    using CFTime: DateTime360Day
    using Dates: Dates, Second, Day

    tomlpath = joinpath(@__DIR__, "sbm_config.toml")

    # 30 days in each month
    starttime = DateTime360Day(2000, 2, 29)
    dt = Day(1)
    clock = Wflow.Clock(starttime, 0, Second(dt))

    Wflow.advance!(clock)
    Wflow.advance!(clock)
    @test clock.time == DateTime360Day(2000, 3, 1)
    @test clock.iteration == 2
    @test clock.dt == dt

    Wflow.rewind!(clock)
    @test clock.time == DateTime360Day(2000, 2, 30)
    @test clock.iteration == 1
    @test clock.dt == dt

    config = Wflow.Config(tomlpath)
    config.time.starttime = "2020-02-29"
    config.time.calendar = "360_day"
    config.time.timestepsecs = Dates.value(Second(dt))

    Wflow.reset_clock!(clock, config)
    @test clock.time isa DateTime360Day
    @test string(clock.time) == "2020-02-29T00:00:00"
    @test clock.iteration == 0
    @test clock.dt == dt
end

@testitem "Calendar noleap (DateTimeNoLeap) for time and clock" begin
    using CFTime: DateTimeNoLeap
    using Dates: Second, Day

    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    # test calendar setting `noleap` in forcing netCDF file (including `_FillValue` in time
    # dimension) and in TOML file (Clock{DateTimeNoLeap}).
    config = Wflow.Config(tomlpath)
    config.input.path_forcing = "forcing-calendar-noleap.nc"
    config.time.calendar = "noleap"

    # with `_FillValue` in time dimension Wflow throws a warning
    reader = @test_logs (
        :warn,
        "Time dimension contains `_FillValue` attribute, this is not in line with CF conventions.",
    ) match_mode = :any Wflow.NCReader(config)
    @test eltype(reader.dataset_times) == DateTimeNoLeap
    @test ismissing(reader.dataset_times) == false # missing in time dimension is not allowed
    @test reader.dataset_times ==
          collect(DateTimeNoLeap(2000, 1, 2):Day(1):DateTimeNoLeap(2000, 1, 6))

    # test Clock{DateTimeNoLeap}
    clock = Wflow.Clock(config, reader)
    @test clock.time isa DateTimeNoLeap
    @test clock.time == DateTimeNoLeap(2000, 1, 1)

    starttime = DateTimeNoLeap(2000, 2, 28)
    dt = Day(1)
    clock = Wflow.Clock(starttime, 0, Second(dt))
    Wflow.advance!(clock)
    @test clock.time == DateTimeNoLeap(2000, 3, 1)
end

@testitem "CFTime" begin
    using CFTime: DateTimeStandard, DateTimeProlepticGregorian, DateTime360Day
    using Dates: DateTime, Date
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

@testitem "timecycles" begin
    using Dates: Date, Day, monthday
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

@testitem "reducer" begin
    V = [6, 5, 4, 1]
    @test Wflow.function_map[Wflow.ReducerType.maximum](V) == 6
    @test Wflow.function_map[Wflow.ReducerType.minimum](V) == 1
    @test Wflow.function_map[Wflow.ReducerType.mean](V) == 4
    @test Wflow.function_map[Wflow.ReducerType.median](V) == 4.5
    @test Wflow.function_map[Wflow.ReducerType.first](V) == 6
    @test Wflow.function_map[Wflow.ReducerType.last](V) == 1
    @test Wflow.function_map[Wflow.ReducerType.sum](V) == 16
end

@testitem "State checking" begin
    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)

    # Extracting required states and test if some are covered (not all are tested!)
    required_states = Wflow.extract_required_states(config)
    @test "soil_water_saturated_zone__depth" in required_states
    @test "soil_layer_water_unsaturated_zone__depth" in required_states
    @test "vegetation_canopy_water__depth" in required_states
    @test "subsurface_water__volume_flow_rate" in required_states
    @test "river_water__instantaneous_volume_flow_rate" in required_states
    @test "reservoir_water_surface__elevation" in required_states
    @test "snowpack_liquid_water__depth" in required_states
    @test !("reservoir_water_level__elevation" in required_states)

    # Adding an unused state the see if the right warning message is thrown
    config.state.variables.additional_state = "additional_state"
    @test_logs (
        :warn,
        string(
            "State variable `additional_state` provided, but is not used in ",
            "model setup, skipping.",
        ),
    ) Wflow.check_states(config)

    # Removing the unused and required state, to test the exception being thrown
    delete!(config.state.variables, "additional_state")
    delete!(config.state.variables, "snowpack_dry_snow__leq_depth")
    @test_throws ArgumentError Wflow.check_states(config)

    # Extracting required states for model type sbm_gwf and test if some are covered
    tomlpath = joinpath(@__DIR__, "sbm_gwf_config.toml")
    config = Wflow.Config(tomlpath)
    required_states = Wflow.extract_required_states(config)
    @test "soil_water_saturated_zone__depth" in required_states
    @test "soil_layer_water_unsaturated_zone__depth" in required_states
    @test "vegetation_canopy_water__depth" in required_states
    @test "subsurface_water__hydraulic_head" in required_states
    @test "river_water__instantaneous_volume_flow_rate" in required_states
    @test "river_water__depth" in required_states
    @test "land_surface_water__depth" in required_states
end

@testitem "Model initialization" begin
    using NCDatasets: dimnames

    tomlpath = joinpath(@__DIR__, "sbm_config.toml")
    config = Wflow.Config(tomlpath)
    @testset "config settings" begin
        @test config.model.cold_start__flag
        config.model.cold_start__flag = false
        @test !config.model.cold_start__flag
        # test using an absolute path for the forcing
        @test !isabspath(config.input.path_forcing)
        abs_path_forcing = Wflow.input_path(config, config.input.path_forcing)
        config.input.path_forcing = abs_path_forcing
        @test isabspath(config.input.path_forcing)
    end

    model = Wflow.Model(config)
    Wflow.advance!(model.clock)
    Wflow.load_dynamic_input!(model)

    (; clock, reader, writer) = model

    @testset "output and state names" begin
        ncdims = ("lon", "lat", "layer", "time")
        @test dimnames(writer.dataset["ustorelayerdepth"]) == ncdims
        ncvars = [k for k in keys(writer.dataset) if !in(k, ncdims)]
        @test "snow" in ncvars
        @test "q_av_river" in ncvars
        @test "q_av_land" in ncvars
        @test length(writer.state_parameters) == 12
    end

    @testset "warm states" begin
        (; land) = model
        nt = Wflow.standard_name_map(model.land)
        lens = Wflow.get_lens("reservoir_water_surface__elevation", land)
        @test lens(model)[1] ≈ 3.6172022486284856
        lens = Wflow.get_lens("soil_water_saturated_zone__depth", land)
        @test lens(model)[9115] ≈ 477.13548089422125
        lens = Wflow.get_lens("snowpack_dry_snow__leq_depth", land)
        @test lens(model)[5] ≈ 11.019233179897599
        lens = Wflow.get_lens("soil_surface__temperature", land)
        @test lens(model)[5] ≈ 0.21814478119608938
        lens = Wflow.get_lens("soil_layer_water_unsaturated_zone__depth", land)
        @test lens(model)[50063][1] ≈ 9.969116007201725
        lens = Wflow.get_lens("snowpack_liquid_water__depth", land)
        @test lens(model)[5] ≈ 0.0
        lens = Wflow.get_lens("vegetation_canopy_water__depth", land)
        @test lens(model)[50063] ≈ 0.0
        lens = Wflow.get_lens("soil_water_saturated_zone__depth", land)
        @test lens(model)[50063] ≈ 558.8578304603327
        lens = Wflow.get_lens("subsurface_water__volume_flow_rate", land)
        @test lens(model)[10606] ≈ 39.972334552895816
        lens = Wflow.get_lens("river_water__instantaneous_volume_flow_rate", land)
        @test lens(model)[149] ≈ 53.48673634956338
        lens = Wflow.get_lens("river_water__depth", land)
        @test lens(model)[149] ≈ 1.167635369628945
        @test model.routing.river_flow.variables.storage[149] ≈ 63854.60119358985
        lens = Wflow.get_lens("land_surface_water__instantaneous_volume_flow_rate", land)
        @test lens(model)[2075] ≈ 3.285909284322251
        lens = Wflow.get_lens("land_surface_water__depth", land)
        @test lens(model)[2075] ≈ 0.052076262033771775
        @test model.routing.overland_flow.variables.storage[2075] ≈ 29920.754983235012
    end

    @testset "domain" begin
        (; domain) = model
        (; indices, reverse_indices) = domain.land.network
        # test if the reverse index reverses the index
        linear_index = 100
        cartesian_index = indices[linear_index]
        @test cartesian_index === CartesianIndex(168, 8)
        @test reverse_indices[cartesian_index] === linear_index
        # test active indices of different domains
        floodplain_inds = Wflow.active_indices(domain, "floodplain_water__volume")
        river_inds = Wflow.active_indices(domain, "river_water__volume")
        reservoir_inds = Wflow.active_indices(domain, "reservoir_water__volume")
        land_inds = Wflow.active_indices(domain, "soil_surface_water__runoff_volume_flux")
        @test floodplain_inds == river_inds
        @test length(land_inds) == 50063
        @test land_inds[100] == CartesianIndex(168, 8)
        @test reservoir_inds == [CartesianIndex(41, 134), CartesianIndex(60, 252)]
    end

    @testset "initial parameter values" begin
        (; land) = model
        @test land.snow.parameters.cfmax[1] ≈ 3.7565300464630127
        @test land.soil.parameters.soilthickness[1] ≈ 2000.0
        @test land.atmospheric_forcing.precipitation[49951] ≈ 2.2100000381469727
        @test land.soil.parameters.c[1] ≈
              [9.152995289601465, 8.919674421902961, 8.70537452585209, 8.690681062890977]
    end

    @testset "changing parameter values" begin
        config.input.static["snowpack__degree_day_coefficient"] = 2.0
        config.input.static["soil__thickness"] = Wflow.InputEntry(;
            scale = [3.0],
            offset = [100.0],
            netcdf_variable_name = "SoilThickness",
        )
        config.input.forcing["atmosphere_water__precipitation_volume_flux"] =
            Wflow.InputEntry(; scale = [1.5], netcdf_variable_name = "precip")
        config.input.static["soil_layer_water__brooks_corey_exponent"] = Wflow.InputEntry(;
            scale = [2.0, 3.0],
            offset = [0.0, 0.0],
            layer = [1, 3],
            netcdf_variable_name = "c",
        )

        model = Wflow.Model(config)
        Wflow.advance!(model.clock)
        Wflow.load_dynamic_input!(model)

        (; land) = model
        @test land.snow.parameters.cfmax[1] == 2.0
        @test land.soil.parameters.soilthickness[1] ≈ 2000.0 * 3.0 + 100.0
        @test land.atmospheric_forcing.precipitation[49951] ≈ 1.5 * 2.2100000381469727
        @test land.soil.parameters.c[1] ≈ [
            2.0 * 9.152995289601465,
            8.919674421902961,
            3.0 * 8.70537452585209,
            8.690681062890977,
        ]
    end
    Wflow.close_files(model; delete_output = false)
end

@testitem "NetCDF files and logging" begin
    using NCDatasets: NCDataset
    using LoggingExtras
    using TOML

    @testset "NetCDF creation" begin
        path = Base.Filesystem.tempname()
        _ = Wflow.create_tracked_netcdf(path)
        # safe to open the same path twice
        ds = Wflow.create_tracked_netcdf(path)
        close(ds)  # path is removed on process exit
    end

    @testset "NetCDF read variants" begin
        NCDataset(normpath("data/input/staticmaps-moselle.nc")) do ds
            @test Wflow.is_increasing(ds[:lon])
            @test !Wflow.is_increasing(ds[:lat])

            @test Wflow.nc_dim_name(ds, :time) == :time
            @test Wflow.nc_dim_name([:longitude], :x) == :longitude
            @test Wflow.nc_dim_name([:lat], :y) == :lat
            @test_throws ErrorException Wflow.nc_dim_name(ds, :not_present)

            x = collect(Wflow.nc_dim(ds, :lon))
            @test length(x) == 291
            @test x isa Vector{Union{Missing, Float64}}

            @test Wflow.internal_dim_name(:lon) == :x
            @test Wflow.internal_dim_name(:latitude) == :y
            @test Wflow.internal_dim_name(:time) == :time

            @test_throws ArgumentError Wflow.read_dims(ds["c"], (x = :, y = :))
            @test_throws ArgumentError Wflow.read_dims(ds["LAI"], (x = :, y = :))
            data, data_dim_order = Wflow.read_dims(ds["wflow_dem"], (x = :, y = :))
            @test data isa Matrix{Union{Float32, Missing}}
            @test data[end, end] === missing
            @test data[125, 1] ≈ 647.187
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

    # test logging and copy of TOML file to output
    @testset "Logging and copy TOML file" begin
        @test Wflow.convert_value(LogLevel, "InfO") == Logging.Info
        @test Wflow.convert_value(LogLevel, 0) == Logging.Info

        tomlpath = joinpath(@__DIR__, "sbm_simple.toml")
        Wflow.run(tomlpath; silent = true)

        config = Wflow.Config(tomlpath)
        output = normpath(abspath(abspath(config.dir_output)))
        toml_archive = Wflow.output_path(config, "sbm_simple.toml")
        path_log = Wflow.output_path(config, "log.txt")
        @test isfile(toml_archive)
        @test isfile(path_log)
        lines = readlines(path_log)
        @test count(startswith(line, "[ Info: ") for line in lines) == 58
        @test count(startswith(line, "┌ Debug: ") for line in lines) == 0

        # Another run with debug log level and a non-default path_log.
        # Must write the modified config to disk first, since the logging
        # applies only to the `Wflow.run(tomlpath)` method.
        # This also add an error to the config.
        tomlpath_debug = joinpath(@__DIR__, "sbm_simple-debug.toml")
        config.logging.loglevel = "debug"
        config.logging.path_log = "log-debug.txt"
        config.fews_run__flag = true
        config.logging.silent = true
        config.input.path_forcing = "doesnt-exist.nc"
        open(tomlpath_debug, "w") do io
            TOML.print(io, Wflow.to_dict(config))
        end
        @test_throws ErrorException Wflow.run(tomlpath_debug)
        rm(tomlpath_debug)
        path_log = Wflow.output_path(config, "log-debug.txt")
        @test isfile(path_log)
        lines = readlines(path_log)
        @test count(contains(line, " | Wflow | [Info] ") for line in lines) > 3
        @test count(contains(line, " | Wflow | [Debug] ") for line in lines) > 0
        msg = " | Wflow | [Error] Wflow simulation failed |exception = ErrorException(\"No files found with name 'doesnt-exist.nc' in '"
        @test contains(lines[end], msg)

        # Final run to test error handling during simulation
        tomlpath_error = joinpath(@__DIR__, "sbm_simple-error.toml")
        config.input.static["river__width"] = Wflow.InputEntry(;
            scale = [0.0],
            offset = [0.0],
            netcdf_variable_name = "wflow_riverwidth",
        )
        open(tomlpath_error, "w") do io
            TOML.print(io, Wflow.to_dict(config))
        end
        @test_throws ErrorException Wflow.run(tomlpath_error)
        rm(tomlpath_error)
    end
end
