@testitem "unit: tosecond" begin
    using Dates:
        Month, Week, Day, Hour, Minute, Second, Millisecond, Microsecond, Nanosecond

    @test_throws MethodError Wflow.tosecond(Month(2))
    @test Wflow.tosecond(Week(2)) == 86400 * 14
    @test Wflow.tosecond(Day(2)) == 86400 * 2
    @test Wflow.tosecond(Hour(2)) == 3600 * 2
    @test Wflow.tosecond(Minute(2)) == 60 * 2
    @test Wflow.tosecond(Second(2)) == 2
    @test Wflow.tosecond(Millisecond(2)) == 2e-3
    @test Wflow.tosecond(Microsecond(2)) == 2e-6
    @test Wflow.tosecond(Nanosecond(2)) == 2e-9
end

@testitem "unit: julian_day (leap days are not counted)" begin
    using Dates: DateTime
    using CFTime
    @test Wflow.julian_day(DateTime(2000, 1, 1)) == 1
    @test Wflow.julian_day(DateTime(2000, 2, 28)) == 59
    @test Wflow.julian_day(DateTime(2000, 2, 29)) == 60
    @test Wflow.julian_day(DateTime(2000, 3, 1)) == 60
    @test Wflow.julian_day(DateTime(2000, 12, 31)) == 365
    @test Wflow.julian_day(DateTime(2001, 1, 1)) == 1
    @test Wflow.julian_day(DateTime(2001, 2, 28)) == 59
    @test Wflow.julian_day(DateTime(2001, 3, 1)) == 60
    @test Wflow.julian_day(DateTime(2001, 12, 31)) == 365
    @test Wflow.julian_day(CFTime.DateTimeStandard(2001, 12, 31)) == 365
    @test Wflow.julian_day(CFTime.DateTime360Day(2001, 12, 30)) == 360
end

@testitem "unit: bounded_divide" begin
    @test Wflow.bounded_divide(1.0, 0.0) == 0.0
    @test Wflow.bounded_divide(1.0, 0.0; default = 0.5) == 0.5
    @test Wflow.bounded_divide(1.0, 0.5) == 1.0
    @test Wflow.bounded_divide(1.0, 0.5; max = 0.75) == 0.75
    @test Wflow.bounded_divide(1.0, 2.0) == 0.5
end

@testitem "unit: scurve" begin
    a = 0.0
    b = 3.0
    c = 2.5

    x = 2.0
    out = Wflow.scurve(x, a, b, c)
    @test out ≈ 0.3325863502664285

    f = π
    @test f * Wflow.scurve(x, a + log(f) / c, f * b, c) ≈ out
end

@testitem "unit: compute_mass_balance_error" begin
    total_in = 5.0
    total_out = 5.0
    storage_rate = 0.0
    error, relative_error =
        Wflow.compute_mass_balance_error(total_in, total_out, storage_rate)
    @test iszero(error)
    @test iszero(relative_error)

    total_out = 6.0
    error, relative_error =
        Wflow.compute_mass_balance_error(total_in, total_out, storage_rate)
    @test error == -1.0
    @test relative_error ≈ -2 / 11
end

@testitem "Lenses" begin
    models = Wflow.Model[]
    # Initialize the first model with mass balance
    do_mass_balance = true
    for file_name in [
        "sbm_config.toml",
        "sbm_gwf_config.toml",
        "sbm_river-floodplain-local-inertial_config.toml",
        "sbm_river-land-local-inertial_config.toml",
        "sbm_gwf_piave_demand_config.toml",
        "sediment_config.toml",
        "sediment_eurosem_engelund_config.toml",
    ]
        config = Wflow.Config(normpath(@__DIR__, file_name))
        config.dir_output = mktempdir()
        if do_mass_balance
            config.model.water_mass_balance__flag = true
            global do_mass_balance = false
        end
        push!(models, Wflow.Model(config))
    end

    for (map_name, standard_name_map) in (
        ("sbm", Wflow.sbm_standard_name_map),
        ("sediment", Wflow.sediment_standard_name_map),
    )
        @testset "Test lenses: $map_name" begin
            dict_data = collect(standard_name_map)
            n = length(dict_data)
            invalids = zeros(Bool, n)
            Wflow.threaded_foreach(1:n; basesize = 25) do i
                (name, data) = dict_data[i]
                (; lens) = data
                invalid = true
                for model in models
                    try
                        lens(model)
                        invalid = false
                        break
                    catch
                        nothing
                    end
                end
                invalids[i] = invalid
            end
            invalid = [dict_data[i][1] for i in findall(invalids)]
            @test isempty(invalid)
        end
    end
end
