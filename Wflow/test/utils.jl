@testitem "tosecond" begin
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

@testitem "Julian day (leap days are not counted)" begin
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

@testitem "Bounded divide" begin
    @test Wflow.bounded_divide(1.0, 0.0) == 0.0
    @test Wflow.bounded_divide(1.0, 0.0; default = 0.5) == 0.5
    @test Wflow.bounded_divide(1.0, 0.5) == 1.0
    @test Wflow.bounded_divide(1.0, 0.5; max = 0.75) == 0.75
    @test Wflow.bounded_divide(1.0, 2.0) == 0.5
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
            invalid = String[]
            for (name, data) in standard_name_map
                (; lens) = data
                valid = false
                for model in models
                    try
                        lens(model)
                        valid = true
                        break
                    catch
                        nothing
                    end
                end
                valid || push!(invalid, name)
            end
            @test isempty(invalid)
        end
    end
end
