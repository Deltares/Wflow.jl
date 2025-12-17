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

@testitem "unit: Units" begin
    using Wflow: Unit, to_SI_factor, to_SI, to_string, ABSOLUTE_DEGREES

    @test to_SI_factor(ABSOLUTE_DEGREES) == 1.0
    @test to_SI(0.0, ABSOLUTE_DEGREES) == 273.15
    @test string(ABSOLUTE_DEGREES) == "°C"

    unit = Unit(; degC = 1) # relative degrees
    @test to_SI_factor(unit) == 1.0
    @test to_SI(0.0, unit) == 0.0

    dt = 86400.0
    unit = Unit(; m = 1, dt = -1)
    @test_throws Exception to_SI(unit, 1.0)
    @test to_SI_factor(unit; dt_val = dt) == inv(dt)
    @test string(unit) == "m dt⁻¹"
    @test to_string(unit; BMI_standard = true) == "m dt-1"

    unit = Unit(; s = 1, m = -1 // 3)
    @test to_SI_factor(unit) == 1.0
    @test string(unit) == "s m⁻¹ᐟ³"
    @test to_string(unit; BMI_standard = true) == "s m-1/3"

    @test_throws Exception Unit(; foo = 42)
end

@testitem "Lenses" begin
    configs = Wflow.Config[]

    # Initialize the first model with mass balance
    config = Wflow.Config(normpath(@__DIR__, "sbm_config.toml"))
    config.model.water_mass_balance__flag = true
    push!(configs, config)

    for file_name in [
        "sbm_gwf_config.toml",
        "sbm_river-floodplain-local-inertial_config.toml",
        "sbm_river-land-local-inertial_config.toml",
        "sbm_gwf_piave_demand_config.toml",
        "sediment_config.toml",
        "sediment_eurosem_engelund_config.toml",
    ]
        push!(configs, Wflow.Config(normpath(@__DIR__, file_name)))
    end

    for transport_method in ("kodatie", "govers", "yalin")
        config = Wflow.Config(normpath(@__DIR__, "sediment_eurosem_engelund_config.toml"))
        config.dir_output = normpath(@__DIR__, "data", "output", transport_method) # Avoid file permission problems
        if transport_method == "kodatie"
            config.model.river_transport = transport_method
        else
            config.model.run_river_model__flag = false
            config.model.land_transport = transport_method
        end
        push!(configs, config)
    end

    models = Wflow.Model.(configs)

    for (map_name, standard_name_map) in (
        ("sbm", Wflow.sbm_standard_name_map),
        ("sediment", Wflow.sediment_standard_name_map),
        ("domain", Wflow.domain_standard_name_map),
        ("routing", Wflow.routing_standard_name_map),
    )
        @testset "Test lenses: $map_name" begin
            invalid = String[]
            for (name, data) in standard_name_map
                (; lens) = data
                isnothing(lens) && continue
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
