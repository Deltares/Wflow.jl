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
    @test string(unit) == "m Δt⁻¹"
    @test to_string(unit; BMI_standard = true) == "m Δt-1"

    unit = Unit(; s = 1, m = -1 // 3)
    @test to_SI_factor(unit) == 1.0
    @test string(unit) == "s m⁻¹ᐟ³"
    @test to_string(unit; BMI_standard = true) == "s m-1/3"

    @test_throws Exception Unit(; foo = 42)
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

@testitem "Variable tags" begin
    for (map_name, map) in Wflow.standard_name_maps
        @testset "Check that each $map_name variable has at least one tag" begin
            vars_without_tags = String[]
            for (name, metadata) in map
                isempty(metadata.tags) && push!(vars_without_tags, name)
            end
            @test isempty(vars_without_tags)
        end
    end
end

@testitem "unit: lenses" begin
    using Accessors: @optic

    configs_sbm = Wflow.Config[]
    configs_sediment = Wflow.Config[]

    for file_name in [
        "sbm_gwf_config.toml",
        "sbm_river-floodplain-local-inertial_config.toml",
        "sbm_river-land-local-inertial_config.toml",
        "sbm_gwf_piave_demand_config.toml",
    ]
        config = Wflow.Config(normpath(@__DIR__, file_name))
        config.model.water_mass_balance__flag = true
        config.dir_output = mktempdir()
        push!(configs_sbm, config)
    end

    for transport_method in ("kodatie", "govers", "yalin", "bagnold")
        config = Wflow.Config(normpath(@__DIR__, "sediment_eurosem_engelund_config.toml"))
        config.dir_output = mktempdir()
        if transport_method in ("kodatie", "bagnold")
            config.model.river_transport = transport_method
            config.model.rainfall_erosion = "answers"
        else
            config.model.run_river_model__flag = false
            config.model.land_transport = transport_method
            config.model.rainfall_erosion = "eurosem"
        end
        push!(configs_sediment, config)
    end

    models_sbm = Wflow.Model.(configs_sbm)
    models_sediment = Wflow.Model.(configs_sediment)

    for (map_name, standard_name_map, models) in (
        ("sbm", Wflow.sbm_standard_name_map, models_sbm),
        ("sediment", Wflow.sediment_standard_name_map, models_sediment),
    )
        @testset "Test lenses: $map_name" begin
            invalids = String[]
            for (name, data) in standard_name_map
                (; lens) = data
                invalid = true
                isnothing(lens) && continue
                for model in models
                    try
                        lens(model)
                        invalid = false
                        break
                    catch
                        nothing
                    end
                end
                invalid && push!(invalids, name)
            end
            @test isempty(invalids)
        end
    end

    # Find duplicate lenses
    lenses = vcat(
        [
            getfield.(values(standard_name_map), :lens) for
            (_, standard_name_map) in Wflow.standard_name_maps
        ]...,
    )
    filter!(!isnothing, lenses)
    duplicates = Set()
    for unique_lens in unique(lenses)
        if count(==(unique_lens), lenses) > 1
            push!(duplicates, unique_lens)
        end
    end
    @test duplicates == Set([
        @optic(_.land.atmospheric_forcing.precipitation),
        @optic(_.land.glacier.variables.glacier_store)
    ])
end
