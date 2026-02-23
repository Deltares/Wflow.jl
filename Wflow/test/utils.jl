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

@testitem "unit: variable tags" begin
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
    configs = Wflow.Config[]

    for file_name in [
        "sbm_gwf_config.toml",
        "sbm_river-floodplain-local-inertial_config.toml",
        "sbm_river-land-local-inertial_config.toml",
        "sbm_gwf_piave_demand_config.toml",
        "sediment_config.toml",
        "sediment_eurosem_engelund_config.toml",
    ]
        config = Wflow.Config(normpath(@__DIR__, file_name))
        config.dir_output = mktempdir()
        config.model.water_mass_balance__flag = true
        push!(configs, config)
    end

    for transport_method in ("kodatie", "govers", "yalin")
        config = Wflow.Config(normpath(@__DIR__, "sediment_eurosem_engelund_config.toml"))
        config.dir_output = mktempdir()
        if transport_method == "kodatie"
            config.model.river_transport = transport_method
        else
            config.model.run_river_model__flag = false
            config.model.land_transport = transport_method
        end
        push!(configs, config)
    end

    models = Wflow.Model.(configs)
    for (map_name, standard_name_map) in Wflow.standard_name_maps
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

@testitem "unit: water_table_change" begin
    using Wflow: to_SI, MM, Unit
    M_PER_DAY = Unit(; m = 1, d = -1)
    using StaticArrays: SVector
    include("testing_utils.jl")
    dt = 86400.0
    n = 1
    N = 5
    specific_yield = 0.43
    i = 1

    soil = init_sbm_soil_model(
        n,
        N;
        ustorelayerdepth = [to_SI.(SVector(100.0, 125.0, 150.0, 175.0, 200.0), Ref(MM))],
        ustorelayerthickness = [
            to_SI.(SVector(110.0, 145.0, 170.0, 205.0, 240.0), Ref(MM)),
        ],
        maxlayers = 5,
        theta_s = [0.98],
        theta_r = [0.02],
        n_unsatlayers = [5],
    )

    # Case: (net_flux ≤ 0.0)
    net_flux = to_SI(-1.0, M_PER_DAY)
    dh, exfilt = Wflow.water_table_change(soil, net_flux, specific_yield, i, dt)
    @test dh ≈ -2.3255813953488373
    @test iszero(exfilt)

    # Case: !(net_flux ≤ 0.0)
    # Hitting both (capacity <= net_flux) and !(capacity <= net_flux)
    net_flux = to_SI(0.05, M_PER_DAY)
    dh, exfilt = Wflow.water_table_change(soil, net_flux, specific_yield, i, dt)
    @test dh ≈ 0.4243119266055048
    @test iszero(exfilt)
end
