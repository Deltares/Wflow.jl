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
    using Accessors: PropertyLens
    using InteractiveUtils: subtypes

    get_fieldname(::PropertyLens{T}) where {T} = T

    function valid_lens(lens; type::Type = Wflow.Model, verbose::Bool = false)
        if isabstracttype(type) || (type isa Union)
            # If the type is abstract, search all subtypes
            sub_types = isabstracttype(type) ? subtypes(type) : Base.uniontypes(type)
            valid = false
            for subtype in sub_types
                if valid_lens(lens; type = subtype, verbose)
                    valid = true
                    break
                end
            end
            return valid
        else
            if lens isa ComposedFunction
                # If the lens is nested
                (; inner, outer) = lens
            else
                # If we are at a leaf
                inner = lens
                outer = nothing
            end

            # Find the field with the expected name
            fieldname = get_fieldname(inner)
            field_index = findfirst(==(fieldname), fieldnames(type))
            return if isnothing(field_index)
                # If the field does not exist, the lens is invalid
                if verbose
                    println("type $type has no field $fieldname.")
                end
                false
            else
                if isnothing(outer)
                    true
                else
                    # Recursion
                    field_types = fieldtypes(type)
                    if Any in field_types
                        error("$type has an unbound type parameter.")
                    end
                    valid_lens(outer; type = field_types[field_index], verbose)
                end
            end
        end
    end

    for (map_name, standard_name_map) in (
        ("sbm", Wflow.sbm_standard_name_map),
        ("sediment", Wflow.sediment_standard_name_map),
    )
        @testset "Test lenses: $map_name" begin
            invalid = String[]
            for (name, value) in standard_name_map
                (; lens) = value
                if !valid_lens(lens)
                    push!(invalid, name)
                    @error "Invalid lens" lens name
                end
            end
            @test isempty(invalid)
        end
    end
end
