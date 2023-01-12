@testset "tosecond" begin
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

@testset "Julian day (leap days are not counted)" begin
    @test Wflow.julian_day(DateTime(2000,1,1)) == 1
    @test Wflow.julian_day(DateTime(2000,2,28)) == 59
    @test Wflow.julian_day(DateTime(2000,2,29)) == 60
    @test Wflow.julian_day(DateTime(2000,3,1)) == 60
    @test Wflow.julian_day(DateTime(2000,12,31)) == 365
    @test Wflow.julian_day(DateTime(2001,1,1)) == 1
    @test Wflow.julian_day(DateTime(2001,2,28)) == 59
    @test Wflow.julian_day(DateTime(2001,3,1)) == 60
    @test Wflow.julian_day(DateTime(2001,12,31)) == 365
end