@testitem "Aqua" begin
    import Aqua
    Aqua.test_all(
        Wflow;
        persistent_tasks = false,
        piracies = (treat_as_own = [Wflow.SbmSoilVariables],),
    )
end
