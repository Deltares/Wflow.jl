@testitem "Aqua" begin
    import Aqua
    import SharedHydrology
    Aqua.test_all(
        WflowServer;
        persistent_tasks = false,
        piracies = (treat_as_own = [SharedHydrology],),
    )
end
