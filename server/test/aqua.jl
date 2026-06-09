@testitem "Aqua" begin
    import Aqua
    import SharedHydrology
    Aqua.test_all(WflowServer; persistent_tasks = false)
end
