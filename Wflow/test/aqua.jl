@testitem "Aqua" begin
    import Aqua
    Aqua.test_all(Wflow; ambiguities = false, persistent_tasks = false)
end